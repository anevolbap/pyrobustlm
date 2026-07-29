# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Cython kernel for one fast-S resampling iteration. Family-generic:
# bisquare, hampel, optimal, lqq, and ggw dispatch on a small family enum.
# ggw uses the same polynomial chi tables as _psi.pyx (duplicated here so
# the kernel stays self-contained and nogil).

cimport cython
from cpython.pycapsule cimport PyCapsule_GetPointer
from libc.math cimport fabs, sqrt, exp, pow as cpow
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free

from numpy.random cimport bitgen_t
from scipy.linalg.cython_lapack cimport dgesv, dgels

import numpy as np
cimport numpy as cnp

cnp.import_array()
# ---------------------------------------------------------------------------
# Family enum, ggw tables, chi/wgt/M-scale kernels and the bounded RNG draw
# live in _psi_kernels.pxi, shared verbatim with the other engine. See the
# header of that file for why the two modules used to carry separate copies.
# ---------------------------------------------------------------------------
include "_psi_kernels.pxi"


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_resample_iter(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[cnp.int64_t, ndim=1, mode="c"] idx,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
    double b0,
    int k_fast_s,
    int max_iter_scale,
    double scale_tol,
    cnp.ndarray[double, ndim=1, mode="c"] beta_out,
):
    """Run one fast-S resampling iteration with the given psi family.

    Parameters
    ----------
    X, y : design and response (C-contiguous).
    idx : p row indices selected for the initial p-subset.
    family : 0=bisquare, 1=hampel, 2=optimal, 3=lqq.
    tuning : up to 3 family-specific tuning constants (k, k1+k2+k3, etc.).
    b0, k_fast_s, max_iter_scale, scale_tol : as in FastSConfig.
    beta_out : preallocated (p,) buffer for the returned coefficient vector.

    Returns
    -------
    (scale, status) where status is 0 (ok), 1 (singular subset), 2 (exact
    fit), 3 (LAPACK error during dgels).
    """
    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    # Initialise ggw tables once (with the GIL). Cheap if already done.
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef double* sub_X = <double*>malloc(p * p * sizeof(double))
    cdef double* sub_y = <double*>malloc(p * sizeof(double))
    cdef double* X_w = <double*>malloc(n * p * sizeof(double))
    cdef double* y_w = <double*>malloc(n * sizeof(double))
    cdef double* r = <double*>malloc(n * sizeof(double))
    cdef double* w = <double*>malloc(n * sizeof(double))
    cdef int* ipiv = <int*>malloc(p * sizeof(int))

    cdef int lwork = max(1, n_int * p_int + 64 * (n_int + p_int))
    cdef double* work = <double*>malloc(lwork * sizeof(double))

    if (sub_X == NULL or sub_y == NULL or X_w == NULL or y_w == NULL
            or r == NULL or w == NULL or ipiv == NULL or work == NULL):
        free(sub_X); free(sub_y); free(X_w); free(y_w)
        free(r); free(w); free(ipiv); free(work)
        raise MemoryError("cy_resample_iter: out of memory")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef cnp.int64_t* idx_data = <cnp.int64_t*>cnp.PyArray_DATA(idx)
    cdef double* beta_data = <double*>cnp.PyArray_DATA(beta_out)
    cdef double* tuning_data = <double*>cnp.PyArray_DATA(tuning)

    cdef Py_ssize_t i, j, kk, row
    cdef double s, sw, dot
    cdef int status = 0
    cdef double scale = 0.0

    # First-tuning value is the "k" for one-parameter families. For
    # hampel/lqq this is just tuning[0]; the others are tuning[1] and
    # tuning[2].
    cdef double k0 = tuning_data[0]

    with nogil:
        for i in range(p):
            row = idx_data[i]
            sub_y[i] = y_data[row]
            for j in range(p):
                sub_X[i + j * p] = X_data[row * p + j]

        dgesv(&p_int, &one, sub_X, &p_int, ipiv, sub_y, &p_int, &info)
        if info != 0:
            status = 1
        else:
            for i in range(p):
                beta_data[i] = sub_y[i]

            for i in range(n):
                dot = 0.0
                for j in range(p):
                    dot += X_data[i * p + j] * beta_data[j]
                r[i] = y_data[i] - dot

            # Initial scale: max |r| / k0 (a crude starting value; the M-scale
            # iteration converges from anywhere positive).
            s = 0.0
            for i in range(n):
                if r[i] >= 0:
                    if r[i] > s:
                        s = r[i]
                else:
                    if -r[i] > s:
                        s = -r[i]
            if s <= 0.0:
                status = 2
                scale = 0.0
            else:
                s = s / k0 if k0 > 0 else s
                if s <= 0.0:
                    s = 1.0

                for kk in range(k_fast_s):
                    s = _mscale_generic(r, n, s, p, family, tuning_data,
                                        b0, max_iter_scale, scale_tol)
                    if s == 0.0:
                        status = 2
                        scale = 0.0
                        break

                    _wgt_zinv(r, w, n, s, family, tuning_data)
                    for i in range(n):
                        sw = sqrt(w[i]) if w[i] > 0 else 0.0
                        y_w[i] = y_data[i] * sw
                        for j in range(p):
                            X_w[i + j * n] = X_data[i * p + j] * sw

                    dgels(b'N', &n_int, &p_int, &one,
                          X_w, &n_int,
                          y_w, &n_int,
                          work, &lwork, &info)
                    if info != 0:
                        status = 3
                        break
                    for i in range(p):
                        beta_data[i] = y_w[i]
                    for i in range(n):
                        dot = 0.0
                        for j in range(p):
                            dot += X_data[i * p + j] * beta_data[j]
                        r[i] = y_data[i] - dot
                else:
                    s = _mscale_generic(r, n, s, p, family, tuning_data,
                                        b0, max_iter_scale, scale_tol)
                    scale = s

    free(sub_X); free(sub_y); free(X_w); free(y_w)
    free(r); free(w); free(ipiv); free(work)

    return scale, status


# ---------------------------------------------------------------------------
# Survivor refinement to convergence. Same iteration as cy_resample_iter
# (m_scale -> IRWLS -> m_scale -> ...) but with a relative-tolerance
# stopping criterion on beta, used after the resampling loop picks the
# best candidates.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_refine_to_convergence(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[double, ndim=1, mode="c"] beta,
    double sigma_init,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
    double b0,
    int max_it,
    double refine_tol,
    int max_iter_scale,
    double scale_tol,
):
    """Iterate (m_scale, IRWLS) on a single survivor candidate until
    ``||beta_{k+1} - beta_k|| / ||beta_k||`` falls below ``refine_tol``
    or ``max_it`` is exhausted.

    The ``beta`` buffer is updated in place and also returned. Returns
    ``(scale, converged_int, n_iter, status)`` where ``status`` is 0 on
    success, 2 on exact fit, 3 on LAPACK error.
    """
    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef double* X_w = <double*>malloc(n * p * sizeof(double))
    cdef double* y_w = <double*>malloc(n * sizeof(double))
    cdef double* r = <double*>malloc(n * sizeof(double))
    cdef double* w = <double*>malloc(n * sizeof(double))
    cdef double* beta_prev = <double*>malloc(p * sizeof(double))
    cdef int lwork = max(1, n_int * p_int + 64 * (n_int + p_int))
    cdef double* work = <double*>malloc(lwork * sizeof(double))

    if (X_w == NULL or y_w == NULL or r == NULL or w == NULL
            or beta_prev == NULL or work == NULL):
        free(X_w); free(y_w); free(r); free(w); free(beta_prev); free(work)
        raise MemoryError("cy_refine_to_convergence: out of memory")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef double* beta_data = <double*>cnp.PyArray_DATA(beta)
    cdef double* tuning_data = <double*>cnp.PyArray_DATA(tuning)

    cdef Py_ssize_t i, j, it
    cdef double s = sigma_init
    cdef double sw, dot, delta, denom, diff
    cdef double scale = 0.0
    cdef int converged = 0
    cdef int status = 0
    cdef int iters = 0

    with nogil:
        for it in range(max_it):
            iters = it + 1
            # residuals from current beta
            for i in range(n):
                dot = 0.0
                for j in range(p):
                    dot += X_data[i * p + j] * beta_data[j]
                r[i] = y_data[i] - dot

            s = _mscale_generic(r, n, s, p, family, tuning_data,
                                b0, max_iter_scale, scale_tol)
            if s == 0.0:
                status = 2
                scale = 0.0
                converged = 1
                break

            # IRWLS weights
            _wgt_zinv(r, w, n, s, family, tuning_data)
            for i in range(n):
                sw = sqrt(w[i]) if w[i] > 0 else 0.0
                y_w[i] = y_data[i] * sw
                for j in range(p):
                    X_w[i + j * n] = X_data[i * p + j] * sw

            # Save previous beta for the convergence check.
            for j in range(p):
                beta_prev[j] = beta_data[j]

            dgels(b'N', &n_int, &p_int, &one,
                  X_w, &n_int,
                  y_w, &n_int,
                  work, &lwork, &info)
            if info != 0:
                status = 3
                scale = s
                break
            for j in range(p):
                beta_data[j] = y_w[j]

            # Relative L2 change on beta.
            delta = 0.0
            denom = 0.0
            for j in range(p):
                diff = beta_data[j] - beta_prev[j]
                delta += diff * diff
                denom += beta_prev[j] * beta_prev[j]
            delta = sqrt(delta)
            if denom < 1e-300:
                denom = 1e-300
            else:
                denom = sqrt(denom)
            if delta / denom < refine_tol:
                converged = 1
                scale = s
                break
            scale = s

    free(X_w); free(y_w); free(r); free(w); free(beta_prev); free(work)
    return scale, converged, iters, status


# ---------------------------------------------------------------------------
# R-faithful kernels for Control(rng="R"). Match robustbase's
# find_scale and refine_fast_s in lmrob.c line-for-line so that the
# rng="R" path can land bit-identical fits.
# ---------------------------------------------------------------------------


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_find_scale_r(
    cnp.ndarray[double, ndim=1, mode="c"] r,
    double b0,
    double initial_scale,
    int p,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
    int max_iter,
    double scale_tol,
):
    """Port of robustbase's ``find_scale`` in ``lmrob.c``.

    Fixed-point M-scale iteration::

        s_{k+1} = s_k * sqrt(sum_rho(r/s_k) / ((n-p) * b0))

    Returns ``(scale, n_iter)``. ``scale`` is the converged M-scale;
    ``n_iter`` is the number of iterations used (0 on immediate
    convergence, ``max_iter`` if the loop ran out).
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()
    cdef Py_ssize_t n = r.shape[0]
    if initial_scale <= 0.0:
        return 0.0, 0
    cdef double* r_ptr = <double*>cnp.PyArray_DATA(r)
    cdef double* t_ptr = <double*>cnp.PyArray_DATA(tuning)
    cdef double s = initial_scale
    cdef double prev = initial_scale
    cdef double inv_npmp = 1.0 / (<double>(n - p))
    cdef double sum_chi, diff
    cdef int it
    cdef int iters_used = max_iter
    with nogil:
        for it in range(max_iter):
            sum_chi = _chi_sum(r_ptr, n, prev, family, t_ptr)
            s = prev * sqrt(sum_chi * inv_npmp / b0)
            diff = s - prev
            if diff < 0:
                diff = -diff
            if diff <= scale_tol * prev:
                iters_used = it
                break
            prev = s
    return s, iters_used


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_refine_fast_s_r(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[double, ndim=1, mode="c"] beta,
    double initial_scale,
    int kk,
    int conv_flag,
    int max_k,
    double rel_tol,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
    double b0,
):
    """Port of robustbase's ``refine_fast_s`` in ``lmrob.c``.

    Two modes (per the ``conv_flag`` argument):

    - ``conv_flag == 0``: do exactly ``kk`` refining steps without a
      convergence check. Used inside the resample loop with ``kk =
      k_fast_s``.
    - ``conv_flag == 1``: iterate up to ``max_k`` steps with R's
      ``del <= rel_tol * fmax(rel_tol, ||beta_cand||_2)`` test (uses
      the OLD beta's norm). Used for survivor refinement.

    Each step:
    - One Newton-step on the M-scale: ``s = s * sqrt(sum_rho(r/s) /
      ((n-p) * b0))``.
    - IRWLS via LAPACK ``dgels`` (QR-based, same as R).

    ``initial_scale``: if ``<= 0``, recomputed via ``MAD(res, center=0)``
    (matching robustbase's path).

    ``beta`` is updated in place.

    Returns ``(scale, n_iter, converged, status)``. Status: 0 ok, 2
    exact fit (perfect zeros), 3 LAPACK error.
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    cdef double* X_w = <double*>malloc(n * p * sizeof(double))
    cdef double* y_w = <double*>malloc(n * sizeof(double))
    cdef double* r_buf = <double*>malloc(n * sizeof(double))
    cdef double* w_buf = <double*>malloc(n * sizeof(double))
    cdef double* beta_prev = <double*>malloc(p * sizeof(double))
    cdef double* aux = <double*>malloc(n * sizeof(double))
    cdef int lwork = max(1, n_int * p_int + 64 * (n_int + p_int))
    cdef double* work = <double*>malloc(lwork * sizeof(double))
    if (X_w == NULL or y_w == NULL or r_buf == NULL or w_buf == NULL
            or beta_prev == NULL or aux == NULL or work == NULL):
        free(X_w); free(y_w); free(r_buf); free(w_buf)
        free(beta_prev); free(aux); free(work)
        raise MemoryError("cy_refine_fast_s_r: out of memory")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef double* beta_data = <double*>cnp.PyArray_DATA(beta)
    cdef double* t_data = <double*>cnp.PyArray_DATA(tuning)

    cdef Py_ssize_t i, j, it
    cdef double s = initial_scale
    cdef double sw, dot, delta, nrmB, diff, sum_chi
    cdef double inv_npmp = 1.0 / (<double>(n - p))
    cdef int converged = 0
    cdef int status = 0
    cdef int iters = 0
    cdef int kmax
    cdef double half_med, q
    cdef Py_ssize_t mid

    if conv_flag != 0:
        kmax = max_k
    else:
        kmax = kk

    with nogil:
        # Residuals at the input beta.
        for i in range(n):
            dot = 0.0
            for j in range(p):
                dot += X_data[i * p + j] * beta_data[j]
            r_buf[i] = y_data[i] - dot

        # If initial_scale <= 0, recompute via MAD around zero (R's
        # ``MAD(res, center=0, ...)``).
        if s <= 0.0:
            for i in range(n):
                aux[i] = r_buf[i] if r_buf[i] >= 0 else -r_buf[i]
            # Inline median of aux[0..n-1] via a copy sort.
            # Use std-style insertion of a sort routine here; for n
            # small (n in 21..few-thousand) a partial sort would do.
            # For simplicity, do a full sort in O(n log n). Allocates
            # nothing new because we sort in-place over aux.
            # Bubble-ish insertion sort is fine for our sizes; switch
            # to qsort if hot.
            for i in range(1, n):
                # insertion sort
                dot = aux[i]
                j = i - 1
                while j >= 0 and aux[j] > dot:
                    aux[j + 1] = aux[j]
                    j -= 1
                aux[j + 1] = dot
            mid = n // 2
            if (n & 1) == 1:
                half_med = aux[mid]
            else:
                half_med = 0.5 * (aux[mid - 1] + aux[mid])
            s = 1.4826 * half_med

        # Main loop.
        for it in range(kmax):
            iters = it + 1
            if s == 0.0:
                status = 2
                converged = 1
                break

            # One Newton step on the M-scale.
            sum_chi = _chi_sum(r_buf, n, s, family, t_data)
            s = s * sqrt(sum_chi * inv_npmp / b0)
            if s == 0.0:
                status = 2
                converged = 1
                break

            # IRWLS weights.
            _wgt_zinv(r_buf, w_buf, n, s, family, t_data)

            # Build weighted system: y_w = sqrt(w) * y, X_w = sqrt(w) * X
            # in column-major layout for dgels.
            for i in range(n):
                sw = sqrt(w_buf[i]) if w_buf[i] > 0 else 0.0
                y_w[i] = y_data[i] * sw
                for j in range(p):
                    X_w[i + j * n] = X_data[i * p + j] * sw

            # Save old beta for the convergence test.
            for j in range(p):
                beta_prev[j] = beta_data[j]

            dgels(b'N', &n_int, &p_int, &one,
                  X_w, &n_int,
                  y_w, &n_int,
                  work, &lwork, &info)
            if info != 0:
                status = 3
                break
            for j in range(p):
                beta_data[j] = y_w[j]

            # Convergence test only when conv_flag != 0.
            if conv_flag != 0:
                delta = 0.0
                nrmB = 0.0
                for j in range(p):
                    diff = beta_data[j] - beta_prev[j]
                    delta += diff * diff
                    nrmB += beta_prev[j] * beta_prev[j]
                delta = sqrt(delta)
                nrmB = sqrt(nrmB)
                # R: del <= rel_tol * fmax2(rel_tol, ||beta_cand||_2)
                q = rel_tol if rel_tol > nrmB else nrmB
                if delta <= rel_tol * q:
                    converged = 1
                    break

            # Update residuals for the next iter.
            for i in range(n):
                dot = 0.0
                for j in range(p):
                    dot += X_data[i * p + j] * beta_data[j]
                r_buf[i] = y_data[i] - dot

    free(X_w); free(y_w); free(r_buf); free(w_buf)
    free(beta_prev); free(aux); free(work)
    return s, iters, converged, status


# ---------------------------------------------------------------------------
# Full per-iteration kernel: draw a p-subset from a bitgen, retry up to
# mts times if the LU factorisation reports a singular pivot, then run the
# same k-step refinement as cy_resample_iter.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_draw_and_iter(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    object bitgen_capsule,
    int mts,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
    double b0,
    int k_fast_s,
    int max_iter_scale,
    double scale_tol,
    cnp.ndarray[double, ndim=1, mode="c"] beta_out,
):
    """End-to-end nogil iter: subset draw + initial solve + k-step refinement.

    ``bitgen_capsule`` is the ``BitGenerator.capsule`` of an
    ``np.random.Generator`` whose state is owned by this call. The caller
    must not use the same Generator from another thread concurrently.

    Returns ``(scale, status)``. Status codes match :func:`cy_resample_iter`:
    0 (ok), 1 (no non-singular subset in mts tries), 2 (exact fit),
    3 (LAPACK error during dgels).
    """
    cdef bitgen_t* bg = <bitgen_t*>PyCapsule_GetPointer(
        bitgen_capsule, "BitGenerator"
    )

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef double* sub_X = <double*>malloc(p * p * sizeof(double))
    cdef double* sub_y = <double*>malloc(p * sizeof(double))
    cdef double* X_w = <double*>malloc(n * p * sizeof(double))
    cdef double* y_w = <double*>malloc(n * sizeof(double))
    cdef double* r_buf = <double*>malloc(n * sizeof(double))
    cdef double* w = <double*>malloc(n * sizeof(double))
    cdef int* ipiv = <int*>malloc(p * sizeof(int))
    cdef Py_ssize_t* perm = <Py_ssize_t*>malloc(n * sizeof(Py_ssize_t))

    cdef int lwork = max(1, n_int * p_int + 64 * (n_int + p_int))
    cdef double* work = <double*>malloc(lwork * sizeof(double))

    if (sub_X == NULL or sub_y == NULL or X_w == NULL or y_w == NULL
            or r_buf == NULL or w == NULL or ipiv == NULL or perm == NULL
            or work == NULL):
        free(sub_X); free(sub_y); free(X_w); free(y_w)
        free(r_buf); free(w); free(ipiv); free(perm); free(work)
        raise MemoryError("cy_draw_and_iter: out of memory")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef double* beta_data = <double*>cnp.PyArray_DATA(beta_out)
    cdef double* tuning_data = <double*>cnp.PyArray_DATA(tuning)

    cdef Py_ssize_t i, j, kk, row, swap
    cdef uint64_t r
    cdef double s, sw, dot
    cdef int status = 0
    cdef double scale = 0.0
    cdef int got_subset = 0
    cdef int try_i
    cdef double k0 = tuning_data[0]

    with nogil:
        # Try up to mts times to draw a non-singular subset and solve it.
        for try_i in range(mts):
            # Floyd's combination algorithm, matching numpy's
            # ``rng.choice(n, p, replace=False)`` for shared PCG64 state:
            # for j in {n-p, n-p+1, ..., n-1} draw idx in [0, j+1); if
            # idx already in result, substitute j. Result is insertion-
            # ordered.
            for i in range(p):
                j = n - p + i
                r = _bounded_uint64(bg, <uint64_t>(j + 1))
                swap = <Py_ssize_t>r
                # Linear scan over the current sample; p is small.
                row = 0
                for kk in range(i):
                    if perm[kk] == swap:
                        row = 1
                        break
                if row:
                    perm[i] = j
                else:
                    perm[i] = swap

            # Build column-major sub_X (p×p) and sub_y (p).
            for i in range(p):
                row = perm[i]
                sub_y[i] = y_data[row]
                for j in range(p):
                    sub_X[i + j * p] = X_data[row * p + j]

            dgesv(&p_int, &one, sub_X, &p_int, ipiv, sub_y, &p_int, &info)
            if info == 0:
                got_subset = 1
                break

        if not got_subset:
            status = 1
        else:
            for i in range(p):
                beta_data[i] = sub_y[i]

            for i in range(n):
                dot = 0.0
                for j in range(p):
                    dot += X_data[i * p + j] * beta_data[j]
                r_buf[i] = y_data[i] - dot

            s = 0.0
            for i in range(n):
                if r_buf[i] >= 0:
                    if r_buf[i] > s:
                        s = r_buf[i]
                else:
                    if -r_buf[i] > s:
                        s = -r_buf[i]
            if s <= 0.0:
                status = 2
                scale = 0.0
            else:
                s = s / k0 if k0 > 0 else s
                if s <= 0.0:
                    s = 1.0

                for kk in range(k_fast_s):
                    s = _mscale_generic(r_buf, n, s, p, family, tuning_data,
                                        b0, max_iter_scale, scale_tol)
                    if s == 0.0:
                        status = 2
                        scale = 0.0
                        break

                    _wgt_zinv(r_buf, w, n, s, family, tuning_data)
                    for i in range(n):
                        sw = sqrt(w[i]) if w[i] > 0 else 0.0
                        y_w[i] = y_data[i] * sw
                        for j in range(p):
                            X_w[i + j * n] = X_data[i * p + j] * sw

                    dgels(b'N', &n_int, &p_int, &one,
                          X_w, &n_int,
                          y_w, &n_int,
                          work, &lwork, &info)
                    if info != 0:
                        status = 3
                        break
                    for i in range(p):
                        beta_data[i] = y_w[i]
                    for i in range(n):
                        dot = 0.0
                        for j in range(p):
                            dot += X_data[i * p + j] * beta_data[j]
                        r_buf[i] = y_data[i] - dot
                else:
                    s = _mscale_generic(r_buf, n, s, p, family, tuning_data,
                                        b0, max_iter_scale, scale_tol)
                    scale = s

    free(sub_X); free(sub_y); free(X_w); free(y_w)
    free(r_buf); free(w); free(ipiv); free(perm); free(work)

    return scale, status


# Backwards-compatible alias kept for the v0.3.0 API.
def cy_resample_iter_bisquare(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[cnp.int64_t, ndim=1, mode="c"] idx,
    double k_chi,
    double b0,
    int k_fast_s,
    int max_iter_scale,
    double scale_tol,
    cnp.ndarray[double, ndim=1, mode="c"] beta_out,
):
    """Backwards-compatible alias: bisquare with single tuning parameter."""
    cdef cnp.ndarray[double, ndim=1, mode="c"] tuning = np.zeros(3, dtype=np.float64)
    tuning[0] = k_chi
    return cy_resample_iter(
        X, y, idx, 0, tuning, b0, k_fast_s, max_iter_scale, scale_tol, beta_out
    )
