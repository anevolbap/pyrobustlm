# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Cython kernel for one fast-S resampling iteration. The hot path of fast_s
# is now: pre-draw a p-subset in Python (small), then call this kernel
# without the GIL. All LAPACK calls go through scipy.linalg.cython_lapack;
# m_scale and the bisquare wgt come from the existing _psi.pyx kernels via
# the equivalent inline implementations here (kept identical to lmrob.c).

cimport cython
from libc.math cimport fabs, sqrt
from libc.stdlib cimport malloc, free

from scipy.linalg.cython_lapack cimport dgesv, dgels

import numpy as np
cimport numpy as cnp

cnp.import_array()


# ---------------------------------------------------------------------------
# Local copies of the bisquare M-scale and weight kernels, so the whole
# inner loop runs nogil. These mirror src/pyrobustlm/_core/_psi.pyx and
# robustbase's lmrob.c::find_scale.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _bisquare_wgt(const double* x, double* out, Py_ssize_t n, double k) nogil:
    """w = (1 - (x/k)^2)^2 if |x/k| <= 1 else 0."""
    cdef Py_ssize_t i
    cdef double a, u
    for i in range(n):
        a = x[i] / k
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            u = 1.0 - a * a
            out[i] = u * u


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _bisquare_chi_sum(const double* r, Py_ssize_t n, double s, double k) nogil:
    """sum of bisquare chi(r/s), with chi normalised so chi(inf)=1."""
    cdef Py_ssize_t i
    cdef double t, ax, total = 0.0
    for i in range(n):
        t = r[i] / s
        ax = t if t >= 0 else -t
        if ax >= k:
            total += 1.0
        else:
            t = ax / k
            t = 1.0 - t * t
            total += 1.0 - t * t * t
    return total


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _bisquare_mscale(
    const double* r,
    Py_ssize_t n,
    double init_scale,
    Py_ssize_t p,
    double k,
    double b0,
    int max_iter,
    double tol,
) nogil:
    """Iterate ``s_{i+1} = s_i * sqrt( sum(chi(r/s_i)) / ((n-p) * b0) )``."""
    cdef double s = init_scale
    cdef double prev = init_scale
    cdef double inv_npmp = 1.0 / (<double>(n - p))
    cdef double sum_chi, diff
    cdef int it
    if s <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = _bisquare_chi_sum(r, n, s, k)
        s = prev * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev
        if diff < 0:
            diff = -diff
        if diff <= tol * prev:
            return s
        prev = s
    return s


# ---------------------------------------------------------------------------
# The resampling iteration.
# ---------------------------------------------------------------------------
# ``X`` is C-contiguous (n, p) row-major. LAPACK's dgesv/dgels expect
# column-major. We allocate scratch buffers in column-major layout inside
# this function (per-thread, freed at the end). Allocation goes through
# malloc/free so the function can stay nogil end-to-end.

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_resample_iter_bisquare(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[long, ndim=1, mode="c"] idx,
    double k_chi,
    double b0,
    int k_fast_s,
    int max_iter_scale,
    double scale_tol,
    cnp.ndarray[double, ndim=1, mode="c"] beta_out,
):
    """Run one fast-S resampling iteration with the bisquare psi/chi.

    Returns a tuple ``(scale, status)`` where ``status`` is 0 on success,
    1 if the initial subset proved singular, 2 if an exact fit was found
    (scale == 0, beta_out holds the candidate), 3 on a LAPACK error.
    """
    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    # Workspaces. Sub_X is column-major (p, p); X_w (n, p) column-major; etc.
    cdef double* sub_X = <double*>malloc(p * p * sizeof(double))
    cdef double* sub_y = <double*>malloc(p * sizeof(double))
    cdef double* X_w = <double*>malloc(n * p * sizeof(double))
    cdef double* y_w = <double*>malloc(n * sizeof(double))
    cdef double* r = <double*>malloc(n * sizeof(double))
    cdef double* w = <double*>malloc(n * sizeof(double))
    cdef int* ipiv = <int*>malloc(p * sizeof(int))

    # dgels workspace. The optimal lwork can be queried with lwork=-1, but
    # a safe over-allocation is fine for the typical (n=20-5000, p=2-50)
    # range we care about.
    cdef int lwork = max(1, n_int * p_int + 64 * (n_int + p_int))
    cdef double* work = <double*>malloc(lwork * sizeof(double))

    if (sub_X == NULL or sub_y == NULL or X_w == NULL or y_w == NULL
            or r == NULL or w == NULL or ipiv == NULL or work == NULL):
        free(sub_X); free(sub_y); free(X_w); free(y_w)
        free(r); free(w); free(ipiv); free(work)
        raise MemoryError("cy_resample_iter_bisquare: out of memory")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef long* idx_data = <long*>cnp.PyArray_DATA(idx)
    cdef double* beta_data = <double*>cnp.PyArray_DATA(beta_out)

    cdef Py_ssize_t i, j, kk, row
    cdef double s, sw, dot

    cdef int status = 0
    cdef double scale = 0.0

    with nogil:
        # ----- Build column-major sub_X (p×p) and sub_y (p) -----
        for i in range(p):
            row = idx_data[i]
            sub_y[i] = y_data[row]
            # sub_X column-major: sub_X[i + j*p] = X[row, j]
            for j in range(p):
                sub_X[i + j * p] = X_data[row * p + j]

        # dgesv: solve sub_X · x = sub_y in place, x replaces sub_y.
        dgesv(&p_int, &one, sub_X, &p_int, ipiv, sub_y, &p_int, &info)
        if info != 0:
            status = 1  # singular
        else:
            # beta_out <- sub_y (the solution)
            for i in range(p):
                beta_data[i] = sub_y[i]

            # ----- residuals r = y - X @ beta -----
            for i in range(n):
                dot = 0.0
                for j in range(p):
                    dot += X_data[i * p + j] * beta_data[j]
                r[i] = y_data[i] - dot

            # ----- initial scale via MAD(r) -----
            # We piggyback on the m-scale iteration with init_scale=MAD; for
            # speed and simplicity we skip MAD and let the bisquare iteration
            # converge from a quick estimate. R uses MAD; using
            # max(|r|)/k_chi as a coarse upper bound also works.
            s = 0.0
            for i in range(n):
                if r[i] >= 0:
                    if r[i] > s:
                        s = r[i]
                else:
                    if -r[i] > s:
                        s = -r[i]
            if s <= 0.0:
                # Exact fit.
                status = 2
                scale = 0.0
            else:
                s = s / k_chi
                if s <= 0.0:
                    s = 1.0

                # ----- k_fast_s rounds of (m_scale, IRWLS step) -----
                for kk in range(k_fast_s):
                    s = _bisquare_mscale(r, n, s, p, k_chi, b0,
                                         max_iter_scale, scale_tol)
                    if s == 0.0:
                        status = 2
                        scale = 0.0
                        break
                    # Build sqrt(weights) in w; weighted X in X_w (col-major);
                    # weighted y in y_w.
                    _bisquare_wgt_zinv(r, w, n, s, k_chi)
                    for i in range(n):
                        sw = sqrt(w[i])
                        y_w[i] = y_data[i] * sw
                        for j in range(p):
                            X_w[i + j * n] = X_data[i * p + j] * sw

                    # dgels: min ||X_w · beta - y_w||_2 in place.
                    dgels(b'N', &n_int, &p_int, &one,
                          X_w, &n_int,
                          y_w, &n_int,
                          work, &lwork, &info)
                    if info != 0:
                        status = 3
                        break
                    # Solution lives in y_w[0:p].
                    for i in range(p):
                        beta_data[i] = y_w[i]
                    # Recompute residuals.
                    for i in range(n):
                        dot = 0.0
                        for j in range(p):
                            dot += X_data[i * p + j] * beta_data[j]
                        r[i] = y_data[i] - dot
                else:
                    # Final m_scale at the resulting beta.
                    s = _bisquare_mscale(r, n, s, p, k_chi, b0,
                                         max_iter_scale, scale_tol)
                    scale = s

    free(sub_X); free(sub_y); free(X_w); free(y_w)
    free(r); free(w); free(ipiv); free(work)

    return scale, status


# bisquare wgt operating on residual r and producing the IRWLS weight
# w = psi(r/s) / (r/s) = (1 - (r/(s*k))^2)^2 (clip when |r|>s*k).
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _bisquare_wgt_zinv(
    const double* r, double* out, Py_ssize_t n, double s, double k
) nogil:
    cdef Py_ssize_t i
    cdef double a, u
    cdef double inv_sk = 1.0 / (s * k)
    for i in range(n):
        a = r[i] * inv_sk
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            u = 1.0 - a * a
            out[i] = u * u
