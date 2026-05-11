# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Monolithic Cython kernel for ``lmrob``. Mirrors robustbase/src/lmrob.c
# in structure: one top-level entry point that owns the whole fast-S +
# survivor-refinement pipeline. All scratch buffers are allocated once at
# the top and reused. LAPACK via scipy.linalg.cython_lapack; RNG via
# numpy's bitgen_t capsule.
#
# Stage 1: bisquare-only fast-S + refinement, returning (beta, scale,
# residuals, rweights, status). MM, D-scale, and vcov will land in later
# stages (the existing Python paths still handle them for now).

cimport cython
from cpython.pycapsule cimport PyCapsule_GetPointer
from libc.math cimport fabs, sqrt
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy

from numpy.random cimport bitgen_t
from scipy.linalg.cython_lapack cimport dgels, dgesv

import numpy as np
cimport numpy as cnp

cnp.import_array()


# ---------------------------------------------------------------------------
# Family enum. Mirrors pyrobustlm._fast_s._FAMILY_IDS for now; will be
# the single source of truth once the monolithic kernel is the default.
# ---------------------------------------------------------------------------
cdef enum:
    FAM_BISQUARE = 0
    FAM_HAMPEL = 1
    FAM_OPTIMAL = 2
    FAM_LQQ = 3
    FAM_GGW = 4


# ---------------------------------------------------------------------------
# Bisquare chi (normalised so chi(inf) = 1) and bisquare-derived weight.
# Inlined nogil so the inner loops never cross the Python boundary.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _bisquare_chi_sum(
    const double* r, Py_ssize_t n, double s, double k,
) nogil:
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
cdef inline void _bisquare_wgt(
    const double* r, double* out, Py_ssize_t n, double s, double k,
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
# Random integer in [0, bound). Lemire's debiased modulo via numpy bitgen.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline uint64_t _bounded_uint64(bitgen_t* bg, uint64_t bound) nogil:
    if bound <= 1:
        return 0
    cdef uint64_t threshold = (-bound) % bound
    cdef uint64_t r
    while True:
        r = bg.next_uint64(bg.state)
        if r >= threshold:
            return r % bound


# ---------------------------------------------------------------------------
# Floyd's combination algorithm. Fills perm[0..p-1] with p distinct indices
# from {0..n-1}, insertion-ordered.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _draw_subset(
    bitgen_t* bg, Py_ssize_t n, Py_ssize_t p, Py_ssize_t* perm,
) nogil:
    cdef Py_ssize_t i, j, idx
    cdef int found
    for i in range(p):
        j = n - p + i
        idx = <Py_ssize_t>_bounded_uint64(bg, <uint64_t>(j + 1))
        found = 0
        for k in range(i):
            if perm[k] == idx:
                found = 1
                break
        if found:
            perm[i] = j
        else:
            perm[i] = idx


# Forward-declared variable for the loop above. Cython demands the `k` be
# typed; using a separate kludge avoids leaking it into other scopes.
ctypedef Py_ssize_t _ssz


# ---------------------------------------------------------------------------
# Per-thread scratch: all the workspaces a fast-S iteration needs. Allocate
# once at the top of cy_lmrob_fit, reuse for every iteration.
# ---------------------------------------------------------------------------

cdef struct _Scratch:
    double* sub_X       # column-major (p, p)
    double* sub_y       # (p,)
    double* X_w         # column-major (n, p), weighted design
    double* y_w         # (n,), weighted response
    double* r           # (n,), residuals
    double* w           # (n,), IRWLS weights
    Py_ssize_t* perm    # (p,), subset indices
    double* beta        # (p,), current beta
    double* beta_prev   # (p,), previous beta (for convergence check)
    double* work        # LAPACK workspace
    int* ipiv           # (p,), dgesv pivots
    int lwork


cdef int _alloc_scratch(_Scratch* s, Py_ssize_t n, Py_ssize_t p) nogil:
    cdef int lwork = max(1, <int>(n * p) + 64 * (<int>n + <int>p))
    s.sub_X = <double*>malloc(p * p * sizeof(double))
    s.sub_y = <double*>malloc(p * sizeof(double))
    s.X_w = <double*>malloc(n * p * sizeof(double))
    s.y_w = <double*>malloc(n * sizeof(double))
    s.r = <double*>malloc(n * sizeof(double))
    s.w = <double*>malloc(n * sizeof(double))
    s.perm = <Py_ssize_t*>malloc(p * sizeof(Py_ssize_t))
    s.beta = <double*>malloc(p * sizeof(double))
    s.beta_prev = <double*>malloc(p * sizeof(double))
    s.work = <double*>malloc(lwork * sizeof(double))
    s.ipiv = <int*>malloc(p * sizeof(int))
    s.lwork = lwork
    if (s.sub_X == NULL or s.sub_y == NULL or s.X_w == NULL or s.y_w == NULL
            or s.r == NULL or s.w == NULL or s.perm == NULL or s.beta == NULL
            or s.beta_prev == NULL or s.work == NULL or s.ipiv == NULL):
        return -1
    return 0


cdef void _free_scratch(_Scratch* s) nogil:
    free(s.sub_X); free(s.sub_y); free(s.X_w); free(s.y_w)
    free(s.r); free(s.w); free(s.perm); free(s.beta); free(s.beta_prev)
    free(s.work); free(s.ipiv)


# ---------------------------------------------------------------------------
# Residuals = y - X @ beta. Both X (n*p row-major) and beta (p) are in C
# layout; the dot product is straightforward.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void _residuals(
    const double* X, const double* y, const double* beta,
    double* out, Py_ssize_t n, Py_ssize_t p,
) nogil:
    cdef Py_ssize_t i, j
    cdef double dot
    for i in range(n):
        dot = 0.0
        for j in range(p):
            dot += X[i * p + j] * beta[j]
        out[i] = y[i] - dot


# ---------------------------------------------------------------------------
# One IRWLS step: build weighted (X_w, y_w) in column-major, then dgels.
# beta is updated in place.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int _irwls_step(
    const double* X, const double* y, const double* r_in,
    double s, double k,
    double* beta, _Scratch* scr,
    Py_ssize_t n, Py_ssize_t p,
) nogil:
    """Returns 0 on success, 3 on LAPACK error."""
    cdef Py_ssize_t i, j
    cdef double sw
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    _bisquare_wgt(r_in, scr.w, n, s, k)
    for i in range(n):
        sw = sqrt(scr.w[i]) if scr.w[i] > 0 else 0.0
        scr.y_w[i] = y[i] * sw
        for j in range(p):
            scr.X_w[i + j * n] = X[i * p + j] * sw

    dgels(b'N', &n_int, &p_int, &one,
          scr.X_w, &n_int,
          scr.y_w, &n_int,
          scr.work, &scr.lwork, &info)
    if info != 0:
        return 3
    for j in range(p):
        beta[j] = scr.y_w[j]
    return 0


# ---------------------------------------------------------------------------
# K-step refinement inside the resampling loop. Body of one candidate.
# Returns final scale; updates scr.beta in place. status is 0 on success,
# 2 on exact fit, 3 on LAPACK error.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _k_step_refine(
    const double* X, const double* y,
    double init_s, double k, double b0,
    int k_fast_s, int max_iter_scale, double scale_tol,
    _Scratch* scr,
    Py_ssize_t n, Py_ssize_t p,
    int* status_out,
) nogil:
    cdef double s = init_s
    cdef int kk, irwls_status
    cdef int status = 0

    for kk in range(k_fast_s):
        s = _bisquare_mscale(scr.r, n, s, p, k, b0, max_iter_scale, scale_tol)
        if s == 0.0:
            status = 2
            status_out[0] = status
            return 0.0
        irwls_status = _irwls_step(X, y, scr.r, s, k, scr.beta, scr, n, p)
        if irwls_status != 0:
            status = irwls_status
            status_out[0] = status
            return s
        _residuals(X, y, scr.beta, scr.r, n, p)
    # Final M-scale at the resulting beta.
    s = _bisquare_mscale(scr.r, n, s, p, k, b0, max_iter_scale, scale_tol)
    status_out[0] = 0
    return s


# ---------------------------------------------------------------------------
# Survivor refinement to convergence on a single candidate.
# Returns final scale; updates beta in place. Sets converged_out.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _refine_to_convergence(
    const double* X, const double* y,
    double init_s, double k, double b0,
    int max_it, double refine_tol,
    int max_iter_scale, double scale_tol,
    double* beta,
    _Scratch* scr,
    Py_ssize_t n, Py_ssize_t p,
    int* converged_out, int* n_iter_out,
) nogil:
    cdef double s = init_s
    cdef double delta, denom, diff
    cdef Py_ssize_t j
    cdef int it
    cdef int irwls_status

    for it in range(max_it):
        _residuals(X, y, beta, scr.r, n, p)
        s = _bisquare_mscale(scr.r, n, s, p, k, b0, max_iter_scale, scale_tol)
        if s == 0.0:
            converged_out[0] = 1
            n_iter_out[0] = it + 1
            return 0.0
        # Save previous beta.
        for j in range(p):
            scr.beta_prev[j] = beta[j]
        irwls_status = _irwls_step(X, y, scr.r, s, k, beta, scr, n, p)
        if irwls_status != 0:
            converged_out[0] = 0
            n_iter_out[0] = it + 1
            return s
        # Relative L2 change.
        delta = 0.0
        denom = 0.0
        for j in range(p):
            diff = beta[j] - scr.beta_prev[j]
            delta += diff * diff
            denom += scr.beta_prev[j] * scr.beta_prev[j]
        delta = sqrt(delta)
        denom = sqrt(denom) if denom > 1e-300 else 1e-150
        if delta / denom < refine_tol:
            converged_out[0] = 1
            n_iter_out[0] = it + 1
            return s
    converged_out[0] = 0
    n_iter_out[0] = max_it
    return s


# ---------------------------------------------------------------------------
# Top-level fast-S entry point.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_lmrob_fast_s(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    object bitgen_capsule,
    double k_chi,
    double b0,
    int nResample,
    int mts,
    int k_fast_s,
    int best_r,
    int max_it,
    double refine_tol,
    int max_iter_scale,
    double scale_tol,
    cnp.ndarray[double, ndim=1, mode="c"] beta_out,
):
    """Bisquare fast-S, single nogil call.

    Runs the full pipeline: ``nResample`` candidates × (subset draw +
    initial solve + ``k_fast_s`` K-step refinements), best-of-``best_r``
    selection, then survivor refinement to convergence. All in one C
    block with one workspace allocation.

    Returns ``(scale, status, n_iter, converged)`` where status is
    0 (ok), 1 (no non-singular subset found), 2 (exact fit), 3 (LAPACK
    error), 4 (alloc failed).
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

    # Per-thread scratch.
    cdef _Scratch scr
    if _alloc_scratch(&scr, n, p) != 0:
        _free_scratch(&scr)
        raise MemoryError("cy_lmrob_fast_s: out of memory")

    # Best-of-best_r heap kept in two parallel C arrays. Size at most best_r.
    cdef double* best_scales = <double*>malloc(best_r * sizeof(double))
    cdef double* best_betas = <double*>malloc(best_r * p * sizeof(double))
    if best_scales == NULL or best_betas == NULL:
        free(best_scales); free(best_betas)
        _free_scratch(&scr)
        raise MemoryError("cy_lmrob_fast_s: out of memory for best-r heap")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef double* beta_out_data = <double*>cnp.PyArray_DATA(beta_out)

    cdef Py_ssize_t i, j, row, swap, try_i, kept
    cdef double s, scale, max_abs, candidate_scale
    cdef int got_subset, status, k_status, converged
    cdef int worst_i, total_iters
    cdef int found
    cdef uint64_t r_u

    status = 0
    scale = 0.0
    kept = 0
    converged = 0
    total_iters = 0

    with nogil:
        # Resampling loop ----------------------------------------------------
        for try_i in range(nResample):
            # Draw a p-subset.
            got_subset = 0
            for _ in range(mts):
                for i in range(p):
                    j = n - p + i
                    r_u = _bounded_uint64(bg, <uint64_t>(j + 1))
                    swap = <Py_ssize_t>r_u
                    found = 0
                    for row in range(i):
                        if scr.perm[row] == swap:
                            found = 1
                            break
                    if found:
                        scr.perm[i] = j
                    else:
                        scr.perm[i] = swap

                # Build sub_X (col-major), sub_y.
                for i in range(p):
                    row = scr.perm[i]
                    scr.sub_y[i] = y_data[row]
                    for j in range(p):
                        scr.sub_X[i + j * p] = X_data[row * p + j]

                dgesv(&p_int, &one, scr.sub_X, &p_int, scr.ipiv,
                      scr.sub_y, &p_int, &info)
                if info == 0:
                    got_subset = 1
                    break

            if not got_subset:
                continue  # try another resample

            for j in range(p):
                scr.beta[j] = scr.sub_y[j]
            _residuals(X_data, y_data, scr.beta, scr.r, n, p)

            # Coarse initial scale.
            max_abs = 0.0
            for i in range(n):
                if scr.r[i] >= 0:
                    if scr.r[i] > max_abs:
                        max_abs = scr.r[i]
                else:
                    if -scr.r[i] > max_abs:
                        max_abs = -scr.r[i]
            if max_abs == 0.0:
                # Exact fit on this subset.
                for j in range(p):
                    beta_out_data[j] = scr.beta[j]
                scale = 0.0
                status = 2
                break
            s = max_abs / k_chi if k_chi > 0 else max_abs
            if s <= 0.0:
                s = 1.0

            candidate_scale = _k_step_refine(
                X_data, y_data, s, k_chi, b0,
                k_fast_s, max_iter_scale, scale_tol,
                &scr, n, p, &k_status,
            )
            if k_status == 2:
                # Exact fit during refinement.
                for j in range(p):
                    beta_out_data[j] = scr.beta[j]
                scale = 0.0
                status = 2
                break
            if k_status != 0:
                continue  # LAPACK error on this candidate

            # Insert into the best-r heap.
            if kept < best_r:
                best_scales[kept] = candidate_scale
                for j in range(p):
                    best_betas[kept * p + j] = scr.beta[j]
                kept += 1
            else:
                worst_i = 0
                for i in range(1, best_r):
                    if best_scales[i] > best_scales[worst_i]:
                        worst_i = i
                if candidate_scale < best_scales[worst_i]:
                    best_scales[worst_i] = candidate_scale
                    for j in range(p):
                        best_betas[worst_i * p + j] = scr.beta[j]

        # Survivor refinement -----------------------------------------------
        if status == 0:
            if kept == 0:
                status = 1
            else:
                # Refine every survivor; keep the one with smallest scale.
                worst_i = -1  # repurposed: index of best after refine
                scale = 1e300
                for i in range(kept):
                    for j in range(p):
                        scr.beta[j] = best_betas[i * p + j]
                    candidate_scale = _refine_to_convergence(
                        X_data, y_data,
                        best_scales[i], k_chi, b0,
                        max_it, refine_tol,
                        max_iter_scale, scale_tol,
                        scr.beta, &scr, n, p,
                        &converged, &total_iters,
                    )
                    if candidate_scale < scale:
                        scale = candidate_scale
                        for j in range(p):
                            beta_out_data[j] = scr.beta[j]
                        worst_i = i

    free(best_scales); free(best_betas)
    _free_scratch(&scr)
    return scale, status, total_iters, converged
