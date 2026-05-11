# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Cython kernel for one fast-S resampling iteration. Family-generic:
# bisquare, hampel, optimal, and lqq dispatch on a small family enum.
# ggw is not supported here yet (the polynomial-chi tables would need to
# be brought over from _psi.pyx); ggw fits go through the NumPy path.

cimport cython
from libc.math cimport fabs, sqrt
from libc.stdlib cimport malloc, free

from scipy.linalg.cython_lapack cimport dgesv, dgels

import numpy as np
cimport numpy as cnp

cnp.import_array()


# Family enum. Keep these in sync with the dispatch table in
# pyrobustlm/_fast_s.py (`_FAMILY_IDS`).
cdef enum:
    FAM_BISQUARE = 0
    FAM_HAMPEL = 1
    FAM_OPTIMAL = 2
    FAM_LQQ = 3


# ---------------------------------------------------------------------------
# Per-family chi-sum and IRWLS weight kernels. All take a ``tuning``
# buffer of up to 3 doubles; layout is family-specific.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _chi_sum(
    const double* r,
    Py_ssize_t n,
    double s,
    int family,
    const double* tuning,
) nogil:
    cdef Py_ssize_t i
    cdef double total = 0.0
    cdef double a, ax, t, xi, u
    cdef double k, a_t, b_t, r_t, c, b_l, s_l, k01, denom, s5, s6, k01_2, end3, s0, dx
    cdef double R1h, R2h, R3h, R4h, ac, a2, nc, inv_nc, bma_inv_half

    if family == FAM_BISQUARE:
        k = tuning[0]
        for i in range(n):
            t = r[i] / s
            ax = t if t >= 0 else -t
            if ax >= k:
                total += 1.0
            else:
                t = ax / k
                t = 1.0 - t * t
                total += 1.0 - t * t * t
    elif family == FAM_HAMPEL:
        a_t = tuning[0]; b_t = tuning[1]; r_t = tuning[2]
        nc = a_t * (b_t + r_t - a_t) * 0.5
        inv_nc = 1.0 / nc
        bma_inv_half = 0.5 / (r_t - b_t)
        for i in range(n):
            xi = r[i] / s
            u = xi if xi >= 0 else -xi
            if u <= a_t:
                total += (xi * xi * 0.5) * inv_nc
            elif u <= b_t:
                total += (u - 0.5 * a_t) * a_t * inv_nc
            elif u <= r_t:
                total += (b_t - 0.5 * a_t + (u - b_t) * (1.0 - (u - b_t) * bma_inv_half)) * a_t * inv_nc
            else:
                total += 1.0
    elif family == FAM_OPTIMAL:
        k = tuning[0]
        R1h = -1.944 * 0.5
        R2h = 1.728 * 0.25
        R3h = -0.312 / 6.0
        R4h = 0.016 / 8.0
        for i in range(n):
            xi = r[i] / s
            ac = xi / k
            ax = ac if ac >= 0 else -ac
            if ax > 3.0:
                total += 1.0
            elif ax > 2.0:
                a2 = ax * ax
                total += (a2 * (R1h + a2 * (R2h + a2 * (R3h + a2 * R4h))) + 1.792) / 3.25
            else:
                total += (ac * ac) / 6.5
    else:  # FAM_LQQ
        b_l = tuning[0]; c = tuning[1]; s_l = tuning[2]
        k01 = b_l + c
        s5 = s_l - 1.0
        s6 = -2.0 * k01 + b_l * s_l
        k01_2 = k01 * k01
        denom = s_l * c * (3.0 * c + 2.0 * b_l) + k01_2
        if s5 == 0.0:
            end3 = k01
        else:
            end3 = k01 - s6 / s5
        for i in range(n):
            xi = r[i] / s
            ax = xi if xi >= 0 else -xi
            if ax <= c:
                total += (3.0 * s_l - 3.0) / denom * xi * xi
            elif ax <= k01:
                s0 = ax - c
                total += (6.0 * s_l - 6.0) / denom * (xi * xi * 0.5 - s_l / b_l * s0 * s0 * s0 / 6.0)
            elif ax < end3:
                dx = ax - k01
                total += (6.0 * s5) / denom * (
                    k01_2 * 0.5 - s_l * b_l * b_l / 6.0
                    - dx * 0.5 * (s6 + dx * (s5 + dx * s5 * s5 / 3.0 / s6))
                )
            else:
                total += 1.0
    return total


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _wgt_zinv(
    const double* r,
    double* out,
    Py_ssize_t n,
    double s,
    int family,
    const double* tuning,
) nogil:
    cdef Py_ssize_t i
    cdef double a, u, ax
    cdef double k, inv_sk, ac, a2
    cdef double a_t, b_t, r_t, xi
    cdef double R1, R2, R3, R4
    cdef double b_l, c, s_l, k01, denom, s5, s6, k01_2, end3, s0, dx
    cdef double rho_p

    if family == FAM_BISQUARE:
        k = tuning[0]
        inv_sk = 1.0 / (s * k)
        for i in range(n):
            a = r[i] * inv_sk
            if a < -1.0 or a > 1.0:
                out[i] = 0.0
            else:
                u = 1.0 - a * a
                out[i] = u * u
    elif family == FAM_HAMPEL:
        a_t = tuning[0]; b_t = tuning[1]; r_t = tuning[2]
        for i in range(n):
            xi = r[i] / s
            u = xi if xi >= 0 else -xi
            if u <= a_t:
                out[i] = 1.0
            elif u <= b_t:
                out[i] = a_t / u if u > 0 else 1.0
            elif u <= r_t:
                out[i] = (a_t * (r_t - u) / (r_t - b_t)) / u
            else:
                out[i] = 0.0
    elif family == FAM_OPTIMAL:
        k = tuning[0]
        R1 = -1.944; R2 = 1.728; R3 = -0.312; R4 = 0.016
        for i in range(n):
            xi = r[i] / s
            ac = xi / k
            ax = ac if ac >= 0 else -ac
            if ax > 3.0:
                out[i] = 0.0
            elif ax > 2.0:
                a2 = ax * ax
                # psi(x) = ac * (R1 + a2*(R2 + a2*(R3 + a2*R4)))
                # wgt = psi(x)/x = (1/k) * (R1 + a2*(R2 + a2*(R3 + a2*R4)))
                # IRWLS weight = wgt(r/s) / ((r/s)/c) where c absorbs into wgt
                rho_p = R1 + a2 * (R2 + a2 * (R3 + a2 * R4))
                out[i] = rho_p / 3.25
            else:
                out[i] = 1.0 / 3.25
    else:  # FAM_LQQ
        b_l = tuning[0]; c = tuning[1]; s_l = tuning[2]
        k01 = b_l + c
        s5 = s_l - 1.0
        s6 = -2.0 * k01 + b_l * s_l
        k01_2 = k01 * k01
        denom = s_l * c * (3.0 * c + 2.0 * b_l) + k01_2
        if s5 == 0.0:
            end3 = k01
        else:
            end3 = k01 - s6 / s5
        for i in range(n):
            xi = r[i] / s
            ax = xi if xi >= 0 else -xi
            if ax <= c:
                out[i] = 1.0
            elif ax <= k01:
                s0 = ax - c
                # psi_lqq(x)/x = 1 - (s_l/b_l) * ((ax-c)^2/2) * sign(x)/x simplified
                # Equivalent: wgt = 1 - s_l * s0^2 / (b_l * ax) for ax in (c, k01]
                if ax > 0:
                    out[i] = 1.0 - s_l * s0 * s0 / (b_l * ax)
                else:
                    out[i] = 1.0
            elif ax < end3:
                dx = ax - k01
                # psi_lqq on the outer plateau-decay branch
                if ax > 0:
                    out[i] = (k01 - dx * (s5 + dx * s5 * s5 / s6) - c) / ax
                    if out[i] < 0:
                        out[i] = 0.0
                else:
                    out[i] = 0.0
            else:
                out[i] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _mscale_generic(
    const double* r,
    Py_ssize_t n,
    double init_scale,
    Py_ssize_t p,
    int family,
    const double* tuning,
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
        sum_chi = _chi_sum(r, n, s, family, tuning)
        s = prev * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev
        if diff < 0:
            diff = -diff
        if diff <= tol * prev:
            return s
        prev = s
    return s


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_resample_iter(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[long, ndim=1, mode="c"] idx,
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
    cdef long* idx_data = <long*>cnp.PyArray_DATA(idx)
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


# Backwards-compatible alias kept for the v0.3.0 API.
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
    """Backwards-compatible alias: bisquare with single tuning parameter."""
    cdef cnp.ndarray[double, ndim=1, mode="c"] tuning = np.zeros(3, dtype=np.float64)
    tuning[0] = k_chi
    return cy_resample_iter(
        X, y, idx, 0, tuning, b0, k_fast_s, max_iter_scale, scale_tol, beta_out
    )
