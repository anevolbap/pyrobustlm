# cython: language_level=3
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Cython-accelerated bisquare psi/wgt/rho/psi_prime + the m_scale iteration.
# Used by the resampling loop in fast_s for the (extremely common) bisquare
# default. Pure-NumPy implementations in pyrobustlm._psifuns remain canonical.

cimport cython
from libc.math cimport fabs, sqrt, exp, pow as cpow
import numpy as np
cimport numpy as cnp

cnp.import_array()

# Largest x such that exp(-x^2/2) does not underflow (matches lmrob.c:945).
cdef double _MAX_EX2_SQR_HALF = 37.7 * 37.7 / 2.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_psi(double[::1] x, double k, double[::1] out):
    """Tukey bisquare psi: x * (1 - (x/k)^2)^2 for |x| <= k, else 0."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double a, u
    for i in range(n):
        a = x[i] / k
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            u = 1.0 - a * a
            out[i] = x[i] * u * u


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_wgt(double[::1] x, double k, double[::1] out):
    """Bisquare weight: psi(x)/x = (1 - (x/k)^2)^2 for |x| <= k, else 0."""
    cdef Py_ssize_t i, n = x.shape[0]
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
def bisquare_rho(double[::1] x, double k, double[::1] out):
    """Normalised chi: t * (3 - 3t + t^2) for t = (x/k)^2 <= 1, else 1."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double t, ax
    for i in range(n):
        ax = x[i] if x[i] >= 0 else -x[i]
        if ax > k:
            out[i] = 1.0
        else:
            t = (x[i] / k) * (x[i] / k)
            out[i] = t * (3.0 + t * (-3.0 + t))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bisquare_psi_prime(double[::1] x, double k, double[::1] out):
    """psi'(x; k) = (1 - (x/k)^2)(1 - 5(x/k)^2)  for |x| <= k, else 0."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double a, a2
    for i in range(n):
        a = x[i] / k
        if a < -1.0 or a > 1.0:
            out[i] = 0.0
        else:
            a2 = a * a
            out[i] = (1.0 - a2) * (1.0 - 5.0 * a2)


# ---------------------------------------------------------------------------
# M-scale inner iteration for bisquare. Returns the converged scale.
# ---------------------------------------------------------------------------
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def m_scale_bisquare(
    double[::1] r,
    double k,
    double b0,
    double init_scale,
    int max_iter,
    double tol,
    int p,
):
    """One M-scale iteration loop for bisquare chi.

    Returns the converged sigma satisfying ``mean_{n-p} chi(r/sigma) = b0``.
    """
    cdef Py_ssize_t i, it, n = r.shape[0]
    cdef double s = init_scale
    cdef double prev = init_scale
    cdef double a, ax, t, sum_chi, diff
    cdef double inv_npmp = 1.0 / (n - p)
    if s <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = 0.0
        for i in range(n):
            a = r[i] / s
            ax = a if a >= 0.0 else -a
            if ax > k:
                sum_chi += 1.0
            else:
                t = (a / k) * (a / k)
                sum_chi += t * (3.0 + t * (-3.0 + t))
        s = s * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev if s >= prev else prev - s
        if diff <= tol * prev:
            return s
        prev = s
    return s


# ===========================================================================
# GGW polynomial rho (table-driven). The 6 hardcoded coefficient tables
# from lmrob.c::rho_ggw, plus an inlined m_scale.
# ===========================================================================
cdef double _GGW_C0[6]
cdef double _GGW_END[6]
# 6 cases x 20 polynomial coefficients each.
cdef double _GGW_POLY[6][20]

cdef int _ggw_tables_init = 0


cdef _init_ggw_tables():
    global _ggw_tables_init
    cdef double[20] case0 = [
        0.094164571656733, -0.168937372816728, 0.00427612218326869,
        0.336876420549802, -0.166472338873754, 0.0436904383670537,
        -0.00732077121233756, 0.000792550423837942, -5.08385693557726e-05,
        1.46908724988936e-06, -0.837547853001024, 0.876392734183528,
        -0.184600387321924, 0.0219985685280105, -0.00156403138825785,
        6.16243137719362e-05, -7.478979895101e-07, -3.99563057938975e-08,
        1.78125589532002e-09, -2.22317669250326e-11
    ]
    cdef double[20] case1 = [
        0.174505224068561, -0.168853188892986, 0.00579250806463694,
        0.624193375180937, -0.419882092234336, 0.150011303015251,
        -0.0342185249354937, 0.00504325944243195, -0.0004404209084091,
        1.73268448820236e-05, -0.842160072154898, 1.19912623576069,
        -0.345595777445623, 0.0566407000764478, -0.00560501531439071,
        0.000319084704541442, -7.4279004383686e-06, -2.02063746721802e-07,
        1.65716101809839e-08, -2.97536178313245e-10
    ]
    cdef double[20] case2 = [
        1.41117142330711, -0.168853741371095, 0.0164713906344165,
        5.04767833986545, -9.65574752971554, 9.80999125035463,
        -6.36344090274658, 2.667031271863, -0.662324374141645,
        0.0740982983873332, -0.84794906554363, 3.4315790970352,
        -2.82958670601597, 1.33442885893807, -0.384812004961396,
        0.0661359078129487, -0.00557221619221031, -5.42574872792348e-05,
        4.92564168111658e-05, -2.80432020951381e-06
    ]
    cdef double[20] case3 = [
        0.104604570079252, 0.0626649856211545, -0.220058184826331,
        0.403388189975896, -0.213020713708997, 0.102623342948069,
        -0.0392618698058543, 0.00937878752829234, -0.00122303709506374,
        6.70669880352453e-05, 0.632651530179424, -1.14744323908043,
        0.981941598165897, -0.341211275272191, 0.0671272892644464,
        -0.00826237596187364, 0.0006529134641922, -3.23468516804340e-05,
        9.17904701930209e-07, -1.14119059405971e-08
    ]
    cdef double[20] case4 = [
        0.205026436642222, 0.0627464477520301, -0.308483319391091,
        0.791480474953874, -0.585521414631968, 0.394979618040607,
        -0.211512515412973, 0.0707208739858416, -0.0129092527174621,
        0.000990938134086886, 0.629919019245325, -1.60049136444912,
        1.91903069049618, -0.933285960363159, 0.256861783311473,
        -0.0442133943831343, 0.00488402902512139, -0.000338084604725483,
        1.33974565571893e-05, -2.32450916247553e-07
    ]
    cdef double[20] case5 = [
        1.35010856132000, 0.0627465630782482, -0.791613168488525,
        5.21196700244212, -9.89433796586115, 17.1277266427962,
        -23.5364159883776, 20.1943966645350, -9.4593988142692,
        1.86332355622445, 0.62986381140768, -4.10676399816156,
        12.6361433997327, -15.7697199271455, 11.1373468568838,
        -4.91933095295458, 1.39443093325178, -0.247689078940725,
        0.0251861553415515, -0.00112130382664914
    ]
    cdef Py_ssize_t i
    for i in range(20):
        _GGW_POLY[0][i] = case0[i]
        _GGW_POLY[1][i] = case1[i]
        _GGW_POLY[2][i] = case2[i]
        _GGW_POLY[3][i] = case3[i]
        _GGW_POLY[4][i] = case4[i]
        _GGW_POLY[5][i] = case5[i]
    _GGW_C0[0] = 1.694
    _GGW_C0[1] = 1.2442567
    _GGW_C0[2] = 0.4375470
    _GGW_C0[3] = 1.063
    _GGW_C0[4] = 0.7593544
    _GGW_C0[5] = 0.2959132
    _GGW_END[0] = 18.5527638190955
    _GGW_END[1] = 13.7587939698492
    _GGW_END[2] = 4.89447236180905
    _GGW_END[3] = 11.4974874371859
    _GGW_END[4] = 8.15075376884422
    _GGW_END[5] = 3.17587939698492
    _ggw_tables_init = 1


cdef inline double _ggw_rho_one(double x, int j) nogil:
    """Polynomial chi for ggw case j (0-based). Mirrors lmrob.c::rho_ggw."""
    cdef double ax = x if x >= 0 else -x
    cdef double c = _GGW_C0[j]
    cdef double end = _GGW_END[j]
    cdef double res
    if ax <= c:
        return _GGW_POLY[j][0] * ax * ax
    if ax <= 3 * c:
        # Horner on coeffs [1..9]
        res = _GGW_POLY[j][9]
        res = res * ax + _GGW_POLY[j][8]
        res = res * ax + _GGW_POLY[j][7]
        res = res * ax + _GGW_POLY[j][6]
        res = res * ax + _GGW_POLY[j][5]
        res = res * ax + _GGW_POLY[j][4]
        res = res * ax + _GGW_POLY[j][3]
        res = res * ax + _GGW_POLY[j][2]
        res = res * ax + _GGW_POLY[j][1]
        return res
    if ax <= end:
        res = _GGW_POLY[j][19]
        res = res * ax + _GGW_POLY[j][18]
        res = res * ax + _GGW_POLY[j][17]
        res = res * ax + _GGW_POLY[j][16]
        res = res * ax + _GGW_POLY[j][15]
        res = res * ax + _GGW_POLY[j][14]
        res = res * ax + _GGW_POLY[j][13]
        res = res * ax + _GGW_POLY[j][12]
        res = res * ax + _GGW_POLY[j][11]
        res = res * ax + _GGW_POLY[j][10]
        return res
    return 1.0


def ggw_rho_case(double[::1] x, int case_idx, double[::1] out):
    """Vectorised polynomial chi for ggw case ``case_idx`` (1..6)."""
    if not _ggw_tables_init:
        _init_ggw_tables()
    cdef int j = case_idx - 1
    cdef Py_ssize_t i, n = x.shape[0]
    for i in range(n):
        out[i] = _ggw_rho_one(x[i], j)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def m_scale_ggw_case(
    double[::1] r,
    int case_idx,
    double b0,
    double init_scale,
    int max_iter,
    double tol,
    int p,
):
    """M-scale loop with ggw polynomial chi inlined for one of the 6 cases."""
    if not _ggw_tables_init:
        _init_ggw_tables()
    cdef int j = case_idx - 1
    cdef Py_ssize_t i, it, n = r.shape[0]
    cdef double s = init_scale, prev = init_scale
    cdef double sum_chi, diff
    cdef double inv_npmp = 1.0 / (n - p)
    if s <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = 0.0
        for i in range(n):
            sum_chi += _ggw_rho_one(r[i] / s, j)
        s = s * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev if s >= prev else prev - s
        if diff <= tol * prev:
            return s
        prev = s
    return s


# ===========================================================================
# Inlined m_scale variants for non-bisquare families. Each fully inlines
# its chi function so the iteration loop stays in C without per-iteration
# Python dispatch.
# ===========================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def m_scale_hampel(
    double[::1] r,
    double a,
    double b_,
    double rr,
    double b0,
    double init_scale,
    int max_iter,
    double tol,
    int p,
):
    """M-scale loop with hampel chi inlined."""
    cdef Py_ssize_t i, it, n = r.shape[0]
    cdef double s = init_scale, prev = init_scale
    cdef double xi, u, sum_chi, diff
    cdef double inv_npmp = 1.0 / (n - p)
    cdef double nc = a * (b_ + rr - a) * 0.5
    cdef double inv_nc = 1.0 / nc
    cdef double bma_inv_half = 0.5 / (rr - b_)
    if s <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = 0.0
        for i in range(n):
            xi = r[i] / s
            u = xi if xi >= 0.0 else -xi
            if u <= a:
                sum_chi += (xi * xi * 0.5) * inv_nc
            elif u <= b_:
                sum_chi += (u - 0.5 * a) * a * inv_nc
            elif u <= rr:
                sum_chi += (b_ - 0.5 * a + (u - b_) * (1.0 - (u - b_) * bma_inv_half)) * a * inv_nc
            else:
                sum_chi += 1.0
        s = s * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev if s >= prev else prev - s
        if diff <= tol * prev:
            return s
        prev = s
    return s


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def m_scale_optimal(
    double[::1] r,
    double k,
    double b0,
    double init_scale,
    int max_iter,
    double tol,
    int p,
):
    """M-scale loop with optimal chi inlined."""
    cdef Py_ssize_t i, it, n = r.shape[0]
    cdef double s = init_scale, prev = init_scale
    cdef double R1h = -1.944 / 2.0, R2h = 1.728 / 4.0, R3h = -0.312 / 6.0, R4h = 0.016 / 8.0
    cdef double xi, ac, ax, a2, sum_chi, diff
    cdef double inv_npmp = 1.0 / (n - p)
    if s <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = 0.0
        for i in range(n):
            xi = r[i] / s
            ac = xi / k
            ax = ac if ac >= 0 else -ac
            if ax > 3.0:
                sum_chi += 1.0
            elif ax > 2.0:
                a2 = ax * ax
                sum_chi += (a2 * (R1h + a2 * (R2h + a2 * (R3h + a2 * R4h))) + 1.792) / 3.25
            else:
                sum_chi += (ac * ac) / 6.5
        s = s * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev if s >= prev else prev - s
        if diff <= tol * prev:
            return s
        prev = s
    return s


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def m_scale_lqq(
    double[::1] r,
    double b_,
    double c,
    double s_,
    double b0,
    double init_scale,
    int max_iter,
    double tol,
    int p,
):
    """M-scale loop with lqq chi inlined."""
    cdef Py_ssize_t i, it, n = r.shape[0]
    cdef double s = init_scale, prev = init_scale
    cdef double xi, ax, dx, s0, sum_chi, diff
    cdef double inv_npmp = 1.0 / (n - p)
    cdef double k01 = b_ + c
    cdef double s5 = s_ - 1.0
    cdef double s6 = -2.0 * k01 + b_ * s_
    cdef double k01_2 = k01 * k01
    cdef double denom = s_ * c * (3.0 * c + 2.0 * b_) + k01_2
    cdef double end3
    if s5 == 0.0:
        end3 = k01
    else:
        end3 = k01 - s6 / s5
    if init_scale <= 0.0:
        return 0.0
    for it in range(max_iter):
        sum_chi = 0.0
        for i in range(n):
            xi = r[i] / s
            ax = xi if xi >= 0 else -xi
            if ax <= c:
                sum_chi += (3.0 * s_ - 3.0) / denom * xi * xi
            elif ax <= k01:
                s0 = ax - c
                sum_chi += (6.0 * s_ - 6.0) / denom * (xi * xi * 0.5 - s_ / b_ * s0 * s0 * s0 / 6.0)
            elif ax < end3:
                dx = ax - k01
                sum_chi += (6.0 * s5) / denom * (
                    k01_2 * 0.5 - s_ * b_ * b_ / 6.0
                    - dx * 0.5 * (s6 + dx * (s5 + dx * s5 * s5 / 3.0 / s6))
                )
            else:
                sum_chi += 1.0
        s = s * sqrt(sum_chi * inv_npmp / b0)
        diff = s - prev if s >= prev else prev - s
        if diff <= tol * prev:
            return s
        prev = s
    return s


# ===========================================================================
# Huber
# ===========================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def huber_psi(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double xi
    for i in range(n):
        xi = x[i]
        if xi <= -k:
            out[i] = -k
        elif xi < k:
            out[i] = xi
        else:
            out[i] = k


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def huber_wgt(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double ax
    for i in range(n):
        ax = fabs(x[i])
        if ax >= k:
            # Avoid div-by-zero; if ax==0 we'd be on the inner branch.
            out[i] = k / (ax if ax > 0 else 1e-300)
        else:
            out[i] = 1.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def huber_rho(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double ax
    for i in range(n):
        ax = fabs(x[i])
        if ax <= k:
            out[i] = 0.5 * x[i] * x[i]
        else:
            out[i] = k * (ax - 0.5 * k)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def huber_psi_prime(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    for i in range(n):
        out[i] = 0.0 if fabs(x[i]) >= k else 1.0


# ===========================================================================
# Hampel  (k = (a, b, r))
# ===========================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def hampel_psi(double[::1] x, double a, double b, double r, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double xi, u, sx
    cdef double bma_inv = 1.0 / (r - b)
    for i in range(n):
        xi = x[i]
        u = xi if xi >= 0 else -xi
        sx = 1.0 if xi >= 0 else -1.0
        if u <= a:
            out[i] = xi
        elif u <= b:
            out[i] = sx * a
        elif u <= r:
            out[i] = sx * a * (r - u) * bma_inv
        else:
            out[i] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def hampel_wgt(double[::1] x, double a, double b, double r, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double u
    cdef double bma_inv = 1.0 / (r - b)
    for i in range(n):
        u = fabs(x[i])
        if u <= a:
            out[i] = 1.0
        elif u <= b:
            out[i] = a / (u if u > 0 else 1e-300)
        elif u <= r:
            out[i] = a * (r - u) * bma_inv / (u if u > 0 else 1e-300)
        else:
            out[i] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def hampel_rho(double[::1] x, double a, double b, double r, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double xi, u
    cdef double nc = a * (b + r - a) * 0.5
    cdef double inv_nc = 1.0 / nc
    cdef double bma_inv_half = 0.5 / (r - b)
    for i in range(n):
        xi = x[i]
        u = xi if xi >= 0 else -xi
        if u <= a:
            out[i] = (xi * xi * 0.5) * inv_nc
        elif u <= b:
            out[i] = (u - 0.5 * a) * a * inv_nc
        elif u <= r:
            out[i] = (b - 0.5 * a + (u - b) * (1.0 - (u - b) * bma_inv_half)) * a * inv_nc
        else:
            out[i] = 1.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def hampel_psi_prime(double[::1] x, double a, double b, double r, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double u
    cdef double slope = a / (b - r)
    for i in range(n):
        u = fabs(x[i])
        if u <= a:
            out[i] = 1.0
        elif u <= b:
            out[i] = 0.0
        elif u <= r:
            out[i] = slope
        else:
            out[i] = 0.0


# ===========================================================================
# Optimal (Yohai). Single tuning constant.
# ===========================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def optimal_psi(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016
    cdef double ac, ax, a2, poly, val
    for i in range(n):
        ac = x[i] / k
        ax = ac if ac >= 0 else -ac
        if ax > 3.0:
            out[i] = 0.0
        elif ax > 2.0:
            a2 = ac * ac
            poly = ((R4 * a2 + R3) * a2 + R2) * a2 + R1
            val = k * poly * ac
            if ac > 0.0:
                out[i] = val if val > 0.0 else 0.0
            else:
                out[i] = -fabs(val)
        else:
            out[i] = x[i]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def optimal_wgt(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016
    cdef double ax, a2, val
    for i in range(n):
        ax = x[i] / k
        if ax < 0:
            ax = -ax
        if ax > 3.0:
            out[i] = 0.0
        elif ax > 2.0:
            a2 = ax * ax
            val = R1 + a2 * (R2 + a2 * (R3 + a2 * R4))
            out[i] = val if val > 0.0 else 0.0
        else:
            out[i] = 1.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def optimal_rho(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double R1h = -1.944 / 2.0, R2h = 1.728 / 4.0, R3h = -0.312 / 6.0, R4h = 0.016 / 8.0
    cdef double ac, ax, a2
    for i in range(n):
        ac = x[i] / k
        ax = ac if ac >= 0 else -ac
        if ax > 3.0:
            out[i] = 1.0
        elif ax > 2.0:
            a2 = ax * ax
            out[i] = (a2 * (R1h + a2 * (R2h + a2 * (R3h + a2 * R4h))) + 1.792) / 3.25
        else:
            out[i] = (ac * ac) / 6.5


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def optimal_psi_prime(double[::1] x, double k, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016
    cdef double ax, a2
    for i in range(n):
        ax = x[i] / k
        if ax < 0:
            ax = -ax
        if ax > 3.0:
            out[i] = 0.0
        elif ax > 2.0:
            a2 = ax * ax
            out[i] = R1 + a2 * (3.0 * R2 + a2 * (5.0 * R3 + a2 * 7.0 * R4))
        else:
            out[i] = 1.0


# ===========================================================================
# LQQ (Linear-Quadratic-Quadratic). k = (b, c, s)
# ===========================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lqq_psi(double[::1] x, double b, double c, double s, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double xi, ax, sx
    cdef double k01 = b + c
    cdef double s5 = s - 1.0
    cdef double s6 = -2.0 * k01 + b * s
    cdef double end3
    cdef double dx
    if s5 == 0.0:
        end3 = k01
    else:
        end3 = k01 - s6 / s5
    for i in range(n):
        xi = x[i]
        ax = xi if xi >= 0 else -xi
        sx = 1.0 if xi >= 0 else -1.0
        if ax <= c:
            out[i] = xi
        elif ax <= k01:
            out[i] = sx * (ax - s * (ax - c) * (ax - c) / b * 0.5)
        elif ax < end3:
            dx = ax - k01
            out[i] = sx * (
                -s6 * 0.5
                - (s5 * s5) / s6 * (dx * dx * 0.5 + s6 / s5 * dx)
            )
        else:
            out[i] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lqq_wgt(double[::1] x, double b, double c, double s, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double ax, s0, dx
    cdef double k01 = b + c
    cdef double s5 = s - 1.0
    cdef double s6 = -2.0 * k01 + b * s
    cdef double end3
    if s5 == 0.0:
        end3 = k01
    else:
        end3 = k01 - s6 / s5
    for i in range(n):
        ax = fabs(x[i])
        if ax <= c:
            out[i] = 1.0
        elif ax <= k01:
            s0 = ax - c
            out[i] = 1.0 - s * s0 * s0 / (2.0 * ax * b)
        elif ax < end3:
            dx = ax - k01
            out[i] = -(
                s6 * 0.5 + s5 * s5 / s6 * dx * (dx * 0.5 + s6 / s5)
            ) / (ax if ax > 0 else 1e-300)
        else:
            out[i] = 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lqq_rho(double[::1] x, double b, double c, double s, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double xi, ax, s0, dx
    cdef double k01 = b + c
    cdef double s5 = s - 1.0
    cdef double s6 = -2.0 * k01 + b * s
    cdef double k01_2 = k01 * k01
    cdef double denom = s * c * (3.0 * c + 2.0 * b) + k01_2
    cdef double end3
    if s5 == 0.0:
        end3 = k01
    else:
        end3 = k01 - s6 / s5
    for i in range(n):
        xi = x[i]
        ax = xi if xi >= 0 else -xi
        if ax <= c:
            out[i] = (3.0 * s - 3.0) / denom * xi * xi
        elif ax <= k01:
            s0 = ax - c
            out[i] = (6.0 * s - 6.0) / denom * (
                xi * xi * 0.5 - s / b * s0 * s0 * s0 / 6.0
            )
        elif ax < end3:
            dx = ax - k01
            out[i] = (6.0 * s5) / denom * (
                k01_2 * 0.5
                - s * b * b / 6.0
                - dx * 0.5 * (
                    s6 + dx * (s5 + dx * s5 * s5 / 3.0 / s6)
                )
            )
        else:
            out[i] = 1.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lqq_psi_prime(double[::1] x, double b, double c, double s, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double ax, dx
    cdef double k01 = b + c
    cdef double s5p = 1.0 - s
    cdef double a_param
    if s5p == 0.0:
        for i in range(n):
            ax = fabs(x[i])
            if ax <= c:
                out[i] = 1.0
            elif ax <= k01:
                out[i] = 1.0 - s * (ax - c) / b
            else:
                out[i] = 0.0
        return
    a_param = (b * s - 2.0 * k01) / s5p
    for i in range(n):
        ax = fabs(x[i])
        if ax <= c:
            out[i] = 1.0
        elif ax <= k01:
            out[i] = 1.0 - s * (ax - c) / b
        elif ax < k01 + a_param:
            dx = ax - k01
            out[i] = -s5p * (dx / a_param - 1.0)
        else:
            out[i] = 0.0


# ===========================================================================
# GGW (Generalised Gauss-Weight) - psi/wgt closed form, given (a, b, c).
# Caller resolves the (case_idx -> a, b, c) mapping.
# ===========================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def ggw_psi(double[::1] x, double a, double b, double c, double[::1] out):
    """psi(x) = x * exp(-(|x|-c)+^b / (2 a))"""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double ax, diff, arg
    for i in range(n):
        ax = fabs(x[i])
        if ax <= c:
            out[i] = x[i]
        else:
            diff = ax - c
            arg = cpow(diff, b) / (2.0 * a)
            if arg > _MAX_EX2_SQR_HALF:
                arg = _MAX_EX2_SQR_HALF
            out[i] = x[i] * exp(-arg)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def ggw_wgt(double[::1] x, double a, double b, double c, double[::1] out):
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double ax, diff, arg
    for i in range(n):
        ax = fabs(x[i])
        if ax <= c:
            out[i] = 1.0
        else:
            diff = ax - c
            arg = cpow(diff, b) / (2.0 * a)
            if arg > _MAX_EX2_SQR_HALF:
                arg = _MAX_EX2_SQR_HALF
            out[i] = exp(-arg)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def ggw_psi_prime(double[::1] x, double a, double b, double c, double[::1] out):
    """d/dx of ggw_psi (lmrob.c:1304-1317)."""
    cdef Py_ssize_t i, n = x.shape[0]
    cdef double xi, ax, diff, arg, e, bracket
    cdef double inv_2a = 1.0 / (2.0 * a)
    for i in range(n):
        xi = x[i]
        ax = fabs(xi)
        if ax < c:
            out[i] = 1.0
            continue
        diff = ax - c
        arg = cpow(diff, b) * inv_2a
        if arg > _MAX_EX2_SQR_HALF:
            out[i] = 0.0
            continue
        e = exp(-arg)
        bracket = 1.0 - (b * inv_2a) * ax * cpow(diff, b - 1.0)
        out[i] = e * bracket

