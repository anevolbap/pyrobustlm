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
from libc.math cimport fabs, sqrt, exp, pow as cpow
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free

from numpy.random cimport bitgen_t
from scipy.linalg.cython_lapack cimport dgels, dgesv

import numpy as np
cimport numpy as cnp

cnp.import_array()


# Largest x such that exp(-x^2/2) does not underflow (matches lmrob.c:945).
cdef double _MAX_EX2_SQR_HALF = 37.7 * 37.7 / 2.0


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
# GGW polynomial chi tables (same data as _psi.pyx / _fast_s.pyx). Lazily
# initialised on first ggw call.
# ---------------------------------------------------------------------------
cdef double _GGW_C0[6]
cdef double _GGW_END[6]
cdef double _GGW_POLY[6][20]
cdef double _GGW_ABC_A[7]
cdef double _GGW_ABC_B[7]
cdef double _GGW_ABC_C[7]
cdef int _ggw_tables_init = 0


cdef _init_ggw_tables():
    global _ggw_tables_init
    if _ggw_tables_init:
        return
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
    _GGW_ABC_A[0] = 0.0; _GGW_ABC_B[0] = 0.0; _GGW_ABC_C[0] = 0.0
    _GGW_ABC_A[1] = 0.648;     _GGW_ABC_B[1] = 1.0; _GGW_ABC_C[1] = 1.694
    _GGW_ABC_A[2] = 0.4760508; _GGW_ABC_B[2] = 1.0; _GGW_ABC_C[2] = 1.2442567
    _GGW_ABC_A[3] = 0.1674046; _GGW_ABC_B[3] = 1.0; _GGW_ABC_C[3] = 0.4375470
    _GGW_ABC_A[4] = 1.387;     _GGW_ABC_B[4] = 1.5; _GGW_ABC_C[4] = 1.063
    _GGW_ABC_A[5] = 0.8372485; _GGW_ABC_B[5] = 1.5; _GGW_ABC_C[5] = 0.7593544
    _GGW_ABC_A[6] = 0.2036741; _GGW_ABC_B[6] = 1.5; _GGW_ABC_C[6] = 0.2959132
    _ggw_tables_init = 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _ggw_rho_one(double x, int j) nogil:
    cdef double ax = x if x >= 0 else -x
    cdef double c = _GGW_C0[j]
    cdef double end = _GGW_END[j]
    cdef double res
    if ax <= c:
        return _GGW_POLY[j][0] * ax * ax
    if ax <= 3 * c:
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


# ---------------------------------------------------------------------------
# Per-family chi-sum and IRWLS weight. Mirror the kernels in _fast_s.pyx.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _chi_sum(
    const double* r, Py_ssize_t n, double s, int family, const double* tuning,
) nogil:
    cdef Py_ssize_t i
    cdef double total = 0.0
    cdef double a, ax, t, xi, u, dx, s0
    cdef double k, a_t, b_t, r_t, c, b_l, s_l, k01, denom, s5, s6, k01_2, end3
    cdef double R1h, R2h, R3h, R4h, ac, a2, nc, inv_nc, bma_inv_half
    cdef int j

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
    elif family == FAM_LQQ:
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
    else:  # FAM_GGW
        j = <int>(tuning[0]) - 1
        if j < 0:
            j = 0
        elif j > 5:
            j = 5
        for i in range(n):
            total += _ggw_rho_one(r[i] / s, j)
    return total


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _wgt_zinv(
    const double* r, double* out, Py_ssize_t n, double s,
    int family, const double* tuning,
) nogil:
    cdef Py_ssize_t i
    cdef double a, u, ax, xi, ac, a2, rho_p, dx, s0
    cdef double k, inv_sk
    cdef double a_t, b_t, r_t
    cdef double R1, R2, R3, R4
    cdef double b_l, c, s_l, k01, s5, s6, end3
    cdef int j

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
            ac = (r[i] / s) / k
            ax = ac if ac >= 0 else -ac
            if ax > 3.0:
                out[i] = 0.0
            elif ax > 2.0:
                a2 = ax * ax
                rho_p = R1 + a2 * (R2 + a2 * (R3 + a2 * R4))
                out[i] = rho_p if rho_p > 0.0 else 0.0
            else:
                out[i] = 1.0
    elif family == FAM_LQQ:
        b_l = tuning[0]; c = tuning[1]; s_l = tuning[2]
        k01 = b_l + c
        s5 = s_l - 1.0
        s6 = -2.0 * k01 + b_l * s_l
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
                if ax > 0:
                    out[i] = 1.0 - s_l * s0 * s0 / (2.0 * ax * b_l)
                else:
                    out[i] = 1.0
            elif ax < end3:
                dx = ax - k01
                if ax > 0:
                    out[i] = -(
                        s6 * 0.5 + s5 * s5 / s6 * dx * (dx * 0.5 + s6 / s5)
                    ) / ax
                else:
                    out[i] = 0.0
            else:
                out[i] = 0.0
    else:  # FAM_GGW
        j = <int>(tuning[0])
        if j < 1:
            j = 1
        elif j > 6:
            j = 6
        a_t = _GGW_ABC_A[j]
        b_t = _GGW_ABC_B[j]
        r_t = _GGW_ABC_C[j]  # c
        for i in range(n):
            xi = r[i] / s
            ax = xi if xi >= 0 else -xi
            if ax <= r_t:
                out[i] = 1.0
            else:
                dx = ax - r_t
                ac = cpow(dx, b_t) / (2.0 * a_t)
                if ac > _MAX_EX2_SQR_HALF:
                    ac = _MAX_EX2_SQR_HALF
                out[i] = exp(-ac)


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
    double s, int family, const double* tuning,
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

    _wgt_zinv(r_in, scr.w, n, s, family, tuning)
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
    double init_s, int family, const double* tuning, double b0,
    int k_fast_s, int max_iter_scale, double scale_tol,
    _Scratch* scr,
    Py_ssize_t n, Py_ssize_t p,
    int* status_out,
) nogil:
    cdef double s = init_s
    cdef int kk, irwls_status
    cdef int status = 0

    for kk in range(k_fast_s):
        s = _mscale_generic(scr.r, n, s, p, family, tuning, b0,
                            max_iter_scale, scale_tol)
        if s == 0.0:
            status = 2
            status_out[0] = status
            return 0.0
        irwls_status = _irwls_step(X, y, scr.r, s, family, tuning,
                                   scr.beta, scr, n, p)
        if irwls_status != 0:
            status = irwls_status
            status_out[0] = status
            return s
        _residuals(X, y, scr.beta, scr.r, n, p)
    s = _mscale_generic(scr.r, n, s, p, family, tuning, b0,
                        max_iter_scale, scale_tol)
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
    double init_s, int family, const double* tuning, double b0,
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
        s = _mscale_generic(scr.r, n, s, p, family, tuning, b0,
                            max_iter_scale, scale_tol)
        if s == 0.0:
            converged_out[0] = 1
            n_iter_out[0] = it + 1
            return 0.0
        # Save previous beta.
        for j in range(p):
            scr.beta_prev[j] = beta[j]
        irwls_status = _irwls_step(X, y, scr.r, s, family, tuning,
                                   beta, scr, n, p)
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
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
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
    """Fast-S in a single nogil call. Family-generic.

    ``family``: 0=bisquare, 1=hampel, 2=optimal, 3=lqq, 4=ggw.
    ``tuning``: family-specific tuning constants (len 3 is enough for all
    supported families; ggw uses tuning[0] as case_idx).

    Runs the full pipeline: ``nResample`` candidates × (subset draw +
    initial solve + ``k_fast_s`` K-step refinements), best-of-``best_r``
    selection, then survivor refinement to convergence. All in one C
    block with one workspace allocation.

    Returns ``(scale, status, n_iter, converged)`` where status is
    0 (ok), 1 (no non-singular subset found), 2 (exact fit), 3 (LAPACK
    error), 4 (alloc failed).
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()
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
    cdef double* tuning_data = <double*>cnp.PyArray_DATA(tuning)

    cdef Py_ssize_t i, j, row, swap, try_i, kept
    cdef double s, scale, max_abs, candidate_scale
    cdef int got_subset, status, k_status, converged
    cdef int worst_i, total_iters
    cdef int found
    cdef uint64_t r_u
    # Initial-scale divisor: use tuning[0] for bisquare/optimal; for
    # ggw we use _GGW_C0[case-1]; for hampel/lqq tuning[2]/tuning[1]
    # serve as the natural scale. Default to 1.0 if the chosen divisor
    # is non-positive.
    cdef double k0 = tuning_data[0]
    if family == FAM_GGW:
        if k0 < 1: k0 = 1
        elif k0 > 6: k0 = 6
        k0 = _GGW_C0[<int>k0 - 1]
    if k0 <= 0.0:
        k0 = 1.0

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
            s = max_abs / k0
            if s <= 0.0:
                s = 1.0

            candidate_scale = _k_step_refine(
                X_data, y_data, s, family, tuning_data, b0,
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
                        best_scales[i], family, tuning_data, b0,
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


# ---------------------------------------------------------------------------
# Per-family psi'(x) and psi(x) per-element kernels for vcov.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _psi_prime_eval(
    const double* x, double* out, Py_ssize_t n,
    int family, const double* tuning,
) nogil:
    cdef Py_ssize_t i
    cdef double a, a2, xi, ax, ac, slope
    cdef double k, a_t, b_t, r_t
    cdef double R1, R2, R3, R4
    cdef double b_l, c, s_l, k01, s5p, a_param, dx
    cdef int j
    cdef double diff, arg, e, bracket, inv_2a, ax_abs, b_g
    if family == FAM_BISQUARE:
        k = tuning[0]
        for i in range(n):
            a = x[i] / k
            if a < -1.0 or a > 1.0:
                out[i] = 0.0
            else:
                a2 = a * a
                out[i] = (1.0 - a2) * (1.0 - 5.0 * a2)
    elif family == FAM_HAMPEL:
        a_t = tuning[0]; b_t = tuning[1]; r_t = tuning[2]
        slope = a_t / (b_t - r_t)
        for i in range(n):
            ax = x[i] if x[i] >= 0 else -x[i]
            if ax <= a_t:
                out[i] = 1.0
            elif ax <= b_t:
                out[i] = 0.0
            elif ax <= r_t:
                out[i] = slope
            else:
                out[i] = 0.0
    elif family == FAM_OPTIMAL:
        k = tuning[0]
        R1 = -1.944; R2 = 1.728; R3 = -0.312; R4 = 0.016
        for i in range(n):
            ax = x[i] / k
            if ax < 0: ax = -ax
            if ax > 3.0:
                out[i] = 0.0
            elif ax > 2.0:
                a2 = ax * ax
                out[i] = R1 + a2 * (3.0 * R2 + a2 * (5.0 * R3 + a2 * 7.0 * R4))
            else:
                out[i] = 1.0
    elif family == FAM_LQQ:
        b_l = tuning[0]; c = tuning[1]; s_l = tuning[2]
        k01 = b_l + c
        s5p = 1.0 - s_l
        if s5p == 0.0:
            for i in range(n):
                ax = x[i] if x[i] >= 0 else -x[i]
                if ax <= c:
                    out[i] = 1.0
                elif ax <= k01:
                    out[i] = 1.0 - s_l * (ax - c) / b_l
                else:
                    out[i] = 0.0
        else:
            a_param = (b_l * s_l - 2.0 * k01) / s5p
            for i in range(n):
                ax = x[i] if x[i] >= 0 else -x[i]
                if ax <= c:
                    out[i] = 1.0
                elif ax <= k01:
                    out[i] = 1.0 - s_l * (ax - c) / b_l
                elif ax < k01 + a_param:
                    dx = ax - k01
                    out[i] = -s5p * (dx / a_param - 1.0)
                else:
                    out[i] = 0.0
    else:  # FAM_GGW
        j = <int>(tuning[0])
        if j < 1: j = 1
        elif j > 6: j = 6
        a_t = _GGW_ABC_A[j]
        b_g = _GGW_ABC_B[j]
        r_t = _GGW_ABC_C[j]  # c
        inv_2a = 1.0 / (2.0 * a_t)
        for i in range(n):
            xi = x[i]
            ax_abs = xi if xi >= 0 else -xi
            if ax_abs < r_t:
                out[i] = 1.0
                continue
            diff = ax_abs - r_t
            arg = cpow(diff, b_g) * inv_2a
            if arg > _MAX_EX2_SQR_HALF:
                out[i] = 0.0
                continue
            e = exp(-arg)
            bracket = 1.0 - (b_g * inv_2a) * ax_abs * cpow(diff, b_g - 1.0)
            out[i] = e * bracket


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _psi_eval(
    const double* x, double* out, Py_ssize_t n,
    int family, const double* tuning,
) nogil:
    """psi(x) per element. For all families psi(x) = x * wgt(x) (which we
    already have inlined via _wgt_zinv). Compute as x * wgt."""
    cdef Py_ssize_t i
    cdef double* tmp = <double*>malloc(n * sizeof(double))
    if tmp == NULL:
        return
    _wgt_zinv(x, tmp, n, 1.0, family, tuning)
    for i in range(n):
        out[i] = x[i] * tmp[i]
    free(tmp)


cdef inline double _chi_prime_factor(int family, const double* tuning) nogil:
    """chi'(x)/psi(x) = 1/rho_unnorm(inf). Family-specific."""
    cdef double c, b_l, s_l, k01_2, denom
    cdef double a_t, b_t, r_t, nc
    cdef int j
    if family == FAM_BISQUARE:
        c = tuning[0]
        return 6.0 / (c * c)
    if family == FAM_HAMPEL:
        a_t = tuning[0]; b_t = tuning[1]; r_t = tuning[2]
        nc = a_t * (b_t + r_t - a_t) * 0.5
        return 1.0 / nc
    if family == FAM_OPTIMAL:
        c = tuning[0]
        return 1.0 / (3.25 * c * c)
    if family == FAM_LQQ:
        b_l = tuning[0]; c = tuning[1]; s_l = tuning[2]
        k01_2 = (b_l + c) * (b_l + c)
        denom = s_l * c * (3.0 * c + 2.0 * b_l) + k01_2
        return 6.0 * (s_l - 1.0) / denom
    # FAM_GGW: factor depends on case. Same constants as our tabulated
    # asympt_corrfact would imply; for now we use the same as bisquare
    # which is approximate. The four "fast" ggw cases have specific
    # tabulated values in pyrobustlm.inference._asympt_corrfact; we use
    # those.
    j = <int>(tuning[0])
    if j == 1: return 1.0 / 1.6047  # case 1: b=1, 95% eff
    if j == 4: return 1.0 / 1.6047  # case 4: b=1.5, 95% eff
    return 6.0 / (tuning[0] * tuning[0])  # fallback


# ---------------------------------------------------------------------------
# vcov_avar1: ports pyrobustlm.inference.vcov_avar1 into a single nogil
# kernel. Per-element psi/chi via the helpers above; matrix ops via
# dgesv (inversion) and dsyev (posdefify).
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef int _compute_vcov_avar1_body(
    const double* X_data,
    const double* r_data,
    const double* r0_data,
    double sigma,
    int family,
    const double* tpsi_data,
    const double* tchi_data,
    double bb,
    Py_ssize_t n,
    Py_ssize_t p,
    double* cov_data,
) nogil except 3:
    """Shared implementation of vcov_avar1. Allocates its own scratch,
    fills ``cov_data`` (C-contiguous p*p), returns status: 0 ok,
    3 LAPACK error / out of memory.
    """
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    cdef double sgma = sigma
    if sgma < 1e-300:
        sgma = 1e-300

    cdef double* r_s = <double*>malloc(n * sizeof(double))
    cdef double* r0_s = <double*>malloc(n * sizeof(double))
    cdef double* w_pp = <double*>malloc(n * sizeof(double))
    cdef double* w0_pp = <double*>malloc(n * sizeof(double))
    cdef double* psi_rs = <double*>malloc(n * sizeof(double))
    cdef double* chi_r0s = <double*>malloc(n * sizeof(double))
    cdef double* A = <double*>malloc(p * p * sizeof(double))
    cdef double* rhs = <double*>malloc(p * p * sizeof(double))
    cdef double* tmp_pv = <double*>malloc(p * sizeof(double))
    cdef double* a_vec = <double*>malloc(p * sizeof(double))
    cdef double* Xww = <double*>malloc(p * sizeof(double))
    cdef double* u_mat = <double*>malloc(p * p * sizeof(double))
    cdef double* outer_buf = <double*>malloc(p * p * sizeof(double))
    cdef int* ipiv = <int*>malloc(p * sizeof(int))
    if (r_s == NULL or r0_s == NULL or w_pp == NULL or w0_pp == NULL or
            psi_rs == NULL or chi_r0s == NULL or A == NULL or rhs == NULL or
            tmp_pv == NULL or a_vec == NULL or Xww == NULL or
            u_mat == NULL or outer_buf == NULL or ipiv == NULL):
        free(r_s); free(r0_s); free(w_pp); free(w0_pp)
        free(psi_rs); free(chi_r0s); free(A); free(rhs)
        free(tmp_pv); free(a_vec); free(Xww)
        free(u_mat); free(outer_buf); free(ipiv)
        return 3

    cdef double cp_factor = _chi_prime_factor(family, tchi_data)
    cdef Py_ssize_t i, j, k
    cdef double denom, dot, scalar
    cdef int status = 0

    # Scaled residuals.
    for i in range(n):
        r_s[i] = r_data[i] / sgma
        r0_s[i] = r0_data[i] / sgma

    # w_pp = psi'(r_s) (using psi-tuning), psi_rs = psi(r_s).
    _psi_prime_eval(r_s, w_pp, n, family, tpsi_data)
    _psi_eval(r_s, psi_rs, n, family, tpsi_data)

    # chi(r0_s) (chi-tuning); chi'(r0_s) = chi_prime_factor * psi(r0_s).
    # Use _chi_sum logic per-element by exploiting _chi_eval indirectly.
    # Instead, compute psi(r0_s) into chi_r0s as a scratch then scale
    # for w0_pp; compute chi via a separate loop.
    _psi_eval(r0_s, w0_pp, n, family, tchi_data)  # w0_pp = psi(r0_s)
    for i in range(n):
        w0_pp[i] = w0_pp[i] * cp_factor  # chi' = factor * psi

    # chi(r0_s) per-element: ax/k via _chi_sum equivalent. The
    # cleanest is to inline the bisquare/optimal/lqq/hampel/ggw rho.
    # Use the fact that chi(x) = rho_unnormalised(x) / rho_inf,
    # and we have ``cp_factor = 1/rho_inf``, so chi(x) = cp_factor *
    # integral psi from 0 to |x|. For ggw and the others, we'd need
    # the closed-form rho. For now compute via inlining of the
    # family-specific rho.

    # Build chi(r0_s) by re-using _chi_sum logic per element. Easiest:
    # call _chi_sum on each single-element slice. That's awkward; we
    # write a per-family chi_eval below in tandem.
    _chi_eval(r0_s, chi_r0s, n, family, tchi_data)

    # XwX = X^T diag(w_pp) X, column-major (we use it for dgesv).
    for j in range(p):
        for k in range(p):
            dot = 0.0
            for i in range(n):
                dot += X_data[i * p + j] * w_pp[i] * X_data[i * p + k]
            A[j + k * p] = dot
            rhs[j + k * p] = 1.0 if j == k else 0.0

    dgesv(&p_int, &p_int, A, &p_int, ipiv, rhs, &p_int, &info)
    if info != 0:
        status = 3
    else:
        # A = sigma * (X'WX)^{-1}. rhs holds the inverse; multiply by sgma.
        for i in range(p * p):
            A[i] = rhs[i] * sgma

        # denom = mean(w0_pp * r0_s)
        denom = 0.0
        for i in range(n):
            denom += w0_pp[i] * r0_s[i]
        denom /= <double>n
        if denom == 0.0:
            status = 3
        else:
            # tmp_pv = X^T (w_pp * r_s)
            for j in range(p):
                dot = 0.0
                for i in range(n):
                    dot += X_data[i * p + j] * (w_pp[i] * r_s[i])
                tmp_pv[j] = dot

            # a_vec = A @ tmp_pv / denom (A is column-major p×p)
            for j in range(p):
                dot = 0.0
                for k in range(p):
                    dot += A[j + k * p] * tmp_pv[k]
                a_vec[j] = dot / denom

            # Xww = X^T (psi_rs * chi_r0s)
            for j in range(p):
                dot = 0.0
                for i in range(n):
                    dot += X_data[i * p + j] * psi_rs[i] * chi_r0s[i]
                Xww[j] = dot

            # u1 = A @ (X^T diag(psi_rs^2) X) @ (n * A)
            # Step 1: u_mat = X^T diag(psi_rs^2) X (column-major)
            for j in range(p):
                for k in range(p):
                    dot = 0.0
                    for i in range(n):
                        dot += X_data[i * p + j] * (psi_rs[i] * psi_rs[i]) * X_data[i * p + k]
                    u_mat[j + k * p] = dot
            # Step 2: outer_buf = A @ u_mat (column-major p×p)
            for j in range(p):
                for k in range(p):
                    dot = 0.0
                    for i in range(p):
                        dot += A[j + i * p] * u_mat[i + k * p]
                    outer_buf[j + k * p] = dot
            # Step 3: u_mat = outer_buf @ (n * A) = n * outer_buf @ A
            for j in range(p):
                for k in range(p):
                    dot = 0.0
                    for i in range(p):
                        dot += outer_buf[j + i * p] * A[i + k * p]
                    u_mat[j + k * p] = <double>n * dot

            # u2 = outer(a_vec, Xww) @ A. outer(a, Xww)[j,k] = a[j] * Xww[k].
            # (outer_a_Xww @ A)[j, k] = sum_i a[j] * Xww[i] * A[i + k*p]
            #                        = a[j] * (Xww^T @ A_col_k)
            # Pre-compute Xww^T @ A → row vector of length p.
            # tmp_pv reused: tmp_pv[k] = sum_i Xww[i] * A[i + k*p]
            for k in range(p):
                dot = 0.0
                for i in range(p):
                    dot += Xww[i] * A[i + k * p]
                tmp_pv[k] = dot
            # Subtract u2 = a_vec[j] * tmp_pv[k] from u_mat in row j, col k.
            for j in range(p):
                for k in range(p):
                    u_mat[j + k * p] -= a_vec[j] * tmp_pv[k]

            # u3 = A @ outer(Xww, a_vec). outer(Xww, a_vec)[j,k] = Xww[j] * a[k].
            # (A @ outer)[j, k] = sum_i A[j + i*p] * Xww[i] * a[k]
            #                  = (A @ Xww)[j] * a[k]
            for j in range(p):
                dot = 0.0
                for i in range(p):
                    dot += A[j + i * p] * Xww[i]
                outer_buf[j] = dot  # store A @ Xww in first p entries
            for j in range(p):
                for k in range(p):
                    u_mat[j + k * p] -= outer_buf[j] * a_vec[k]

            # u4 = mean(chi_r0s^2 - bb^2) * outer(a_vec, a_vec)
            scalar = 0.0
            for i in range(n):
                scalar += chi_r0s[i] * chi_r0s[i] - bb * bb
            scalar /= <double>n
            for j in range(p):
                for k in range(p):
                    u_mat[j + k * p] += scalar * a_vec[j] * a_vec[k]

            # cov = u_mat / n. cov_out is C-contiguous (row-major), while
            # u_mat is column-major. Convert during the write.
            for j in range(p):
                for k in range(p):
                    cov_data[j * p + k] = u_mat[j + k * p] / <double>n

            # Symmetrize (the math should already give a symmetric matrix
            # up to floating-point error).
            for j in range(p):
                for k in range(j + 1, p):
                    scalar = 0.5 * (cov_data[j * p + k] + cov_data[k * p + j])
                    cov_data[j * p + k] = scalar
                    cov_data[k * p + j] = scalar

    free(r_s); free(r0_s); free(w_pp); free(w0_pp)
    free(psi_rs); free(chi_r0s); free(A); free(rhs)
    free(tmp_pv); free(a_vec); free(Xww)
    free(u_mat); free(outer_buf); free(ipiv)
    return status


def cy_lmrob_vcov_avar1(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] residuals,
    cnp.ndarray[double, ndim=1, mode="c"] init_residuals,
    double sigma,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning_psi,
    cnp.ndarray[double, ndim=1, mode="c"] tuning_chi,
    double bb,
    cnp.ndarray[double, ndim=2, mode="c"] cov_out,
):
    """Asymptotic sandwich covariance per robustbase ``.vcov.avar1``.

    All five lmrob-supported families. Mirrors
    pyrobustlm.inference.vcov_avar1 element-wise.

    Returns ``status``: 0 ok, 3 LAPACK error.
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int status
    with nogil:
        status = _compute_vcov_avar1_body(
            <double*>cnp.PyArray_DATA(X),
            <double*>cnp.PyArray_DATA(residuals),
            <double*>cnp.PyArray_DATA(init_residuals),
            sigma, family,
            <double*>cnp.PyArray_DATA(tuning_psi),
            <double*>cnp.PyArray_DATA(tuning_chi),
            bb, n, p,
            <double*>cnp.PyArray_DATA(cov_out),
        )
    return status


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void _chi_eval(
    const double* x, double* out, Py_ssize_t n,
    int family, const double* tuning,
) nogil:
    """chi(x) (normalised: chi(inf)=1) per element. Inlined formulas
    mirror _chi_sum but write per-element."""
    cdef Py_ssize_t i
    cdef double t, ax, xi, u, ac, a2, dx, s0
    cdef double k, a_t, b_t, r_t, nc, inv_nc, bma_inv_half
    cdef double R1h, R2h, R3h, R4h
    cdef double b_l, c, s_l, k01, denom, s5, s6, k01_2, end3
    cdef int j

    if family == FAM_BISQUARE:
        k = tuning[0]
        for i in range(n):
            t = x[i]
            ax = t if t >= 0 else -t
            if ax >= k:
                out[i] = 1.0
            else:
                t = ax / k
                t = 1.0 - t * t
                out[i] = 1.0 - t * t * t
    elif family == FAM_HAMPEL:
        a_t = tuning[0]; b_t = tuning[1]; r_t = tuning[2]
        nc = a_t * (b_t + r_t - a_t) * 0.5
        inv_nc = 1.0 / nc
        bma_inv_half = 0.5 / (r_t - b_t)
        for i in range(n):
            xi = x[i]
            u = xi if xi >= 0 else -xi
            if u <= a_t:
                out[i] = (xi * xi * 0.5) * inv_nc
            elif u <= b_t:
                out[i] = (u - 0.5 * a_t) * a_t * inv_nc
            elif u <= r_t:
                out[i] = (b_t - 0.5 * a_t + (u - b_t) * (1.0 - (u - b_t) * bma_inv_half)) * a_t * inv_nc
            else:
                out[i] = 1.0
    elif family == FAM_OPTIMAL:
        k = tuning[0]
        R1h = -1.944 * 0.5; R2h = 1.728 * 0.25; R3h = -0.312 / 6.0; R4h = 0.016 / 8.0
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
    elif family == FAM_LQQ:
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
            xi = x[i]
            ax = xi if xi >= 0 else -xi
            if ax <= c:
                out[i] = (3.0 * s_l - 3.0) / denom * xi * xi
            elif ax <= k01:
                s0 = ax - c
                out[i] = (6.0 * s_l - 6.0) / denom * (xi * xi * 0.5 - s_l / b_l * s0 * s0 * s0 / 6.0)
            elif ax < end3:
                dx = ax - k01
                out[i] = (6.0 * s5) / denom * (
                    k01_2 * 0.5 - s_l * b_l * b_l / 6.0
                    - dx * 0.5 * (s6 + dx * (s5 + dx * s5 * s5 / 3.0 / s6))
                )
            else:
                out[i] = 1.0
    else:  # FAM_GGW
        j = <int>(tuning[0]) - 1
        if j < 0: j = 0
        elif j > 5: j = 5
        for i in range(n):
            out[i] = _ggw_rho_one(x[i], j)


# ---------------------------------------------------------------------------
# Design-adaptive D-scale (Koller & Stahel 2014). Mirrors
# robustbase/src/lmrob.c::R_find_D_scale and pyrobustlm.d_scale.
# ---------------------------------------------------------------------------
# Kappa and (tfact, tcorr) for tau, per family at default tuning. Mirrors
# the tables in pyrobustlm.d_scale (_TAU_FAST_TABLE) and lmrob.kappa
# tabulation. For non-default tuning the caller falls back to the Python
# path; this kernel only handles the common defaults.

# Computed via scipy.integrate.quad on the default-tuning psi.r/wgt
# integrand (see pyrobustlm.d_scale.kappa) at runtime in Python; the
# values below were captured for the family defaults and re-used here.
cdef double _DSCALE_KAPPA_BISQUARE = 0.8280771566048320
cdef double _DSCALE_KAPPA_HAMPEL = 0.8569775805834327
cdef double _DSCALE_KAPPA_OPTIMAL = 0.9355077953265407
cdef double _DSCALE_KAPPA_LQQ = 0.8626400360440886
# ggw default cases are case 1 (b=1, 95% eff) and case 4 (b=1.5, 95% eff).
# kappa values measured the same way as the others.
cdef double _DSCALE_KAPPA_GGW_C1 = 0.8914986545654882
cdef double _DSCALE_KAPPA_GGW_C4 = 0.8914986545654882


cdef inline int _dscale_tau_factors(
    int family, const double* tuning, double* out_tfact, double* out_tcorr,
) nogil:
    """Look up (tfact, tcorr) for tau = sqrt(1-tfact*h)*(tcorr*h+1).
    Returns 0 on success, 1 if no tabulated factors are available."""
    cdef int case_idx
    if family == FAM_BISQUARE:
        out_tfact[0] = 0.9473684
        out_tcorr[0] = -0.0900833
        return 0
    if family == FAM_HAMPEL:
        out_tfact[0] = 0.94739770
        out_tcorr[0] = -0.04103958
        return 0
    if family == FAM_OPTIMAL:
        out_tfact[0] = 0.94735878
        out_tcorr[0] = -0.09444537
        return 0
    if family == FAM_LQQ:
        out_tfact[0] = 0.94736359
        out_tcorr[0] = -0.08594805
        return 0
    if family == FAM_GGW:
        case_idx = <int>(tuning[0])
        if case_idx == 1:
            out_tfact[0] = 0.9473787
            out_tcorr[0] = -0.1143846
            return 0
        if case_idx == 4:
            out_tfact[0] = 0.94741036
            out_tcorr[0] = -0.08424648
            return 0
    return 1


cdef inline int _dscale_kappa(
    int family, const double* tuning, double* out_kappa,
) nogil:
    """Tabulated kappa per family. Returns 0 on success, 1 if no table."""
    cdef int case_idx
    if family == FAM_BISQUARE:
        out_kappa[0] = _DSCALE_KAPPA_BISQUARE
        return 0
    if family == FAM_HAMPEL:
        out_kappa[0] = _DSCALE_KAPPA_HAMPEL
        return 0
    if family == FAM_OPTIMAL:
        out_kappa[0] = _DSCALE_KAPPA_OPTIMAL
        return 0
    if family == FAM_LQQ:
        out_kappa[0] = _DSCALE_KAPPA_LQQ
        return 0
    if family == FAM_GGW:
        case_idx = <int>(tuning[0])
        if case_idx == 1:
            out_kappa[0] = _DSCALE_KAPPA_GGW_C1
            return 0
        if case_idx == 4:
            out_kappa[0] = _DSCALE_KAPPA_GGW_C4
            return 0
    return 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int _compute_hat_diagonal(
    const double* X, const double* rweights,
    Py_ssize_t n, Py_ssize_t p,
    double* XwX_inv_buf, int* ipiv, double* h_out,
) nogil:
    """Compute h_i = w_i * X[i,:] (X^T W X)^{-1} X[i,:]^T.

    ``XwX_inv_buf`` must be at least p*p doubles. On return it contains
    (X^T W X)^{-1}. Returns 0 on success, 3 on LAPACK error.
    """
    cdef Py_ssize_t i, j, k
    cdef double w_i, dot, hi
    cdef int p_int = <int>p
    cdef int info = 0

    # Build XwX = X^T diag(w) X (column-major in XwX_inv_buf).
    for j in range(p):
        for k in range(p):
            dot = 0.0
            for i in range(n):
                dot += X[i * p + j] * rweights[i] * X[i * p + k]
            XwX_inv_buf[j + k * p] = dot

    # Identity rhs (p x p, column-major).
    # We solve A * Z = I in place; on return XwX_inv_buf is A^{-1}.
    cdef double* rhs = <double*>malloc(p * p * sizeof(double))
    if rhs == NULL:
        return 3
    for j in range(p):
        for k in range(p):
            rhs[j + k * p] = 1.0 if j == k else 0.0
    dgesv(&p_int, &p_int, XwX_inv_buf, &p_int, ipiv, rhs, &p_int, &info)
    if info != 0:
        free(rhs)
        return 3
    # Copy A^{-1} back to XwX_inv_buf.
    for j in range(p * p):
        XwX_inv_buf[j] = rhs[j]
    free(rhs)

    # h_i = w_i * X[i,:] A^{-1} X[i,:]^T. Quadratic form via two-step.
    for i in range(n):
        hi = 0.0
        for j in range(p):
            dot = 0.0
            for k in range(p):
                dot += XwX_inv_buf[j + k * p] * X[i * p + k]
            hi += X[i * p + j] * dot
        h_out[i] = rweights[i] * hi
        if h_out[i] > 1.0:
            h_out[i] = 1.0
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int _d_scale_iterate(
    const double* r, const double* tau_vec, Py_ssize_t n,
    int family, const double* tuning,
    double kappa_val, double init_scale,
    int max_iter, double tol,
    double* out_scale, int* out_converged,
) nogil:
    cdef double sgma = init_scale
    cdef double new_s, tsum1, tsum2, z, wi
    cdef double* w_arr = <double*>malloc(n * sizeof(double))
    if w_arr == NULL:
        return 3
    cdef Py_ssize_t i
    cdef int it
    if sgma <= 0.0:
        out_scale[0] = 0.0
        out_converged[0] = 0
        free(w_arr)
        return 0
    for it in range(max_iter):
        # Compute scaled residuals r_i / (tau_i * sgma) and weights.
        # Use _wgt_zinv: it expects r/s; here we want r / (tau * sgma).
        # Inline the loop to avoid an extra buffer.
        tsum1 = 0.0
        tsum2 = 0.0
        for i in range(n):
            # tmp = r[i] / (tau_vec[i] * sgma); compute wgt(tmp) per family
            # by inlining the bisquare-style branching. Reuse _wgt_zinv via
            # a scratch buffer trick: build z_buf, then call _wgt_zinv.
            w_arr[i] = r[i] / (tau_vec[i] * sgma)  # store z temporarily
        # Now compute weights from z stored in w_arr; reuse via _wgt_zinv
        # treating w_arr as both input residuals (scaled by 1) and output.
        # _wgt_zinv computes w = wgt(r/s); pass s=1.0 to use w_arr directly.
        _wgt_zinv(w_arr, w_arr, n, 1.0, family, tuning)
        for i in range(n):
            wi = w_arr[i]
            tsum1 += r[i] * r[i] * wi
            tsum2 += wi * tau_vec[i] * tau_vec[i]
        if tsum2 == 0.0:
            out_scale[0] = sgma
            out_converged[0] = 0
            free(w_arr)
            return 0
        new_s = sqrt(tsum1 / (tsum2 * kappa_val))
        if fabs(new_s - sgma) < tol * (tol if tol > sgma else sgma):
            sgma = new_s
            out_scale[0] = sgma
            out_converged[0] = 1
            free(w_arr)
            return 0
        sgma = new_s
    out_scale[0] = sgma
    out_converged[0] = 0
    free(w_arr)
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_lmrob_d_scale(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] residuals,
    cnp.ndarray[double, ndim=1, mode="c"] rweights,
    double init_scale,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning,
    int max_iter,
    double tol,
    cnp.ndarray[double, ndim=1, mode="c"] tau_out,
):
    """Design-adaptive D-scale (Koller & Stahel 2014).

    Computes hat values h from rweights-weighted X, then tau from the
    family-tabulated (tfact, tcorr), then iterates the dt1 fixed-point
    until convergence. ``tau_out`` is filled with the per-observation
    tau values so the Python side can stash them for vcov_w.

    Returns ``(scale, converged, status)``. Status 0 = ok, 1 = no
    tabulated tau/kappa (caller should fall back), 3 = LAPACK error.
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]

    cdef double tfact = 0.0, tcorr = 0.0
    cdef double kappa_val = 0.0
    cdef double* tuning_data = <double*>cnp.PyArray_DATA(tuning)
    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* r_data = <double*>cnp.PyArray_DATA(residuals)
    cdef double* w_data = <double*>cnp.PyArray_DATA(rweights)
    cdef double* tau_data = <double*>cnp.PyArray_DATA(tau_out)
    cdef double scale = init_scale
    cdef int converged = 0
    cdef int status = 0
    cdef Py_ssize_t i

    # Lookup tabulated coefficients before going nogil.
    if _dscale_tau_factors(family, tuning_data, &tfact, &tcorr) != 0:
        return scale, 0, 1
    if _dscale_kappa(family, tuning_data, &kappa_val) != 0:
        return scale, 0, 1

    cdef double* h = <double*>malloc(n * sizeof(double))
    cdef double* XwX = <double*>malloc(p * p * sizeof(double))
    cdef int* ipiv = <int*>malloc(p * sizeof(int))
    if h == NULL or XwX == NULL or ipiv == NULL:
        free(h); free(XwX); free(ipiv)
        raise MemoryError("cy_lmrob_d_scale: out of memory")

    cdef int hat_status
    with nogil:
        hat_status = _compute_hat_diagonal(X_data, w_data, n, p, XwX, ipiv, h)
    if hat_status != 0:
        free(h); free(XwX); free(ipiv)
        return scale, 0, hat_status

    cdef double hi
    with nogil:
        # tau_i = sqrt(1 - tfact*h_i) * (tcorr*h_i + 1)
        for i in range(n):
            hi = h[i]
            tau_data[i] = sqrt(1.0 - tfact * hi) * (tcorr * hi + 1.0)

        # Starting value matches d_scale.py: sqrt(sum(w r^2) / kappa / sum(tau^2 w)).
        # If non-finite or non-positive, fall back to init_scale.
        cdef_num = 0.0
        cdef_den = 0.0
        # (No cdef inside nogil-with-block; declare above instead.)
    # Re-do starting value in normal Python territory to avoid Cython grief
    cdef double num = 0.0
    cdef double den = 0.0
    for i in range(n):
        num += w_data[i] * r_data[i] * r_data[i]
        den += w_data[i] * tau_data[i] * tau_data[i]
    cdef double start
    if den == 0.0 or kappa_val == 0.0:
        start = init_scale
    else:
        start = sqrt(num / (kappa_val * den))
        if start <= 0.0 or start != start:  # NaN check
            start = init_scale

    cdef int iter_status
    with nogil:
        iter_status = _d_scale_iterate(
            r_data, tau_data, n, family, tuning_data,
            kappa_val, start, max_iter, tol,
            &scale, &converged,
        )

    free(h); free(XwX); free(ipiv)
    if iter_status != 0:
        return scale, converged, iter_status
    return scale, converged, 0


# ---------------------------------------------------------------------------
# MM iteration kernel. Port of robustbase/src/lmrob.c::rwls and the existing
# pyrobustlm._mm.mm_iterate. IRWLS with fixed scale and the L1 convergence
# test ``d_beta <= rel_tol * max(rel_tol, ||beta||_1)``.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline int _mm_loop(
    const double* X, const double* y,
    double sigma, int family, const double* tuning,
    int max_it, double rel_tol,
    double* beta,
    _Scratch* scr,
    Py_ssize_t n, Py_ssize_t p,
    int* converged_out, int* n_iter_out,
) nogil:
    """Returns 0 on success, 3 on LAPACK error."""
    cdef int it, irwls_status
    cdef Py_ssize_t j
    cdef double d_beta, norm1_new, diff
    if sigma == 0.0:
        converged_out[0] = 1
        n_iter_out[0] = 0
        return 0
    for it in range(max_it):
        _residuals(X, y, beta, scr.r, n, p)
        for j in range(p):
            scr.beta_prev[j] = beta[j]
        irwls_status = _irwls_step(X, y, scr.r, sigma, family, tuning,
                                   beta, scr, n, p)
        if irwls_status != 0:
            converged_out[0] = 0
            n_iter_out[0] = it + 1
            return irwls_status
        d_beta = 0.0
        norm1_new = 0.0
        for j in range(p):
            diff = beta[j] - scr.beta_prev[j]
            d_beta += diff if diff >= 0 else -diff
            norm1_new += beta[j] if beta[j] >= 0 else -beta[j]
        if d_beta <= rel_tol * (rel_tol if rel_tol > norm1_new else norm1_new):
            converged_out[0] = 1
            n_iter_out[0] = it + 1
            return 0
    converged_out[0] = 0
    n_iter_out[0] = max_it
    return 0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_lmrob_mm(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    cnp.ndarray[double, ndim=1, mode="c"] beta,
    double sigma,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning_psi,
    int max_it,
    double rel_tol,
):
    """Standalone MM IRWLS loop at a fixed scale.

    ``beta`` is updated in place. ``tuning_psi`` is the 95%-efficiency
    tuning (vs ``tuning_chi`` used for the S step).

    Returns ``(n_iter, converged, status)``. Status 0 = ok, 3 = LAPACK.
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef _Scratch scr
    if _alloc_scratch(&scr, n, p) != 0:
        _free_scratch(&scr)
        raise MemoryError("cy_lmrob_mm: out of memory")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef double* beta_data = <double*>cnp.PyArray_DATA(beta)
    cdef double* tuning_data = <double*>cnp.PyArray_DATA(tuning_psi)
    cdef int converged = 0
    cdef int n_iter = 0
    cdef int status = 0

    with nogil:
        status = _mm_loop(
            X_data, y_data, sigma, family, tuning_data,
            max_it, rel_tol, beta_data, &scr, n, p,
            &converged, &n_iter,
        )

    _free_scratch(&scr)
    return n_iter, converged, status


# ---------------------------------------------------------------------------
# Combined fast-S + MM kernel. One nogil block, one workspace.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def cy_lmrob_fit(
    cnp.ndarray[double, ndim=2, mode="c"] X,
    cnp.ndarray[double, ndim=1, mode="c"] y,
    object bitgen_capsule,
    int family,
    cnp.ndarray[double, ndim=1, mode="c"] tuning_chi,
    cnp.ndarray[double, ndim=1, mode="c"] tuning_psi,
    double b0,
    int nResample,
    int mts,
    int k_fast_s,
    int best_r,
    int max_it_s,
    double refine_tol,
    int max_it_mm,
    double rel_tol_mm,
    int max_iter_scale,
    double scale_tol,
    cnp.ndarray[double, ndim=1, mode="c"] beta_out,
    cnp.ndarray[double, ndim=1, mode="c"] residuals_out,
    cnp.ndarray[double, ndim=1, mode="c"] rweights_out,
    cnp.ndarray[double, ndim=1, mode="c"] beta_init_out=None,
    cnp.ndarray[double, ndim=2, mode="c"] cov_out=None,
    cnp.ndarray[double, ndim=1, mode="c"] tuning_chi_for_cov=None,
    double bb_for_cov=0.5,
):
    """Full fast-S + MM fit in a single nogil block.

    Pipeline:
      1. fast-S resampling (chi-tuning)
      2. survivor refinement to convergence
      3. MM IRWLS (psi-tuning) with the resulting scale held fixed
      4. Compute final residuals and IRWLS weights (Mwgt at psi tuning)

    If ``beta_init_out`` is provided, the post-S, pre-MM coefficients are
    written into it (caller uses this to compute the initial residuals
    that vcov_avar1 needs).

    Returns ``(scale, status, n_iter_s, conv_s, n_iter_mm, conv_mm)``.
    """
    if family == FAM_GGW and not _ggw_tables_init:
        _init_ggw_tables()

    cdef bitgen_t* bg = <bitgen_t*>PyCapsule_GetPointer(
        bitgen_capsule, "BitGenerator"
    )

    cdef Py_ssize_t n = X.shape[0]
    cdef Py_ssize_t p = X.shape[1]
    cdef int n_int = <int>n
    cdef int p_int = <int>p
    cdef int one = 1
    cdef int info = 0

    cdef _Scratch scr
    if _alloc_scratch(&scr, n, p) != 0:
        _free_scratch(&scr)
        raise MemoryError("cy_lmrob_fit: out of memory")
    cdef double* best_scales = <double*>malloc(best_r * sizeof(double))
    cdef double* best_betas = <double*>malloc(best_r * p * sizeof(double))
    if best_scales == NULL or best_betas == NULL:
        free(best_scales); free(best_betas)
        _free_scratch(&scr)
        raise MemoryError("cy_lmrob_fit: out of memory for best-r heap")

    cdef double* X_data = <double*>cnp.PyArray_DATA(X)
    cdef double* y_data = <double*>cnp.PyArray_DATA(y)
    cdef double* beta_out_data = <double*>cnp.PyArray_DATA(beta_out)
    cdef double* r_out_data = <double*>cnp.PyArray_DATA(residuals_out)
    cdef double* w_out_data = <double*>cnp.PyArray_DATA(rweights_out)
    cdef double* tuning_chi_data = <double*>cnp.PyArray_DATA(tuning_chi)
    cdef double* tuning_psi_data = <double*>cnp.PyArray_DATA(tuning_psi)

    cdef Py_ssize_t i, j, row, swap, try_i, kept
    cdef double s, scale, max_abs, candidate_scale
    cdef int got_subset, status, k_status
    cdef int conv_s, n_iter_s
    cdef int conv_mm, n_iter_mm
    cdef int found
    cdef uint64_t r_u

    cdef double k0 = tuning_chi_data[0]
    if family == FAM_GGW:
        if k0 < 1: k0 = 1
        elif k0 > 6: k0 = 6
        k0 = _GGW_C0[<int>k0 - 1]
    if k0 <= 0.0:
        k0 = 1.0

    status = 0
    scale = 0.0
    kept = 0
    conv_s = 0
    n_iter_s = 0
    conv_mm = 0
    n_iter_mm = 0

    with nogil:
        # fast-S resampling -------------------------------------------------
        for try_i in range(nResample):
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
                continue

            for j in range(p):
                scr.beta[j] = scr.sub_y[j]
            _residuals(X_data, y_data, scr.beta, scr.r, n, p)

            max_abs = 0.0
            for i in range(n):
                if scr.r[i] >= 0:
                    if scr.r[i] > max_abs:
                        max_abs = scr.r[i]
                else:
                    if -scr.r[i] > max_abs:
                        max_abs = -scr.r[i]
            if max_abs == 0.0:
                for j in range(p):
                    beta_out_data[j] = scr.beta[j]
                scale = 0.0
                status = 2
                break
            s = max_abs / k0
            if s <= 0.0:
                s = 1.0

            candidate_scale = _k_step_refine(
                X_data, y_data, s, family, tuning_chi_data, b0,
                k_fast_s, max_iter_scale, scale_tol,
                &scr, n, p, &k_status,
            )
            if k_status == 2:
                for j in range(p):
                    beta_out_data[j] = scr.beta[j]
                scale = 0.0
                status = 2
                break
            if k_status != 0:
                continue

            if kept < best_r:
                best_scales[kept] = candidate_scale
                for j in range(p):
                    best_betas[kept * p + j] = scr.beta[j]
                kept += 1
            else:
                row = 0  # repurposed as worst_i
                for i in range(1, best_r):
                    if best_scales[i] > best_scales[row]:
                        row = i
                if candidate_scale < best_scales[row]:
                    best_scales[row] = candidate_scale
                    for j in range(p):
                        best_betas[row * p + j] = scr.beta[j]

        # Survivor refinement and pick best ---------------------------------
        if status == 0:
            if kept == 0:
                status = 1
            else:
                scale = 1e300
                for i in range(kept):
                    for j in range(p):
                        scr.beta[j] = best_betas[i * p + j]
                    candidate_scale = _refine_to_convergence(
                        X_data, y_data,
                        best_scales[i], family, tuning_chi_data, b0,
                        max_it_s, refine_tol,
                        max_iter_scale, scale_tol,
                        scr.beta, &scr, n, p,
                        &conv_s, &n_iter_s,
                    )
                    if candidate_scale < scale:
                        scale = candidate_scale
                        for j in range(p):
                            beta_out_data[j] = scr.beta[j]

        # MM step -----------------------------------------------------------
        if status == 0 and scale > 0.0:
            # Save the post-S beta if caller wants init residuals.
            if beta_init_out is not None:
                for j in range(p):
                    (<double*>cnp.PyArray_DATA(beta_init_out))[j] = beta_out_data[j]
            _mm_loop(
                X_data, y_data, scale, family, tuning_psi_data,
                max_it_mm, rel_tol_mm,
                beta_out_data, &scr, n, p,
                &conv_mm, &n_iter_mm,
            )

        # Final residuals and IRWLS weights at psi tuning -------------------
        if status == 0 or status == 2:
            _residuals(X_data, y_data, beta_out_data, r_out_data, n, p)
            if scale > 0.0:
                _wgt_zinv(r_out_data, w_out_data, n, scale, family,
                          tuning_psi_data)
            else:
                for i in range(n):
                    w_out_data[i] = 1.0

    # ------------------------------------------------------------------
    # Optional vcov_avar1 computed in the same call. Needs init_residuals
    # = y - X @ beta_init (reconstructed from beta_init_out). Uses the
    # caller-provided tuning_chi_for_cov; if not provided we fall back
    # to the same tuning the S-step used.
    # ------------------------------------------------------------------
    cdef double* tuning_chi_ptr
    cdef double* init_res = NULL
    cdef double dot2
    cdef int vcov_status = 0
    if (
        cov_out is not None and beta_init_out is not None
        and status == 0 and scale > 0.0
    ):
        init_res = <double*>malloc(n * sizeof(double))
        if init_res == NULL:
            free(best_scales); free(best_betas)
            _free_scratch(&scr)
            raise MemoryError("cy_lmrob_fit: out of memory for init_residuals")
        if tuning_chi_for_cov is not None:
            tuning_chi_ptr = <double*>cnp.PyArray_DATA(tuning_chi_for_cov)
        else:
            tuning_chi_ptr = tuning_chi_data
        with nogil:
            for i in range(n):
                dot2 = 0.0
                for j in range(p):
                    dot2 += X_data[i * p + j] * (<double*>cnp.PyArray_DATA(beta_init_out))[j]
                init_res[i] = y_data[i] - dot2
            vcov_status = _compute_vcov_avar1_body(
                X_data, r_out_data, init_res, scale, family,
                tuning_psi_data, tuning_chi_ptr, bb_for_cov,
                n, p, <double*>cnp.PyArray_DATA(cov_out),
            )
        free(init_res)
    else:
        vcov_status = -1  # not computed

    free(best_scales); free(best_betas)
    _free_scratch(&scr)
    return scale, status, n_iter_s, conv_s, n_iter_mm, conv_mm, vcov_status
