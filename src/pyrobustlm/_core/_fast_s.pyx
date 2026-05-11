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


# Family enum. Keep these in sync with the dispatch table in
# pyrobustlm/_fast_s.py (`_FAMILY_IDS`).
cdef enum:
    FAM_BISQUARE = 0
    FAM_HAMPEL = 1
    FAM_OPTIMAL = 2
    FAM_LQQ = 3
    FAM_GGW = 4


# Largest x such that exp(-x^2/2) does not underflow (matches lmrob.c:945).
cdef double _MAX_EX2_SQR_HALF = 37.7 * 37.7 / 2.0


# ---------------------------------------------------------------------------
# GGW polynomial chi tables. Mirror the tables in _psi.pyx exactly; we
# duplicate them here so the resampling loop stays nogil without cross-
# module cdef linking. Initialised lazily on first ggw call.
# ---------------------------------------------------------------------------
cdef double _GGW_C0[6]
cdef double _GGW_END[6]
cdef double _GGW_POLY[6][20]
cdef double _GGW_ABC_A[7]
cdef double _GGW_ABC_B[7]
cdef double _GGW_ABC_C[7]
cdef int _ggw_tables_init = 0


cdef _init_ggw_tables():
    """Populate the ggw polynomial coefficient tables and the (a, b, c)
    lookup for case_idx in [1, 6]. Idempotent."""
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
    # (a, b, c) lookup; mirrors _GGW_ABC in pyrobustlm.scale.
    _GGW_ABC_A[0] = 0.0; _GGW_ABC_B[0] = 0.0; _GGW_ABC_C[0] = 0.0  # unused
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
    """Polynomial chi for ggw case j (0-based). Mirrors lmrob.c::rho_ggw."""
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
    cdef int j  # ggw polynomial case index

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
        # tuning[0] holds the case index (1..6). Tables initialised at the
        # top-level entry point before we drop the GIL.
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
    cdef int j  # ggw case index

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
        # Mirrors optimal_wgt in _psi.pyx (line 616).
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
        # Mirrors lqq_wgt in _psi.pyx (line 712) exactly.
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
        # tuning[0] = case index (1..6); we look up (a, b, c) and apply
        # the analytical wgt = exp(-((|x|-c)^b)/(2a)).
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
                # cpow(diff, b_t) / (2 * a_t)
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
# Bounded random integer in [0, bound), using numpy's bitgen. Lemire's debiased
# method to avoid modulo bias for moderately-sized bounds.
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline uint64_t _bounded_uint64(bitgen_t* bg, uint64_t bound) nogil:
    if bound <= 1:
        return 0
    cdef uint64_t threshold = (-bound) % bound  # (2^64 - bound) % bound
    cdef uint64_t r
    while True:
        r = bg.next_uint64(bg.state)
        if r >= threshold:
            return r % bound


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
