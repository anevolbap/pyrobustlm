# SPDX-License-Identifier: GPL-3.0-or-later
#
# Shared psi/chi/M-scale kernels, textually included by both
# ``_fast_s.pyx`` and ``_lmrob.pyx``.
#
# These functions were copy-pasted into both modules and drifted.
# ``_lmrob.pyx`` picked up reciprocal hoisting (``r[i] * inv_sk`` for
# ``r[i] / s``, ``ax >= 1.0`` for ``ax >= k``, precomputed lqq
# constants) that ``_fast_s.pyx`` never received, so the two engines
# rounded differently on identical input. The same habit produced an
# outright bug next door: a second hardcoded ``kappa`` table whose ggw
# case-4 entry was a copy of case 1, worth 3.8% on every ggw D-scale.
# This file is the ``_lmrob.pyx`` version, so ``_fast_s.pyx`` gains the
# optimisation it was missing.
#
# ``include`` rather than ``cimport`` on purpose: these are ``cdef
# inline`` helpers plus module-static lookup tables. Textual inclusion
# keeps each extension self-contained (its own copy of the read-only
# tables, no cross-module linkage) while leaving one place to edit.
#
# Requires the including module to have already cimported:
#     cython, libc.math (fabs, sqrt, exp, pow as cpow),
#     libc.stdint.uint64_t, numpy.random.bitgen_t
#
cdef double _MAX_EX2_SQR_HALF = 37.7 * 37.7 / 2.0


# ---------------------------------------------------------------------------
# Family enum. Mirrors pylmrob._fast_s._FAMILY_IDS for now; will be
# the single source of truth once the monolithic kernel is the default.
# ---------------------------------------------------------------------------
cdef enum:
    FAM_BISQUARE = 0
    FAM_HAMPEL = 1
    FAM_OPTIMAL = 2
    FAM_LQQ = 3
    FAM_GGW = 4
    FAM_WELSH = 5


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
    cdef double inv_s, inv_k, inv_sk, inv_denom
    cdef double c0_lqq, c1_lqq, c2_lqq, c3_lqq, s_l_over_b_l, s5_sq_over_3s6
    cdef double inv_3p25, inv_6p5
    cdef int j

    if family == FAM_BISQUARE:
        k = tuning[0]
        inv_sk = 1.0 / (s * k)
        for i in range(n):
            t = r[i] * inv_sk
            ax = t if t >= 0 else -t
            if ax >= 1.0:
                total += 1.0
            else:
                t = 1.0 - ax * ax
                total += 1.0 - t * t * t
    elif family == FAM_HAMPEL:
        a_t = tuning[0]; b_t = tuning[1]; r_t = tuning[2]
        nc = a_t * (b_t + r_t - a_t) * 0.5
        inv_nc = 1.0 / nc
        inv_s = 1.0 / s
        bma_inv_half = 0.5 / (r_t - b_t)
        for i in range(n):
            xi = r[i] * inv_s
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
        inv_sk = 1.0 / (s * k)
        R1h = -1.944 * 0.5
        R2h = 1.728 * 0.25
        R3h = -0.312 / 6.0
        R4h = 0.016 / 8.0
        inv_3p25 = 1.0 / 3.25
        inv_6p5 = 1.0 / 6.5
        for i in range(n):
            ac = r[i] * inv_sk
            ax = ac if ac >= 0 else -ac
            if ax > 3.0:
                total += 1.0
            elif ax > 2.0:
                a2 = ax * ax
                total += (a2 * (R1h + a2 * (R2h + a2 * (R3h + a2 * R4h))) + 1.792) * inv_3p25
            else:
                total += (ac * ac) * inv_6p5
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
        inv_s = 1.0 / s
        inv_denom = 1.0 / denom
        c0_lqq = (3.0 * s_l - 3.0) * inv_denom
        c1_lqq = (6.0 * s_l - 6.0) * inv_denom
        c2_lqq = (6.0 * s5) * inv_denom
        c3_lqq = k01_2 * 0.5 - s_l * b_l * b_l / 6.0
        s_l_over_b_l = s_l / b_l
        s5_sq_over_3s6 = s5 * s5 / (3.0 * s6) if s6 != 0.0 else 0.0
        for i in range(n):
            xi = r[i] * inv_s
            ax = xi if xi >= 0 else -xi
            if ax <= c:
                total += c0_lqq * xi * xi
            elif ax <= k01:
                s0 = ax - c
                total += c1_lqq * (xi * xi * 0.5 - s_l_over_b_l * s0 * s0 * s0 / 6.0)
            elif ax < end3:
                dx = ax - k01
                total += c2_lqq * (
                    c3_lqq - dx * 0.5 * (s6 + dx * (s5 + dx * s5_sq_over_3s6))
                )
            else:
                total += 1.0
    elif family == FAM_WELSH:
        # rho(x; c) = 1 - exp(-(x/(c*s))^2 / 2)
        k = tuning[0]
        inv_sk = 1.0 / (s * k)
        for i in range(n):
            t = r[i] * inv_sk
            total += 1.0 - exp(-0.5 * t * t)
    else:  # FAM_GGW
        j = <int>(tuning[0]) - 1
        if j < 0:
            j = 0
        elif j > 5:
            j = 5
        total += _ggw_chi_sum(r, n, 1.0 / s, j)
    return total


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double _ggw_chi_sum(
    const double* r, Py_ssize_t n, double inv_s, int j,
) nogil:
    """Sum of GGW chi(r[i]/s) over i, with per-case constants hoisted out
    of the inner loop. Same numerical result as repeated _ggw_rho_one;
    measurable speedup on the M-scale inner loop (~30% on n=2000)."""
    cdef double c0 = _GGW_C0[j]
    cdef double end = _GGW_END[j]
    cdef double p0  = _GGW_POLY[j][0]
    cdef double p1  = _GGW_POLY[j][1]
    cdef double p2  = _GGW_POLY[j][2]
    cdef double p3  = _GGW_POLY[j][3]
    cdef double p4  = _GGW_POLY[j][4]
    cdef double p5  = _GGW_POLY[j][5]
    cdef double p6  = _GGW_POLY[j][6]
    cdef double p7  = _GGW_POLY[j][7]
    cdef double p8  = _GGW_POLY[j][8]
    cdef double p9  = _GGW_POLY[j][9]
    cdef double p10 = _GGW_POLY[j][10]
    cdef double p11 = _GGW_POLY[j][11]
    cdef double p12 = _GGW_POLY[j][12]
    cdef double p13 = _GGW_POLY[j][13]
    cdef double p14 = _GGW_POLY[j][14]
    cdef double p15 = _GGW_POLY[j][15]
    cdef double p16 = _GGW_POLY[j][16]
    cdef double p17 = _GGW_POLY[j][17]
    cdef double p18 = _GGW_POLY[j][18]
    cdef double p19 = _GGW_POLY[j][19]
    cdef double total = 0.0
    cdef double x, ax, res
    cdef double three_c0 = 3.0 * c0
    cdef Py_ssize_t i
    for i in range(n):
        x = r[i] * inv_s
        ax = x if x >= 0 else -x
        if ax <= c0:
            total += p0 * ax * ax
        elif ax <= three_c0:
            res = p9
            res = res * ax + p8
            res = res * ax + p7
            res = res * ax + p6
            res = res * ax + p5
            res = res * ax + p4
            res = res * ax + p3
            res = res * ax + p2
            res = res * ax + p1
            total += res
        elif ax <= end:
            res = p19
            res = res * ax + p18
            res = res * ax + p17
            res = res * ax + p16
            res = res * ax + p15
            res = res * ax + p14
            res = res * ax + p13
            res = res * ax + p12
            res = res * ax + p11
            res = res * ax + p10
            total += res
        else:
            total += 1.0
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
    cdef double k, inv_sk, inv_s, a_t_inv_rmb, inv_2a
    cdef double a_t, b_t, r_t
    cdef double R1, R2, R3, R4
    cdef double b_l, c, s_l, k01, s5, s6, end3
    cdef double s_l_over_2bl, s5_sq_over_s6, s6_over_s5, s6_half
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
        inv_s = 1.0 / s
        a_t_inv_rmb = a_t / (r_t - b_t)
        for i in range(n):
            xi = r[i] * inv_s
            u = xi if xi >= 0 else -xi
            if u <= a_t:
                out[i] = 1.0
            elif u <= b_t:
                out[i] = a_t / u if u > 0 else 1.0
            elif u <= r_t:
                out[i] = a_t_inv_rmb * (r_t - u) / u
            else:
                out[i] = 0.0
    elif family == FAM_OPTIMAL:
        k = tuning[0]
        inv_sk = 1.0 / (s * k)
        R1 = -1.944; R2 = 1.728; R3 = -0.312; R4 = 0.016
        for i in range(n):
            ac = r[i] * inv_sk
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
        inv_s = 1.0 / s
        s_l_over_2bl = s_l / (2.0 * b_l)
        s5_sq_over_s6 = s5 * s5 / s6 if s6 != 0.0 else 0.0
        s6_over_s5 = s6 / s5 if s5 != 0.0 else 0.0
        s6_half = s6 * 0.5
        for i in range(n):
            xi = r[i] * inv_s
            ax = xi if xi >= 0 else -xi
            if ax <= c:
                out[i] = 1.0
            elif ax <= k01:
                s0 = ax - c
                if ax > 0:
                    out[i] = 1.0 - s_l_over_2bl * s0 * s0 / ax
                else:
                    out[i] = 1.0
            elif ax < end3:
                dx = ax - k01
                if ax > 0:
                    out[i] = -(
                        s6_half + s5_sq_over_s6 * dx * (dx * 0.5 + s6_over_s5)
                    ) / ax
                else:
                    out[i] = 0.0
            else:
                out[i] = 0.0
    elif family == FAM_WELSH:
        # wgt(x; c) = exp(-(x/(c*s))^2/2)
        k = tuning[0]
        inv_sk = 1.0 / (s * k)
        for i in range(n):
            a = r[i] * inv_sk
            out[i] = exp(-0.5 * a * a)
    else:  # FAM_GGW
        # The b parameter is 1.0 (cases 1..3) or 1.5 (cases 4..6) in the
        # standard table. Hardcoding those two saves one libm pow() call
        # per element. User-supplied b falls back to cpow.
        j = <int>(tuning[0])
        if j < 1:
            j = 1
        elif j > 6:
            j = 6
        a_t = _GGW_ABC_A[j]
        b_t = _GGW_ABC_B[j]
        r_t = _GGW_ABC_C[j]  # c
        inv_2a = 1.0 / (2.0 * a_t)
        if b_t == 1.0:
            for i in range(n):
                xi = r[i] / s
                ax = xi if xi >= 0 else -xi
                if ax <= r_t:
                    out[i] = 1.0
                else:
                    dx = ax - r_t
                    ac = dx * inv_2a
                    if ac > _MAX_EX2_SQR_HALF:
                        ac = _MAX_EX2_SQR_HALF
                    out[i] = exp(-ac)
        elif b_t == 1.5:
            for i in range(n):
                xi = r[i] / s
                ax = xi if xi >= 0 else -xi
                if ax <= r_t:
                    out[i] = 1.0
                else:
                    dx = ax - r_t
                    ac = dx * sqrt(dx) * inv_2a
                    if ac > _MAX_EX2_SQR_HALF:
                        ac = _MAX_EX2_SQR_HALF
                    out[i] = exp(-ac)
        else:
            for i in range(n):
                xi = r[i] / s
                ax = xi if xi >= 0 else -xi
                if ax <= r_t:
                    out[i] = 1.0
                else:
                    dx = ax - r_t
                    ac = cpow(dx, b_t) * inv_2a
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
