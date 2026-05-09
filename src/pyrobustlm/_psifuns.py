# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-NumPy reference implementations of the six robustbase psi families.

Direct port of ``robustbase/src/lmrob.c`` (rho_/psi_/psip_/wgt_ huber, biwgt,
hmpl, opt, ggw, lqq). The Cython kernels in ``pyrobustlm._core._psi`` will
mirror this module element-for-element; the Python version keeps the math
testable without the Cython build.

For each family, six functions of ``x`` and tuning constants ``k``:

- ``rho``       (a.k.a. R's ``Mchi`` for these chi-shaped families)
- ``psi``       (a.k.a. R's ``Mpsi``; psi = rho')
- ``psi_prime`` (psi' = rho'')
- ``wgt``       (psi(x)/x; a.k.a. R's ``Mwgt``)

All accept scalar or NumPy array ``x``; ``k`` is a 1-D float array (length 1
for huber/biwgt/opt, length 3 for hmpl, length 4 for ggw, length 3 for lqq).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

# ---------------------------------------------------------------------------
# Numerical safety constants matching the C source.
# ---------------------------------------------------------------------------
# Largest x such that exp(-x^2/2) does not underflow (lmrob.c:945).
_MAX_EX2 = 37.7


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _to_k(k: ArrayLike, n_required: int) -> NDArray[np.float64]:
    arr = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()
    if arr.size < n_required:
        raise ValueError(f"psi family requires {n_required} tuning constants, got {arr.size}")
    return arr


def _xa(x: ArrayLike) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return (x, |x|) as float64 arrays."""
    x = np.asarray(x, dtype=np.float64)
    return x, np.abs(x)


def _safe_div(a: NDArray, b: NDArray) -> NDArray:
    """Element-wise a/b, with 0/0 -> 0."""
    out = np.zeros_like(a, dtype=np.float64)
    mask = b != 0
    out[mask] = a[mask] / b[mask]
    return out


# ===========================================================================
# Huber
# ===========================================================================
def rho_huber(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, ax = _xa(x)
    inner = ax <= kk
    return np.where(inner, 0.5 * x * x, kk * (ax - 0.5 * kk))


def psi_huber(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, ax = _xa(x)
    return np.where(ax < kk, x, np.sign(x) * kk)


def psi_prime_huber(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    _, ax = _xa(x)
    return np.where(ax >= kk, 0.0, 1.0)


def wgt_huber(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    _, ax = _xa(x)
    return np.where(ax >= kk, kk / np.maximum(ax, 1e-300), 1.0)


# ===========================================================================
# Bisquare / Tukey
# ===========================================================================
def rho_biwgt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, ax = _xa(x)
    out = np.ones_like(x)
    inner = ax <= kk
    t = (x[inner] / kk) ** 2
    out[inner] = t * (3.0 + t * (-3.0 + t))
    return out


def psi_biwgt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, ax = _xa(x)
    out = np.zeros_like(x)
    inner = ax <= kk
    a = x[inner] / kk
    u = 1.0 - a * a
    out[inner] = x[inner] * u * u
    return out


def psi_prime_biwgt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, ax = _xa(x)
    out = np.zeros_like(x)
    inner = ax <= kk
    a = x[inner] / kk
    a2 = a * a
    out[inner] = (1.0 - a2) * (1.0 - 5.0 * a2)
    return out


def wgt_biwgt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, ax = _xa(x)
    out = np.zeros_like(x)
    inner = ax <= kk
    a = x[inner] / kk
    out[inner] = (1.0 - a) ** 2 * (1.0 + a) ** 2
    return out


# ===========================================================================
# Hampel (3 tuning constants a, b, r)
# ===========================================================================
def rho_hmpl(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    a, b, r = kk[0], kk[1], kk[2]
    x, u = _xa(x)
    nc = a * (b + r - a) / 2.0
    out = np.ones_like(x)
    m1 = u <= a
    m2 = (~m1) & (u <= b)
    m3 = (~m1) & (~m2) & (u <= r)
    out[m1] = (x[m1] ** 2) / 2.0 / nc
    out[m2] = (u[m2] - a / 2.0) * a / nc
    out[m3] = (b - a / 2.0 + (u[m3] - b) * (1.0 - (u[m3] - b) / (r - b) / 2.0)) * a / nc
    return out


def psi_hmpl(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    a, b, r = kk[0], kk[1], kk[2]
    x, u = _xa(x)
    sx = np.sign(x)
    out = np.zeros_like(x)
    m1 = u <= a
    m2 = (~m1) & (u <= b)
    m3 = (~m1) & (~m2) & (u <= r)
    out[m1] = x[m1]
    out[m2] = sx[m2] * a
    out[m3] = sx[m3] * a * (r - u[m3]) / (r - b)
    return out


def psi_prime_hmpl(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    a, b, r = kk[0], kk[1], kk[2]
    _, u = _xa(x)
    out = np.zeros_like(u)
    m1 = u <= a
    m3 = (u > b) & (u <= r)
    out[m1] = 1.0
    out[m3] = a / (b - r)
    return out


def wgt_hmpl(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    a, b, r = kk[0], kk[1], kk[2]
    _, u = _xa(x)
    out = np.zeros_like(u)
    m1 = u <= a
    m2 = (~m1) & (u <= b)
    m3 = (~m1) & (~m2) & (u <= r)
    out[m1] = 1.0
    out[m2] = a / np.maximum(u[m2], 1e-300)
    out[m3] = a * (r - u[m3]) / (r - b) / np.maximum(u[m3], 1e-300)
    return out


# ===========================================================================
# Optimal (Yohai / robust package). Single tuning constant.
# ===========================================================================
def _opt_polyR(t: NDArray[np.float64]) -> NDArray[np.float64]:
    # Coefficients in psi_opt branch: R1 a + R2 a^3 + R3 a^5 + R4 a^7 with a=ax
    R1, R2, R3, R4 = -1.944, 1.728, -0.312, 0.016
    a2 = t * t
    return ((R4 * a2 + R3) * a2 + R2) * a2 + R1


def rho_opt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    R1, R2, R3, R4 = -1.944 / 2.0, 1.728 / 4.0, -0.312 / 6.0, 0.016 / 8.0
    x, ax = _xa(x)
    ac = x / kk
    a = np.abs(ac)
    out = np.empty_like(x)
    m_far = a > 3.0
    m_mid = (~m_far) & (a > 2.0)
    m_in = (~m_far) & (~m_mid)
    out[m_far] = 1.0
    a2 = a[m_mid] ** 2
    out[m_mid] = (a2 * (R1 + a2 * (R2 + a2 * (R3 + a2 * R4))) + 1.792) / 3.25
    out[m_in] = (ac[m_in] ** 2) / 6.5
    return out


def psi_opt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    x, _ = _xa(x)
    ac = x / kk
    a = np.abs(ac)
    out = np.empty_like(x)
    m_far = a > 3.0
    m_mid = (~m_far) & (a > 2.0)
    m_in = (~m_far) & (~m_mid)
    out[m_far] = 0.0
    poly = kk * _opt_polyR(a[m_mid]) * ac[m_mid]
    out_mid = np.where(ac[m_mid] > 0, np.maximum(0.0, poly), -np.abs(poly))
    out[m_mid] = out_mid
    out[m_in] = x[m_in]
    return out


def psi_prime_opt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    R1, R2, R3, R4 = -1.944, 1.728, -0.312, 0.016
    x, _ = _xa(x)
    a = np.abs(x / kk)
    out = np.empty_like(x)
    m_far = a > 3.0
    m_mid = (~m_far) & (a > 2.0)
    m_in = (~m_far) & (~m_mid)
    out[m_far] = 0.0
    a2 = a[m_mid] ** 2
    out[m_mid] = R1 + a2 * (3 * R2 + a2 * (5 * R3 + a2 * 7 * R4))
    out[m_in] = 1.0
    return out


def wgt_opt(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 1)[0]
    R1, R2, R3, R4 = -1.944, 1.728, -0.312, 0.016
    x, _ = _xa(x)
    a = np.abs(x / kk)
    out = np.empty_like(x)
    m_far = a > 3.0
    m_mid = (~m_far) & (a > 2.0)
    m_in = (~m_far) & (~m_mid)
    out[m_far] = 0.0
    a2 = a[m_mid] ** 2
    out[m_mid] = np.maximum(0.0, R1 + a2 * (R2 + a2 * (R3 + a2 * R4)))
    out[m_in] = 1.0
    return out


# ===========================================================================
# LQQ (Linear-Quadratic-Quadratic). 3 tuning constants (b, c, s).
# ===========================================================================
def psi_lqq(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    b, c, s = kk[0], kk[1], kk[2]
    x, ax = _xa(x)
    out = np.empty_like(x)
    sx = np.sign(x)

    k01 = b + c
    s5 = s - 1.0
    s6 = -2.0 * k01 + b * s

    m_lin = ax <= c
    m_quad = (~m_lin) & (ax <= k01)
    end3 = k01 - s6 / s5 if s5 != 0 else k01  # protect div-by-zero
    m_quad2 = (~m_lin) & (~m_quad) & (ax < end3)
    m_zero = ~(m_lin | m_quad | m_quad2)

    out[m_lin] = x[m_lin]
    out[m_quad] = sx[m_quad] * (ax[m_quad] - s * (ax[m_quad] - c) ** 2 / b / 2.0)
    out[m_quad2] = sx[m_quad2] * (
        -s6 / 2.0 - (s5**2) / s6 * ((ax[m_quad2] - k01) ** 2 / 2.0 + s6 / s5 * (ax[m_quad2] - k01))
    )
    out[m_zero] = 0.0
    return out


def rho_lqq(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    b, c, s = kk[0], kk[1], kk[2]
    x, ax = _xa(x)
    out = np.ones_like(x)

    k01 = b + c
    s5 = s - 1.0
    s6 = -2.0 * k01 + b * s
    denom = s * c * (3.0 * c + 2.0 * b) + k01**2

    m_lin = ax <= c
    m_quad = (~m_lin) & (ax <= k01)
    end3 = k01 - s6 / s5 if s5 != 0 else k01
    m_quad2 = (~m_lin) & (~m_quad) & (ax < end3)

    out[m_lin] = (3.0 * s - 3.0) / denom * x[m_lin] ** 2
    s0 = ax[m_quad] - c
    out[m_quad] = (6.0 * s - 6.0) / denom * (x[m_quad] ** 2 / 2.0 - s / b * s0**3 / 6.0)
    s7 = ax[m_quad2] - k01
    k01_2 = k01**2
    out[m_quad2] = (
        (6.0 * s5)
        / denom
        * (k01_2 / 2.0 - s * b * b / 6.0 - s7 / 2.0 * (s6 + s7 * (s5 + s7 * s5 * s5 / 3.0 / s6)))
    )
    return out


def psi_prime_lqq(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    """psi'_lqq matching robustbase's ``psip_lqq`` (lmrob.c:1339).

    Uses the lmrob.c convention internally: ``s5 = 1 - s``, ``a = (b*s - 2*k01)/s5``.
    """
    kk = _to_k(k, 3)
    b, c, s = kk[0], kk[1], kk[2]
    x, ax = _xa(x)
    out = np.zeros_like(x)
    k01 = b + c
    s5p = 1.0 - s  # convention used in psip_lqq, opposite sign to psi_lqq
    if s5p == 0.0:
        # degenerate; psi' undefined past k01
        m_lin = ax <= c
        m_quad = (~m_lin) & (ax <= k01)
        out[m_lin] = 1.0
        out[m_quad] = 1.0 - s * (ax[m_quad] - c) / b
        return out
    a = (b * s - 2.0 * k01) / s5p

    m_lin = ax <= c
    m_quad = (~m_lin) & (ax <= k01)
    m_quad2 = (~m_lin) & (~m_quad) & (ax < k01 + a)

    out[m_lin] = 1.0
    out[m_quad] = 1.0 - s * (ax[m_quad] - c) / b
    out[m_quad2] = -s5p * ((ax[m_quad2] - k01) / a - 1.0)
    return out


def wgt_lqq(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    kk = _to_k(k, 3)
    b, c, s = kk[0], kk[1], kk[2]
    _, ax = _xa(x)
    out = np.zeros_like(ax)
    k01 = b + c
    s5 = s - 1.0
    s6 = -2.0 * k01 + b * s
    end3 = k01 - s6 / s5 if s5 != 0 else k01

    m_lin = ax <= c
    m_quad = (~m_lin) & (ax <= k01)
    m_quad2 = (~m_lin) & (~m_quad) & (ax < end3)

    out[m_lin] = 1.0
    s0 = ax[m_quad] - c
    out[m_quad] = 1.0 - s * s0**2 / (2.0 * ax[m_quad] * b)
    s7 = ax[m_quad2] - k01
    out[m_quad2] = -(s6 / 2.0 + s5 * s5 / s6 * s7 * (s7 / 2.0 + s6 / s5)) / np.maximum(
        ax[m_quad2], 1e-300
    )
    return out


# ===========================================================================
# GGW (Generalized Gauss-Weight) — Koller & Stahel (2011)
# ===========================================================================
# When k[0] > 0, robustbase uses precomputed polynomial approximations for one
# of six (b, eff/bp) cases. We mirror that table verbatim. When k[0] == 0 the
# C code falls back to numerical integration of psi(t)*psi'(t); we use SciPy
# QUADPACK in the same role.

# Hardcoded cases (case index = int(k[0]) - 1)
_GGW_C = np.array(
    [
        # case 0: b=1, 95% efficiency
        [
            0.094164571656733,
            -0.168937372816728,
            0.00427612218326869,
            0.336876420549802,
            -0.166472338873754,
            0.0436904383670537,
            -0.00732077121233756,
            0.000792550423837942,
            -5.08385693557726e-05,
            1.46908724988936e-06,
            -0.837547853001024,
            0.876392734183528,
            -0.184600387321924,
            0.0219985685280105,
            -0.00156403138825785,
            6.16243137719362e-05,
            -7.478979895101e-07,
            -3.99563057938975e-08,
            1.78125589532002e-09,
            -2.22317669250326e-11,
        ],
        # case 1: b=1, 85% efficiency
        [
            0.174505224068561,
            -0.168853188892986,
            0.00579250806463694,
            0.624193375180937,
            -0.419882092234336,
            0.150011303015251,
            -0.0342185249354937,
            0.00504325944243195,
            -0.0004404209084091,
            1.73268448820236e-05,
            -0.842160072154898,
            1.19912623576069,
            -0.345595777445623,
            0.0566407000764478,
            -0.00560501531439071,
            0.000319084704541442,
            -7.4279004383686e-06,
            -2.02063746721802e-07,
            1.65716101809839e-08,
            -2.97536178313245e-10,
        ],
        # case 2: b=1, bp 0.5
        [
            1.41117142330711,
            -0.168853741371095,
            0.0164713906344165,
            5.04767833986545,
            -9.65574752971554,
            9.80999125035463,
            -6.36344090274658,
            2.667031271863,
            -0.662324374141645,
            0.0740982983873332,
            -0.84794906554363,
            3.4315790970352,
            -2.82958670601597,
            1.33442885893807,
            -0.384812004961396,
            0.0661359078129487,
            -0.00557221619221031,
            -5.42574872792348e-05,
            4.92564168111658e-05,
            -2.80432020951381e-06,
        ],
        # case 3: b=1.5, 95% efficiency
        [
            0.104604570079252,
            0.0626649856211545,
            -0.220058184826331,
            0.403388189975896,
            -0.213020713708997,
            0.102623342948069,
            -0.0392618698058543,
            0.00937878752829234,
            -0.00122303709506374,
            6.70669880352453e-05,
            0.632651530179424,
            -1.14744323908043,
            0.981941598165897,
            -0.341211275272191,
            0.0671272892644464,
            -0.00826237596187364,
            0.0006529134641922,
            -3.23468516804340e-05,
            9.17904701930209e-07,
            -1.14119059405971e-08,
        ],
        # case 4: b=1.5, 85% efficiency
        [
            0.205026436642222,
            0.0627464477520301,
            -0.308483319391091,
            0.791480474953874,
            -0.585521414631968,
            0.394979618040607,
            -0.211512515412973,
            0.0707208739858416,
            -0.0129092527174621,
            0.000990938134086886,
            0.629919019245325,
            -1.60049136444912,
            1.91903069049618,
            -0.933285960363159,
            0.256861783311473,
            -0.0442133943831343,
            0.00488402902512139,
            -0.000338084604725483,
            1.33974565571893e-05,
            -2.32450916247553e-07,
        ],
        # case 5: b=1.5, bp 0.5
        [
            1.35010856132000,
            0.0627465630782482,
            -0.791613168488525,
            5.21196700244212,
            -9.89433796586115,
            17.1277266427962,
            -23.5364159883776,
            20.1943966645350,
            -9.4593988142692,
            1.86332355622445,
            0.62986381140768,
            -4.10676399816156,
            12.6361433997327,
            -15.7697199271455,
            11.1373468568838,
            -4.91933095295458,
            1.39443093325178,
            -0.247689078940725,
            0.0251861553415515,
            -0.00112130382664914,
        ],
    ],
    dtype=np.float64,
)
_GGW_END = np.array(
    [
        18.5527638190955,
        13.7587939698492,
        4.89447236180905,
        11.4974874371859,
        8.15075376884422,
        3.17587939698492,
    ],
    dtype=np.float64,
)
# SET_ABC_GGW table from robustbase/src/lmrob.c (lines 1279-1293).
# (a, b, c) for cases 1..6. Format used by psi/wgt/psip:
#   psi(x; a, b, c) = x * exp(-(|x|-c)+^b / (2 a))
# Note: this `b` is the *exponent*, not the K&S "b" — it matches the C k[2].
_GGW_ABC: dict[int, tuple[float, float, float]] = {
    1: (0.648, 1.0, 1.694),
    2: (0.4760508, 1.0, 1.2442567),
    3: (0.1674046, 1.0, 0.4375470),
    4: (1.387, 1.5, 1.063),
    5: (0.8372485, 1.5, 0.7593544),
    6: (0.2036741, 1.5, 0.2959132),
}
# c-only table (still used by rho_ggw which has its own c lookup).
_GGW_C0 = np.array(
    [1.694, 1.2442567, 0.4375470, 1.063, 0.7593544, 0.2959132],
    dtype=np.float64,
)


def _polyval_horner(coeffs: NDArray[np.float64], x: NDArray[np.float64]) -> NDArray[np.float64]:
    """Horner evaluation of a polynomial whose constant term is coeffs[0]."""
    res = np.zeros_like(x)
    # Match the C code's right-to-left nesting: coeffs[0] + x*(coeffs[1] + x*(...))
    for c in coeffs[::-1]:
        res = res * x + c
    return res


def _ggw_unpack(k: NDArray[np.float64]) -> tuple[float, float, float, int]:
    """Return (a, b, c, case) from the internal-form ggw tuning vector.

    Internal form mirrors the C kernel's view of ``k``:

    - Length-1: ``(case_idx,)`` selecting one of the 6 hardcoded cases (1..6).
    - Length-4 starting with 0: ``(0, a, b, c)`` for user-specified parameters.
    - Length-5 starting with 0: ``(0, a, b, c, rho_inf)`` (rho_inf only used by
      rho via numerical integration).

    Note: this is the *internal* form. For the user-facing form
    ``c(min_slope, b, eff, bdp)``, see :func:`pyrobustlm.psi.psi`.
    """
    k = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()
    if k.size == 0:
        raise ValueError("ggw: empty tuning vector")
    case = int(k[0])
    if case >= 1:
        if case > 6:
            raise ValueError(f"ggw: case index {case} out of range 1..6")
        a, b, c = _GGW_ABC[case]
        return a, b, c, case
    if case == 0:
        if k.size < 4:
            raise ValueError("ggw user-specified case needs at least (0, a, b, c)")
        return float(k[1]), float(k[2]), float(k[3]), 0
    raise ValueError(f"ggw: invalid case selector {case}")


def rho_ggw(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    a, b, c, case = _ggw_unpack(np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel())
    if case >= 1:
        # Hardcoded polynomial approximation
        j = case - 1
        ax = np.abs(np.asarray(x, dtype=np.float64))
        c_local = _GGW_C0[j]
        end = _GGW_END[j]
        out = np.ones_like(ax)
        m1 = ax <= c_local
        m2 = (~m1) & (ax <= 3 * c_local)
        m3 = (~m1) & (~m2) & (ax <= end)
        out[m1] = _GGW_C[j, 0] * ax[m1] * ax[m1]
        out[m2] = _polyval_horner(_GGW_C[j, 1:10], ax[m2])
        out[m3] = _polyval_horner(_GGW_C[j, 10:20], ax[m3])
        return out
    return _ggw_general_rho(x, np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel())


def _ggw_psi_kernel(x: NDArray[np.float64], a: float, b: float, c: float) -> NDArray[np.float64]:
    """Closed form psi(x; a, b, c) = x * exp(-(|x|-c)+^b / (2 a))."""
    ax = np.abs(x)
    inner = ax <= c
    diff = np.where(inner, 0.0, ax - c)
    arg = np.where(inner, 0.0, np.power(diff, b) / (2.0 * a))
    arg = np.minimum(arg, _MAX_EX2 * _MAX_EX2 / 2.0)
    return x * np.exp(-arg)


def psi_ggw(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    a, b, c, _ = _ggw_unpack(np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel())
    return _ggw_psi_kernel(np.asarray(x, dtype=np.float64), a, b, c)


def psi_prime_ggw(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    """psi'_ggw mirroring lmrob.c:1304-1317.

    For |x|<c: returns 1. For |x|>=c:
      let arg = (|x|-c)^b / (2a)    (in C: ea = -arg with a*=2 first)
      psi' = exp(-arg) * (1 - b/(2a) * |x| * (|x|-c)^(b-1))
    """
    a, b, c, _ = _ggw_unpack(np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel())
    x = np.asarray(x, dtype=np.float64)
    ax = np.abs(x)
    out = np.ones_like(x)
    outer = ax >= c
    if not outer.any():
        return out
    diff = np.where(outer, ax - c, 0.0)
    arg = np.where(outer, np.power(diff, b) / (2.0 * a), 0.0)
    arg = np.minimum(arg, _MAX_EX2 * _MAX_EX2 / 2.0)
    e = np.exp(-arg)
    # Note: the C code uses ax = |x| (not x) inside the bracket.
    bracket = 1.0 - (b / (2.0 * a)) * ax * np.power(np.where(outer, diff, 1.0), b - 1.0)
    out[outer] = (e * bracket)[outer]
    return out


def wgt_ggw(x: ArrayLike, k: ArrayLike) -> NDArray[np.float64]:
    a, b, c, _ = _ggw_unpack(np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel())
    x = np.asarray(x, dtype=np.float64)
    ax = np.abs(x)
    inner = ax <= c
    diff = np.where(inner, 0.0, ax - c)
    arg = np.where(inner, 0.0, np.power(diff, b) / (2.0 * a))
    arg = np.minimum(arg, _MAX_EX2 * _MAX_EX2 / 2.0)
    return np.exp(-arg)


def _ggw_general_rho(x: ArrayLike, k: NDArray[np.float64]) -> NDArray[np.float64]:
    """Numerical-integration fallback for ggw when k[0] == 0.

    rho(x) = (1 / rho_inf) * Integral_{0}^{|x|} psi(t) dt
    where rho_inf = k[3] is the supplied normalising constant.
    """
    from scipy import integrate

    a, b, c = float(k[1]), float(k[2]), float(k[3])

    def _psi_scalar(t: float) -> float:
        diff = max(0.0, abs(t) - c)
        arg = (diff**a) / (2.0 * b)
        return t * np.exp(-arg) if arg < (_MAX_EX2 * _MAX_EX2 / 2.0) else 0.0

    rho_inf = max(k[3], 1.0)  # caller normalises; safety floor
    x_arr = np.asarray(x, dtype=np.float64).ravel()
    out = np.empty_like(x_arr)
    for i, xi in enumerate(x_arr):
        if xi == 0.0:
            out[i] = 0.0
        else:
            val, _ = integrate.quad(_psi_scalar, 0.0, abs(xi), epsabs=1e-10)
            out[i] = val / rho_inf
    return out.reshape(np.asarray(x).shape)


# ===========================================================================
# Family dispatch table
# ===========================================================================
PsiFn = Callable[[ArrayLike, ArrayLike], NDArray[np.float64]]


# Normalisation: for the chi-shaped families used in lmrob, ``rho``/``Mchi``
# is normalised so that ``rho(infinity) = 1``, while ``psi``/``Mpsi`` returns
# the unnormalised derivative. Therefore ``chi'(x) = (factor) * psi(x)``
# where ``factor`` depends on the family and tuning constants.
def _chi_prime_factor(family: str, k_arr: NDArray[np.float64]) -> float:
    """Constant ``chi'(x) / psi(x)`` for each family.

    Verified against ``robustbase::Mchi(x, cc, fam, deriv=1) /
    robustbase::Mpsi(x, cc, fam)``.
    """
    fam = family.lower()
    if fam in ("bisquare", "biweight"):
        # rho_unnorm(inf) = c^2/6  =>  chi'(x) = (6/c^2) * psi(x)
        c = float(k_arr[0])
        return 6.0 / (c * c)
    if fam == "huber":
        # rho_unnorm(inf) = +inf for huber; Mchi normalises differently
        # but we don't use vcov_avar1 with huber today. Leave as 1.
        return 1.0
    if fam == "hampel":
        a, b, r = float(k_arr[0]), float(k_arr[1]), float(k_arr[2])
        nc = a * (b + r - a) / 2.0
        return 1.0 / nc
    if fam == "optimal":
        # rho_opt(x) = (x/c)^2 / 6.5 in the inner branch; rho'(x)/psi(x) =
        # (2x/(6.5 c^2)) / x = 1/(3.25 c^2)
        c = float(k_arr[0])
        return 1.0 / (3.25 * c * c)
    if fam == "lqq":
        # leading inner branch: rho = (3s-3)/denom * x^2  =>  rho'/psi = 6(s-1)/denom
        b, c, s = float(k_arr[0]), float(k_arr[1]), float(k_arr[2])
        denom = s * c * (3.0 * c + 2.0 * b) + (b + c) ** 2
        return 6.0 * (s - 1.0) / denom
    if fam == "ggw":
        # Case-dependent; sampled from R's Mchi(deriv=1)/Mpsi at default tuning.
        # Cases 1..6 mirror SET_ABC_GGW.
        ggw_factors = {
            1: 0.1883291308,
            2: 0.3565452618,
            3: 2.6680355468,
            4: 0.2092091351,
            5: 0.4087348267,
            6: 2.4955990111,
        }
        case_idx = int(k_arr[0])
        if 1 <= case_idx <= 6:
            return ggw_factors[case_idx]
        return 1.0  # user-specified ggw: would need numerical integration
    return 1.0


def chi_prime(x: ArrayLike, family: str, k: ArrayLike) -> NDArray[np.float64]:
    """``Mchi(x, c, family, deriv=1)`` — derivative of the normalised chi."""
    k_arr = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()
    psi_fn = _dispatch(family, "psi")
    factor = _chi_prime_factor(family, k_arr)
    return factor * psi_fn(x, k_arr)


_FAMILY_FNS: dict[str, dict[str, PsiFn]] = {
    "huber": {
        "rho": rho_huber,
        "psi": psi_huber,
        "psi_prime": psi_prime_huber,
        "wgt": wgt_huber,
    },
    "bisquare": {
        "rho": rho_biwgt,
        "psi": psi_biwgt,
        "psi_prime": psi_prime_biwgt,
        "wgt": wgt_biwgt,
    },
    "biweight": {  # alias
        "rho": rho_biwgt,
        "psi": psi_biwgt,
        "psi_prime": psi_prime_biwgt,
        "wgt": wgt_biwgt,
    },
    "hampel": {
        "rho": rho_hmpl,
        "psi": psi_hmpl,
        "psi_prime": psi_prime_hmpl,
        "wgt": wgt_hmpl,
    },
    "optimal": {
        "rho": rho_opt,
        "psi": psi_opt,
        "psi_prime": psi_prime_opt,
        "wgt": wgt_opt,
    },
    "lqq": {
        "rho": rho_lqq,
        "psi": psi_lqq,
        "psi_prime": psi_prime_lqq,
        "wgt": wgt_lqq,
    },
    "ggw": {
        "rho": rho_ggw,
        "psi": psi_ggw,
        "psi_prime": psi_prime_ggw,
        "wgt": wgt_ggw,
    },
}


def _dispatch(family: str, kind: str) -> PsiFn:
    fam = family.lower()
    if fam not in _FAMILY_FNS:
        raise ValueError(f"unknown psi family {family!r}; expected one of {sorted(_FAMILY_FNS)}")
    if kind not in _FAMILY_FNS[fam]:
        raise ValueError(f"unknown kind {kind!r}; expected one of {sorted(_FAMILY_FNS[fam])}")
    return _FAMILY_FNS[fam][kind]
