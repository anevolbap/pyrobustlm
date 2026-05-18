# SPDX-License-Identifier: GPL-3.0-or-later
"""Robust M-scale estimator.

Direct port of robustbase/src/lmrob.c::find_scale (lmrob.c:2611). Iteration:

    scale_{k+1} = scale_k * sqrt(mean(chi(r/scale_k)) / b0)

where ``mean`` averages over n - p (not n), matching the C source's
``sum_rho_sc`` (lmrob.c:2656). MAD is used as the initial scale when none is
supplied.
"""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np

from pylmrob import _psifuns as _pf


def _try_import_cpsi() -> Any | None:
    """Return the compiled ``pylmrob._core._psi`` module, or ``None``.

    Imported lazily through ``importlib`` so static type checkers don't
    flag the Cython extension as missing.
    """
    try:
        return importlib.import_module("pylmrob._core._psi")
    except ImportError:
        return None


# GGW (case_idx 1..6) -> (a, b, c) used by both Cython and NumPy paths.
# Mirrors SET_ABC_GGW in robustbase/src/lmrob.c (lines 1279-1293).
_GGW_ABC: dict[int, tuple[float, float, float]] = {
    1: (0.648, 1.0, 1.694),
    2: (0.4760508, 1.0, 1.2442567),
    3: (0.1674046, 1.0, 0.4375470),
    4: (1.387, 1.5, 1.063),
    5: (0.8372485, 1.5, 0.7593544),
    6: (0.2036741, 1.5, 0.2959132),
}


def _cython_wgt(
    z: np.ndarray,
    family: str,
    k: tuple[float, ...] | float | np.ndarray,
) -> np.ndarray | None:
    """Compute weights via the Cython kernel for ``family``, or return None.

    Returns ``None`` when the compiled kernel is unavailable or the family
    isn't supported on the fast path; callers fall back to NumPy.
    """
    cpsi = _try_import_cpsi()
    if cpsi is None:
        return None
    fam = family.lower()
    z_buf = np.ascontiguousarray(z, dtype=np.float64)
    out = np.empty_like(z_buf)
    k_arr = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()

    if fam in ("bisquare", "biweight"):
        cpsi.bisquare_wgt(z_buf, float(k_arr[0]), out)
        return out
    if fam == "huber":
        cpsi.huber_wgt(z_buf, float(k_arr[0]), out)
        return out
    if fam == "hampel":
        cpsi.hampel_wgt(z_buf, float(k_arr[0]), float(k_arr[1]), float(k_arr[2]), out)
        return out
    if fam == "optimal":
        cpsi.optimal_wgt(z_buf, float(k_arr[0]), out)
        return out
    if fam == "lqq":
        cpsi.lqq_wgt(z_buf, float(k_arr[0]), float(k_arr[1]), float(k_arr[2]), out)
        return out
    if fam == "ggw":
        case_idx = int(k_arr[0])
        if 1 <= case_idx <= 6:
            a, b, c = _GGW_ABC[case_idx]
            cpsi.ggw_wgt(z_buf, a, b, c, out)
            return out
        if case_idx == 0 and k_arr.size >= 4:
            cpsi.ggw_wgt(z_buf, float(k_arr[1]), float(k_arr[2]), float(k_arr[3]), out)
            return out
    return None


def _cython_rho(
    z: np.ndarray,
    family: str,
    k: tuple[float, ...] | float | np.ndarray,
) -> np.ndarray | None:
    """Same shape as ``_cython_wgt`` but for the chi/rho function.

    GGW rho still needs the polynomial table from _psifuns; we don't have
    a Cython port yet, so we return None for ggw to fall back.
    """
    cpsi = _try_import_cpsi()
    if cpsi is None:
        return None
    fam = family.lower()
    z_buf = np.ascontiguousarray(z, dtype=np.float64)
    out = np.empty_like(z_buf)
    k_arr = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()

    if fam in ("bisquare", "biweight"):
        cpsi.bisquare_rho(z_buf, float(k_arr[0]), out)
        return out
    if fam == "huber":
        cpsi.huber_rho(z_buf, float(k_arr[0]), out)
        return out
    if fam == "hampel":
        cpsi.hampel_rho(z_buf, float(k_arr[0]), float(k_arr[1]), float(k_arr[2]), out)
        return out
    if fam == "optimal":
        cpsi.optimal_rho(z_buf, float(k_arr[0]), out)
        return out
    if fam == "lqq":
        cpsi.lqq_rho(z_buf, float(k_arr[0]), float(k_arr[1]), float(k_arr[2]), out)
        return out
    # ggw rho uses tabulated polynomials in _psifuns; no Cython port yet.
    return None


_DEFAULT_K_CHI: dict[str, tuple[float, ...]] = {
    "huber": (1.345,),
    "bisquare": (1.547645,),
    "biweight": (1.547645,),
    "hampel": (1.5 * 0.2119163, 3.5 * 0.2119163, 8.0 * 0.2119163),
    "optimal": (0.4047,),
    "lqq": (0.4015457, 0.2676971, 1.5),
    "ggw": (3,),  # case 3: b=1, bp=0.5
}


def _mad(x: np.ndarray) -> float:
    """Median Absolute Deviation, scaled to be a consistent estimator of sd
    under the Gaussian. Matches R's ``mad(x, constant = 1.4826)`` with
    ``center=median(x)``.
    """
    med = np.median(x)
    return float(1.4826 * np.median(np.abs(x - med)))


def m_scale(
    r: np.ndarray,
    family: str = "bisquare",
    k: float | tuple[float, ...] | None = None,
    b0: float = 0.5,
    max_iter: int = 200,
    tol: float = 1e-10,
    init_scale: float | None = None,
    p: int = 0,
) -> float:
    """Compute the robust M-scale of residuals ``r``.

    Parameters
    ----------
    r :
        Residual vector (1-D, finite).
    family :
        Chi family. Defaults to ``"bisquare"``.
    k :
        Tuning constant(s). When ``None``, uses the family's default 50%
        breakdown constants (matching R).
    b0 :
        Consistency constant (target value of ``E[chi(Z)]`` under N(0,1)).
        Defaults to 0.5 for 50% breakdown.
    max_iter, tol :
        Iteration limit and relative-change convergence threshold.
    init_scale :
        Starting scale. ``None`` -> use MAD.
    p :
        Number of regression parameters; the inner mean is over ``n - p``,
        not ``n``. Default 0 reduces to the plain mean.

    Returns
    -------
    float
        Estimated scale. May raise ``RuntimeError`` if non-finite residuals
        are encountered; emits a ``RuntimeWarning`` (not exception) on
        non-convergence to match R's behaviour.
    """
    r = np.asarray(r, dtype=np.float64).ravel()
    if not np.isfinite(r).all():
        raise ValueError("m_scale: residuals must be finite")

    n = r.size
    if n - p <= 0:
        raise ValueError(f"m_scale: n - p = {n - p} <= 0")

    if k is None:
        fam_lower = family.lower()
        if fam_lower not in _DEFAULT_K_CHI:
            raise ValueError(f"m_scale: unknown family {family!r}")
        k = _DEFAULT_K_CHI[fam_lower]

    if init_scale is None:
        s = _mad(r)
        # If MAD is zero (perfect fit), match R's zero-tol behaviour.
        if s == 0.0:
            return 0.0
    else:
        if init_scale <= 0:
            import warnings

            warnings.warn(
                f"find_scale(*, initial_scale = {init_scale} <= 0) -> final scale = 0",
                RuntimeWarning,
                stacklevel=2,
            )
            return 0.0
        s = float(init_scale)

    fam_lower = family.lower()
    # Cython fast path: fully-inlined m_scale loops per family.
    cpsi = _try_import_cpsi()
    if cpsi is not None:
        k_arr = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()
        r_buf = np.ascontiguousarray(r, dtype=np.float64)
        if fam_lower in ("bisquare", "biweight"):
            return float(cpsi.m_scale_bisquare(r_buf, float(k_arr[0]), b0, s, max_iter, tol, p))
        if fam_lower == "hampel":
            return float(
                cpsi.m_scale_hampel(
                    r_buf,
                    float(k_arr[0]),
                    float(k_arr[1]),
                    float(k_arr[2]),
                    b0,
                    s,
                    max_iter,
                    tol,
                    p,
                )
            )
        if fam_lower == "optimal":
            return float(cpsi.m_scale_optimal(r_buf, float(k_arr[0]), b0, s, max_iter, tol, p))
        if fam_lower == "lqq":
            return float(
                cpsi.m_scale_lqq(
                    r_buf,
                    float(k_arr[0]),
                    float(k_arr[1]),
                    float(k_arr[2]),
                    b0,
                    s,
                    max_iter,
                    tol,
                    p,
                )
            )
        if fam_lower == "ggw":
            case_idx = int(k_arr[0])
            if 1 <= case_idx <= 6:
                return float(cpsi.m_scale_ggw_case(r_buf, case_idx, b0, s, max_iter, tol, p))

    # Generic loop. Use the Cython chi when available, else NumPy.
    rho_fn_np = _pf._dispatch(family, "rho")
    k_arr = np.asarray(k, dtype=np.float64)

    prev = s
    for _it in range(max_iter):
        z = r / s
        chi_cy = _cython_rho(z, family, k_arr)
        chi_vals = chi_cy if chi_cy is not None else rho_fn_np(z, k_arr)
        mean_chi = float(np.sum(chi_vals)) / (n - p)
        s = s * np.sqrt(mean_chi / b0)
        if abs(s - prev) <= tol * prev:
            return float(s)
        prev = s

    import warnings

    warnings.warn(
        f"m_scale did not converge in {max_iter} iterations (tol={tol})",
        RuntimeWarning,
        stacklevel=2,
    )
    return float(s)
