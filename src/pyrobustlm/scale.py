# SPDX-License-Identifier: GPL-3.0-or-later
"""Robust M-scale estimator.

Direct port of robustbase/src/lmrob.c::find_scale (lmrob.c:2611). Iteration:

    scale_{k+1} = scale_k * sqrt(mean(chi(r/scale_k)) / b0)

where ``mean`` averages over n - p (not n), matching the C source's
``sum_rho_sc`` (lmrob.c:2656). MAD is used as the initial scale when none is
supplied.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pyrobustlm import psi as _psi

if TYPE_CHECKING:
    from pyrobustlm.control import PsiFamily


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
    family: PsiFamily = "bisquare",
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

    prev = s
    for _it in range(max_iter):
        # mean(chi(r/s)) over (n - p) -- matches sum_rho_sc in the C source.
        chi_vals = _psi.rho(r / s, family, k)
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
