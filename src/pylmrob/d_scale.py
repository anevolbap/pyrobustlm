# SPDX-License-Identifier: GPL-3.0-or-later
"""Design-adaptive D-scale estimator (Koller & Stahel 2014).

Direct port of robustbase's ``lmrob..D..fit`` (R/lmrob.MM.R:817-862) and
the ``R_find_D_scale`` C iteration (src/lmrob.c:2762, type ``dt1`` only,
which is the sole option remaining in robustbase 0.99-7).

Algorithm (after MM with residuals ``r``, robustness weights ``w``,
hat values ``h``, scale ``sigma``)::

    kappa = E[ psi(Z) Z - kappa wgt(Z) ] = 0  under Z ~ N(0,1)
    tau_i = sqrt(1 - tfact * h_i) * (tcorr * h_i + 1)        # fast form
    iterate: w_i  <- wgt(r_i / (tau_i sigma))
             sigma_new = sqrt(sum(r_i^2 w_i) / (sum(w_i tau_i^2) kappa))

The D-scale replaces the S-scale in the final result; the SMDM method
re-fits MM with this new scale.
"""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np

from pylmrob import psi as _psi


def _try_cpsi() -> Any | None:
    try:
        return importlib.import_module("pylmrob._core._psi")
    except ImportError:
        return None


# ---------------------------------------------------------------------------
# kappa: solution of  E[psi(Z) Z - kappa wgt(Z)] = 0  under Z ~ N(0,1)
# Equivalently:       kappa = E[psi(Z) Z] / E[wgt(Z)]
# ---------------------------------------------------------------------------
def kappa(family: str, c_psi: float | tuple[float, ...] | np.ndarray) -> float:
    """Compute the consistency kappa under the standard normal.

    Mirrors ``robustbase::lmrob.kappa``. The integrand for the numerator
    vanishes at zero, which trips up ``scipy.integrate.quad`` defaults; we
    explicitly hand it break points so the adaptive grid samples both lobes.
    """
    from scipy import integrate

    def _psi_at(t: float) -> float:
        return float(_psi.psi(np.array([t]), family, c_psi)[0])

    def _wgt_at(t: float) -> float:
        return float(_psi.wgt(np.array([t]), family, c_psi)[0])

    def num(t: float) -> float:
        return _psi_at(t) * t * np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

    def den(t: float) -> float:
        return _wgt_at(t) * np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

    # Break points at zero and around the typical psi-tuning bandwidths so
    # quad's adaptive grid hits the bell curve's lobes.
    breaks = [0.0, -1.0, 1.0, -3.0, 3.0]
    a, _ = integrate.quad(num, -50.0, 50.0, points=breaks, epsabs=1e-12)
    b, _ = integrate.quad(den, -50.0, 50.0, points=breaks, epsabs=1e-12)
    if b == 0.0:
        raise FloatingPointError("kappa: denominator E[wgt(Z)] = 0")
    return a / b


# ---------------------------------------------------------------------------
# tau: per-observation design-adaptive scale factor.
# robustbase has hardcoded (tfact, tcorr) for the default tunings of each
# family; for other tunings R falls back to numerical integration. We
# replicate the hardcoded table.
# ---------------------------------------------------------------------------
# (psi family, tuning constants converted to R's user-facing form for matching)
# robustbase R/lmrob.MM.R:920-952. Keys are (family, k_tuple) using the
# **internal** tuning shape.
_TAU_FAST_TABLE = {
    ("bisquare", (4.685061,)): (0.9473684, -0.0900833),
    ("biweight", (4.685061,)): (0.9473684, -0.0900833),
    ("optimal", (1.060158,)): (0.94735878, -0.09444537),
    (
        "hampel",
        (
            1.5 * 0.9016085,
            3.5 * 0.9016085,
            8.0 * 0.9016085,
        ),
    ): (0.94739770, -0.04103958),
    # lqq internal form (b, c, s) for default user-facing (-0.5, 1.5, 0.95, NA)
    ("lqq", (1.4734061, 0.9822707, 1.5)): (0.94736359, -0.08594805),
    # ggw cases 1..6 -> internal cc (case_idx,)
    ("ggw", (1,)): (0.9473787, -0.1143846),  # b=1, 95% eff
    ("ggw", (4,)): (0.94741036, -0.08424648),  # b=1.5, 95% eff
}


def tau(
    h: np.ndarray,
    family: str,
    c_psi: float | tuple[float, ...] | np.ndarray,
    fast: bool = True,
) -> np.ndarray:
    """Compute per-observation tau_i values.

    Fast form uses the precomputed ``(tfact, tcorr)`` table for the
    default tuning of each family. Non-fast (numerical-integration)
    fallback is not implemented yet; callers should stick with default
    tuning if using KS2011.
    """
    fam = family.lower()
    k_arr = np.atleast_1d(np.asarray(c_psi, dtype=np.float64)).ravel()
    key = (fam, tuple(k_arr.tolist()))
    if fast and key in _TAU_FAST_TABLE:
        tfact, tcorr = _TAU_FAST_TABLE[key]
        return np.sqrt(1.0 - tfact * h) * (tcorr * h + 1.0)

    # Best-effort numerical fallback: use the bisquare-like values as a
    # starting point and warn. In practice users would compute the table
    # via lmrob.tau.fast.coefs() in R; for non-default tunings this branch
    # is approximate.
    import warnings

    warnings.warn(
        f"tau: no fast coefficients for family={family!r} c_psi={c_psi!r}; "
        "using bisquare-default approximation",
        RuntimeWarning,
        stacklevel=2,
    )
    tfact, tcorr = 0.9473684, -0.0900833
    return np.sqrt(1.0 - tfact * h) * (tcorr * h + 1.0)


# ---------------------------------------------------------------------------
# D-scale iteration (type "dt1" matching robustbase's only enabled type)
# ---------------------------------------------------------------------------
def find_d_scale(
    r: np.ndarray,
    tau_vec: np.ndarray,
    kappa_val: float,
    family: str,
    c_psi: float | tuple[float, ...] | np.ndarray,
    init_scale: float,
    max_iter: int = 200,
    tol: float = 1e-7,
) -> tuple[float, bool]:
    """Iterate the D-scale to convergence.

    Returns ``(scale, converged)``.

    The iteration matches ``R_find_D_scale`` with ``type=3`` (dt1):

    .. math::
        \\sigma_{k+1} = \\sqrt{
            \\frac{\\sum_i r_i^2 \\, w_i}{\\kappa \\, \\sum_i w_i \\, \\tau_i^2}
        }

    where ``w_i = wgt(r_i / (tau_i sigma_k))``.
    """
    sgma = float(init_scale)
    if sgma <= 0.0:
        return 0.0, False
    converged = False
    for _ in range(max_iter):
        z = r / (tau_vec * sgma)
        w = _psi.wgt(z, family, c_psi)
        tsum1 = float(np.sum(r * r * w))
        tsum2 = float(np.sum(w * tau_vec * tau_vec))
        if tsum2 == 0.0:
            return sgma, False
        new = np.sqrt(tsum1 / (tsum2 * kappa_val))
        if abs(new - sgma) < tol * max(tol, sgma):
            sgma = float(new)
            converged = True
            break
        sgma = float(new)
    return sgma, converged


def d_scale(
    X: np.ndarray,
    residuals: np.ndarray,
    rweights: np.ndarray,
    init_scale: float,
    family: str,
    c_psi: float | tuple[float, ...] | np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-7,
) -> tuple[float, bool, np.ndarray, np.ndarray]:
    """Compute the design-adaptive D-scale.

    Convenience wrapper that computes hat values, tau, kappa, and runs
    :func:`find_d_scale`. Returns ``(scale, converged, tau, h)`` so callers
    can reuse the intermediate quantities.
    """
    # Hat diagonal weighted by robustness weights, matching R's
    # ``.lmrob.hat(x, obj$rweights)``.
    sw = np.sqrt(np.maximum(rweights, 0.0))
    Xw = X * sw[:, None]
    Q, _ = np.linalg.qr(Xw, mode="reduced")
    h = np.minimum(1.0, np.sum(Q * Q, axis=1))

    kp = kappa(family, c_psi)
    tau_vec = tau(h, family, c_psi)
    # Starting value (R/lmrob.MM.R:848): sqrt(sum(w r^2) / kappa / sum(tau^2 w))
    num = float(np.sum(rweights * residuals * residuals))
    den = float(np.sum(rweights * tau_vec * tau_vec))
    if den == 0.0 or kp == 0.0:
        return init_scale, False, tau_vec, h
    start = np.sqrt(num / (kp * den))
    if not np.isfinite(start) or start <= 0:
        start = init_scale
    sgma, converged = find_d_scale(residuals, tau_vec, kp, family, c_psi, start, max_iter, tol)
    return sgma, converged, tau_vec, h
