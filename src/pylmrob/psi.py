# SPDX-License-Identifier: GPL-3.0-or-later
"""Public interface to psi/chi/weight functions.

The numerical kernels live in :mod:`pylmrob._psifuns` (NumPy reference)
and will be mirrored by :mod:`pylmrob._core._psi` (Cython) in a future
performance pass.

R cross-reference:

- ``psi.psi(x, family, k)``       <-> ``robustbase::Mpsi(x, k, family)``
- ``psi.rho(x, family, k)``       <-> ``robustbase::Mchi(x, k, family)``
  (for the chi-shaped families used in lmrob; identical here)
- ``psi.psi_prime(x, family, k)`` <-> ``robustbase::Mpsi(x, k, family, deriv=1)``
- ``psi.wgt(x, family, k)``       <-> ``robustbase::Mwgt(x, k, family)``
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from pylmrob import _psifuns

# Public type aliases. ``PsiFamily`` is the constrained Literal users should
# pass; functions also accept plain ``str`` for runtime convenience and
# validate at call time inside the dispatch table.
PsiFamily = Literal["bisquare", "huber", "hampel", "optimal", "ggw", "lqq"]
TuningK = float | tuple[float, ...] | ArrayLike


def _eval(kind: str, x: ArrayLike, family: str, k: TuningK) -> np.ndarray:
    fn = _psifuns._dispatch(family, kind)
    out = fn(x, k)
    if np.isscalar(x):
        return out.reshape(())
    return out


def rho(x: ArrayLike, family: str, k: TuningK) -> np.ndarray:
    return _eval("rho", x, family, k)


def psi(x: ArrayLike, family: str, k: TuningK) -> np.ndarray:
    return _eval("psi", x, family, k)


def psi_prime(x: ArrayLike, family: str, k: TuningK) -> np.ndarray:
    return _eval("psi_prime", x, family, k)


def wgt(x: ArrayLike, family: str, k: TuningK) -> np.ndarray:
    return _eval("wgt", x, family, k)


# ---------------------------------------------------------------------------
# E[psi^2] and E[psi'] under the standard normal (Phase 7 needs these for
# .vcov.avar1). For Phase 2 we expose them as numerical-integration based
# helpers; closed-form versions for huber/bisquare can land later.
# ---------------------------------------------------------------------------
def Epsi2(family: str, k: TuningK) -> float:
    from scipy import integrate

    def integrand(t: float) -> float:
        psi_val = float(psi(np.array([t]), family, k)[0])
        return psi_val * psi_val * np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

    val, _ = integrate.quad(integrand, -50.0, 50.0, epsabs=1e-12)
    return float(val)


def EDpsi(family: str, k: TuningK) -> float:
    from scipy import integrate

    def integrand(t: float) -> float:
        d_val = float(psi_prime(np.array([t]), family, k)[0])
        return d_val * np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

    val, _ = integrate.quad(integrand, -50.0, 50.0, epsabs=1e-12)
    return float(val)


# ---------------------------------------------------------------------------
# Tuning constants
# ---------------------------------------------------------------------------
# Default tuning constants matching R's lmrob.control() for each family.
# These are the values used at 95% efficiency for the M-step ("psi") and at
# 50% breakdown for the S-step ("chi"). Source: robustbase R/lmrob.R
# .Mpsi.tuning.default() and .Mchi.tuning.default() tables.
# Values here are taken directly from those R tables.
_PSI_TUNING_DEFAULT_PSI: dict[str, tuple[float, ...]] = {
    "huber": (1.345,),
    "bisquare": (4.685061,),
    "biweight": (4.685061,),
    "hampel": (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085),
    "optimal": (1.060158,),
    "ggw": (1.0, 1.5, 0.5, 1.694),  # case 0 (b=1, 95% eff)
    "lqq": (1.4734061, 0.9826779, 1.5),
}

_PSI_TUNING_DEFAULT_CHI: dict[str, tuple[float, ...]] = {
    "huber": (0.6745,),  # for ~50% bdp via huber proposal-2 scale
    "bisquare": (1.547645,),
    "biweight": (1.547645,),
    "hampel": (1.5 * 0.2119163, 3.5 * 0.2119163, 8.0 * 0.2119163),
    "optimal": (0.4047,),
    "ggw": (1.0, 1.5, 0.5, 0.4375470),  # case 2 (b=1, bp 0.5)
    "lqq": (0.4015457, 0.2676971, 1.5),
}


def tuning_for_efficiency(family: str, efficiency: float = 0.95) -> tuple[float, ...]:
    """Return tuning constants giving the requested asymptotic efficiency.

    For Phase 2 we only support the canonical 95%-efficiency tuning that R
    uses by default. Other efficiency targets will land in Phase 8 alongside
    the full ``Control`` machinery.
    """
    if abs(efficiency - 0.95) > 1e-6:
        raise NotImplementedError("tuning_for_efficiency only supports efficiency=0.95 in Phase 2.")
    fam = family.lower()
    if fam not in _PSI_TUNING_DEFAULT_PSI:
        raise ValueError(f"unknown psi family: {family!r}")
    return _PSI_TUNING_DEFAULT_PSI[fam]


def tuning_for_breakdown(family: str, breakdown: float = 0.5) -> tuple[float, ...]:
    """Return tuning constants giving the requested breakdown (chi side)."""
    if abs(breakdown - 0.5) > 1e-6:
        raise NotImplementedError("tuning_for_breakdown only supports breakdown=0.5 in Phase 2.")
    fam = family.lower()
    if fam not in _PSI_TUNING_DEFAULT_CHI:
        raise ValueError(f"unknown psi family: {family!r}")
    return _PSI_TUNING_DEFAULT_CHI[fam]


__all__ = [
    "EDpsi",
    "Epsi2",
    "PsiFamily",
    "psi",
    "psi_prime",
    "rho",
    "tuning_for_breakdown",
    "tuning_for_efficiency",
    "wgt",
]
