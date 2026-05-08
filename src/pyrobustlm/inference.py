# SPDX-License-Identifier: GPL-3.0-or-later
"""Covariance estimators for lmrob.

Phase 7. Implements ``vcov_avar1`` matching robustbase's ``.vcov.avar1``
(default for ``setting="KS2014"`` and the legacy MM setting).

References
----------
- Yohai (1987), eqs (3.1)-(3.4)
- Maronna, Martin & Yohai (2006), Chapter 4

The asymptotic-variance matrix is

    V = sigma^2 * E[psi'(z)]^-2 * E[psi(z)^2] * (X' X / n)^{-1}

with ``z = r / sigma``. We estimate the expectations from the finite
sample (sandwich form):

    A = sum(psi'(r_i / sigma)) / n
    B = sum(psi(r_i / sigma)^2) / n
    V_hat = sigma^2 * (B / A^2) * (X'X / n)^{-1}

This is the same expression robustbase uses inside ``.vcov.avar1`` modulo a
``n / (n - p)`` finite-sample correction we apply explicitly.

``vcov_w`` (KS2011) and ``vcov_asymp`` are placeholders for now.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from pyrobustlm import psi as _psi


def vcov_avar1(
    X: NDArray[np.float64],
    residuals: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: float | tuple[float, ...],
) -> NDArray[np.float64]:
    """Asymptotic sandwich covariance matrix.

    Parameters
    ----------
    X :
        Design matrix used in the fit, shape (n, p).
    residuals :
        Residual vector ``y - X @ beta_MM``.
    sigma :
        Final scale (from the S step).
    psi_family, psi_k :
        Efficiency-tuned psi family and constants used in the M step.
    """
    n, p = X.shape
    z = residuals / sigma if sigma != 0 else residuals
    psi_vals = _psi.psi(z, psi_family, psi_k)
    psi_prime_vals = _psi.psi_prime(z, psi_family, psi_k)

    A = float(np.mean(psi_prime_vals))
    B = float(np.mean(psi_vals**2))
    if A == 0.0:
        raise FloatingPointError("vcov_avar1: mean(psi'(r/sigma)) = 0")

    XtX = X.T @ X
    XtX_inv = np.linalg.inv(XtX)
    factor = (sigma**2) * B / (A * A)
    # Apply the n/(n-p) bias correction R uses for small n.
    factor *= n / max(n - p, 1)
    return factor * (XtX_inv * n) / n  # left factored so units stay clear


def vcov_w(
    X: NDArray[np.float64],
    residuals: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: float | tuple[float, ...],
) -> NDArray[np.float64]:
    """Koller & Stahel (2011) finite-sample-corrected sandwich.

    Phase 7+ (deferred). Falls back to ``vcov_avar1`` for now with a
    warning.
    """
    import warnings

    warnings.warn(
        "vcov_w (KS2011 finite-sample correction) not yet implemented; falling back to vcov_avar1.",
        RuntimeWarning,
        stacklevel=2,
    )
    return vcov_avar1(X, residuals, sigma, psi_family, psi_k)


def vcov_asymp(
    X: NDArray[np.float64],
    residuals: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: float | tuple[float, ...],
) -> NDArray[np.float64]:
    """Legacy ``Asymp`` covariance estimator.

    For the families used in lmrob, this is identical to ``vcov_avar1`` up
    to the small-sample bias correction. Provided as an alias.
    """
    return vcov_avar1(X, residuals, sigma, psi_family, psi_k)
