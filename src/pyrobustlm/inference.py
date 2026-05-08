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
    init_residuals: NDArray[np.float64] | None = None,
    chi_family: str | None = None,
    chi_k: float | tuple[float, ...] | None = None,
    bb: float = 0.5,
) -> NDArray[np.float64]:
    """Asymptotic sandwich covariance per robustbase ``.vcov.avar1``
    (R/lmrob.MM.R:510-577).

    Formula::

        r.s  = r / scale   (final scaled residuals)
        r0.s = r0 / scale  (initial S scaled residuals)
        w  = Mpsi'(r.s, c.psi)        # psi'
        w0 = Mchi'(r0.s, c.chi)       # chi'
        A  = (X' diag(w) X)^{-1} * scale
        a  = A * X' (Mpsi(r.s) * r.s) / mean(w0 * r0.s)
        Xww = X' (Mpsi(r.s) * Mchi(r0.s))
        u1 = A * X' diag(Mpsi(r.s)^2) X * (n A)
        u2 = a * Xww' * A
        u3 = A * Xww * a'
        u4 = mean(Mchi(r0.s)^2 - bb^2) * a a'
        cov = (u1 - u2 - u3 + u4) / n

    Parameters
    ----------
    X :
        Design matrix used in the fit, shape (n, p).
    residuals :
        Final MM residuals.
    sigma :
        Final scale (from the S step).
    psi_family, psi_k :
        Efficiency-tuned psi family and constants used in the M step.
    init_residuals :
        Initial S-step residuals. Defaults to ``residuals`` (which
        approximates the formula but is less accurate).
    chi_family, chi_k :
        Chi family and tuning. Default to ``psi_family`` / ``psi_k`` if not
        supplied (for default lmrob, ``chi_family == psi_family``).
    bb :
        Consistency constant from the S-step (default 0.5).
    """
    n, p = X.shape
    if init_residuals is None:
        init_residuals = residuals
    if chi_family is None:
        chi_family = psi_family
    if chi_k is None:
        chi_k = psi_k

    sgma = max(sigma, 1e-300)
    r_s = residuals / sgma
    r0_s = init_residuals / sgma

    # In robustbase: Mpsi(., deriv=1) = psi' (the unnormalised psi-prime).
    # Mchi is the normalised rho such that chi(infinity) = 1, so its
    # derivative chi' = (1/rho_unnorm(inf)) * psi.
    from pyrobustlm._psifuns import chi_prime as _chi_prime

    w_pp = _psi.psi_prime(r_s, psi_family, psi_k)  # Mpsi(r.s, deriv=1) = psi'
    w0_pp = _chi_prime(r0_s, chi_family, chi_k)  # Mchi(r0.s, deriv=1) = chi'
    psi_rs = _psi.psi(r_s, psi_family, psi_k)  # Mpsi(r.s, deriv=0) = psi
    chi_r0s = _psi.rho(r0_s, chi_family, chi_k)  # Mchi(r0.s, deriv=0) = chi

    # x.wx = X' diag(w_pp) X
    XwX = X.T @ (X * w_pp[:, None])
    try:
        A = np.linalg.solve(XwX, np.eye(p)) * sgma
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError("vcov_avar1: X' W X is singular; consider cov='.vcov.w'") from exc

    denom = float(np.mean(w0_pp * r0_s))
    if denom == 0.0:
        raise FloatingPointError("vcov_avar1: mean(chi'(r0/s) * r0/s) = 0")

    # ``a`` uses w = psi' (R/lmrob.MM.R:547,563), not w = psi.
    a_vec = A @ (X.T @ (w_pp * r_s)) / denom

    Xww = X.T @ (psi_rs * chi_r0s)

    # u1: A * X' diag(psi^2) X * (n A)
    u1 = A @ (X.T @ (X * (psi_rs * psi_rs)[:, None])) @ (n * A)
    # u2: a Xww' A   (a column outer with Xww row -> p x p)
    u2 = np.outer(a_vec, Xww) @ A
    # u3: A Xww a'
    u3 = A @ np.outer(Xww, a_vec)
    # u4: scalar * a a'
    u4 = float(np.mean(chi_r0s * chi_r0s - bb * bb)) * np.outer(a_vec, a_vec)

    cov = (u1 - u2 - u3 + u4) / n

    # Force symmetry; project to PSD if needed (matches R's "posdefify").
    cov = 0.5 * (cov + cov.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    if (eigvals < 0).any():
        eigvals = np.where(eigvals < 0, 0.0, eigvals)
        cov = (eigvecs * eigvals) @ eigvecs.T
        cov = 0.5 * (cov + cov.T)

    return cov


def _asympt_corrfact(psi_family: str, psi_k: float | tuple[float, ...]) -> float:
    """Asymptotic correction factor used by ``.vcov.w``
    (robustbase R/lmrob.MM.R:415-437).

    For the default tuning of each family, robustbase has hardcoded values.
    Otherwise we compute via numerical integration of psi^2 / (r psi)^2
    under the standard normal.
    """
    fam = psi_family.lower()
    k_arr = np.atleast_1d(np.asarray(psi_k, dtype=float)).ravel()

    DEFAULT_CORRFACT = {
        "bisquare": 1.0526317574,
        "welsh": 1.0526704649,
        "optimal": 1.0526419204,
        "hampel": 1.0526016980,
        "lqq": 1.0526365291,
    }
    DEFAULT_TUNING = {
        "bisquare": np.array([4.685061]),
        "huber": np.array([1.345]),
        "hampel": np.array([1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085]),
        "optimal": np.array([1.060158]),
        "lqq": np.array([1.4734061, 0.9826779, 1.5]),
    }
    if fam in DEFAULT_CORRFACT:
        d = DEFAULT_TUNING.get(fam)
        if d is not None and d.shape == k_arr.shape and np.allclose(d, k_arr):
            return DEFAULT_CORRFACT[fam]

    # Fallback: numerical integration
    from scipy import integrate

    def _phi(t: float) -> float:
        return float(np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi))

    def num_inner(t: float) -> float:
        v = float(_psi.psi(np.array([t]), fam, psi_k)[0])
        return v * v * _phi(t)

    def den_inner(t: float) -> float:
        v = float(_psi.psi(np.array([t]), fam, psi_k)[0])
        return t * v * _phi(t)

    num, _ = integrate.quad(num_inner, -50.0, 50.0, epsabs=1e-12)
    den, _ = integrate.quad(den_inner, -50.0, 50.0, epsabs=1e-12)
    return num / (den * den)


def vcov_w(
    X: NDArray[np.float64],
    residuals: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: float | tuple[float, ...],
    rweights: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Koller & Stahel (2011) sandwich, ``cov.corrfact = "asympt"`` branch.

    Implements the asymptotic-correction-factor branch of robustbase's
    ``.vcov.w`` (R/lmrob.MM.R:354):

        V = sigma^2 * corrfact * (X' W X)^{-1}

    where ``W = diag(rweights)`` and ``corrfact`` is the per-family value
    tabulated at default tuning, falling back to numerical integration
    otherwise.

    The Huber finite-sample correction (``cov.hubercorr=TRUE``) and the
    other ``cov.corrfact`` branches (``empirical``, ``hybrid``, ``tau``,
    ``tauold``) are not yet implemented; they primarily affect inference
    in small samples and do not change point estimates.
    """
    if rweights is None:
        z = residuals / sigma if sigma != 0 else residuals
        rweights = _psi.wgt(z, psi_family, psi_k)
    sw = np.sqrt(np.maximum(rweights, 0.0))
    Xw = X * sw[:, None]
    XtWX = Xw.T @ Xw
    XtWX_inv = np.linalg.inv(XtWX)
    corrfact = _asympt_corrfact(psi_family, psi_k)
    return (sigma**2) * corrfact * XtWX_inv


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
