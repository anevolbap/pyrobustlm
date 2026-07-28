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

from pylmrob import psi as _psi


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
    from pylmrob._psifuns import chi_prime as _chi_prime

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
        # Internal form from robustbase:::.psi.conv.cc("lqq", <default>).
        # This was 0.9826779, which no longer matched Control's value, so
        # np.allclose below always failed and lqq silently took the
        # numerical-integration fallback instead of this table.
        "lqq": np.array([1.4734061, 0.9822707, 1.5]),
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


_VALID_CORRFACT = ("tau", "empirical", "asympt", "hybrid", "tauold")
_VALID_DFCORR = ("mean", "none", "mn.vc", "varc", "mn.df")
_VALID_RESID = ("final", "initial", "trick")


def _hybrid_mpp2(psi_family: str, psi_k: float | tuple[float, ...]) -> float:
    """Asymptotic ``E[r psi(r)]^2`` under N(0,1) for the hybrid correction.

    Hardcoded for the default tuning of each family (R/lmrob.MM.R:476-487).
    Falls back to numerical integration otherwise.
    """
    fam = psi_family.lower()
    k_arr = np.atleast_1d(np.asarray(psi_k, dtype=float)).ravel()

    DEFAULT_MPP2 = {
        "bisquare": 0.5742327,
        "welsh": 0.5445068,
        "optimal": 0.8598825,
        "hampel": 0.6775217,
        "lqq": 0.6883393,
    }
    DEFAULT_TUNING = {
        "bisquare": np.array([4.685061]),
        "huber": np.array([1.345]),
        "hampel": np.array([1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085]),
        "optimal": np.array([1.060158]),
        "lqq": np.array([1.4734061, 0.9822707, 1.5]),
    }
    if fam in DEFAULT_MPP2:
        d = DEFAULT_TUNING.get(fam)
        if d is not None and d.shape == k_arr.shape and np.allclose(d, k_arr):
            return DEFAULT_MPP2[fam]

    from scipy import integrate

    def integrand(t: float) -> float:
        v = float(_psi.psi(np.array([t]), fam, psi_k)[0])
        return t * v * np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi)

    val, _ = integrate.quad(
        integrand,
        -50.0,
        50.0,
        points=[0.0, -1.0, 1.0, -3.0, 3.0],
        epsabs=1e-12,
    )
    return val * val


def vcov_w(
    X: NDArray[np.float64],
    residuals: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: float | tuple[float, ...],
    rweights: NDArray[np.float64] | None = None,
    *,
    init_residuals: NDArray[np.float64] | None = None,
    init_scale: float | None = None,
    chi_k: float | tuple[float, ...] | None = None,
    method: str = "MM",
    cov_corrfact: str | None = None,
    cov_hubercorr: bool | None = None,
    cov_dfcorr: str | None = None,
    cov_resid: str = "final",
    cov_xwx: bool = True,
    tau: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Koller & Stahel (2011) sandwich, all five ``cov.corrfact`` branches.

    Direct port of robustbase's ``.vcov.w`` (R/lmrob.MM.R:354-507). The
    formula is::

        V = scale^2 * sscorr * corrfact * Cinv

    where ``Cinv = (R^{-1})(R^{-1})'`` from the QR of ``sqrt(w) X`` (with
    ``w = rweights`` when ``cov_xwx=True``), ``corrfact`` and ``sscorr``
    depend on the chosen branches, and ``scale`` is the final scale.

    Parameters
    ----------
    X, residuals, sigma, psi_family, psi_k :
        Required model inputs.
    rweights :
        Robustness weights from the MM fit. Recomputed from residuals
        and ``psi_k`` if omitted.
    init_residuals, init_scale, chi_k :
        S-step residuals, S-step scale, and chi tuning. Needed for
        ``cov_resid in {"initial", "trick"}``.
    method :
        Fit method (``"MM"``, ``"SMD"``, ``"SMDM"``, ...). The default
        ``cov_hubercorr`` is ``True`` when ``"D"`` is not in the method.
    cov_corrfact :
        One of ``{"asympt", "empirical", "tau", "hybrid", "tauold"}``.
        Default depends on ``cov_hubercorr``.
    cov_hubercorr :
        Apply Huber's finite-sample correction. Defaults to
        ``"D" not in method``.
    cov_dfcorr :
        One of ``{"mean", "mn.vc", "none", "varc", "mn.df"}``.
    cov_resid :
        ``"final"`` (default), ``"initial"``, or ``"trick"`` -- selects
        which residuals/scale standardise the inputs.
    cov_xwx :
        Whether ``Cinv`` uses the rweights-weighted QR (default True).
    tau :
        Per-observation design factors. Required for
        ``cov_corrfact in {"tau", "hybrid", "tauold"}``.
    """
    if cov_corrfact is not None and cov_corrfact not in _VALID_CORRFACT:
        raise ValueError(f"cov_corrfact must be one of {_VALID_CORRFACT}")
    if cov_dfcorr is not None and cov_dfcorr not in _VALID_DFCORR:
        raise ValueError(f"cov_dfcorr must be one of {_VALID_DFCORR}")
    if cov_resid not in _VALID_RESID:
        raise ValueError(f"cov_resid must be one of {_VALID_RESID}")

    # Defaults (mirror lmrob.MM.R:362-377)
    if cov_hubercorr is None:
        cov_hubercorr = "D" not in method
    if cov_corrfact is None:
        cov_corrfact = "empirical" if cov_hubercorr else "tau"
    if cov_dfcorr is None:
        cov_dfcorr = "mn.vc" if (cov_hubercorr or cov_corrfact in ("tau", "hybrid")) else "mean"

    # Tuning constant for the residuals' standardisation. R picks tuning.chi
    # for cov_resid='initial' or method in (S, SD); else tuning.psi.
    if cov_resid == "initial" or method in ("S", "SD"):
        c_psi = chi_k if chi_k is not None else psi_k
    else:
        c_psi = psi_k

    if rweights is None:
        z = residuals / sigma if sigma != 0 else residuals
        rweights = _psi.wgt(z, psi_family, psi_k)

    n = X.shape[0]
    w = rweights if cov_xwx else np.ones(n)
    sw = np.sqrt(np.maximum(w, 0.0))
    Xw = X * sw[:, None]
    # Cinv = (X' W X)^{-1} via QR of sqrt(w) * X.
    _, R = np.linalg.qr(Xw, mode="reduced")
    try:
        Rinv = np.linalg.solve(R, np.eye(R.shape[0]))
    except np.linalg.LinAlgError:
        Rinv = np.linalg.pinv(R)
    cinv = Rinv @ Rinv.T
    p = R.shape[0]

    # Correction factor
    if cov_corrfact == "asympt":
        corrfact = _asympt_corrfact(psi_family, c_psi)
        varcorr = 1.0
    else:
        # Standardised residuals
        if cov_resid == "initial":
            if init_residuals is None or init_scale is None:
                raise ValueError("cov_resid='initial' needs init_residuals and init_scale")
            rstand = np.asarray(init_residuals) / float(init_scale)
        elif cov_resid == "trick":
            if init_residuals is None or init_scale is None:
                raise ValueError("cov_resid='trick' needs init_residuals and init_scale")
            rstand = np.asarray(init_residuals) / float(init_scale)
        else:  # "final"
            rstand = np.asarray(residuals) / float(sigma)

        # tau per-observation factor
        if cov_corrfact in ("tau", "hybrid", "tauold"):
            if tau is None:
                raise ValueError(f"cov_corrfact={cov_corrfact!r} needs the tau vector")
            tau_vec = np.asarray(tau, dtype=np.float64)
        else:
            tau_vec = np.ones(n)

        rstand = rstand / tau_vec
        r_psi = _psi.psi(rstand, psi_family, c_psi)
        r_psipr = _psi.psi_prime(rstand, psi_family, c_psi)

        mpp = float(np.mean(r_psipr))
        mpp2 = mpp * mpp

        # Huber correction
        if cov_hubercorr:
            vpp = float(np.sum((r_psipr - mpp) ** 2)) / n
            hcorr = (1.0 + p / n * vpp / mpp2) ** 2
        else:
            hcorr = 1.0

        # varcorr (R/lmrob.MM.R:472-473)
        if cov_corrfact == "tau" and not np.allclose(tau_vec, 1.0):
            varcorr = 1.0 / float(np.mean(tau_vec * tau_vec))
        else:
            varcorr = n / max(n - p, 1)

        # Hybrid: replace mpp2 by tabulated asymptotic E[r psi(r)]^2
        if cov_corrfact == "hybrid":
            mpp2 = _hybrid_mpp2(psi_family, c_psi)

        tau_factor = np.ones(n) if cov_corrfact == "tauold" else tau_vec * tau_vec
        corrfact = float(np.mean(tau_factor * r_psi * r_psi)) / mpp2 * hcorr

    # Sample-size correction (sscorr)
    mean_w = float(np.mean(w))
    if cov_dfcorr == "mean":
        sscorr = mean_w
    elif cov_dfcorr == "mn.vc":
        sscorr = mean_w * varcorr
    elif cov_dfcorr == "none":
        sscorr = 1.0
    elif cov_dfcorr == "varc":
        sscorr = varcorr
    elif cov_dfcorr == "mn.df":
        sw_sum = float(np.sum(w))
        if sw_sum <= p:
            raise FloatingPointError("vcov_w: sum(w) <= p; cov_dfcorr='mn.df' undefined")
        sscorr = mean_w * mean_w / (1.0 - p / sw_sum)
    else:  # pragma: no cover  -- already validated above
        raise ValueError(f"unknown cov_dfcorr: {cov_dfcorr}")

    cov = (sigma**2) * sscorr * corrfact * cinv
    cov = 0.5 * (cov + cov.T)
    return cov


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
