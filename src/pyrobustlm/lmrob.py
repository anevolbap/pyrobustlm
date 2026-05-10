# SPDX-License-Identifier: GPL-3.0-or-later
"""Top-level ``lmrob`` entry points.

End-to-end pipeline: formula -> design matrix -> fast-S init -> MM
iteration -> covariance -> :class:`pyrobustlm.results.LmRobResults`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from pyrobustlm import psi as _psi
from pyrobustlm._fast_s import FastSConfig, fast_s
from pyrobustlm._mm import mm_iterate
from pyrobustlm.control import Control
from pyrobustlm.d_scale import d_scale
from pyrobustlm.formula import model_matrix
from pyrobustlm.inference import vcov_avar1, vcov_w
from pyrobustlm.ms_estimator import m_s_fit
from pyrobustlm.results import LmRobResults

if TYPE_CHECKING:
    import pandas as pd


def _to_chi_psi_family(psi_family: str) -> str:
    """The chi family used in S iteration is the same as the user's psi family."""
    return psi_family


def lmrob(
    formula: str,
    data: pd.DataFrame,
    control: Control | None = None,
    weights: np.ndarray | None = None,
    na_action: str = "drop",
    seed: int | np.random.Generator | None = None,
    **kwargs: Any,
) -> LmRobResults:
    """Fit a robust MM linear regression.

    Parameters
    ----------
    formula :
        R-style formula, e.g. ``"y ~ x1 + x2 + x3"``. Parsed by
        :mod:`formulaic`.
    data :
        DataFrame containing the columns referenced by ``formula``.
    control :
        Algorithm parameters; defaults to ``Control()`` (KS2014 preset).
    weights :
        Optional case weights. Currently raises if non-None (Phase 8+).
    na_action :
        ``"drop"`` (default) drops rows with any NA before fitting.
    seed :
        Seed for the resampling RNG.

    Returns
    -------
    LmRobResults
    """
    if weights is not None:
        raise NotImplementedError("Per-case weights are not yet supported.")
    if kwargs:
        raise TypeError(f"unexpected keyword arguments: {sorted(kwargs)!r}")
    if control is None:
        control = Control()
    # Control.__post_init__ guarantees these are populated; narrow for the
    # type checker.
    assert control.psi is not None
    assert control.method is not None
    assert control.cov is not None
    psi_family: str = control.psi
    cov_kind: str = control.cov

    # ------------------------------------------------------------------
    # Design matrix
    # ------------------------------------------------------------------
    if na_action == "drop":
        data = data.dropna()
    elif na_action == "raise":
        if data.isna().any().any():
            raise ValueError("data contains NA values; pass na_action='drop' to skip them")
    else:
        raise ValueError(f"unknown na_action: {na_action!r}")

    design = model_matrix(formula, data)
    y, X, term_names = design.y, design.X, design.term_names
    n, p = X.shape
    if n <= p:
        raise ValueError(f"need n > p; got n={n}, p={p}")

    # ------------------------------------------------------------------
    # Pick init method. ``init="auto"`` chooses M-S when factor columns
    # exist and S otherwise. Anything else honors the explicit user request.
    # ------------------------------------------------------------------
    init_method = control.init
    if init_method == "auto":
        init_method = "M-S" if design.is_factor_col.any() else "S"

    s_seed = seed if seed is not None else control.seed
    k_chi_tuple = tuple(np.atleast_1d(np.asarray(control.tuning_chi, dtype=float)).ravel())

    if init_method == "S":
        # ------------------------------------------------------------------
        # Initial S estimate via fast-S resampling
        # ------------------------------------------------------------------
        cfg = FastSConfig(
            psi_chi=psi_family,
            k_chi=k_chi_tuple,
            b0=control.bb,
            nResample=control.nResample,
            k_fast_s=control.k_fast_s,
            best_r=control.best_r_s,
            max_it=control.max_it,
            refine_tol=control.refine_tol,
            scale_tol=control.scale_tol,
            max_iter_scale=control.k_max,
            mts=control.mts,
            n_workers=control.n_workers,
        )
        s_result = fast_s(X, y, cfg=cfg, seed=s_seed)
        beta_init = s_result.coef
        sigma_init = s_result.scale
        init_info: dict[str, object] = {
            "coef": s_result.coef.copy(),
            "scale": s_result.scale,
            "n_iter": s_result.n_iter,
            "method": "S",
        }
    elif init_method == "M-S":
        if not design.is_factor_col.any():
            # No factor cols; M-S degenerates to S.
            cfg = FastSConfig(
                psi_chi=psi_family,
                k_chi=k_chi_tuple,
                b0=control.bb,
                nResample=control.nResample,
                k_fast_s=control.k_fast_s,
                best_r=control.best_r_s,
                max_it=control.max_it,
                refine_tol=control.refine_tol,
                scale_tol=control.scale_tol,
                max_iter_scale=control.k_max,
                mts=control.mts,
            )
            s_result = fast_s(X, y, cfg=cfg, seed=s_seed)
            beta_init = s_result.coef
            sigma_init = s_result.scale
            init_info = {
                "coef": s_result.coef.copy(),
                "scale": s_result.scale,
                "n_iter": s_result.n_iter,
                "method": "S",
            }
        else:
            cat_cols = design.is_factor_col
            X_cat = X[:, cat_cols]
            X_cont = X[:, ~cat_cols]
            ms_result = m_s_fit(
                X_cat=X_cat,
                X_cont=X_cont,
                y=y,
                psi_chi=psi_family,
                k_chi=k_chi_tuple,
                b0=control.bb,
                k_m_s=control.k_m_s,
                nResample=control.nResample,
                max_it=control.max_it,
                rel_tol=control.rel_tol,
                seed=s_seed,
            )
            # Re-stitch into the original column order.
            beta_init = np.empty(p, dtype=np.float64)
            beta_init[cat_cols] = ms_result.coef_cat
            beta_init[~cat_cols] = ms_result.coef_cont
            sigma_init = ms_result.scale
            init_info = {
                "coef": beta_init.copy(),
                "scale": ms_result.scale,
                "n_iter": ms_result.n_iter,
                "method": "M-S",
            }
    else:
        raise NotImplementedError(f"init={init_method!r} not implemented")

    # ------------------------------------------------------------------
    # MM step holding sigma fixed.
    # ------------------------------------------------------------------
    psi_k_eff = tuple(np.atleast_1d(np.asarray(control.tuning_psi, dtype=float)).ravel())
    mm = mm_iterate(
        X=X,
        y=y,
        beta_init=beta_init,
        sigma=sigma_init,
        psi_family=psi_family,
        psi_k=psi_k_eff,
        max_it=control.max_it,
        rel_tol=control.rel_tol,
    )

    coef = mm.coef
    sigma = sigma_init
    residuals = y - X @ coef
    fitted = X @ coef
    z = residuals / sigma if sigma != 0 else residuals
    rweights = _psi.wgt(z, psi_family, psi_k_eff)

    # ------------------------------------------------------------------
    # D-step (used by methods SMD and SMDM): refines scale via design-
    # adaptive weighting. SMDM also re-fits MM with the new scale.
    # ``tau_vec`` is also stashed for ``vcov_w(corrfact='tau' | 'hybrid' | 'tauold')``.
    # ------------------------------------------------------------------
    method = control.method or "MM"
    tau_vec: np.ndarray | None = None
    if "D" in method:
        sigma_d, d_converged, tau_vec, _h = d_scale(
            X=X,
            residuals=residuals,
            rweights=rweights,
            init_scale=sigma,
            family=psi_family,
            c_psi=psi_k_eff,
            max_iter=control.k_max,
            tol=control.rel_tol,
        )
        if d_converged and sigma_d > 0:
            sigma = sigma_d
            init_info["d_scale"] = sigma_d
            init_info["d_converged"] = True
            # The trailing M in SMDM: re-fit MM with the new D-scale.
            if method.endswith("M"):
                mm2 = mm_iterate(
                    X=X,
                    y=y,
                    beta_init=coef,
                    sigma=sigma,
                    psi_family=psi_family,
                    psi_k=psi_k_eff,
                    max_it=control.max_it,
                    rel_tol=control.rel_tol,
                )
                coef = mm2.coef
                residuals = y - X @ coef
                fitted = X @ coef
                z = residuals / sigma if sigma != 0 else residuals
                rweights = _psi.wgt(z, psi_family, psi_k_eff)

    # ------------------------------------------------------------------
    # Covariance
    # ------------------------------------------------------------------
    init_residuals = y - X @ beta_init
    if cov_kind == ".vcov.avar1":
        cov = vcov_avar1(
            X=X,
            residuals=residuals,
            sigma=sigma,
            psi_family=psi_family,
            psi_k=psi_k_eff,
            init_residuals=init_residuals,
            chi_family=psi_family,
            chi_k=k_chi_tuple,
            bb=control.bb,
        )
    elif cov_kind == ".vcov.w":
        # Pre-compute tau if it isn't already populated by the D-step.
        # vcov_w with corrfact="tau"/"hybrid"/"tauold" needs it.
        if tau_vec is None:
            from pyrobustlm.d_scale import tau as _compute_tau

            sw = np.sqrt(np.maximum(rweights, 0.0))
            Xw = X * sw[:, None]
            Q, _ = np.linalg.qr(Xw, mode="reduced")
            h = np.minimum(1.0, np.sum(Q * Q, axis=1))
            tau_vec = _compute_tau(h, psi_family, psi_k_eff)
        s_init_scale_raw = init_info.get("scale", sigma)
        s_init_scale = (
            float(s_init_scale_raw) if isinstance(s_init_scale_raw, (int, float)) else float(sigma)
        )  # type: ignore[arg-type]
        cov = vcov_w(
            X=X,
            residuals=residuals,
            sigma=sigma,
            psi_family=psi_family,
            psi_k=psi_k_eff,
            rweights=rweights,
            init_residuals=init_residuals,
            init_scale=s_init_scale,
            chi_k=k_chi_tuple,
            method=method,
            tau=tau_vec,
        )
    else:
        # Fallback: legacy / unknown
        cov = vcov_avar1(
            X=X,
            residuals=residuals,
            sigma=sigma,
            psi_family=psi_family,
            psi_k=psi_k_eff,
        )

    return LmRobResults(
        coef_=np.asarray(coef, dtype=np.float64),
        scale_=float(sigma),
        weights_=np.asarray(rweights, dtype=np.float64),
        rweights_=np.asarray(rweights, dtype=np.float64),
        residuals_=np.asarray(residuals, dtype=np.float64),
        fitted_=np.asarray(fitted, dtype=np.float64),
        cov_=np.asarray(cov, dtype=np.float64),
        df_residual_=n - p,
        converged_=mm.converged,
        n_iter_=mm.n_iter,
        nobs_=n,
        term_names_=term_names,
        control=control,
        init_=init_info,
        rhs_spec_=design.rhs_spec,
        design_x_=np.asarray(X, dtype=np.float64),
        design_y_=np.asarray(y, dtype=np.float64),
    )


class LmRob:
    """scikit-learn-style estimator wrapper around :func:`lmrob`."""

    def __init__(self, control: Control | None = None) -> None:
        self.control = control or Control()
        self._result: LmRobResults | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, seed: int | None = None) -> LmRob:
        # Build a tiny pandas DataFrame so we can reuse the formula path.
        import pandas as pd

        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if X.ndim != 2:
            raise ValueError("X must be 2-D")
        cols = [f"x{i}" for i in range(X.shape[1])]
        df = pd.DataFrame(X, columns=cols)
        df["y"] = y
        formula = "y ~ " + " + ".join(cols)
        self._result = lmrob(formula, df, control=self.control, seed=seed)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict on a new design matrix (raw, without intercept column).

        ``LmRob`` always fits with an intercept, so we wrap ``X`` in a
        DataFrame with the same column names used at fit time and let the
        stored formula spec re-add the intercept.
        """
        if self._result is None:
            raise RuntimeError("call .fit before .predict")
        import pandas as pd

        X = np.asarray(X, dtype=np.float64)
        cols = [f"x{i}" for i in range(X.shape[1])]
        return self._result.predict(pd.DataFrame(X, columns=cols))

    @property
    def result_(self) -> LmRobResults:
        if self._result is None:
            raise RuntimeError("call .fit before accessing result_")
        return self._result
