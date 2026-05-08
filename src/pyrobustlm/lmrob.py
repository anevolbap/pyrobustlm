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
from pyrobustlm.formula import model_matrix
from pyrobustlm.inference import vcov_avar1, vcov_w
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

    y, X, term_names = model_matrix(formula, data)
    n, p = X.shape
    if n <= p:
        raise ValueError(f"need n > p; got n={n}, p={p}")

    # ------------------------------------------------------------------
    # Initial S estimate via fast-S resampling
    # ------------------------------------------------------------------
    cfg = FastSConfig(
        psi_chi=control.psi,
        k_chi=tuple(np.atleast_1d(np.asarray(control.tuning_chi, dtype=float)).ravel()),
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
    s_seed = seed if seed is not None else control.seed
    s_result = fast_s(X, y, cfg=cfg, seed=s_seed)

    init_info: dict[str, object] = {
        "coef": s_result.coef.copy(),
        "scale": s_result.scale,
        "n_iter": s_result.n_iter,
        "method": "S",
    }

    # ------------------------------------------------------------------
    # MM step holding sigma fixed.
    # ------------------------------------------------------------------
    psi_k_eff = tuple(np.atleast_1d(np.asarray(control.tuning_psi, dtype=float)).ravel())
    mm = mm_iterate(
        X=X,
        y=y,
        beta_init=s_result.coef,
        sigma=s_result.scale,
        psi_family=control.psi,
        psi_k=psi_k_eff,
        max_it=control.max_it,
        rel_tol=control.rel_tol,
    )

    coef = mm.coef
    sigma = s_result.scale
    residuals = y - X @ coef
    fitted = X @ coef
    z = residuals / sigma if sigma != 0 else residuals
    rweights = _psi.wgt(z, control.psi, psi_k_eff)

    # ------------------------------------------------------------------
    # Covariance
    # ------------------------------------------------------------------
    cov_fn = {
        ".vcov.avar1": vcov_avar1,
        ".vcov.w": vcov_w,
        "Asymp": vcov_avar1,
    }.get(control.cov, vcov_avar1)
    cov = cov_fn(X, residuals, sigma, control.psi, psi_k_eff)

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
        if self._result is None:
            raise RuntimeError("call .fit before .predict")
        # Re-attach intercept column.
        X = np.asarray(X, dtype=np.float64)
        # Add intercept if the fit included one (always true in our default formula).
        ones = np.ones((X.shape[0], 1), dtype=np.float64)
        return self._result.predict(np.hstack([ones, X]))

    @property
    def result_(self) -> LmRobResults:
        if self._result is None:
            raise RuntimeError("call .fit before accessing result_")
        return self._result
