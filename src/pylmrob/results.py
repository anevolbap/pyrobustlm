# SPDX-License-Identifier: GPL-3.0-or-later
"""Result object returned by ``lmrob`` fits."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    # NOTE: do *not* import pylmrob.bootstrap here. bootstrap.py imports
    # LmRobResults from this module (also under TYPE_CHECKING), and
    # CodeQL's py/unsafe-cyclic-import flags the bidirectional pair
    # even though both are TYPE_CHECKING-only. The bootstrap() method
    # below uses a forward-reference string annotation instead.
    from pylmrob.control import Control
    from pylmrob.summary import SummaryLmRob


@dataclass
class LmRobResults:
    """Output of an ``lmrob`` fit.

    Attributes mirror R's ``lmrob`` object where practical.
    """

    coef_: np.ndarray
    scale_: float
    weights_: np.ndarray  # final IRWLS weights = wgt(r/sigma)
    rweights_: np.ndarray  # alias for R's $rweights (same as weights_ here)
    residuals_: np.ndarray
    fitted_: np.ndarray
    cov_: np.ndarray
    df_residual_: int
    converged_: bool
    n_iter_: int
    nobs_: int
    term_names_: list[str]
    control: Control
    init_: dict[str, object] = field(default_factory=dict)
    # formulaic ModelSpec for the RHS, used by ``predict`` to re-apply the
    # design transformation (factor encoding, I(x**2), etc.) to new data.
    rhs_spec_: object | None = None
    # Design matrix and response stashed for downstream operations
    # (e.g. anova(test="Deviance") refits reduced models on the full scale).
    design_x_: np.ndarray | None = None
    design_y_: np.ndarray | None = None

    # Backwards-compatible alias to mirror R's $coefficients.
    @property
    def coefficients_(self) -> np.ndarray:
        return self.coef_

    @property
    def standard_errors_(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov_))

    def confint(
        self,
        level: float = 0.95,
        method: str = "wald",
        *,
        n_boot: int = 1000,
        seed: int | None = None,
        n_workers: int = 1,
        kind: str = "percentile",
    ) -> np.ndarray:
        """Confidence intervals for the regression coefficients.

        Two methods:

        - ``"wald"`` (default): asymptotic normal CIs from the sandwich
          covariance. ``z * se`` where ``z`` is the standard-normal
          quantile at ``(1 + level) / 2``.
        - ``"bootstrap"``: case-resampling bootstrap; runs
          :meth:`bootstrap` internally and returns the requested ``kind``
          (``"percentile"`` or ``"basic"``) CIs.

        Parameters
        ----------
        level
            Coverage level, e.g. ``0.95``.
        method
            ``"wald"`` (default) or ``"bootstrap"``.
        n_boot, seed, n_workers
            Forwarded to :meth:`bootstrap` when ``method="bootstrap"``.
        kind
            Bootstrap CI kind: ``"percentile"`` or ``"basic"``. Ignored
            for ``method="wald"``.
        """
        if method == "wald":
            from scipy.stats import norm

            z = norm.ppf((1 + level) / 2)
            se = self.standard_errors_
            lo = self.coef_ - z * se
            hi = self.coef_ + z * se
            return np.column_stack([lo, hi])
        if method == "bootstrap":
            boot = self.bootstrap(n_boot=n_boot, level=level, seed=seed, n_workers=n_workers)
            ci = boot.percentile_ci if kind == "percentile" else boot.basic_ci
            return np.asarray(ci, dtype=np.float64)
        raise ValueError(f"method must be 'wald' or 'bootstrap', got {method!r}")

    # ------------------------------------------------------------------
    # statsmodels-style attribute aliases. Let pylmrob fits drop into
    # statsmodels.regression-shaped code without adapters.
    # ------------------------------------------------------------------

    @property
    def params(self) -> np.ndarray:
        return self.coef_

    @property
    def bse(self) -> np.ndarray:
        return self.standard_errors_

    @property
    def tvalues(self) -> np.ndarray:
        se = self.standard_errors_
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(se > 0, self.coef_ / np.where(se > 0, se, 1.0), np.nan)

    @property
    def pvalues(self) -> np.ndarray:
        from scipy.stats import t as t_dist

        return 2.0 * t_dist.sf(np.abs(self.tvalues), df=self.df_residual_)

    def conf_int(self, alpha: float = 0.05) -> np.ndarray:
        """``statsmodels`` spelling of :meth:`confint`. Uses ``1 - alpha``."""
        return self.confint(level=1.0 - alpha)

    def predict(
        self,
        new_data: object,
        *,
        interval: str = "none",
        level: float = 0.95,
    ) -> np.ndarray:
        """Predict on new data, optionally with confidence/prediction bands.

        Accepts either:

        - a pandas ``DataFrame`` with the columns referenced by the original
          formula. The fit's stored formulaic ``ModelSpec`` re-applies any
          factor encoding, interactions, ``I(x**2)`` transforms, etc.
        - a 2-D NumPy array already shaped ``(n, p)``, matching the original
          design (intercept column included if the formula had one).

        Parameters
        ----------
        interval :
            ``"none"`` (default) returns the point predictions, shape ``(n,)``.
            ``"confidence"`` returns ``(n, 3)`` columns ``(fit, lwr, upr)``
            with the confidence interval for the *mean* response at each
            new observation (Var = X^T cov X).
            ``"prediction"`` returns ``(n, 3)`` with the prediction interval
            for a single *new* observation (Var = sigma^2 + X^T cov X).
        level :
            Confidence level for the interval. Default 0.95.

        Bands use the t-distribution with ``df_residual_`` degrees of
        freedom, mirroring R's ``predict.lm`` / ``predict.lmrob`` convention.
        """
        # pandas is a hard dependency, but only required when ``new_data`` is
        # a DataFrame. We dispatch by duck-typing to avoid a top-level import.
        is_dataframe = type(new_data).__name__ == "DataFrame" and hasattr(new_data, "columns")

        arr: np.ndarray
        if is_dataframe:
            if self.rhs_spec_ is None:
                # Fit used the simple-formula fast path; rebuild the design
                # by selecting term columns from new_data directly. Cast to
                # ``Any`` because we duck-typed ``new_data`` as a DataFrame
                # without importing pandas.
                df: Any = new_data
                cols = [c for c in self.term_names_ if c not in ("Intercept", "(Intercept)")]
                missing = [c for c in cols if c not in df.columns]
                if missing:
                    raise ValueError(
                        f"predict(DataFrame): missing columns in new_data: {missing!r}"
                    )
                feat = np.asarray(df.loc[:, cols].to_numpy(), dtype=np.float64)
                has_intercept = bool(
                    self.term_names_ and self.term_names_[0] in ("Intercept", "(Intercept)")
                )
                arr = np.column_stack([np.ones(feat.shape[0]), feat]) if has_intercept else feat
            else:
                from pylmrob.formula import apply_spec

                arr = apply_spec(self.rhs_spec_, new_data)
        else:
            arr = np.asarray(new_data, dtype=np.float64)

        if arr.ndim != 2:
            raise ValueError(f"predict expected a 2-D matrix; got shape {arr.shape}")
        if arr.shape[1] != self.coef_.size:
            raise ValueError(
                f"predict: design has {arr.shape[1]} columns but the fit "
                f"has {self.coef_.size} coefficients"
            )
        point = arr @ self.coef_
        if interval == "none":
            return point
        if interval not in ("confidence", "prediction"):
            raise ValueError(
                f"interval must be one of 'none', 'confidence', 'prediction'; got {interval!r}"
            )
        # Var(X beta) per row = sum((X cov) * X, axis=1) = diag(X cov X^T).
        var_fit = np.einsum("ij,jk,ik->i", arr, self.cov_, arr)
        var_fit = np.maximum(var_fit, 0.0)  # numerical safety
        var = var_fit + (self.scale_**2 if interval == "prediction" else 0.0)
        from scipy.stats import t as t_dist

        q = t_dist.ppf((1.0 + level) / 2.0, df=self.df_residual_)
        se = np.sqrt(var)
        return np.column_stack([point, point - q * se, point + q * se])

    def predict_std(self, new_data: object, *, kind: str = "confidence") -> np.ndarray:
        """Standard deviation of the prediction at each new observation.

        Returns ``sqrt(Var(X^T beta_hat))`` (``kind="confidence"``, default)
        or ``sqrt(sigma^2 + Var(X^T beta_hat))`` (``kind="prediction"``).
        Use this when you want the raw SE and intend to build your own
        intervals (e.g. with a non-Gaussian distribution or for a Bayesian
        update); :meth:`predict` with ``interval="confidence"`` already
        returns t-quantile-scaled bands.
        """
        # Build the design matrix the same way predict() does, then return
        # only the standard error column.
        # Re-use predict() by extracting the (upr - point) / q value, but
        # that's wasteful. Cheaper to inline the design build.
        is_dataframe = type(new_data).__name__ == "DataFrame" and hasattr(new_data, "columns")
        arr: np.ndarray
        if is_dataframe:
            if self.rhs_spec_ is None:
                df: Any = new_data
                cols = [c for c in self.term_names_ if c not in ("Intercept", "(Intercept)")]
                missing = [c for c in cols if c not in df.columns]
                if missing:
                    raise ValueError(f"predict_std(DataFrame): missing columns: {missing!r}")
                feat = np.asarray(df.loc[:, cols].to_numpy(), dtype=np.float64)
                has_intercept = bool(
                    self.term_names_ and self.term_names_[0] in ("Intercept", "(Intercept)")
                )
                arr = np.column_stack([np.ones(feat.shape[0]), feat]) if has_intercept else feat
            else:
                from pylmrob.formula import apply_spec

                arr = apply_spec(self.rhs_spec_, new_data)
        else:
            arr = np.asarray(new_data, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != self.coef_.size:
            raise ValueError(
                f"predict_std: design shape {arr.shape} incompatible with "
                f"{self.coef_.size} coefficients"
            )
        var_fit = np.einsum("ij,jk,ik->i", arr, self.cov_, arr)
        var_fit = np.maximum(var_fit, 0.0)
        if kind == "confidence":
            return np.sqrt(var_fit)
        if kind == "prediction":
            return np.sqrt(var_fit + self.scale_ * self.scale_)
        raise ValueError(f"kind must be 'confidence' or 'prediction'; got {kind!r}")

    def diagnostics(self, outlier_threshold: float = 2.5) -> object:
        """Per-observation diagnostic statistics.

        Returns a :class:`pylmrob.diagnostics.DiagnosticsTable` with
        leverage, robust Cook's distance, standardized residuals, the
        robust weights, and a boolean outlier flag
        (``|std_residuals| > outlier_threshold``).

        Requires the fit to have a stashed design matrix
        (``design_x_``); the default ``lmrob()`` call always stashes it.
        """
        if self.design_x_ is None:
            raise RuntimeError(
                "diagnostics() needs the design matrix; fit was created without design_x_"
            )
        from pylmrob.diagnostics import DiagnosticsTable, cooks_distance, hatvalues

        h = hatvalues(self, self.design_x_)
        cd = cooks_distance(self, self.design_x_)
        sigma = max(self.scale_, 1e-300)
        z = self.residuals_ / sigma
        outliers = np.abs(z) > outlier_threshold
        # Masked outliers: rows that the robust fit rejected (rweights ~ 0)
        # that ALSO look far from the clean data's design centre. The "clean
        # data" is everything with rweights > 0.5. We compute leverage of
        # every row against the clean-data X'X; the contamination ends up
        # with high values because they sit outside the clean centroid in
        # X-space. Pure OLS leverage doesn't work here because the
        # contamination's own weight pulls the OLS centroid toward itself
        # (this is the masking phenomenon).
        X = self.design_x_
        n = X.shape[0]
        p = self.coef_.size
        wt_eps = 1e-3
        clean = self.rweights_ > 0.5
        masked_outliers = np.zeros(n, dtype=bool)
        if int(clean.sum()) > p:
            X_clean = X[clean]
            try:
                xtx_inv_clean = np.linalg.inv(X_clean.T @ X_clean)
                clean_lev = np.einsum("ij,jk,ik->i", X, xtx_inv_clean, X)
                # Threshold: median of the clean-rows' own leverages plus
                # 3 * MAD. Robust to the contamination block.
                clean_lev_in = clean_lev[clean]
                med = float(np.median(clean_lev_in))
                mad = float(np.median(np.abs(clean_lev_in - med))) * 1.4826
                lev_thresh = med + 3.0 * max(mad, 1e-12)
                masked_outliers = (clean_lev > lev_thresh) & (self.rweights_ < wt_eps)
            except np.linalg.LinAlgError:
                pass
        return DiagnosticsTable(
            leverage=h,
            cooks_distance=cd,
            std_residuals=z,
            rweights=self.rweights_,
            outliers=outliers,
            masked_outliers=masked_outliers,
        )

    def bootstrap(
        self,
        n_boot: int = 1000,
        level: float = 0.95,
        seed: int | np.random.Generator | None = None,
        n_workers: int = 1,
    ) -> Any:  # BootstrapResult; see note at top of file re cyclic import
        """Method-style spelling of :func:`pylmrob.bootstrap`.

        Equivalent to ``pylmrob.bootstrap(self, n_boot=..., ...)``;
        matches the ``fit.anova()`` / ``fit.diagnostics()`` style.
        """
        from pylmrob.bootstrap import bootstrap as _bootstrap

        return _bootstrap(self, n_boot=n_boot, level=level, seed=seed, n_workers=n_workers)

    def anova(self, *others: LmRobResults, test: str = "Wald") -> object:
        """Method-style spelling of :func:`pylmrob.anova`.

        Equivalent to ``pylmrob.anova(self, *others, test=test)``; lets
        you write ``full.anova(reduced)`` instead of the free-function
        form, matching R's idiom.
        """
        from pylmrob.anova import anova as _anova

        return _anova(self, *others, test=test)

    def summary(self, style: str = "r") -> SummaryLmRob:
        """Return a ``SummaryLmRob`` matching R's ``summary.lmrob``.

        Parameters
        ----------
        style
            ``"r"`` (default): R-style ``summary.lmrob`` output, matching
            ``robustbase`` line-for-line where practical.
            ``"statsmodels"``: a fixed-width table matching the
            ``statsmodels.iolib.summary.Summary`` layout. Use this when
            piping pylmrob fits into statsmodels-shaped reporting code.

        The returned object stringifies to the chosen style via
        ``str()`` / ``print()``; ``.render(style=...)`` overrides the
        choice without rebuilding the object.
        """
        from pylmrob.summary import make_summary

        out = make_summary(self)
        out.style = style
        return out

    def __repr__(self) -> str:
        """Compact R-style ``print.lmrob`` output.

        For the full coefficient table with std errors, t / p values and
        R-squared, call ``print(self.summary())``.
        """
        method = self.control.method or "MM"
        lines: list[str] = [
            f'lmrob(method="{method}", psi="{self.control.psi}")',
            "",
            "Coefficients:",
        ]
        # Column-aligned coefficient row matching R's print.
        name_w = max(11, max((len(n) for n in self.term_names_), default=11))
        header = "  ".join(f"{n:>{name_w}}" for n in self.term_names_)
        values = "  ".join(f"{v:>{name_w}.4g}" for v in self.coef_)
        lines.append(header)
        lines.append(values)
        return "\n".join(lines)
