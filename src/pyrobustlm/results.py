# SPDX-License-Identifier: GPL-3.0-or-later
"""Result object returned by ``lmrob`` fits."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pyrobustlm.control import Control
    from pyrobustlm.summary import SummaryLmRob


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

    def confint(self, level: float = 0.95) -> np.ndarray:
        from scipy.stats import norm

        z = norm.ppf((1 + level) / 2)
        se = self.standard_errors_
        lo = self.coef_ - z * se
        hi = self.coef_ + z * se
        return np.column_stack([lo, hi])

    def predict(self, new_data: object) -> np.ndarray:
        """Predict on new data.

        Accepts either:

        - a pandas ``DataFrame`` with the columns referenced by the original
          formula. The fit's stored formulaic ``ModelSpec`` re-applies any
          factor encoding, interactions, ``I(x**2)`` transforms, etc.
        - a 2-D NumPy array already shaped ``(n, p)``, matching the original
          design (intercept column included if the formula had one).
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
                from typing import Any

                df: Any = new_data
                cols = [c for c in self.term_names_ if c not in ("Intercept", "(Intercept)")]
                missing = [c for c in cols if c not in df.columns]
                if missing:
                    raise ValueError(
                        f"predict(DataFrame): missing columns in new_data: {missing!r}"
                    )
                feat = np.asarray(df.loc[:, cols].to_numpy(), dtype=np.float64)
                has_intercept = bool(
                    self.term_names_
                    and self.term_names_[0] in ("Intercept", "(Intercept)")
                )
                arr = (
                    np.column_stack([np.ones(feat.shape[0]), feat])
                    if has_intercept
                    else feat
                )
            else:
                from pyrobustlm.formula import apply_spec

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
        return arr @ self.coef_

    def summary(self) -> SummaryLmRob:
        """Return a ``SummaryLmRob`` matching R's ``summary.lmrob``.

        The returned object stringifies to an R-style printout (use
        ``str(fit.summary())`` or just ``print(fit.summary())``) and exposes
        the underlying coefficient table and R-squared as attributes.
        """
        from pyrobustlm.summary import make_summary

        return make_summary(self)

    def __repr__(self) -> str:
        coefs = ", ".join(f"{n}={v:.4g}" for n, v in zip(self.term_names_, self.coef_, strict=True))
        return f"LmRobResults({coefs}; scale={self.scale_:.4g})"
