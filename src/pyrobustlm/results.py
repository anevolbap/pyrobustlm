# SPDX-License-Identifier: GPL-3.0-or-later
"""Result object returned by ``lmrob`` fits."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pyrobustlm.control import Control


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

    def predict(self, new_X: np.ndarray) -> np.ndarray:
        return np.asarray(new_X, dtype=np.float64) @ self.coef_

    def summary(self) -> str:
        from scipy.stats import norm

        se = self.standard_errors_
        z = self.coef_ / np.where(se > 0, se, 1)
        pvals = 2 * (1 - norm.cdf(np.abs(z)))
        rows = []
        rows.append(
            f"lmrob fit (setting={self.control.setting}, psi={self.control.psi}, init={self.control.init})"
        )
        rows.append(
            f"n = {self.nobs_}, df residual = {self.df_residual_}, scale = {self.scale_:.6g}"
        )
        rows.append(f"converged = {self.converged_} (iters: {self.n_iter_})")
        rows.append("")
        header = f"{'term':<20} {'estimate':>12} {'std.err':>10} {'z':>8} {'p':>10}"
        rows.append(header)
        rows.append("-" * len(header))
        for name, b, s, zi, p in zip(self.term_names_, self.coef_, se, z, pvals, strict=True):
            rows.append(f"{name:<20} {b:12.6g} {s:10.4g} {zi:8.3g} {p:10.4g}")
        return "\n".join(rows)

    def __repr__(self) -> str:
        coefs = ", ".join(f"{n}={v:.4g}" for n, v in zip(self.term_names_, self.coef_, strict=True))
        return f"LmRobResults({coefs}; scale={self.scale_:.4g})"
