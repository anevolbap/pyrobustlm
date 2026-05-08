# SPDX-License-Identifier: GPL-3.0-or-later
"""Result object returned by ``lmrob`` fits.

Phase 8 will fill this in. Defined as a minimal dataclass for now so that
type hints in other modules resolve.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

    from pyrobustlm.control import Control


@dataclass
class LmRobResults:
    """Output of an ``lmrob`` fit.

    Attribute names match R's ``lmrob`` object where practical, with
    ``snake_case`` versions for Pythonic access.
    """

    coef_: np.ndarray
    scale_: float
    weights_: np.ndarray
    residuals_: np.ndarray
    fitted_: np.ndarray
    cov_: np.ndarray
    df_residual_: int
    converged_: bool
    rweights_: np.ndarray
    nobs_: int
    control: Control
    init_: dict[str, object] = field(default_factory=dict)

    def summary(self) -> str:
        raise NotImplementedError(
            "LmRobResults.summary is not implemented yet; see plan.md Phase 8."
        )

    def predict(
        self,
        new_X: np.ndarray,
        interval: str | None = None,
    ) -> np.ndarray:
        raise NotImplementedError(
            "LmRobResults.predict is not implemented yet; see plan.md Phase 8."
        )

    def confint(self, level: float = 0.95) -> np.ndarray:
        raise NotImplementedError(
            "LmRobResults.confint is not implemented yet; see plan.md Phase 8."
        )
