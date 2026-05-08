# SPDX-License-Identifier: GPL-3.0-or-later
"""Top-level ``lmrob`` entry points.

Phase 8 will implement the full fit pipeline. This file currently exposes
the public names with stub bodies so callers can import them.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

    from pyrobustlm.control import Control
    from pyrobustlm.results import LmRobResults


def lmrob(
    formula: str,
    data: pd.DataFrame,
    control: Control | None = None,
    weights: np.ndarray | None = None,
    na_action: str = "drop",
    **kwargs: Any,
) -> LmRobResults:
    """Fit a robust MM linear regression.

    Phase 8 deliverable. Currently raises ``NotImplementedError``.
    """

    raise NotImplementedError("lmrob() is not implemented yet; see plan.md Phase 8.")


class LmRob:
    """scikit-learn-style estimator wrapper around :func:`lmrob`."""

    def __init__(self, control: Control | None = None) -> None:
        self.control = control

    def fit(self, X: np.ndarray, y: np.ndarray) -> LmRob:
        raise NotImplementedError("LmRob.fit is not implemented yet; see plan.md Phase 8.")

    def predict(self, X: np.ndarray) -> np.ndarray:
        raise NotImplementedError("LmRob.predict is not implemented yet; see plan.md Phase 8.")
