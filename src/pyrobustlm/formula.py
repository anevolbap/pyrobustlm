# SPDX-License-Identifier: GPL-3.0-or-later
"""Formula handling.

Wraps :mod:`formulaic` to produce design matrices that mirror R's
``model.matrix``. Phase 8.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd


def model_matrix(
    formula: str,
    data: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return ``(y, X, term_names)`` from a formula and a DataFrame."""

    raise NotImplementedError("formula.model_matrix — Phase 8")
