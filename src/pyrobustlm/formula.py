# SPDX-License-Identifier: GPL-3.0-or-later
"""Formula handling.

Wraps ``formulaic`` to produce design matrices that mirror R's
``model.matrix``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd


def model_matrix(
    formula: str,
    data: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return ``(y, X, term_names)`` from a formula and a DataFrame."""
    from formulaic import Formula

    parsed = Formula(formula)
    # formulaic returns ModelMatrix(es) with .lhs/.rhs by default; we pass
    # .get_model_matrix(data) which returns a DataFrame for each side.
    mm = parsed.get_model_matrix(data)
    if not hasattr(mm, "lhs") or not hasattr(mm, "rhs"):
        raise ValueError(f"formula {formula!r} must be two-sided (have a y ~ ... LHS)")
    y_df = mm.lhs
    X_df = mm.rhs
    if y_df.shape[1] != 1:
        raise ValueError(f"formula LHS must produce one column; got {y_df.shape[1]}")
    y = np.asarray(y_df.iloc[:, 0].to_numpy(), dtype=np.float64)
    X = np.asarray(X_df.to_numpy(), dtype=np.float64)
    term_names = [str(c) for c in X_df.columns]
    return y, X, term_names
