# SPDX-License-Identifier: GPL-3.0-or-later
"""Formula handling.

Wraps ``formulaic`` to produce design matrices that mirror R's
``model.matrix``. We also tag the categorical columns so the caller can
dispatch to the M-S estimator when needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd


@dataclass
class DesignMatrix:
    y: np.ndarray
    X: np.ndarray
    term_names: list[str]
    # Boolean mask: True where the column is from a categorical encoding.
    is_factor_col: np.ndarray
    # formulaic ModelSpec for the RHS, so :meth:`LmRobResults.predict` can
    # re-apply the same factor-encoding / interaction / I(...) transforms
    # to a new DataFrame.
    rhs_spec: object | None = None


def model_matrix(
    formula: str,
    data: pd.DataFrame,
) -> DesignMatrix:
    """Return a :class:`DesignMatrix` from a formula and a DataFrame.

    formulaic encodes factor levels as ``Var[T.level]``; we use that
    bracket pattern to flag categorical columns.
    """
    from formulaic import Formula

    parsed = Formula(formula)
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
    is_factor = np.array(
        ["[T." in c or c.startswith("C(") for c in term_names],
        dtype=bool,
    )
    return DesignMatrix(
        y=y,
        X=X,
        term_names=term_names,
        is_factor_col=is_factor,
        rhs_spec=X_df.model_spec,
    )


def apply_spec(rhs_spec: object, data: object) -> np.ndarray:
    """Re-apply a stored ``formulaic.ModelSpec`` to a new DataFrame.

    ``rhs_spec`` is typed as ``object`` because we do not want to depend
    on ``formulaic`` at type-check time, but at runtime it must be a
    ``formulaic.ModelSpec``. Returns the design matrix as a NumPy array.
    """
    get_mm = getattr(rhs_spec, "get_model_matrix", None)
    if get_mm is None:
        raise TypeError("rhs_spec is not a formulaic ModelSpec (no get_model_matrix)")
    new_mm = get_mm(data)
    return np.asarray(new_mm.to_numpy(), dtype=np.float64)
