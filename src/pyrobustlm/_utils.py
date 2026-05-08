# SPDX-License-Identifier: GPL-3.0-or-later
"""Internal utilities: design-matrix prep, rank checks, dataset loaders.

Most helpers land in their consuming phase. This module currently only
provides a couple of trivial helpers.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


def is_factor_column(col: np.ndarray) -> bool:
    """Phase 5 helper. Returns True if ``col`` is a categorical encoding."""

    raise NotImplementedError("_utils.is_factor_column — Phase 5")


def split_design(
    X: np.ndarray,
    factor_cols: list[int],
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(X_cat, X_cont)`` split by column indices."""

    raise NotImplementedError("_utils.split_design — Phase 5")
