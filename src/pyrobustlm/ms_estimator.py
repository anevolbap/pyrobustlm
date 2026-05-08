# SPDX-License-Identifier: GPL-3.0-or-later
"""M-S estimator for designs with categorical predictors.

Maronna & Yohai (2000). Phase 5.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

    from pyrobustlm.control import Control


def m_s_fit(
    X_cat: np.ndarray,
    X_cont: np.ndarray,
    y: np.ndarray,
    control: Control,
) -> tuple[np.ndarray, float]:
    raise NotImplementedError("ms_estimator.m_s_fit — Phase 5")
