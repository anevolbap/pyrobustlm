# SPDX-License-Identifier: GPL-3.0-or-later
"""M-S estimator for factor designs (Maronna-Yohai 2000).

Status: stub. ``init="M-S"`` is not yet implemented in pyrobustlm. For
designs with categorical predictors that produce frequently singular
random subsamples, set ``Control(mts=...)`` higher (default 1000 already
covers the reference corpus).

Tracking: plan.md Phase 5.
"""

from __future__ import annotations

import numpy as np

from pyrobustlm.control import Control


def m_s_fit(
    X_cat: np.ndarray,
    X_cont: np.ndarray,
    y: np.ndarray,
    control: Control,
) -> tuple[np.ndarray, float]:
    raise NotImplementedError(
        "init='M-S' is not yet implemented (plan.md Phase 5). "
        "Use init='S' (the default) and increase Control(mts=...) "
        "if the resampling cannot find non-singular subsets."
    )
