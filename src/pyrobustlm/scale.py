# SPDX-License-Identifier: GPL-3.0-or-later
"""Robust M-scale estimator.

Phase 3 wraps :mod:`pyrobustlm._core._scale` and exposes a NumPy-friendly
``m_scale`` function that mirrors R's ``robustbase::lmrob.mscale``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

    from pyrobustlm.control import PsiFamily


def m_scale(
    r: np.ndarray,
    family: PsiFamily = "bisquare",
    k: float | tuple[float, ...] | None = None,
    b0: float = 0.5,
    max_iter: int = 200,
    tol: float = 1e-10,
    init_scale: float | None = None,
) -> float:
    raise NotImplementedError("scale.m_scale — Phase 3")
