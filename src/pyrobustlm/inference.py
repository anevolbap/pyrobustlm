# SPDX-License-Identifier: GPL-3.0-or-later
"""Covariance estimators and Wald-type inference.

Phase 7 implements:

- ``vcov_avar1``  (R's ``.vcov.avar1``, default for KS2014)
- ``vcov_w``      (R's ``.vcov.w``, Koller & Stahel 2011)
- ``vcov_asymp``  (legacy MM)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

    from pyrobustlm.control import PsiFamily


def vcov_avar1(
    X: np.ndarray,
    weights: np.ndarray,
    psi_family: PsiFamily,
    psi_k: float | tuple[float, ...],
    sigma: float,
) -> np.ndarray:
    raise NotImplementedError("inference.vcov_avar1 — Phase 7")


def vcov_w(
    X: np.ndarray,
    weights: np.ndarray,
    psi_family: PsiFamily,
    psi_k: float | tuple[float, ...],
    sigma: float,
) -> np.ndarray:
    raise NotImplementedError("inference.vcov_w — Phase 7")


def vcov_asymp(
    X: np.ndarray,
    weights: np.ndarray,
    psi_family: PsiFamily,
    psi_k: float | tuple[float, ...],
    sigma: float,
) -> np.ndarray:
    raise NotImplementedError("inference.vcov_asymp — Phase 7")
