# SPDX-License-Identifier: GPL-3.0-or-later
"""Public interface to psi/chi/weight functions.

The numerical kernels live in :mod:`pyrobustlm._core._psi` (Cython).
This module is a thin Python wrapper that handles dispatch by family
name and broadcasts scalar/array inputs.

Implemented in Phase 2.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

    from pyrobustlm.control import PsiFamily


def rho(
    x: float | np.ndarray,
    family: PsiFamily,
    k: float | tuple[float, ...],
) -> float | np.ndarray:
    raise NotImplementedError("psi.rho — Phase 2")


def psi(
    x: float | np.ndarray,
    family: PsiFamily,
    k: float | tuple[float, ...],
) -> float | np.ndarray:
    raise NotImplementedError("psi.psi — Phase 2")


def psi_prime(
    x: float | np.ndarray,
    family: PsiFamily,
    k: float | tuple[float, ...],
) -> float | np.ndarray:
    raise NotImplementedError("psi.psi_prime — Phase 2")


def wgt(
    x: float | np.ndarray,
    family: PsiFamily,
    k: float | tuple[float, ...],
) -> float | np.ndarray:
    raise NotImplementedError("psi.wgt — Phase 2")


def Epsi2(family: PsiFamily, k: float | tuple[float, ...]) -> float:
    raise NotImplementedError("psi.Epsi2 — Phase 2")


def EDpsi(family: PsiFamily, k: float | tuple[float, ...]) -> float:
    raise NotImplementedError("psi.EDpsi — Phase 2")


def tuning_for_efficiency(family: PsiFamily, efficiency: float) -> tuple[float, ...]:
    """Return the tuning constants giving the requested asymptotic efficiency.

    Mirrors R's ``.psi.const`` / ``.psi.conv.cc`` machinery.
    Phase 2.
    """

    raise NotImplementedError("psi.tuning_for_efficiency — Phase 2")


def tuning_for_breakdown(family: PsiFamily, breakdown: float) -> tuple[float, ...]:
    raise NotImplementedError("psi.tuning_for_breakdown — Phase 2")
