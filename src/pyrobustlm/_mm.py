# SPDX-License-Identifier: GPL-3.0-or-later
"""MM iteration on top of an initial S estimate.

Direct port of robustbase/src/lmrob.c::rwls (lmrob.c:1524). Plain IRWLS
with the L1-norm convergence test::

    d_beta <= epsilon * max(epsilon, ||beta_new||_1)

No step-halving (the R/C source doesn't do it either).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pyrobustlm._fast_s import _irwls_step


@dataclass
class MMResult:
    coef: NDArray[np.float64]
    converged: bool
    n_iter: int


def mm_iterate(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    beta_init: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: float | tuple[float, ...],
    max_it: int = 50,
    rel_tol: float = 1e-7,
) -> MMResult:
    """Iterate IRWLS with the efficiency-tuned psi until beta converges."""
    beta_cur = beta_init.copy().astype(np.float64)
    if sigma == 0.0:
        return MMResult(coef=beta_cur, converged=True, n_iter=0)

    converged = False
    it = 0
    for it in range(max_it):
        beta_new = _irwls_step(X, y, beta_cur, sigma, psi_family, psi_k)
        d_beta = float(np.sum(np.abs(beta_new - beta_cur)))
        norm1_new = float(np.sum(np.abs(beta_new)))
        beta_cur = beta_new
        if d_beta <= rel_tol * max(rel_tol, norm1_new):
            converged = True
            break
    return MMResult(coef=beta_cur, converged=converged, n_iter=it + 1)
