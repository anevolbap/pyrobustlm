# SPDX-License-Identifier: GPL-3.0-or-later
"""L1 (least absolute deviation) regression for the ``init="L1"`` path.

Used as an initial estimate for the MM-step when fast-S resampling is
either unstable or undesirable. L1 is high breakdown (50%) but lower
Gaussian efficiency than M-S or S; for that reason it's not the
default initial estimator. R's ``lmrob`` exposes it as ``init="L1"``.

The L1 problem ``min_beta sum |y - X beta|`` reformulates as an LP:

    min_{beta, u}  sum u_i
    s.t.           u_i >= y_i - X[i] beta
                   u_i >= -(y_i - X[i] beta)
                   u_i >= 0

We solve via ``scipy.optimize.linprog`` (HiGHS backend).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linprog


@dataclass
class L1Result:
    coef: NDArray[np.float64]
    scale: float
    converged: bool


def l1_fit(X: NDArray[np.float64], y: NDArray[np.float64]) -> L1Result:
    """L1 least-absolute-deviation regression.

    Returns a ``L1Result`` with ``coef`` (length p), ``scale``
    (median absolute deviation of residuals around zero, the M-scale
    convention robustbase uses for the L1 init), and a ``converged``
    flag from the LP solver.
    """
    X = np.ascontiguousarray(X, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    n, p = X.shape

    # Variables ordered (beta[0..p-1], u[0..n-1]).
    # Objective: minimise sum(u_i) -> c = [0]*p + [1]*n
    c = np.concatenate([np.zeros(p), np.ones(n)])

    # Constraints A_ub x <= b_ub for x = (beta, u):
    #   X beta - u <= y           ->   A1 = [X | -I_n], b1 = y
    #   -X beta - u <= -y         ->   A2 = [-X | -I_n], b2 = -y
    eye_n = -np.eye(n)
    A_ub = np.vstack(
        [
            np.hstack([X, eye_n]),
            np.hstack([-X, eye_n]),
        ]
    )
    b_ub = np.concatenate([y, -y])

    # Bounds: beta unbounded, u >= 0
    bounds = [(None, None)] * p + [(0.0, None)] * n

    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(f"L1 linprog failed: {res.message}")

    beta = np.asarray(res.x[:p], dtype=np.float64)
    resid = y - X @ beta
    # Scale via MAD around zero (matches robustbase's L1 path; the M-scale
    # step refines this further).
    sigma = float(np.median(np.abs(resid))) / 0.6744897501960817
    return L1Result(coef=beta, scale=sigma, converged=True)
