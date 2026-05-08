# SPDX-License-Identifier: GPL-3.0-or-later
"""M-S estimator for designs with categorical predictors.

Simplified Maronna & Yohai (2000) alternating descent. The full robustbase
implementation lives in C (``src/lmrob.c::m_s_descent``); we mirror its
high-level structure in NumPy:

1. Split design ``X = [X_cat | X_cont]``.
2. Initialise ``beta_cat`` by L1 fit of ``y`` on ``X_cat`` (no continuous
   part). L1 is solved via SciPy's ``linprog``.
3. Iterate ``k_m_s`` times:

   a. ``beta_cont`` <- S fit of ``y - X_cat @ beta_cat`` on ``X_cont``.
   b. ``beta_cat``  <- L1 fit of ``y - X_cont @ beta_cont`` on ``X_cat``.

4. Stitch ``beta = [beta_cat, beta_cont]`` and return it together with the
   final M-scale of the full residual vector.

This is good enough for designs where pure-S resampling fails (random
p-subsets are systematically singular). It is **not** bit-identical with
robustbase's M-S; the C implementation does sophisticated subsample +
descent with multiple restarts. Documented in ``docs/numerical-notes.md``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linprog

from pyrobustlm._fast_s import FastSConfig, fast_s
from pyrobustlm.scale import m_scale


@dataclass
class MSResult:
    coef: NDArray[np.float64]  # full-design beta (cat slots first, then cont)
    coef_cat: NDArray[np.float64]
    coef_cont: NDArray[np.float64]
    scale: float
    converged: bool
    n_iter: int


def _l1_fit(X: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    """Solve ``min_b sum |y - X b|`` via linear programming.

    Reformulation: introduce slack ``u >= |y - X b|``. The constraints
    ``u >= y - X b`` and ``u >= X b - y`` give two linear inequalities per
    observation. In linprog's ``A_ub @ x <= b_ub`` form:

        [-X, -I]  [b; u] <= -y         (u >= y - X b)
        [+X, -I]  [b; u] <= +y         (u >= X b - y)
    """
    n, p = X.shape
    c = np.concatenate([np.zeros(p), np.ones(n)])
    A_ub = np.vstack(
        [
            np.hstack([-X, -np.eye(n)]),
            np.hstack([+X, -np.eye(n)]),
        ]
    )
    b_ub = np.concatenate([-y, +y])
    bounds = [(None, None)] * p + [(0.0, None)] * n
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(f"L1 fit failed: {res.message}")
    return res.x[:p].astype(np.float64)


def m_s_fit(
    X_cat: NDArray[np.float64],
    X_cont: NDArray[np.float64],
    y: NDArray[np.float64],
    psi_chi: str = "bisquare",
    k_chi: tuple[float, ...] = (1.547645,),
    b0: float = 0.5,
    k_m_s: int = 20,
    nResample: int = 200,
    max_it: int = 50,
    rel_tol: float = 1e-7,
    seed: int | np.random.Generator | None = None,
) -> MSResult:
    """Alternating L1/S estimator for designs with factor + continuous parts.

    Parameters
    ----------
    X_cat, X_cont :
        Categorical and continuous parts of the design matrix.
    y :
        Response vector.
    psi_chi, k_chi, b0 :
        Chi family / tuning constants / consistency constant for the S step.
    k_m_s :
        Number of alternating descent iterations.
    nResample, max_it, rel_tol, seed :
        Forwarded to :func:`pyrobustlm._fast_s.fast_s` for the inner S step.
    """
    X_cat = np.ascontiguousarray(X_cat, dtype=np.float64)
    X_cont = np.ascontiguousarray(X_cont, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    if X_cat.ndim != 2 or X_cont.ndim != 2:
        raise ValueError("X_cat and X_cont must be 2-D")
    if X_cat.shape[0] != X_cont.shape[0] or X_cat.shape[0] != y.shape[0]:
        raise ValueError("row counts of X_cat, X_cont, y must agree")

    p_cat = X_cat.shape[1]
    p_cont = X_cont.shape[1]

    # 1. L1 init for beta_cat (no continuous part yet)
    beta_cat = _l1_fit(X_cat, y) if p_cat > 0 else np.empty(0)
    beta_cont = np.zeros(p_cont, dtype=np.float64)

    cfg = FastSConfig(
        psi_chi=psi_chi,
        k_chi=k_chi,
        b0=b0,
        nResample=nResample,
        max_it=max_it,
        refine_tol=rel_tol,
    )

    converged = False
    last_beta_cat = beta_cat.copy()
    last_beta_cont = beta_cont.copy()
    it = 0
    for it in range(k_m_s):
        # (a) S fit of partial residual on X_cont
        partial_resid = y - (X_cat @ beta_cat if p_cat > 0 else 0.0)
        if p_cont > 0:
            s_res = fast_s(X_cont, partial_resid, cfg=cfg, seed=seed)
            beta_cont = s_res.coef
        # (b) L1 fit of partial residual on X_cat
        partial_resid = y - (X_cont @ beta_cont if p_cont > 0 else 0.0)
        if p_cat > 0:
            beta_cat = _l1_fit(X_cat, partial_resid)

        # Convergence: changes in both blocks
        d_cat = float(np.linalg.norm(beta_cat - last_beta_cat))
        d_cont = float(np.linalg.norm(beta_cont - last_beta_cont))
        denom = max(
            float(np.linalg.norm(beta_cat)) + float(np.linalg.norm(beta_cont)),
            1e-300,
        )
        if (d_cat + d_cont) / denom < rel_tol:
            converged = True
            break
        last_beta_cat = beta_cat.copy()
        last_beta_cont = beta_cont.copy()

    # Final scale on the full residual.
    full_resid = (
        y - (X_cat @ beta_cat if p_cat > 0 else 0.0) - (X_cont @ beta_cont if p_cont > 0 else 0.0)
    )
    sigma = m_scale(
        full_resid,
        family=psi_chi,
        k=k_chi,
        b0=b0,
        p=p_cat + p_cont,
    )

    coef = np.concatenate([beta_cat, beta_cont])
    return MSResult(
        coef=coef,
        coef_cat=beta_cat,
        coef_cont=beta_cont,
        scale=sigma,
        converged=converged,
        n_iter=it + 1,
    )
