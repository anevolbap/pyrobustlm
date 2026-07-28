# SPDX-License-Identifier: GPL-3.0-or-later
"""M-S estimator for designs with categorical predictors.

Direct port of robustbase's ``R_lmrob_M_S`` (src/lmrob.c:401), in four
phases:

1. **Orthogonalize** via L1 regression. ``y`` and each column of
   ``X_cont`` are L1-regressed on ``X_cat``; the residuals form the
   "orthogonal" working data.
2. **Subsample** ``nResample`` random ``p2``-subsets of the orthogonal
   ``(X_cont_orth, y_orth)``. Each yields a candidate coefficient pair;
   we keep the best (smallest M-scale) across resamples.
3. **Transform back** the orthogonal-space coefficients to the original
   design space.
4. **Descent**: alternate L1 (cat block) and weighted-LS (cont block)
   updates until convergence or K_m_s steps without improvement.

References: Maronna & Yohai (2000); robustbase ``src/lmrob.c::m_s_subsample``,
``m_s_descent``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import linprog

from pylmrob import psi as _psi
from pylmrob._fast_s import _draw_nonsingular_subset
from pylmrob.scale import m_scale


@dataclass
class MSResult:
    coef: NDArray[np.float64]  # full beta in original column order
    coef_cat: NDArray[np.float64]
    coef_cont: NDArray[np.float64]
    scale: float
    converged: bool
    n_iter: int


def _l1_fit(X: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    """Solve ``min_b sum |y - X b|`` via linprog.

    Reformulation: introduce slack ``u >= |y - X b|`` with two linear
    inequalities per observation. ``A_ub @ x <= b_ub`` form::

        [-X, -I] [b; u] <= -y         (u >= y - X b)
        [+X, -I] [b; u] <= +y         (u >= X b - y)
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


def _l1_fit_residuals(
    X: NDArray[np.float64], y: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """L1 fit returning ``(coef, y - X coef)``."""
    coef = _l1_fit(X, y)
    return coef, y - X @ coef


def _wls_fit(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    w: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Weighted-LS fit: ``min_b sum w_i (y_i - x_i b)^2`` via QR/lstsq."""
    sw = np.sqrt(np.maximum(w, 0.0))
    Xw = X * sw[:, None]
    yw = y * sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return np.ascontiguousarray(beta, dtype=np.float64)


def m_s_fit(
    X_cat: NDArray[np.float64],
    X_cont: NDArray[np.float64],
    y: NDArray[np.float64],
    psi_chi: str = "bisquare",
    k_chi: tuple[float, ...] = (1.54764,),
    b0: float = 0.5,
    nResample: int = 200,
    k_m_s: int = 20,
    max_k: int = 200,
    max_it: int = 50,
    rel_tol: float = 1e-7,
    scale_tol: float = 1e-10,
    zero_tol: float = 1e-10,
    mts: int = 1000,
    seed: int | np.random.Generator | None = None,
) -> MSResult:
    """Maronna-Yohai 2000 M-S estimator.

    Parameters mirror robustbase's ``lmrob.control`` for the M-S subset
    (``k.m_s``, ``max.k``, ``rel.tol``, ``scale.tol``, ``zero.tol``,
    ``mts``).
    """
    X_cat = np.ascontiguousarray(X_cat, dtype=np.float64)
    X_cont = np.ascontiguousarray(X_cont, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    n, p1 = X_cat.shape
    _, p2 = X_cont.shape
    if X_cont.shape[0] != n or y.shape[0] != n:
        raise ValueError("row counts of X_cat, X_cont, y must agree")
    if p1 == 0 or p2 == 0:
        raise ValueError("M-S requires both categorical and continuous columns")

    rng = np.random.default_rng(seed)

    # ------------------------------------------------------------------
    # Phase 1: Orthogonalize via L1
    # ------------------------------------------------------------------
    ot1, y_orth = _l1_fit_residuals(X_cat, y)
    oT2 = np.empty((p1, p2), dtype=np.float64)
    X_cont_orth = np.empty_like(X_cont)
    for j in range(p2):
        coef_j, X_cont_orth[:, j] = _l1_fit_residuals(X_cat, X_cont[:, j])
        oT2[:, j] = coef_j

    # ------------------------------------------------------------------
    # Phase 2: Subsample (in orth space). Track best by M-scale.
    # ------------------------------------------------------------------
    best_b1_orth = np.zeros(p1)
    best_b2 = np.zeros(p2)
    best_scale = float("inf")
    p_total = p1 + p2

    for _ in range(nResample):
        idx = _draw_nonsingular_subset(X_cont_orth, rng, p2, mts)
        if idx is None:
            continue
        try:
            t2 = np.linalg.solve(X_cont_orth[idx], y_orth[idx])
        except np.linalg.LinAlgError:
            continue

        partial = y_orth - X_cont_orth @ t2
        t1 = _l1_fit(X_cat, partial)
        res = partial - X_cat @ t1

        # Cheap "looks promising" gate (lmrob.c:2291): mean(rho(res / best_sc)) < bb
        sc_test = max(_mad(res), 1e-12) if best_scale == float("inf") else best_scale
        rho_vals = _psi.rho(res / sc_test, psi_chi, k_chi)
        if float(np.mean(rho_vals)) * n / max(n - p_total, 1) < b0:
            new_scale = m_scale(
                res,
                family=psi_chi,
                k=k_chi,
                b0=b0,
                max_iter=max_it,
                tol=scale_tol,
                p=p_total,
            )
            if new_scale > 0 and new_scale < best_scale:
                best_scale = new_scale
                best_b1_orth = t1.copy()
                best_b2 = t2.copy()
                if new_scale < zero_tol:
                    # perfect fit; stop early
                    break

    if best_scale == float("inf"):
        # No candidate cleared the gate; fall back to the last-iter L1 + a
        # full M-scale just so we return *something* plausible.
        coef_init = _l1_fit(np.column_stack([X_cat, X_cont]), y)
        b_cat = coef_init[:p1]
        b_cont = coef_init[p1:]
        sc = m_scale(
            y - X_cat @ b_cat - X_cont @ b_cont,
            family=psi_chi,
            k=k_chi,
            b0=b0,
            p=p_total,
        )
        return MSResult(
            coef=np.concatenate([b_cat, b_cont]),
            coef_cat=b_cat,
            coef_cont=b_cont,
            scale=float(sc),
            converged=False,
            n_iter=0,
        )

    # ------------------------------------------------------------------
    # Phase 3: Transform back
    #   b1_orig = ot1 + b1_orth - oT2 @ b2_orth
    # ------------------------------------------------------------------
    b1 = ot1 + best_b1_orth - oT2 @ best_b2
    b2 = best_b2.copy()
    sc = best_scale
    res = y - X_cat @ b1 - X_cont @ b2

    # ------------------------------------------------------------------
    # Phase 4: Descent. Alternate L1 (cat) and WLS (cont) until k_m_s
    # consecutive non-improvements or max_k iterations.
    # ------------------------------------------------------------------
    no_improve = 0
    converged = False
    iter_used = 0
    t1, t2 = b1.copy(), b2.copy()
    res2 = res.copy()
    for k in range(max_k):
        iter_used = k + 1
        # WLS update of t2
        y_tilde = y - X_cat @ t1
        if sc == 0.0:
            break
        w = _psi.wgt(res2 / sc, psi_chi, k_chi)
        t2_new = _wls_fit(X_cont, y_tilde, w)
        res2 = y - X_cont @ t2_new

        # L1 update of t1 against res2
        t1_new, res2 = _l1_fit_residuals(X_cat, res2)

        # New scale on full residuals
        sc_new = m_scale(
            res2,
            family=psi_chi,
            k=k_chi,
            b0=b0,
            max_iter=max_it,
            tol=scale_tol,
            p=p_total,
            init_scale=sc,
        )

        # Convergence: relative change in the joint coef vector
        delta = np.sqrt(float(np.sum((t1_new - t1) ** 2)) + float(np.sum((t2_new - t2) ** 2)))
        nrm_b = np.sqrt(float(np.sum(t1_new**2)) + float(np.sum(t2_new**2)))
        if delta < rel_tol * max(rel_tol, nrm_b):
            converged = True

        t1, t2 = t1_new, t2_new
        if sc_new < sc and sc_new > 0:
            b1, b2, sc = t1.copy(), t2.copy(), sc_new
            no_improve = 0
        else:
            no_improve += 1

        if converged:
            break
        if no_improve >= k_m_s:
            break

    return MSResult(
        coef=np.concatenate([b1, b2]),
        coef_cat=b1,
        coef_cont=b2,
        scale=float(sc),
        converged=converged,
        n_iter=iter_used,
    )


def _mad(x: NDArray[np.float64]) -> float:
    med = float(np.median(x))
    return float(1.4826 * np.median(np.abs(x - med)))
