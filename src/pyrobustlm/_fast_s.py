# SPDX-License-Identifier: GPL-3.0-or-later
"""Fast-S regression estimator (Salibian-Barrera & Yohai 2006).

Pure-NumPy implementation. Deferred Cython/parallel acceleration to a later
phase. Algorithm summary:

    1. Repeat ``nResample`` times:
        a. draw a random p-subset (with non-singular check, retry up to mts)
        b. solve X[I] beta = y[I]
        c. K-step refinement: alternate (M-scale -> IRWLS) ``k_fast_s`` times
        d. keep the candidate in a best-of-``best_r`` heap by sigma
    2. For each survivor, refine to convergence.
    3. Return the best.

This module is consumed by ``pyrobustlm.lmrob`` for the initial S estimate
when ``init="S"``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from pyrobustlm import psi as _psi
from pyrobustlm.scale import _mad, m_scale


@dataclass(frozen=True)
class FastSConfig:
    psi_chi: str = "bisquare"
    k_chi: tuple[float, ...] = (1.547645,)  # 50% bdp default for bisquare
    b0: float = 0.5
    nResample: int = 500
    k_fast_s: int = 1  # number of refinement steps inside the resampling loop
    best_r: int = 2  # number of candidates to refine to convergence
    max_it: int = 50  # max IRWLS iterations during full refinement
    refine_tol: float = 1e-7
    scale_tol: float = 1e-10
    max_iter_scale: int = 200
    mts: int = 1000  # max attempts when drawing a non-singular p-subset


@dataclass
class FastSResult:
    coef: NDArray[np.float64]
    scale: float
    converged: bool
    n_iter: int
    n_candidates_kept: int


def _draw_nonsingular_subset(
    X: NDArray[np.float64],
    rng: np.random.Generator,
    p: int,
    mts: int,
    rcond: float = 1e-10,
) -> NDArray[np.intp] | None:
    """Sample a p-subset of rows whose submatrix is non-singular."""
    n = X.shape[0]
    for _ in range(mts):
        idx = rng.choice(n, size=p, replace=False)
        sub = X[idx]
        # Cheap rank check via the smallest singular value relative to the largest.
        try:
            sv = np.linalg.svd(sub, compute_uv=False)
        except np.linalg.LinAlgError:
            continue
        if sv[-1] > rcond * sv[0]:
            return idx
    return None


def _irwls_step(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    beta: NDArray[np.float64],
    sigma: float,
    psi_family: str,
    psi_k: tuple[float, ...] | float,
) -> NDArray[np.float64]:
    """One IRWLS step: solve diag(sqrt(w)) X * beta_new = diag(sqrt(w)) y.

    Weights ``w = psi(r/sigma) / (r/sigma)`` for ``r = y - X beta``. Where
    r/sigma == 0 we set w := 1 (since psi'(0) = 1 for all the families used).
    """
    r = y - X @ beta
    if sigma == 0.0:
        return beta.copy()
    z = r / sigma
    w = _psi.wgt(z, psi_family, psi_k)
    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    yw = y * sw
    beta_new, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return beta_new


def _refine_to_convergence(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    beta: NDArray[np.float64],
    sigma: float,
    cfg: FastSConfig,
) -> tuple[NDArray[np.float64], float, bool, int]:
    """Iterate (M-scale, IRWLS) until ||beta_{k+1} - beta_k|| / ||beta_k|| < tol."""
    p = X.shape[1]
    beta_cur = beta.copy()
    s_cur = sigma
    converged = False
    it = 0
    for it in range(cfg.max_it):
        s_new = m_scale(
            y - X @ beta_cur,
            family=cfg.psi_chi,
            k=cfg.k_chi,
            b0=cfg.b0,
            max_iter=cfg.max_iter_scale,
            tol=cfg.scale_tol,
            init_scale=s_cur,
            p=p,
        )
        if s_new == 0.0:
            return beta_cur, 0.0, True, it
        beta_new = _irwls_step(X, y, beta_cur, s_new, cfg.psi_chi, cfg.k_chi)
        delta = np.linalg.norm(beta_new - beta_cur)
        denom = max(np.linalg.norm(beta_cur), 1e-300)
        s_cur = s_new
        beta_cur = beta_new
        if delta / denom < cfg.refine_tol:
            converged = True
            break
    return beta_cur, s_cur, converged, it + 1


def fast_s(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    cfg: FastSConfig | None = None,
    seed: int | np.random.Generator | None = None,
) -> FastSResult:
    """Fast-S regression estimator.

    Parameters
    ----------
    X :
        Design matrix, shape (n, p).
    y :
        Response vector, length n.
    cfg :
        Algorithm configuration; defaults to ``FastSConfig()``.
    seed :
        Seed for the per-fit random generator. Reproducible across runs.

    Returns
    -------
    FastSResult

    Notes
    -----
    Bit-identical reproducibility with R's resampling is impossible
    (R uses Mersenne Twister; we use PCG64). Coefficients agree to within
    a basin of attraction tolerance documented in plan.md §5.1.
    """
    cfg = cfg or FastSConfig()
    rng = np.random.default_rng(seed)
    X = np.ascontiguousarray(X, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    n, p = X.shape
    if y.shape != (n,):
        raise ValueError(f"y shape {y.shape} incompatible with X shape {X.shape}")
    if n <= p:
        raise ValueError(f"fast_s requires n > p; got n={n}, p={p}")

    # ------------------------------------------------------------------
    # Resampling loop
    # ------------------------------------------------------------------
    best_betas: list[NDArray[np.float64]] = []
    best_scales: list[float] = []

    for _ in range(cfg.nResample):
        idx = _draw_nonsingular_subset(X, rng, p, cfg.mts)
        if idx is None:
            continue
        try:
            beta_t = np.linalg.solve(X[idx], y[idx])
        except np.linalg.LinAlgError:
            continue

        # K-step refinement
        r = y - X @ beta_t
        s_t = _mad(r)
        if s_t == 0.0:
            return FastSResult(
                coef=beta_t, scale=0.0, converged=True, n_iter=0, n_candidates_kept=1
            )
        for _kk in range(cfg.k_fast_s):
            s_t = m_scale(
                y - X @ beta_t,
                family=cfg.psi_chi,
                k=cfg.k_chi,
                b0=cfg.b0,
                max_iter=cfg.max_iter_scale,
                tol=cfg.scale_tol,
                init_scale=s_t,
                p=p,
            )
            if s_t == 0.0:
                return FastSResult(
                    coef=beta_t, scale=0.0, converged=True, n_iter=0, n_candidates_kept=1
                )
            beta_t = _irwls_step(X, y, beta_t, s_t, cfg.psi_chi, cfg.k_chi)

        # Final scale at this candidate
        s_t = m_scale(
            y - X @ beta_t,
            family=cfg.psi_chi,
            k=cfg.k_chi,
            b0=cfg.b0,
            max_iter=cfg.max_iter_scale,
            tol=cfg.scale_tol,
            init_scale=s_t,
            p=p,
        )

        # Insert into the bounded heap (small best_r; linear scan is fine).
        if len(best_scales) < cfg.best_r:
            best_betas.append(beta_t)
            best_scales.append(s_t)
        else:
            worst_i = int(np.argmax(best_scales))
            if s_t < best_scales[worst_i]:
                best_betas[worst_i] = beta_t
                best_scales[worst_i] = s_t

    if not best_scales:
        raise RuntimeError("fast_s: no non-singular subsamples were found")

    # ------------------------------------------------------------------
    # Refine survivors to convergence
    # ------------------------------------------------------------------
    refined_beta: list[NDArray[np.float64]] = []
    refined_sigma: list[float] = []
    refined_iter: list[int] = []
    for beta_t, s_t in zip(best_betas, best_scales, strict=True):
        bb, ss, _conv, it = _refine_to_convergence(X, y, beta_t, s_t, cfg)
        refined_beta.append(bb)
        refined_sigma.append(ss)
        refined_iter.append(it)

    best_idx = int(np.argmin(refined_sigma))
    return FastSResult(
        coef=refined_beta[best_idx],
        scale=refined_sigma[best_idx],
        converged=True,
        n_iter=refined_iter[best_idx],
        n_candidates_kept=len(best_scales),
    )
