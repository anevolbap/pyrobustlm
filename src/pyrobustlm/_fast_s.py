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

import os
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from numpy.random import SeedSequence
from numpy.typing import NDArray

from pyrobustlm import _psifuns as _pf
from pyrobustlm.scale import _cython_wgt, _mad, m_scale

# Generic Cython kernel signature. Returns ``(scale, status)`` where
# ``status`` is 0 (ok), 1 (singular subset), 2 (exact fit), 3 (LAPACK error).
_CyIterFn = Callable[
    [
        NDArray[np.float64],  # X
        NDArray[np.float64],  # y
        NDArray[np.int64],  # idx
        int,  # family id
        NDArray[np.float64],  # tuning (len 3)
        float,  # b0
        int,  # k_fast_s
        int,  # max_iter_scale
        float,  # scale_tol
        NDArray[np.float64],  # beta_out
    ],
    tuple[float, int],
]


# Maps psi-family name to the Cython family id. Keep in sync with the
# enum in src/pyrobustlm/_core/_fast_s.pyx.
_FAMILY_IDS: dict[str, int] = {
    "bisquare": 0,
    "biweight": 0,
    "hampel": 1,
    "optimal": 2,
    "lqq": 3,
    "ggw": 4,
}


def _try_import_cy_iter() -> _CyIterFn | None:
    try:
        import importlib

        mod = importlib.import_module("pyrobustlm._core._fast_s")
        return mod.cy_resample_iter  # type: ignore[attr-defined,no-any-return]
    except Exception:
        return None


def _try_import_cy_refine() -> Callable[..., tuple[float, int, int, int]] | None:
    try:
        import importlib

        mod = importlib.import_module("pyrobustlm._core._fast_s")
        return mod.cy_refine_to_convergence  # type: ignore[attr-defined,no-any-return]
    except Exception:
        return None


def _try_import_cy_draw_and_iter() -> Callable[..., tuple[float, int]] | None:
    try:
        import importlib

        mod = importlib.import_module("pyrobustlm._core._fast_s")
        return mod.cy_draw_and_iter  # type: ignore[attr-defined,no-any-return]
    except Exception:
        return None


def _try_import_cy_lmrob_fast_s() -> Callable[..., tuple[float, int, int, int]] | None:
    try:
        import importlib

        mod = importlib.import_module("pyrobustlm._core._lmrob")
        return mod.cy_lmrob_fast_s  # type: ignore[attr-defined,no-any-return]
    except Exception:
        return None


_CY_ITER: _CyIterFn | None = _try_import_cy_iter()
_CY_REFINE: Callable[..., tuple[float, int, int, int]] | None = _try_import_cy_refine()
_CY_DRAW_AND_ITER: Callable[..., tuple[float, int]] | None = _try_import_cy_draw_and_iter()
_CY_LMROB_FAST_S: Callable[..., tuple[float, int, int, int]] | None = _try_import_cy_lmrob_fast_s()


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
    # Parallelism:
    #   n_workers=1: serial (default; bit-identical with pre-parallel releases).
    #   n_workers=0: auto = min(os.cpu_count(), max(1, nResample // 32)).
    #   n_workers>1: that many worker threads.
    # Each worker draws from a SeedSequence-spawned PCG64, so results are
    # deterministic for a given (seed, n_workers, nResample).
    n_workers: int = 1
    # When True, route subset draws through the Cython BitGenerator path
    # (Floyd's combination algorithm). Faster at small n but produces a
    # different sequence of subsets than ``np.random.Generator.choice``,
    # so the basin of attraction can shift.
    fast_rng: bool = False
    # Bisquare-only end-to-end Cython fast-S engine
    # (``pyrobustlm._core._lmrob``). When True, the entire resampling +
    # survivor refinement runs in one nogil C block with all workspace
    # pre-allocated. Currently bisquare only; other families fall back.
    # Off by default while stages 2-6 of the monolithic port are in
    # progress.
    engine_c: bool = False


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

    w_cy = _cython_wgt(z, psi_family, psi_k)
    if w_cy is not None:
        w = w_cy
    else:
        wgt_fn = _pf._dispatch(psi_family, "wgt")
        w = wgt_fn(z, np.asarray(psi_k, dtype=np.float64))

    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    yw = y * sw
    beta_new, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    return np.ascontiguousarray(beta_new, dtype=np.float64)


def _resample_chunk_cython(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    cfg: FastSConfig,
    rng: np.random.Generator,
    n_iter: int,
) -> _ChunkResult:
    """Cython fast path. Dispatches on the psi family.

    Prefers ``cy_draw_and_iter`` (subset draw + iter, all in nogil C using
    numpy's BitGenerator capsule) when available; falls back to
    ``cy_resample_iter`` with a Python-side subset draw.
    """
    cy_iter = _CY_ITER
    cy_draw = _CY_DRAW_AND_ITER
    assert cy_iter is not None  # caller already checked
    p = X.shape[1]
    family_id = _FAMILY_IDS[cfg.psi_chi]
    tuning = np.zeros(3, dtype=np.float64)
    for i, v in enumerate(cfg.k_chi[:3]):
        tuning[i] = float(v)
    best_betas: list[NDArray[np.float64]] = []
    best_scales: list[float] = []

    beta_buf = np.empty(p, dtype=np.float64)

    if cy_draw is not None and cfg.fast_rng:
        bg_capsule = rng.bit_generator.capsule
        for _ in range(n_iter):
            scale, status = cy_draw(
                X,
                y,
                bg_capsule,
                cfg.mts,
                family_id,
                tuning,
                cfg.b0,
                cfg.k_fast_s,
                cfg.max_iter_scale,
                cfg.scale_tol,
                beta_buf,
            )
            if status == 1 or status == 3:
                continue
            if status == 2:
                return _ChunkResult(
                    [],
                    [],
                    FastSResult(
                        coef=beta_buf.copy(),
                        scale=0.0,
                        converged=True,
                        n_iter=0,
                        n_candidates_kept=1,
                    ),
                )
            beta_t = beta_buf.copy()
            s_t = float(scale)
            if len(best_scales) < cfg.best_r:
                best_betas.append(beta_t)
                best_scales.append(s_t)
            else:
                worst_i = int(np.argmax(best_scales))
                if s_t < best_scales[worst_i]:
                    best_betas[worst_i] = beta_t
                    best_scales[worst_i] = s_t
        return _ChunkResult(best_betas, best_scales, None)

    for _ in range(n_iter):
        idx = _draw_nonsingular_subset(X, rng, p, cfg.mts)
        if idx is None:
            continue
        idx_long = np.asarray(idx, dtype=np.int64)
        scale, status = cy_iter(
            X,
            y,
            idx_long,
            family_id,
            tuning,
            cfg.b0,
            cfg.k_fast_s,
            cfg.max_iter_scale,
            cfg.scale_tol,
            beta_buf,
        )
        if status == 1 or status == 3:
            # Singular subset or LAPACK error.
            continue
        if status == 2:
            # Exact fit: short-circuit.
            return _ChunkResult(
                [],
                [],
                FastSResult(
                    coef=beta_buf.copy(),
                    scale=0.0,
                    converged=True,
                    n_iter=0,
                    n_candidates_kept=1,
                ),
            )
        beta_t = beta_buf.copy()
        s_t = float(scale)
        if len(best_scales) < cfg.best_r:
            best_betas.append(beta_t)
            best_scales.append(s_t)
        else:
            worst_i = int(np.argmax(best_scales))
            if s_t < best_scales[worst_i]:
                best_betas[worst_i] = beta_t
                best_scales[worst_i] = s_t

    return _ChunkResult(best_betas, best_scales, None)


@dataclass
class _ChunkResult:
    best_betas: list[NDArray[np.float64]]
    best_scales: list[float]
    early_exact_fit: FastSResult | None


def _resolve_n_workers(req: int, n_iter: int) -> int:
    """Map ``cfg.n_workers`` to a concrete worker count.

    Special values: ``1`` = serial; ``0`` = auto.

    Auto picks ``min(os.cpu_count(), max(2, n_iter // 64))``. The cap
    on ``n_iter`` keeps per-worker chunks above ~64 candidates so the
    per-task overhead amortises. The caller is responsible for only
    invoking auto mode when the per-iteration BLAS cost is large enough
    to dominate Python/GIL overhead (see ``fast_s`` for the cost-based
    gate).
    """
    if req == 1:
        return 1
    if req == 0:
        cores = os.cpu_count() or 1
        return max(1, min(cores, max(2, n_iter // 64)))
    return max(1, int(req))


def _auto_use_threads(n: int, p: int, n_iter: int) -> bool:
    """Heuristic for when threading is likely to be a net win.

    For small problems the GIL-locked Python work between BLAS calls
    dominates and threading hurts. The break-even on this hardware
    (16-core Linux x86_64, OpenBLAS) is roughly when ``n * p^2 >= 1e6``
    *and* ``n_iter >= 250``. Tune as needed; users can always force
    ``n_workers`` explicitly.
    """
    # Calibrated on a 16-core OpenBLAS x86_64 box: threading goes from
    # neutral (~1.0x) at n*p^2 ~= 5,000 to clearly profitable (>= 2x) at
    # n*p^2 >= 100,000. We pick a threshold near the neutral point so
    # auto-mode reliably helps anywhere with non-trivial computation.
    return (n * p * p >= 10_000) and (n_iter >= 250)


def _to_seed_sequence(seed: int | np.random.Generator | None) -> SeedSequence:
    """Turn the user-facing ``seed`` argument into a ``SeedSequence``.

    For an existing ``Generator`` we draw 32 bits from its state to seed a
    fresh sequence; that keeps determinism while letting parallel workers
    each get their own PCG64.
    """
    if isinstance(seed, np.random.Generator):
        # Pull a 64-bit seed deterministically out of the generator.
        s = int(seed.integers(0, 2**63 - 1))
        return SeedSequence(s)
    return SeedSequence(seed)


def _split_iters(total: int, k: int) -> list[int]:
    """Split ``total`` iterations across ``k`` workers as evenly as possible."""
    base, rem = divmod(total, k)
    return [base + (1 if i < rem else 0) for i in range(k)]


def _resample_chunk(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    cfg: FastSConfig,
    seed_seq: SeedSequence,
    n_iter: int,
) -> _ChunkResult:
    """Run ``n_iter`` resampling iterations and return the local best-of-best_r heap.

    Each call has its own PCG64; the return value is reduced together with
    other chunk outputs in :func:`fast_s`.
    """
    rng = np.random.default_rng(seed_seq)
    p = X.shape[1]
    best_betas: list[NDArray[np.float64]] = []
    best_scales: list[float] = []

    # Cython fast-path. Keeps the per-iteration body in nogil C, which is
    # what makes the thread pool actually scale on small-n problems. Covers
    # bisquare/hampel/optimal/lqq; ggw still goes through the NumPy path
    # below.
    if _CY_ITER is not None and cfg.psi_chi in _FAMILY_IDS:
        return _resample_chunk_cython(X, y, cfg, rng, n_iter)

    for _ in range(n_iter):
        idx = _draw_nonsingular_subset(X, rng, p, cfg.mts)
        if idx is None:
            continue
        try:
            beta_t = np.linalg.solve(X[idx], y[idx]).astype(np.float64, copy=False)
        except np.linalg.LinAlgError:
            continue

        r = y - X @ beta_t
        s_t = _mad(r)
        if s_t == 0.0:
            exact = FastSResult(
                coef=beta_t,
                scale=0.0,
                converged=True,
                n_iter=0,
                n_candidates_kept=1,
            )
            return _ChunkResult([], [], exact)

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
                exact = FastSResult(
                    coef=beta_t,
                    scale=0.0,
                    converged=True,
                    n_iter=0,
                    n_candidates_kept=1,
                )
                return _ChunkResult([], [], exact)
            beta_t = _irwls_step(X, y, beta_t, s_t, cfg.psi_chi, cfg.k_chi)

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

        if len(best_scales) < cfg.best_r:
            best_betas.append(beta_t)
            best_scales.append(s_t)
        else:
            worst_i = int(np.argmax(best_scales))
            if s_t < best_scales[worst_i]:
                best_betas[worst_i] = beta_t
                best_scales[worst_i] = s_t

    return _ChunkResult(best_betas, best_scales, None)


def _refine_to_convergence(
    X: NDArray[np.float64],
    y: NDArray[np.float64],
    beta: NDArray[np.float64],
    sigma: float,
    cfg: FastSConfig,
) -> tuple[NDArray[np.float64], float, bool, int]:
    """Iterate (M-scale, IRWLS) until ||beta_{k+1} - beta_k|| / ||beta_k|| < tol."""
    # Cython fast path for the supported families.
    if _CY_REFINE is not None and cfg.psi_chi in _FAMILY_IDS:
        family_id = _FAMILY_IDS[cfg.psi_chi]
        tuning = np.zeros(3, dtype=np.float64)
        for i, v in enumerate(cfg.k_chi[:3]):
            tuning[i] = float(v)
        beta_cur = np.ascontiguousarray(beta, dtype=np.float64).copy()
        scale, converged_int, n_iter, _status = _CY_REFINE(
            X,
            y,
            beta_cur,
            float(sigma),
            family_id,
            tuning,
            cfg.b0,
            cfg.max_it,
            cfg.refine_tol,
            cfg.max_iter_scale,
            cfg.scale_tol,
        )
        return beta_cur, float(scale), bool(converged_int), int(n_iter)

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
    X = np.ascontiguousarray(X, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    n, p = X.shape
    if y.shape != (n,):
        raise ValueError(f"y shape {y.shape} incompatible with X shape {X.shape}")
    if n <= p:
        raise ValueError(f"fast_s requires n > p; got n={n}, p={p}")

    # ------------------------------------------------------------------
    # Monolithic Cython engine (opt-in, all lmrob psi families)
    # ------------------------------------------------------------------
    if cfg.engine_c and _CY_LMROB_FAST_S is not None and cfg.psi_chi in _FAMILY_IDS:
        rng_e = np.random.default_rng(seed)
        beta_out = np.empty(p, dtype=np.float64)
        tuning = np.zeros(3, dtype=np.float64)
        for i, v in enumerate(cfg.k_chi[:3]):
            tuning[i] = float(v)
        family_id = _FAMILY_IDS[cfg.psi_chi]
        # Optional OpenMP parallelism via per-thread bitgens.
        ec_workers = _resolve_n_workers(cfg.n_workers, cfg.nResample)
        if ec_workers > 1:
            seedseq = _to_seed_sequence(seed)
            child_seeds = seedseq.spawn(ec_workers)
            capsules = [np.random.default_rng(s).bit_generator.capsule for s in child_seeds]
        else:
            capsules = None
        scale, status, n_iter, converged = _CY_LMROB_FAST_S(
            X,
            y,
            rng_e.bit_generator.capsule,
            family_id,
            tuning,
            cfg.b0,
            cfg.nResample,
            cfg.mts,
            cfg.k_fast_s,
            cfg.best_r,
            cfg.max_it,
            cfg.refine_tol,
            cfg.max_iter_scale,
            cfg.scale_tol,
            beta_out,
            bitgen_capsules=capsules,
            n_workers=ec_workers,
        )
        if status == 1:
            raise RuntimeError("fast_s: no non-singular subsamples were found")
        return FastSResult(
            coef=beta_out,
            scale=float(scale),
            converged=bool(converged),
            n_iter=int(n_iter),
            n_candidates_kept=cfg.best_r,
        )

    # ------------------------------------------------------------------
    # Resampling loop (optionally parallel via a thread pool)
    # ------------------------------------------------------------------
    if cfg.n_workers == 0 and not _auto_use_threads(n, p, cfg.nResample):
        # Auto mode but problem is small enough that threading hurts.
        n_workers = 1
    else:
        n_workers = _resolve_n_workers(cfg.n_workers, cfg.nResample)
    seed_seq = _to_seed_sequence(seed)

    if n_workers == 1:
        result = _resample_chunk(X, y, cfg, seed_seq, cfg.nResample)
        if result.early_exact_fit is not None:
            return result.early_exact_fit
        best_betas = result.best_betas
        best_scales = result.best_scales
    else:
        # Split nResample across n_workers chunks. Each worker uses a
        # SeedSequence-spawned PCG64, so results are deterministic per
        # (seed, n_workers, nResample).
        chunk_sizes = _split_iters(cfg.nResample, n_workers)
        child_seeds = seed_seq.spawn(n_workers)
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            futures = [
                ex.submit(_resample_chunk, X, y, cfg, ss, m)
                for ss, m in zip(child_seeds, chunk_sizes, strict=True)
                if m > 0
            ]
            chunk_results = [f.result() for f in futures]

        # An exact fit (m_scale==0) short-circuits everything.
        for cr in chunk_results:
            if cr.early_exact_fit is not None:
                return cr.early_exact_fit

        best_betas = []
        best_scales = []
        for cr in chunk_results:
            best_betas.extend(cr.best_betas)
            best_scales.extend(cr.best_scales)
        if len(best_scales) > cfg.best_r:
            order = np.argsort(best_scales)[: cfg.best_r]
            best_betas = [best_betas[i] for i in order]
            best_scales = [best_scales[i] for i in order]

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
