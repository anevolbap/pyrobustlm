# SPDX-License-Identifier: GPL-3.0-or-later
"""Bootstrap inference for ``lmrob`` fits.

Case-resampling bootstrap: each replicate samples ``n`` observations
with replacement from the original ``(X, y)`` and refits ``lmrob`` on
the resample. Returns percentile and basic confidence intervals over
the bootstrap distribution of coefficients.

Matches R's ``boot::boot`` driven workflow for ``lmrob.S`` (R has
``robustbase::lmrob`` integration with the ``boot`` package via
``boot()`` -- this is the equivalent without the external dependency).
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pylmrob.control import Control
    from pylmrob.results import LmRobResults


@dataclass
class BootstrapResult:
    """Output of :func:`bootstrap`.

    Attributes
    ----------
    coefs :
        ``(n_boot, p)`` matrix of bootstrap coefficient draws. Replicates
        that failed to converge are dropped from this matrix.
    scales :
        ``(n_boot,)`` vector of bootstrap M-scale draws (same row
        indexing as ``coefs``).
    n_converged :
        Number of replicates that converged. Bootstrap CIs use only
        these.
    n_boot :
        Number of replicates that were requested.
    level :
        Confidence level used to build the intervals (default 0.95).
    percentile_ci :
        ``(p, 2)`` lower/upper from the empirical quantiles of ``coefs``.
        Standard percentile bootstrap, the most common default.
    basic_ci :
        ``(p, 2)`` lower/upper from the basic (reflected) bootstrap.
        ``2 * theta_hat - q_{1-alpha/2}`` and ``2 * theta_hat - q_{alpha/2}``.
    se :
        ``(p,)`` standard deviation of the bootstrap coefficient draws,
        a robust standard-error estimate.
    bias :
        ``(p,)`` ``mean(coefs) - theta_hat``. Bootstrap bias estimate.
    term_names :
        Names of the coefficients, parallel to ``percentile_ci`` rows.
    """

    coefs: np.ndarray
    scales: np.ndarray
    n_converged: int
    n_boot: int
    level: float
    percentile_ci: np.ndarray
    basic_ci: np.ndarray
    se: np.ndarray
    bias: np.ndarray
    term_names: list[str]


def _one_replicate(
    X: np.ndarray,
    y: np.ndarray,
    idx: np.ndarray,
    psi_family: str,
    control: Control,
    inner_seed: int,
) -> tuple[np.ndarray | None, float | None]:
    """Run one bootstrap replicate. Returns ``(coef, scale)`` or
    ``(None, None)`` on failure (singular subsample, exact fit, etc.)."""
    from pylmrob._fast_s import FastSConfig, fast_s
    from pylmrob._mm import mm_iterate

    Xb = X[idx]
    yb = y[idx]
    cfg = FastSConfig(
        psi_chi=psi_family,
        k_chi=tuple(np.atleast_1d(np.asarray(control.tuning_chi, dtype=float)).ravel()),
        b0=control.bb,
        nResample=int(control.nResample),
        k_fast_s=int(control.k_fast_s),
        best_r=int(control.best_r_s),
        max_it=int(control.max_it),
        refine_tol=float(control.refine_tol),
        scale_tol=float(control.scale_tol),
        max_iter_scale=int(control.k_max),
        mts=int(control.mts),
        engine_c=False,  # safer on the bootstrap resamples; basin
        # drift can flip a subset draw to a near-singular X' W X.
    )
    try:
        s = fast_s(Xb, yb, cfg=cfg, seed=inner_seed)
    except Exception:
        return None, None
    if s.scale == 0.0 or not np.isfinite(s.scale):
        return None, None
    psi_k = tuple(np.atleast_1d(np.asarray(control.tuning_psi, dtype=float)).ravel())
    mm = mm_iterate(
        X=Xb,
        y=yb,
        beta_init=s.coef,
        sigma=s.scale,
        psi_family=psi_family,
        psi_k=psi_k,
        max_it=int(control.max_it),
        rel_tol=float(control.rel_tol),
    )
    if not mm.converged:
        return None, None
    return np.asarray(mm.coef, dtype=np.float64), float(s.scale)


def bootstrap(
    fit: LmRobResults,
    n_boot: int = 1000,
    level: float = 0.95,
    seed: int | np.random.Generator | None = None,
    n_workers: int = 1,
) -> BootstrapResult:
    """Case-resampling bootstrap for an ``lmrob`` fit.

    Each replicate draws ``n`` row indices with replacement, refits
    ``lmrob`` on the resampled design, and stores the converged
    coefficient vector. Replicates that fail (singular subsample, MM
    didn't converge) are dropped silently; the converged count is on
    the result.

    Parameters
    ----------
    fit :
        Fitted :class:`LmRobResults` with the design matrix stashed
        (``design_x_`` / ``design_y_``). The default ``lmrob()`` call
        always stashes them.
    n_boot :
        Number of bootstrap replicates. Default 1000.
    level :
        Confidence level for the CIs. Default 0.95.
    seed :
        Seed for the resampling RNG (NumPy ``Generator`` or int).
        Replicates are deterministic for a fixed seed and ``n_workers``.
    n_workers :
        Number of Python threads to use. Each thread refits ``lmrob``
        on its share of replicates. 1 (default) runs serially.

    Returns
    -------
    BootstrapResult
    """
    if fit.design_x_ is None or fit.design_y_ is None:
        raise RuntimeError(
            "bootstrap needs the design matrix; fit was created without design_x_/design_y_"
        )
    X = np.ascontiguousarray(fit.design_x_, dtype=np.float64)
    y = np.ascontiguousarray(fit.design_y_, dtype=np.float64)
    n = X.shape[0]
    psi_family = fit.control.psi or "bisquare"

    rng = np.random.default_rng(seed)
    # Pre-generate all bootstrap indices and per-replicate inner seeds
    # so threads can run in parallel without contending on the RNG.
    # The inner seed is passed to ``fast_s`` to make the entire
    # replicate deterministic given ``seed``.
    all_idx = rng.integers(0, n, size=(n_boot, n))
    inner_seeds = rng.integers(0, 2**31 - 1, size=n_boot)

    coefs: list[np.ndarray] = []
    scales: list[float] = []

    def _run(b: int) -> tuple[np.ndarray | None, float | None]:
        return _one_replicate(X, y, all_idx[b], psi_family, fit.control, int(inner_seeds[b]))

    if n_workers > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            for c, s in pool.map(_run, range(n_boot)):
                if c is not None and s is not None:
                    coefs.append(c)
                    scales.append(s)
    else:
        for b in range(n_boot):
            c, s = _run(b)
            if c is not None and s is not None:
                coefs.append(c)
                scales.append(s)

    if not coefs:
        raise RuntimeError(
            f"bootstrap: 0 of {n_boot} replicates converged; check the data or n_boot"
        )

    coef_arr = np.asarray(coefs, dtype=np.float64)
    scale_arr = np.asarray(scales, dtype=np.float64)

    alpha = (1.0 - level) / 2.0
    q_lo = np.quantile(coef_arr, alpha, axis=0)
    q_hi = np.quantile(coef_arr, 1.0 - alpha, axis=0)
    percentile_ci = np.column_stack([q_lo, q_hi])
    # Basic CI reflects the percentile interval around the point estimate.
    basic_ci = np.column_stack([2.0 * fit.coef_ - q_hi, 2.0 * fit.coef_ - q_lo])
    se = np.std(coef_arr, axis=0, ddof=1)
    bias = np.mean(coef_arr, axis=0) - fit.coef_

    return BootstrapResult(
        coefs=coef_arr,
        scales=scale_arr,
        n_converged=int(coef_arr.shape[0]),
        n_boot=int(n_boot),
        level=float(level),
        percentile_ci=percentile_ci,
        basic_ci=basic_ci,
        se=se,
        bias=bias,
        term_names=list(fit.term_names_),
    )


__all__ = ["BootstrapResult", "bootstrap"]
