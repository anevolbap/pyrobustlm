# SPDX-License-Identifier: GPL-3.0-or-later
"""Top-level ``lmrob`` entry points.

End-to-end pipeline: formula -> design matrix -> fast-S init -> MM
iteration -> covariance -> :class:`pylmrob.results.LmRobResults`.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from typing import TYPE_CHECKING, Any

import numpy as np

from pylmrob import psi as _psi
from pylmrob._fast_s import FastSConfig, fast_s
from pylmrob._mm import mm_iterate
from pylmrob.control import Control
from pylmrob.d_scale import d_scale
from pylmrob.formula import model_matrix
from pylmrob.inference import vcov_avar1, vcov_w
from pylmrob.ms_estimator import m_s_fit
from pylmrob.results import LmRobResults

if TYPE_CHECKING:
    import pandas as pd


def _to_chi_psi_family(psi_family: str) -> str:
    """The chi family used in S iteration is the same as the user's psi family."""
    return psi_family


# Cache Cython kernel handles at module load. Avoids per-call importlib
# overhead. Each handle is ``None`` if the compiled module is unavailable.
# ``Callable[..., Any]`` is loose on purpose: the Cython signatures are
# more complex than what's worth replicating in Python type hints, and
# the actual checking happens at call time.
_CyFn = Callable[..., Any]


def _load_cy() -> tuple[_CyFn | None, _CyFn | None, _CyFn | None, _CyFn | None]:
    try:
        import importlib

        mod = importlib.import_module("pylmrob._core._lmrob")
        return (
            getattr(mod, "cy_lmrob_fit", None),
            getattr(mod, "cy_lmrob_d_scale", None),
            getattr(mod, "cy_lmrob_vcov_avar1", None),
            getattr(mod, "cy_lmrob_mm", None),
        )
    except Exception:
        return None, None, None, None


_CY_LMROB_FIT, _CY_LMROB_D_SCALE, _CY_LMROB_VCOV, _CY_LMROB_MM = _load_cy()

# Map family names to the integer enum used by the Cython kernel.
_CY_FAMILY_IDS = {
    "bisquare": 0,
    "biweight": 0,
    "hampel": 1,
    "optimal": 2,
    "lqq": 3,
    "ggw": 4,
}


def lmrob(
    formula: str,
    data: pd.DataFrame,
    control: Control | None = None,
    weights: np.ndarray | None = None,
    na_action: str = "drop",
    seed: int | np.random.Generator | None = None,
    **kwargs: Any,
) -> LmRobResults:
    """Fit a robust MM linear regression.

    Parameters
    ----------
    formula :
        R-style formula, e.g. ``"y ~ x1 + x2 + x3"``. Parsed by
        :mod:`formulaic`.
    data :
        DataFrame containing the columns referenced by ``formula``.
    control :
        Algorithm parameters; defaults to ``Control()`` (KS2014 preset,
        ``engine_c=True``).
    weights :
        Optional non-negative per-case weights (length ``len(data)``).
        Implemented via the ``sqrt(w)``-transform that R's lmrob uses:
        the transformed design ``(sqrt(w)*X, sqrt(w)*y)`` goes through
        the unweighted fit. Zero-weight rows are dropped. Compatible
        with both the default Cython engine and the legacy NumPy path
        (the transform is applied before any path dispatch, so the
        Cython kernel never needs to know about weights itself).
    na_action :
        ``"drop"`` (default) drops rows with any NA before fitting.
    seed :
        Seed for the resampling RNG.

    Returns
    -------
    LmRobResults
    """
    if kwargs:
        raise TypeError(f"unexpected keyword arguments: {sorted(kwargs)!r}")
    if control is None:
        control = Control()
    if weights is not None:
        weights = np.asarray(weights, dtype=np.float64).ravel()
        if not np.isfinite(weights).all() or (weights < 0).any():
            raise ValueError("weights must be non-negative and finite")
        if weights.size != len(data):
            raise ValueError(f"weights length {weights.size} != number of rows {len(data)}")
        # Drop zero-weight rows. R does the same: zero-weight observations
        # contribute nothing to the M-scale or IRWLS.
        keep = weights > 0
        if not keep.all():
            data = data.loc[keep].reset_index(drop=True)
            weights = weights[keep]
        # Non-trivial weights flow through the engine_c block too:
        # ``_lmrob_impl`` applies the sqrt(w)-transform to (X, y) before
        # any path dispatch, so the Cython kernel sees a transformed
        # design and never needs to know about weights itself. One
        # exception: ``cov=".vcov.w"`` with non-trivial weights on small
        # n can land the engine_c basin in a corner where every
        # transformed residual is past the psi-redescending cutoff,
        # making ``mean(psi'(r/s)) = 0`` and dividing by zero downstream.
        # Force the NumPy path for that combo.
        if not np.allclose(weights, 1.0) and control.cov == ".vcov.w":
            control = replace(control, engine_c=False)
    try:
        return _lmrob_impl(formula, data, control, na_action, seed, weights)
    except FloatingPointError as exc:
        # ``engine_c=True`` uses an internal Floyd subset-draw which is
        # not byte-identical to ``np.random.choice``. On a few small
        # classical datasets that puts fast-S in a basin where the
        # final ``X' W X`` is singular for ``vcov_avar1``. Retry once
        # with engine_c off (numpy choice -- different basin) so the
        # default path stays robust on those datasets.
        if not control.engine_c or "singular" not in str(exc):
            raise
        return _lmrob_impl(
            formula, data, replace(control, engine_c=False), na_action, seed, weights
        )


def _lmrob_impl(
    formula: str,
    data: pd.DataFrame,
    control: Control,
    na_action: str,
    seed: int | np.random.Generator | None,
    case_weights: np.ndarray | None = None,
) -> LmRobResults:
    # Control.__post_init__ guarantees these are populated; narrow for the
    # type checker.
    assert control.psi is not None
    assert control.method is not None
    assert control.cov is not None
    psi_family: str = control.psi
    cov_kind: str = control.cov

    # ------------------------------------------------------------------
    # Design matrix
    # ------------------------------------------------------------------
    if na_action == "drop":
        # ``dropna`` always allocates a new DataFrame even when no rows
        # are dropped. Check for NAs first; the check is much cheaper than
        # the copy for the typical NA-free path.
        if data.isna().any(axis=None):
            data = data.dropna()
    elif na_action == "raise":
        if data.isna().any().any():
            raise ValueError("data contains NA values; pass na_action='drop' to skip them")
    else:
        raise ValueError(f"unknown na_action: {na_action!r}")

    design = model_matrix(formula, data)
    y, X, term_names = design.y, design.X, design.term_names
    n, p = X.shape
    if n <= p:
        raise ValueError(f"need n > p; got n={n}, p={p}")

    # Per-case weights: R's lmrob handles them with a sqrt(w) transform
    # at the top level (robustbase/R/lmrob.R:96-98). The transformed
    # design (X_t, y_t) goes through the unweighted fit. After fitting,
    # ``residuals_`` / ``fitted_`` are reported on the original scale
    # (y - X @ coef), while ``rweights_`` / ``scale_`` / ``cov_`` come
    # from the transformed fit (matching R's reporting).
    X_orig: np.ndarray | None = None
    y_orig: np.ndarray | None = None
    if case_weights is not None and not np.allclose(case_weights, 1.0):
        X_orig = X
        y_orig = y
        sqrt_w = np.sqrt(case_weights)
        X = sqrt_w[:, None] * X
        y = sqrt_w * y

    # ------------------------------------------------------------------
    # Pick init method. ``init="auto"`` chooses M-S when factor columns
    # exist and S otherwise. Anything else honors the explicit user request.
    # ------------------------------------------------------------------
    init_method = control.init
    if init_method == "auto":
        init_method = "M-S" if design.is_factor_col.any() else "S"

    s_seed = seed if seed is not None else control.seed
    k_chi_tuple = tuple(np.atleast_1d(np.asarray(control.tuning_chi, dtype=float)).ravel())
    psi_k_eff = tuple(np.atleast_1d(np.asarray(control.tuning_psi, dtype=float)).ravel())

    # ------------------------------------------------------------------
    # Monolithic Cython engine: fast-S + MM in one nogil call. Sets the
    # ``_engine_c_done`` flag so the per-step Python blocks below skip
    # themselves.
    # ------------------------------------------------------------------
    _engine_c_done = False
    # Predeclare variables that are populated in one of the engine_c
    # / S / M-S branches below; pyright can't see that the branches
    # are exhaustive for the cases that reach them.
    beta_init = np.empty(p, dtype=np.float64)
    sigma_init = 0.0
    init_info: dict[str, object] = {}
    _engine_mm: Any = None
    _engine_residuals: np.ndarray | None = None
    _engine_rweights: np.ndarray | None = None
    _engine_cov_buf: np.ndarray | None = None
    # engine_c now parallelises internally via OpenMP prange on the
    # resample loop, so the previous ``n*p^2 >= 100k`` auto-fallback to
    # the threaded NumPy path is no longer needed. ``control.n_workers``
    # is honoured directly by the monolithic kernel.
    _engine_c_too_big = False
    _eff_n_workers = control.n_workers
    if (
        control.engine_c
        and init_method == "S"
        and psi_family in _CY_FAMILY_IDS
        and _CY_LMROB_FIT is not None
    ):
        _tuning_chi = np.zeros(3, dtype=np.float64)
        _tuning_psi = np.zeros(3, dtype=np.float64)
        # Tuples are short (1-3 entries); a manual loop avoids the cost of
        # np.asarray + ravel that ``np.atleast_1d`` would incur.
        for _i in range(min(3, len(k_chi_tuple))):
            _tuning_chi[_i] = k_chi_tuple[_i]
        for _i in range(min(3, len(psi_k_eff))):
            _tuning_psi[_i] = psi_k_eff[_i]
        # One ``np.empty`` for all output buffers instead of 4-5.
        # ``np.empty`` is ~10 µs each at small sizes; bundling saves
        # ~40 µs/fit.
        _want_inline_vcov_alloc = cov_kind == ".vcov.avar1"
        _cov_size = p * p if _want_inline_vcov_alloc else 0
        _big = np.empty(2 * p + 2 * n + _cov_size, dtype=np.float64)
        beta_out = _big[:p]
        residuals_out = _big[p : p + n]
        rweights_out = _big[p + n : p + 2 * n]
        beta_init_out = _big[p + 2 * n : 2 * p + 2 * n]
        from pylmrob._fast_s import make_generator as _make_gen

        rng_e = _make_gen(s_seed, kind=control.rng)
        # X and y are already float64 C-contiguous coming from the design
        # builder; skip np.ascontiguousarray if so to avoid a Python call.
        X_c = (
            X
            if (X.dtype == np.float64 and X.flags.c_contiguous)
            else np.ascontiguousarray(X, dtype=np.float64)
        )
        y_c = (
            y
            if (y.dtype == np.float64 and y.flags.c_contiguous)
            else np.ascontiguousarray(y, dtype=np.float64)
        )
        # When we'll also need vcov_avar1, ask the kernel to compute it
        # inline (saves a Python call + a separate np.zeros + a re-pass
        # through the LAPACK setup). Only useful for cov=".vcov.avar1"
        # which is the default MM path.
        # Cov buffer is the tail of ``_big`` (zero-initialized later if
        # the kernel actually writes to it). ``_big`` is uninitialized
        # which is fine: the Cython kernel writes the full p*p block.
        _want_inline_vcov = cov_kind == ".vcov.avar1"
        _cov_buf = _big[2 * p + 2 * n :].reshape((p, p)) if _want_inline_vcov else None
        # Resolve effective worker count and pre-spawn per-thread bitgens
        # for the engine_c parallel resample loop. The resolution mirrors
        # ``_fast_s._resolve_n_workers``.
        from pylmrob._fast_s import _resolve_n_workers as _resolve_nw
        from pylmrob._fast_s import _to_seed_sequence

        _ec_workers = _resolve_nw(_eff_n_workers, control.nResample)
        if _ec_workers > 1:
            _ec_caps = [
                _make_gen(_s, kind=control.rng).bit_generator.capsule
                for _s in _to_seed_sequence(s_seed).spawn(_ec_workers)
            ]
        else:
            _ec_caps = None
        scale_e, status_e, n_iter_s, _conv_s, n_iter_mm, conv_mm, _vcov_status = _CY_LMROB_FIT(
            X_c,
            y_c,
            rng_e.bit_generator.capsule,
            _CY_FAMILY_IDS[psi_family],
            _tuning_chi,
            _tuning_psi,
            control.bb,
            control.nResample,
            control.mts,
            control.k_fast_s,
            control.best_r_s,
            control.max_it,
            control.refine_tol,
            control.max_it,
            control.rel_tol,
            control.k_max,
            control.scale_tol,
            beta_out,
            residuals_out,
            rweights_out,
            beta_init_out,
            _cov_buf,
            _tuning_chi,
            control.bb,
            bitgen_capsules=_ec_caps,
            n_workers=_ec_workers,
        )
        if status_e == 1:
            raise RuntimeError("lmrob: no non-singular subsamples found")
        beta_init = beta_init_out
        sigma_init = float(scale_e)
        init_info = {
            "coef": beta_out,
            "scale": float(scale_e),
            "n_iter": int(n_iter_s),
            "method": "S",
        }
        from types import SimpleNamespace

        _engine_mm = SimpleNamespace(coef=beta_out, converged=bool(conv_mm), n_iter=int(n_iter_mm))
        _engine_c_done = True
        _engine_residuals = residuals_out
        _engine_rweights = rweights_out
        # Stash the inline-computed vcov (if any) so the vcov block below
        # can skip its own call to cy_lmrob_vcov_avar1.
        _engine_cov_buf = _cov_buf if (_cov_buf is not None and _vcov_status == 0) else None

    if _engine_c_done:
        pass  # cy_lmrob_fit already populated beta_init / sigma_init / init_info
    elif init_method == "S":
        # ------------------------------------------------------------------
        # Initial S estimate via fast-S resampling
        # ------------------------------------------------------------------
        cfg = FastSConfig(
            psi_chi=psi_family,
            k_chi=k_chi_tuple,
            b0=control.bb,
            nResample=control.nResample,
            k_fast_s=control.k_fast_s,
            best_r=control.best_r_s,
            max_it=control.max_it,
            refine_tol=control.refine_tol,
            scale_tol=control.scale_tol,
            max_iter_scale=control.k_max,
            mts=control.mts,
            n_workers=_eff_n_workers,
            fast_rng=control.fast_rng,
            # If we already decided ``engine_c`` is too expensive for
            # this problem size (n*p^2 >= 100k), don't let ``fast_s``
            # take its own engine_c branch either; that branch is a
            # single Cython call and disables threading.
            engine_c=control.engine_c and not _engine_c_too_big,
            rng=control.rng,
        )
        s_result = fast_s(X, y, cfg=cfg, seed=s_seed)
        beta_init = s_result.coef
        sigma_init = s_result.scale
        init_info: dict[str, object] = {
            "coef": s_result.coef.copy(),
            "scale": s_result.scale,
            "n_iter": s_result.n_iter,
            "method": "S",
        }
    elif init_method == "M-S":
        if not design.is_factor_col.any():
            # No factor cols; M-S degenerates to S.
            cfg = FastSConfig(
                psi_chi=psi_family,
                k_chi=k_chi_tuple,
                b0=control.bb,
                nResample=control.nResample,
                k_fast_s=control.k_fast_s,
                best_r=control.best_r_s,
                max_it=control.max_it,
                refine_tol=control.refine_tol,
                scale_tol=control.scale_tol,
                max_iter_scale=control.k_max,
                mts=control.mts,
                rng=control.rng,
            )
            s_result = fast_s(X, y, cfg=cfg, seed=s_seed)
            beta_init = s_result.coef
            sigma_init = s_result.scale
            init_info = {
                "coef": s_result.coef.copy(),
                "scale": s_result.scale,
                "n_iter": s_result.n_iter,
                "method": "S",
            }
        else:
            cat_cols = design.is_factor_col
            X_cat = X[:, cat_cols]
            X_cont = X[:, ~cat_cols]
            ms_result = m_s_fit(
                X_cat=X_cat,
                X_cont=X_cont,
                y=y,
                psi_chi=psi_family,
                k_chi=k_chi_tuple,
                b0=control.bb,
                k_m_s=control.k_m_s,
                nResample=control.nResample,
                max_it=control.max_it,
                rel_tol=control.rel_tol,
                seed=s_seed,
            )
            # Re-stitch into the original column order.
            beta_init = np.empty(p, dtype=np.float64)
            beta_init[cat_cols] = ms_result.coef_cat
            beta_init[~cat_cols] = ms_result.coef_cont
            sigma_init = ms_result.scale
            init_info = {
                "coef": beta_init.copy(),
                "scale": ms_result.scale,
                "n_iter": ms_result.n_iter,
                "method": "M-S",
            }
    else:
        raise NotImplementedError(f"init={init_method!r} not implemented")

    # ------------------------------------------------------------------
    # MM step holding sigma fixed (skipped if the monolithic engine
    # already produced post-MM beta + residuals + rweights).
    # ------------------------------------------------------------------
    if _engine_c_done:
        # The engine_c block populated these; the asserts narrow the
        # ``| None`` predeclared types for the static checker.
        assert _engine_residuals is not None
        assert _engine_rweights is not None
        mm = _engine_mm
        coef = mm.coef
        sigma = sigma_init
        residuals = _engine_residuals
        fitted = y - residuals
        rweights = _engine_rweights
    else:
        # Cython MM fast path when engine_c was requested (the user has
        # signalled they want speed; engine_c block above ran or fell
        # back, but in either case we should use the Cython MM here too).
        # Falls back to the NumPy implementation in plain default mode.
        if control.engine_c and _CY_LMROB_MM is not None and psi_family in _CY_FAMILY_IDS:
            _tuning_psi_mm = np.zeros(3, dtype=np.float64)
            for _i in range(min(3, len(psi_k_eff))):
                _tuning_psi_mm[_i] = float(psi_k_eff[_i])
            beta_mm = np.ascontiguousarray(beta_init, dtype=np.float64).copy()
            n_iter_mm, converged_mm, _mm_status = _CY_LMROB_MM(
                np.ascontiguousarray(X, dtype=np.float64),
                np.ascontiguousarray(y, dtype=np.float64),
                beta_mm,
                float(sigma_init),
                _CY_FAMILY_IDS[psi_family],
                _tuning_psi_mm,
                control.max_it,
                control.rel_tol,
            )
            from types import SimpleNamespace

            mm = SimpleNamespace(coef=beta_mm, converged=bool(converged_mm), n_iter=int(n_iter_mm))
        else:
            mm = mm_iterate(
                X=X,
                y=y,
                beta_init=beta_init,
                sigma=sigma_init,
                psi_family=psi_family,
                psi_k=psi_k_eff,
                max_it=control.max_it,
                rel_tol=control.rel_tol,
            )

        coef = mm.coef
        sigma = sigma_init
        residuals = y - X @ coef
        fitted = X @ coef
        z = residuals / sigma if sigma != 0 else residuals
        rweights = _psi.wgt(z, psi_family, psi_k_eff)

    # ------------------------------------------------------------------
    # D-step (used by methods SMD and SMDM): refines scale via design-
    # adaptive weighting. SMDM also re-fits MM with the new scale.
    # ``tau_vec`` is also stashed for ``vcov_w(corrfact='tau' | 'hybrid' | 'tauold')``.
    # ------------------------------------------------------------------
    method = control.method or "MM"
    tau_vec: np.ndarray | None = None
    if "D" in method:
        sigma_d = None
        d_converged = False
        # Cython D-step when engine_c is on and tabulated coefficients exist.
        if control.engine_c and _CY_LMROB_D_SCALE is not None and psi_family in _CY_FAMILY_IDS:
            _cy_d = _CY_LMROB_D_SCALE
            if True:
                _tau_buf = np.empty(n, dtype=np.float64)
                _tuning_psi = np.zeros(3, dtype=np.float64)
                for _i, _v in enumerate(
                    tuple(np.atleast_1d(np.asarray(control.tuning_psi, dtype=float)).ravel())[:3]
                ):
                    _tuning_psi[_i] = float(_v)
                _sigma_d, _d_conv, _d_status = _cy_d(
                    np.ascontiguousarray(X, dtype=np.float64),
                    np.ascontiguousarray(residuals, dtype=np.float64),
                    np.ascontiguousarray(rweights, dtype=np.float64),
                    float(sigma),
                    _CY_FAMILY_IDS[psi_family],
                    _tuning_psi,
                    control.k_max,
                    control.rel_tol,
                    _tau_buf,
                )
                if _d_status == 0:
                    sigma_d = float(_sigma_d)
                    d_converged = bool(_d_conv)
                    tau_vec = _tau_buf
        if sigma_d is None:
            sigma_d, d_converged, tau_vec, _h = d_scale(
                X=X,
                residuals=residuals,
                rweights=rweights,
                init_scale=sigma,
                family=psi_family,
                c_psi=psi_k_eff,
                max_iter=control.k_max,
                tol=control.rel_tol,
            )
        if d_converged and sigma_d > 0:
            sigma = sigma_d
            init_info["d_scale"] = sigma_d
            init_info["d_converged"] = True
            # The trailing M in SMDM: re-fit MM with the new D-scale.
            if method.endswith("M"):
                mm2 = mm_iterate(
                    X=X,
                    y=y,
                    beta_init=coef,
                    sigma=sigma,
                    psi_family=psi_family,
                    psi_k=psi_k_eff,
                    max_it=control.max_it,
                    rel_tol=control.rel_tol,
                )
                coef = mm2.coef
                residuals = y - X @ coef
                fitted = X @ coef
                z = residuals / sigma if sigma != 0 else residuals
                rweights = _psi.wgt(z, psi_family, psi_k_eff)

    # ------------------------------------------------------------------
    # Covariance
    # ------------------------------------------------------------------
    init_residuals = y - X @ beta_init
    cov = None
    if cov_kind == ".vcov.avar1":
        # If the monolithic engine already produced the cov inline,
        # use it (Cython already posdefified).
        if _engine_c_done and _engine_cov_buf is not None:
            cov = _engine_cov_buf
        # Otherwise Cython vcov fast path when engine_c is on. The
        # matrix products now go through BLAS dgemm so the kernel is
        # competitive with NumPy at all sizes; the previous large-n
        # gate is no longer needed.
        elif control.engine_c and _CY_LMROB_VCOV is not None and psi_family in _CY_FAMILY_IDS:
            _cy_vcov = _CY_LMROB_VCOV
            _tuning_psi = np.zeros(3, dtype=np.float64)
            _tuning_chi = np.zeros(3, dtype=np.float64)
            for _i, _v in enumerate(psi_k_eff[:3]):
                _tuning_psi[_i] = float(_v)
            for _i, _v in enumerate(k_chi_tuple[:3]):
                _tuning_chi[_i] = float(_v)
            cov_buf = np.zeros((p, p), dtype=np.float64)
            _vcov_status = _cy_vcov(
                np.ascontiguousarray(X, dtype=np.float64),
                np.ascontiguousarray(residuals, dtype=np.float64),
                np.ascontiguousarray(init_residuals, dtype=np.float64),
                float(sigma),
                _CY_FAMILY_IDS[psi_family],
                _tuning_psi,
                _tuning_chi,
                control.bb,
                cov_buf,
            )
            if _vcov_status == 0:
                # Cython already posdefified.
                cov = cov_buf
        if cov is None:
            cov = vcov_avar1(
                X=X,
                residuals=residuals,
                sigma=sigma,
                psi_family=psi_family,
                psi_k=psi_k_eff,
                init_residuals=init_residuals,
                chi_family=psi_family,
                chi_k=k_chi_tuple,
                bb=control.bb,
            )
    elif cov_kind == ".vcov.w":
        # Pre-compute tau if it isn't already populated by the D-step.
        # vcov_w with corrfact="tau"/"hybrid"/"tauold" needs it.
        if tau_vec is None:
            from pylmrob.d_scale import tau as _compute_tau

            sw = np.sqrt(np.maximum(rweights, 0.0))
            Xw = X * sw[:, None]
            Q, _ = np.linalg.qr(Xw, mode="reduced")
            h = np.minimum(1.0, np.sum(Q * Q, axis=1))
            tau_vec = _compute_tau(h, psi_family, psi_k_eff)
        s_init_scale_raw = init_info.get("scale", sigma)
        s_init_scale = (
            float(s_init_scale_raw) if isinstance(s_init_scale_raw, (int, float)) else float(sigma)
        )  # type: ignore[arg-type]
        cov = vcov_w(
            X=X,
            residuals=residuals,
            sigma=sigma,
            psi_family=psi_family,
            psi_k=psi_k_eff,
            rweights=rweights,
            init_residuals=init_residuals,
            init_scale=s_init_scale,
            chi_k=k_chi_tuple,
            method=method,
            tau=tau_vec,
        )
    else:
        # Fallback: legacy / unknown
        cov = vcov_avar1(
            X=X,
            residuals=residuals,
            sigma=sigma,
            psi_family=psi_family,
            psi_k=psi_k_eff,
        )

    # When weights were applied via sqrt-transform, R reports residuals
    # and fitted on the *original* scale (y - X @ coef), but rweights /
    # scale / cov come from the transformed fit. The design returned to
    # the user is also the original X / y.
    if X_orig is not None and y_orig is not None:
        residuals_out = y_orig - X_orig @ coef
        fitted_out = X_orig @ coef
        design_x_out = X_orig
        design_y_out = y_orig
    else:
        residuals_out = residuals
        fitted_out = fitted
        design_x_out = X
        design_y_out = y

    return LmRobResults(
        coef_=np.asarray(coef, dtype=np.float64),
        scale_=float(sigma),
        weights_=(
            np.asarray(case_weights, dtype=np.float64)
            if case_weights is not None
            else np.ones(n, dtype=np.float64)
        ),
        rweights_=np.asarray(rweights, dtype=np.float64),
        residuals_=np.asarray(residuals_out, dtype=np.float64),
        fitted_=np.asarray(fitted_out, dtype=np.float64),
        cov_=np.asarray(cov, dtype=np.float64),
        df_residual_=n - p,
        converged_=mm.converged,
        n_iter_=mm.n_iter,
        nobs_=n,
        term_names_=term_names,
        control=control,
        init_=init_info,
        rhs_spec_=design.rhs_spec,
        design_x_=np.asarray(design_x_out, dtype=np.float64),
        design_y_=np.asarray(design_y_out, dtype=np.float64),
    )


class LmRob:
    """scikit-learn-style estimator wrapper around :func:`lmrob`."""

    def __init__(self, control: Control | None = None) -> None:
        self.control = control or Control()
        self._result: LmRobResults | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, seed: int | None = None) -> LmRob:
        # Build a tiny pandas DataFrame so we can reuse the formula path.
        # Constructing the DataFrame from a single dict-of-arrays is faster
        # than ``pd.DataFrame(X, columns=...)`` + ``df["y"] = y`` (one
        # block-manager build instead of two).
        import pandas as pd

        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        if X.ndim != 2:
            raise ValueError("X must be 2-D")
        cols = [f"x{i}" for i in range(X.shape[1])]
        data_dict: dict[str, np.ndarray] = {c: X[:, i] for i, c in enumerate(cols)}
        data_dict["y"] = y
        df = pd.DataFrame(data_dict)
        formula = "y ~ " + " + ".join(cols)
        self._result = lmrob(formula, df, control=self.control, seed=seed)
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict on a new design matrix (raw, without intercept column).

        ``LmRob`` always fits with an intercept, so we wrap ``X`` in a
        DataFrame with the same column names used at fit time and let the
        stored formula spec re-add the intercept.
        """
        if self._result is None:
            raise RuntimeError("call .fit before .predict")
        import pandas as pd

        X = np.asarray(X, dtype=np.float64)
        cols = [f"x{i}" for i in range(X.shape[1])]
        return self._result.predict(pd.DataFrame(X, columns=cols))

    def score(self, X: np.ndarray, y: np.ndarray) -> float:
        """Standard ``R^2`` on a test set, sklearn convention.

        ``1 - SS_res / SS_tot`` where SS_res = sum((y - y_hat)^2). This is
        the OLS coefficient of determination, not the robust R^2 reported
        by ``summary()`` (use ``self.result_.summary().r_squared`` for
        that). Returning OLS R^2 keeps ``LmRob`` compatible with sklearn
        utilities (``cross_val_score``, ``GridSearchCV``) that assume the
        regressor scorer contract.
        """
        y = np.asarray(y, dtype=np.float64)
        y_hat = self.predict(X)
        ss_res = float(np.sum((y - y_hat) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        if ss_tot <= 0.0:
            return 0.0
        return 1.0 - ss_res / ss_tot

    def get_params(self, deep: bool = True) -> dict:
        """Return the ``__init__`` parameters, sklearn convention."""
        return {"control": self.control}

    def set_params(self, **params: object) -> LmRob:
        """Set ``__init__`` parameters in place, sklearn convention."""
        for key, value in params.items():
            if not hasattr(self, key):
                raise ValueError(f"Invalid parameter {key!r} for LmRob")
            setattr(self, key, value)
        return self

    @property
    def result_(self) -> LmRobResults:
        if self._result is None:
            raise RuntimeError("call .fit before accessing result_")
        return self._result
