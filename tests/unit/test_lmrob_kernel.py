# SPDX-License-Identifier: GPL-3.0-or-later
"""Stage 1: tests for the monolithic Cython lmrob kernel."""

from __future__ import annotations

import numpy as np
import pytest


def _make_problem(n: int, p: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n), rng.standard_normal((n, p - 1))])
    beta = rng.standard_normal(p)
    y = X @ beta + rng.standard_normal(n)
    return np.ascontiguousarray(X), np.ascontiguousarray(y)


def test_kernel_import() -> None:
    """The monolithic kernel module should be importable."""
    from pylmrob._core._lmrob import cy_lmrob_fast_s  # noqa: F401


def test_kernel_matches_fast_rng_path() -> None:
    """The monolithic kernel and the fast_rng=True NumPy path land on the
    same fit (same Cython bitgen draw sequence, identical math)."""
    from pylmrob._core._lmrob import cy_lmrob_fast_s

    from pylmrob._fast_s import FastSConfig, fast_s

    X, y = _make_problem(200, 6, seed=1)
    p = X.shape[1]

    # Existing path with the same Cython bitgen draw.
    cfg = FastSConfig(nResample=500, n_workers=1, fast_rng=True)
    ref = fast_s(X, y, cfg=cfg, seed=42)

    # New monolithic kernel.
    beta_out = np.empty(p, dtype=np.float64)
    rng = np.random.default_rng(42)
    tuning = np.array([1.54764, 0.0, 0.0], dtype=np.float64)
    scale, status, _n_iter, _conv = cy_lmrob_fast_s(
        X,
        y,
        rng.bit_generator.capsule,
        0,  # family: bisquare
        tuning,
        0.5,
        500,
        1000,
        1,
        2,
        50,
        1e-7,
        200,
        1e-10,
        beta_out,
    )
    assert status == 0
    np.testing.assert_allclose(scale, ref.scale, rtol=1e-6)
    np.testing.assert_allclose(beta_out, ref.coef, rtol=1e-6, atol=1e-8)


def test_engine_c_via_lmrob_api() -> None:
    """Control(engine_c=True) end-to-end produces sensible fits."""
    import pandas as pd

    from pylmrob import Control, lmrob

    X, y = _make_problem(200, 6, seed=2)
    df = pd.DataFrame(X[:, 1:], columns=[f"x{i}" for i in range(X.shape[1] - 1)])
    df["y"] = y
    cols = " + ".join(f"x{i}" for i in range(X.shape[1] - 1))

    fit = lmrob(f"y ~ {cols}", df, control=Control(engine_c=True), seed=42)
    assert fit.converged_
    assert fit.scale_ > 0
    assert fit.coef_.shape == (X.shape[1],)


@pytest.mark.skipif(
    bool(__import__("os").environ.get("CI")),
    reason="Wall-clock perf comparison is too noisy on CI runners",
)
def test_engine_c_speedup_at_small_n() -> None:
    """``engine_c=True`` should be faster than the default path at small n.

    Skipped on CI because shared-runner timing is too noisy to make a
    perf assertion that's both meaningful and reliable.
    """
    import time

    import pandas as pd

    from pylmrob import Control, lmrob

    X, y = _make_problem(200, 6, seed=3)
    df = pd.DataFrame(X[:, 1:], columns=[f"x{i}" for i in range(X.shape[1] - 1)])
    df["y"] = y
    cols = " + ".join(f"x{i}" for i in range(X.shape[1] - 1))
    formula = f"y ~ {cols}"

    # Warm up.
    lmrob(formula, df, seed=0)
    lmrob(formula, df, control=Control(engine_c=True), seed=0)

    def bench(ctrl, reps: int = 10) -> float:
        t0 = time.perf_counter()
        for _ in range(reps):
            lmrob(formula, df, control=ctrl, seed=0)
        return (time.perf_counter() - t0) / reps

    t_default = bench(Control(nResample=500, engine_c=False))
    t_engine_c = bench(Control(nResample=500, engine_c=True))

    assert t_engine_c < t_default, (
        f"engine_c not faster: default={t_default * 1000:.1f}ms, engine_c={t_engine_c * 1000:.1f}ms"
    )


def test_engine_c_ignored_for_non_bisquare() -> None:
    """``engine_c=True`` silently falls back for non-bisquare families."""
    import pandas as pd

    from pylmrob import Control, lmrob

    X, y = _make_problem(200, 6, seed=4)
    df = pd.DataFrame(X[:, 1:], columns=[f"x{i}" for i in range(X.shape[1] - 1)])
    df["y"] = y
    cols = " + ".join(f"x{i}" for i in range(X.shape[1] - 1))

    # Should not raise; should produce a valid fit via the NumPy path.
    fit = lmrob(
        f"y ~ {cols}",
        df,
        control=Control(psi="hampel", engine_c=True),
        seed=42,
    )
    assert fit.converged_
