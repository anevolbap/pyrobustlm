# SPDX-License-Identifier: GPL-3.0-or-later
"""Tests for the parallel fast-S resampling path."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pylmrob._fast_s import (
    FastSConfig,
    _auto_use_threads,
    _resolve_n_workers,
    _split_iters,
    fast_s,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def _make_problem(n: int, p: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n), rng.standard_normal((n, p - 1))])
    beta = rng.standard_normal(p)
    y = X @ beta + rng.standard_normal(n)
    return X, y


def test_split_iters_balanced() -> None:
    assert _split_iters(10, 3) == [4, 3, 3]
    assert _split_iters(7, 7) == [1, 1, 1, 1, 1, 1, 1]
    assert _split_iters(0, 4) == [0, 0, 0, 0]


def test_resolve_n_workers_explicit() -> None:
    assert _resolve_n_workers(1, 500) == 1
    assert _resolve_n_workers(4, 500) == 4
    assert _resolve_n_workers(99, 500) == 99


def test_auto_use_threads_small_off() -> None:
    # Tiny problems should not auto-thread.
    assert _auto_use_threads(n=21, p=4, n_iter=500) is False
    assert _auto_use_threads(n=100, p=5, n_iter=500) is False


def test_auto_use_threads_large_on() -> None:
    # Big problems should auto-thread.
    assert _auto_use_threads(n=5000, p=30, n_iter=2000) is True


def test_parallel_matches_serial_within_basin() -> None:
    """The parallel path lands in the same basin as the serial path.

    Bit-identical equivalence is not guaranteed because chunking changes
    the sequence of subsamples each PCG64 sees, but the final scale and
    coefficients should agree to the same basin tolerance we already
    accept (rtol=1e-4 on the test problem).
    """
    X, y = _make_problem(n=2000, p=10, seed=42)

    serial = fast_s(X, y, cfg=FastSConfig(nResample=200, n_workers=1), seed=42)
    parallel = fast_s(X, y, cfg=FastSConfig(nResample=200, n_workers=4), seed=42)

    np.testing.assert_allclose(parallel.scale, serial.scale, rtol=1e-4)
    np.testing.assert_allclose(parallel.coef, serial.coef, rtol=1e-4, atol=1e-6)


def test_parallel_deterministic_repeated() -> None:
    """Same seed + same n_workers => identical fast-S result on every call."""
    X, y = _make_problem(n=500, p=8, seed=0)

    cfg = FastSConfig(nResample=200, n_workers=4)
    a = fast_s(X, y, cfg=cfg, seed=99)
    b = fast_s(X, y, cfg=cfg, seed=99)

    np.testing.assert_array_equal(a.coef, b.coef)
    assert a.scale == b.scale


def test_n_workers_via_control() -> None:
    """Control.n_workers propagates to FastSConfig and runs end-to-end."""
    import pandas as pd

    from pylmrob import Control, lmrob

    X, y = _make_problem(n=100, p=4, seed=1)
    df = pd.DataFrame(X, columns=["intercept", "x1", "x2", "x3"])
    df["y"] = y

    fit_serial = lmrob(
        "y ~ x1 + x2 + x3",
        df,
        control=Control(nResample=100, n_workers=1),
        seed=42,
    )
    fit_parallel = lmrob(
        "y ~ x1 + x2 + x3",
        df,
        control=Control(nResample=100, n_workers=2),
        seed=42,
    )

    np.testing.assert_allclose(fit_parallel.coef_, fit_serial.coef_, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(fit_parallel.scale_, fit_serial.scale_, rtol=1e-4)


@pytest.mark.parametrize("n_workers", [1, 0, 2, 4])
def test_classical_dataset_parity(n_workers: int) -> None:
    """All four worker settings produce the same fit on stackloss
    within basin tolerance."""
    import pandas as pd

    from pylmrob import Control, lmrob

    df = pd.read_csv(REPO_ROOT / "tests" / "data" / "stackloss.csv")
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        df,
        control=Control(nResample=500, n_workers=n_workers),
        seed=42,
    )
    # Regression: this fit should converge for every worker setting.
    assert fit.converged_
    # And land near R's reference scale (rtol=1e-2 to absorb basin drift).
    np.testing.assert_allclose(fit.scale_, 1.91234, rtol=1e-2)
