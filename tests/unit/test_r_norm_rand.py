# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit tests for :func:`pylmrob.r_qnorm` and :func:`pylmrob.r_norm_rand`.

Property-level checks. Byte-identical comparison against R's actual
``rnorm()`` / ``qnorm()`` lives in
``tests/validation/test_r_norm_rand_vs_R.py``.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm

from pylmrob import r_norm_rand, r_qnorm, r_set_seed


def test_qnorm_central() -> None:
    assert r_qnorm(0.5) == 0.0


def test_qnorm_symmetry() -> None:
    """qnorm(1-p) == -qnorm(p) for p in (0, 1)."""
    for p in [0.1, 0.25, 0.4, 0.49, 0.001]:
        assert math.isclose(r_qnorm(1 - p), -r_qnorm(p), abs_tol=1e-15)


def test_qnorm_boundaries() -> None:
    assert r_qnorm(0.0) == float("-inf")
    assert r_qnorm(1.0) == float("inf")
    assert math.isnan(r_qnorm(-0.1))
    assert math.isnan(r_qnorm(1.5))


def test_qnorm_agrees_with_scipy_to_15_digits() -> None:
    """The two algorithms (Wichura vs scipy's) agree to ~1e-15 except in
    the last ULP near the central peak; never worse than 1e-14."""
    for p in [0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999]:
        diff = abs(r_qnorm(p) - norm.ppf(p))
        assert diff < 1e-14


def test_norm_rand_reproducible() -> None:
    rng_a = r_set_seed(42)
    rng_b = r_set_seed(42)
    a = [r_norm_rand(rng_a) for _ in range(20)]
    b = [r_norm_rand(rng_b) for _ in range(20)]
    assert a == b


def test_norm_rand_finite() -> None:
    rng = r_set_seed(42)
    draws = np.array([r_norm_rand(rng) for _ in range(500)])
    assert np.isfinite(draws).all()


def test_norm_rand_advances_state_by_two_per_draw() -> None:
    """``r_norm_rand`` uses two ``unif_rand`` draws (R's Inversion)."""
    rng_a = r_set_seed(42)
    r_norm_rand(rng_a)
    next_after_norm = rng_a.unif_rand()

    rng_b = r_set_seed(42)
    rng_b.unif_rand()
    rng_b.unif_rand()
    next_after_two_unifs = rng_b.unif_rand()

    assert next_after_norm == next_after_two_unifs


def test_norm_rand_distribution_sanity() -> None:
    """1000 draws have approximately mean 0 and std 1 (loose check, only
    to catch a totally broken implementation)."""
    rng = r_set_seed(123)
    draws = np.array([r_norm_rand(rng) for _ in range(2000)])
    assert abs(draws.mean()) < 0.1
    assert abs(draws.std(ddof=1) - 1.0) < 0.1
