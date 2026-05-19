# SPDX-License-Identifier: GPL-3.0-or-later
"""Byte-identical agreement with R's ``rnorm()`` and ``qnorm()``.

R's default ``RNGkind()[2]`` is ``"Inversion"``: each ``rnorm`` draw
consumes two ``unif_rand`` outputs combined into a single >27-bit
uniform, then ``qnorm5()`` (Wichura AS 241). Both are ported in
:mod:`pylmrob.rng`.
"""

from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pytest

from pylmrob import r_norm_rand, r_qnorm, r_set_seed


def _have_r() -> bool:
    return shutil.which("Rscript") is not None


pytestmark = pytest.mark.skipif(not _have_r(), reason="Rscript not on PATH")


def _r_rnorm(seed: int, n: int) -> np.ndarray:
    script = f'set.seed({seed}); cat(sprintf("%.17g\\n", rnorm({n})), sep="")'
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True,
        check=True,
        text=True,
        timeout=20,
    )
    return np.array([float(x) for x in out.stdout.splitlines() if x], dtype=np.float64)


def _r_qnorm_at(ps: list[float]) -> np.ndarray:
    arg = ", ".join(repr(p) for p in ps)
    script = f'cat(sprintf("%.17g\\n", qnorm(c({arg}))), sep="")'
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True,
        check=True,
        text=True,
        timeout=10,
    )
    return np.array([float(x) for x in out.stdout.splitlines() if x], dtype=np.float64)


@pytest.mark.parametrize("seed", [0, 1, 42, 12345, 2**31 - 1])
def test_rnorm_matches_R_short(seed: int) -> None:
    """First 20 draws match R exactly."""
    expected = _r_rnorm(seed, 20)
    rng = r_set_seed(seed)
    actual = np.array([r_norm_rand(rng) for _ in range(20)], dtype=np.float64)
    np.testing.assert_array_equal(actual, expected)


def test_rnorm_matches_R_long() -> None:
    """500 draws (crosses the 624-word regenerate boundary at ~312 draws
    since each rnorm consumes two unifs) match R."""
    expected = _r_rnorm(42, 500)
    rng = r_set_seed(42)
    actual = np.array([r_norm_rand(rng) for _ in range(500)], dtype=np.float64)
    np.testing.assert_array_equal(actual, expected)


def test_qnorm_central_region() -> None:
    """AS 241's |q| <= 0.425 branch matches R exactly."""
    ps = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    expected = _r_qnorm_at(ps)
    actual = np.array([r_qnorm(p) for p in ps], dtype=np.float64)
    np.testing.assert_array_equal(actual, expected)


def test_qnorm_tail_region() -> None:
    """AS 241's r <= 5 branch (tails) matches R exactly."""
    ps = [0.01, 0.001, 0.999, 0.9999, 1e-5, 1 - 1e-5]
    expected = _r_qnorm_at(ps)
    actual = np.array([r_qnorm(p) for p in ps], dtype=np.float64)
    np.testing.assert_array_equal(actual, expected)


def test_qnorm_extreme_tail() -> None:
    """AS 241's r > 5 branch (extreme tails) matches R exactly."""
    ps = [1e-12, 1e-15, 1 - 1e-12, 1 - 1e-15]
    expected = _r_qnorm_at(ps)
    actual = np.array([r_qnorm(p) for p in ps], dtype=np.float64)
    np.testing.assert_array_equal(actual, expected)
