# SPDX-License-Identifier: GPL-3.0-or-later
"""Bit-identical agreement with R's ``unif_rand`` for the same seed.

Runs ``Rscript`` in a subprocess to capture R's actual output and compares
to :func:`pylmrob.r_set_seed`. Each test skips cleanly if ``Rscript`` is
not on ``PATH``.
"""

from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pytest

from pylmrob import r_set_seed


def _have_r() -> bool:
    return shutil.which("Rscript") is not None


pytestmark = pytest.mark.skipif(not _have_r(), reason="Rscript not on PATH")


def _r_unif_rand(seed: int, n: int) -> np.ndarray:
    """Draw ``n`` uniforms from R after ``set.seed(seed)``.

    Uses ``sprintf("%.17g")`` so the printed text round-trips losslessly
    back to ``float64``.
    """
    script = f'set.seed({seed}); cat(sprintf("%.17g\\n", runif({n})), sep="")'
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True,
        check=True,
        text=True,
        timeout=20,
    )
    lines = [line for line in out.stdout.splitlines() if line]
    return np.array([float(line) for line in lines], dtype=np.float64)


@pytest.mark.parametrize("seed", [0, 1, 42, 12345, 2**31 - 1])
def test_unif_rand_matches_R_short(seed: int) -> None:
    """First 20 draws match R exactly."""
    expected = _r_unif_rand(seed, 20)
    actual = r_set_seed(seed).unif_rand_n(20)
    np.testing.assert_array_equal(actual, expected)


def test_unif_rand_matches_R_across_regenerate() -> None:
    """1000 draws (crosses the 624-word regenerate boundary) match R."""
    seed = 42
    expected = _r_unif_rand(seed, 1000)
    actual = r_set_seed(seed).unif_rand_n(1000)
    np.testing.assert_array_equal(actual, expected)


def test_state_matches_R_for_seed_42() -> None:
    """Initial state words after ``set.seed(42)`` match R's ``.Random.seed``."""
    out = subprocess.run(
        ["Rscript", "-e", "set.seed(42); cat(.Random.seed[3:626], sep='\\n')"],
        capture_output=True,
        check=True,
        text=True,
        timeout=10,
    )
    # R stores state as signed Int32; cast to uint32 via two's complement.
    r_state = np.array(
        [int(line) & 0xFFFFFFFF for line in out.stdout.splitlines() if line],
        dtype=np.uint32,
    )
    assert r_state.shape == (624,)
    rng = r_set_seed(42)
    np.testing.assert_array_equal(rng.state, r_state)
