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

from pylmrob import r_sample_noreplace, r_set_seed


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


def _have_robustbase() -> bool:
    if not _have_r():
        return False
    out = subprocess.run(
        ["Rscript", "-e", 'cat(requireNamespace("robustbase", quietly=TRUE))'],
        capture_output=True,
        text=True,
        timeout=10,
    )
    return out.stdout.strip() == "TRUE"


@pytest.mark.skipif(not _have_robustbase(), reason="robustbase not installed in R")
@pytest.mark.parametrize("seed", [1, 42, 12345])
@pytest.mark.parametrize("n", [10, 25, 100])
def test_sample_noreplace_matches_robustbase_c(seed: int, n: int) -> None:
    """``r_sample_noreplace(seed, n, n)`` matches the permutation
    robustbase's C ``sample_noreplace`` produces after ``set.seed(seed)``.

    Driven through ``.C(robustbase:::R_subsample, ..., sample = TRUE)``,
    whose first step is ``sample_noreplace(ind_space, n, n, idc)``.
    """
    script = f"""
        suppressMessages(library(robustbase))
        n <- {n}; p <- 3L
        set.seed({seed})
        X <- matrix(seq_len(n * p) * 1.0, n, p)
        storage.mode(X) <- "double"
        y <- as.double(seq_len(n))
        set.seed({seed})
        res <- .C(robustbase:::R_subsample,
            x = X, y = y, n = as.integer(n), m = p,
            beta = double(p), ind_space = integer(n), idc = integer(n),
            idr = integer(n), lu = matrix(double(1), p, p), v = double(p),
            pivot = integer(p - 1L), Dr = double(n), Dc = double(p),
            rowequ = integer(1L), colequ = integer(1L), status = integer(1L),
            sample = TRUE, mts = 0L, ss = 1L, tolinv = 1e-7, solve = TRUE)
        cat(res$ind_space, sep = "\\n")
    """
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True,
        check=True,
        text=True,
        timeout=20,
    )
    expected = np.array(
        [int(line) for line in out.stdout.splitlines() if line],
        dtype=np.int64,
    )
    actual = r_sample_noreplace(r_set_seed(seed), n, n)
    np.testing.assert_array_equal(actual, expected)
