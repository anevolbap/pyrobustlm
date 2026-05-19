# SPDX-License-Identifier: GPL-3.0-or-later
"""Bit-identical agreement of ``r_subsample_nonsingular`` with robustbase.

Drives robustbase's C ``R_subsample`` with ``sample = TRUE``, ``ss = 1``
(nonsingular) and compares the returned ``idc`` row indices to those
:func:`pylmrob.r_subsample_nonsingular` picks given the same seed and
the same ``X``.
"""

from __future__ import annotations

import shutil
import subprocess

import numpy as np
import pytest

from pylmrob import r_set_seed, r_subsample_nonsingular


def _have_r() -> bool:
    return shutil.which("Rscript") is not None


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


pytestmark = pytest.mark.skipif(not _have_robustbase(), reason="R + robustbase not available")


def _r_subsample_idc(X: np.ndarray, seed: int, p: int) -> np.ndarray:
    """Call ``.C(R_subsample, ..., sample = TRUE, ss = 1)`` after
    ``set.seed(seed)`` and return ``idc[1:p]`` (0-indexed)."""
    n, _m = X.shape
    np.savetxt("/tmp/_pylmrob_test_X.txt", X)
    script = f"""
        suppressMessages(library(robustbase))
        X <- as.matrix(read.table("/tmp/_pylmrob_test_X.txt"))
        n <- {n}L; p <- {p}L
        storage.mode(X) <- "double"
        y <- as.double(seq_len(n))
        set.seed({seed})
        res <- .C(robustbase:::R_subsample,
            x = X, y = y, n = as.integer(n), m = p,
            beta = double(p), ind_space = integer(n), idc = integer(n),
            idr = integer(n), lu = matrix(double(1), p, p), v = double(p),
            pivot = integer(p - 1L), Dr = double(n), Dc = double(p),
            rowequ = integer(1L), colequ = integer(1L), status = integer(1L),
            sample = TRUE, mts = 1000L, ss = 1L, tolinv = 1e-7, solve = TRUE)
        if (res$status != 0) stop(paste("R_subsample failed status=", res$status))
        cat(res$idc[1:p], sep = "\\n")
    """
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True,
        check=True,
        text=True,
        timeout=20,
    )
    return np.array([int(x) for x in out.stdout.splitlines() if x], dtype=np.int64)


@pytest.mark.parametrize("seed", [1, 42, 12345])
def test_well_conditioned_idc_matches_R(seed: int) -> None:
    rng = np.random.default_rng(0)
    X = rng.standard_normal((25, 4))
    r_idc = _r_subsample_idc(X, seed, 4)
    py_idc = r_subsample_nonsingular(r_set_seed(seed), X, 4)
    assert py_idc is not None
    np.testing.assert_array_equal(py_idc, r_idc)


def test_collinear_block_skipped_same_as_R() -> None:
    """Force the LU row-skip path: rows 0..9 are rank-deficient. Both R
    and pylmrob should walk past them in the same way."""
    rng = np.random.default_rng(0)
    X = rng.standard_normal((30, 4))
    X[:10, 3] = 2 * X[:10, 0]
    r_idc = _r_subsample_idc(X, 42, 4)
    py_idc = r_subsample_nonsingular(r_set_seed(42), X, 4)
    assert py_idc is not None
    np.testing.assert_array_equal(py_idc, r_idc)


@pytest.mark.parametrize("p", [3, 4, 5])
def test_varying_p_matches_R(p: int) -> None:
    rng = np.random.default_rng(0)
    X = rng.standard_normal((40, p + 2))
    r_idc = _r_subsample_idc(X, 7, p)
    py_idc = r_subsample_nonsingular(r_set_seed(7), X, p)
    assert py_idc is not None
    np.testing.assert_array_equal(py_idc, r_idc)
