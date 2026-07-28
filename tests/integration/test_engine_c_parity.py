# SPDX-License-Identifier: GPL-3.0-or-later
"""Whole-fit comparison between ``engine_c=True`` and the NumPy path.

The two engines draw their p-subsets from different streams (the Cython
kernel uses an internal Floyd draw, the NumPy path uses
``Generator.choice``), so they can land in different basins of attraction
and a strict element-wise assertion over the whole corpus would be
testing the RNG, not the code. Arithmetic parity is covered separately and
without any RNG by ``test_kernel_parity.py``.

What this file pins is the part that must always hold:

* every fit on the corpus completes, converges, and reports a finite,
  strictly positive scale. This is the regression guard for the kernel
  bug where a degenerate zero-scale candidate won the best-of-``best_r``
  comparison, after which ``cy_lmrob_fit`` returned success with an
  unwritten ``beta_init`` buffer and the caller read uninitialised memory.
* the two engines agree element-wise on the large majority of fits, so
  basin drift stays the exception rather than the rule.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "tests" / "data"

_DATASETS = [
    ("stackloss", "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."),
    ("coleman", "Y ~ salaryP + fatherWc + sstatus + teacherSc + motherLev"),
    ("salinity", "Y ~ X1 + X2 + X3"),
    ("delivery", "delTime ~ n.prod + distance"),
    ("phosphor", "plant ~ inorg + organic"),
    ("aircraft", "Y ~ X1 + X2 + X3 + X4"),
    ("pension", "Reserves ~ Income"),
    ("starsCYG", "log.light ~ log.Te"),
    ("hbk", "Y ~ X1 + X2 + X3"),
    ("wood", "y ~ x1 + x2 + x3 + x4 + x5"),
]
_SEEDS = (1, 6, 15, 42, 123)

# Measured agreement on the reference machine: 98/100 fits agree to
# <1e-8 over 10 seeds. Two (salinity/seed 6, hbk/seed 1) land in a
# different basin, where the NumPy path is the one that matches R.
# Keep the floor a little below the measured rate so ordinary
# platform-to-platform BLAS noise does not turn this red.
_MIN_AGREEMENT = 0.85


def _load(dataset: str) -> pd.DataFrame:
    path = DATA_DIR / f"{dataset}.csv"
    if not path.exists():
        pytest.skip(f"data file missing: {path}")
    return pd.read_csv(path)


@pytest.mark.parametrize("dataset,formula", _DATASETS)
@pytest.mark.parametrize("seed", _SEEDS)
def test_engine_c_produces_a_valid_fit(dataset: str, formula: str, seed: int) -> None:
    """The default engine must never report success with a degenerate fit.

    No ``pytest.skip`` on exception here on purpose: a raise *is* the
    failure this test exists to catch.
    """
    fit = lmrob(formula, _load(dataset), control=Control(nResample=500), seed=seed)

    assert fit.converged_, f"{dataset}/seed={seed}: did not converge"
    assert np.isfinite(fit.scale_), f"{dataset}/seed={seed}: scale={fit.scale_}"
    assert fit.scale_ > 0.0, f"{dataset}/seed={seed}: non-positive scale {fit.scale_}"
    assert np.all(np.isfinite(fit.coef_)), f"{dataset}/seed={seed}: coef={fit.coef_}"
    assert np.all(np.isfinite(fit.cov_)), f"{dataset}/seed={seed}: non-finite cov"
    # The uninitialised-buffer bug showed up as absurd magnitudes
    # (~1e241) in the init-S coefficients that feed vcov_avar1.
    init_coef = np.asarray(fit.init_.get("coef", fit.coef_), dtype=np.float64)
    assert np.all(np.abs(init_coef) < 1e100), (
        f"{dataset}/seed={seed}: init coef looks uninitialised: {init_coef}"
    )


def test_engines_agree_on_most_fits() -> None:
    """Element-wise agreement on the large majority of the corpus."""
    agree = 0
    total = 0
    drifted: list[tuple[str, int, float]] = []

    for dataset, formula in _DATASETS:
        df = _load(dataset)
        for seed in _SEEDS:
            total += 1
            fit_c = lmrob(formula, df, control=Control(nResample=500, engine_c=True), seed=seed)
            fit_np = lmrob(formula, df, control=Control(nResample=500, engine_c=False), seed=seed)
            coef_err = float(
                np.max(np.abs(fit_c.coef_ - fit_np.coef_) / np.maximum(np.abs(fit_np.coef_), 1.0))
            )
            scale_err = abs(fit_c.scale_ - fit_np.scale_) / fit_np.scale_
            err = max(coef_err, scale_err)
            if err < 1e-8:
                agree += 1
            else:
                drifted.append((dataset, seed, err))

    rate = agree / total
    assert rate >= _MIN_AGREEMENT, (
        f"engines agreed on only {agree}/{total} fits ({rate:.0%}); "
        f"drifted cases: {sorted(drifted, key=lambda t: -t[2])[:5]}"
    )


@pytest.mark.parametrize("dataset,formula", _DATASETS)
def test_init_coef_is_the_s_estimate_on_both_engines(dataset: str, formula: str) -> None:
    """``init_["coef"]`` must be the initial S estimate, not the MM one.

    The engine_c branch used to store ``beta_out``, the buffer the Cython
    MM loop mutates in place, so ``init_["coef"]`` held the post-MM
    coefficients there and the post-S coefficients on the NumPy path.
    The tell: on the engine_c path it was equal to ``coef_``.
    """
    df = _load(dataset)
    fit_c = lmrob(formula, df, control=Control(nResample=500, engine_c=True), seed=42)
    fit_np = lmrob(formula, df, control=Control(nResample=500, engine_c=False), seed=42)

    for label, fit in (("engine_c", fit_c), ("numpy", fit_np)):
        init_coef = np.asarray(fit.init_["coef"], dtype=np.float64)
        assert init_coef.shape == fit.coef_.shape
        # MM moves off the S estimate, so the two must not be identical.
        assert not np.array_equal(init_coef, fit.coef_), (
            f"{dataset}/{label}: init_['coef'] is the MM estimate, not the S estimate"
        )
        # and it must not alias the reported coefficients
        assert init_coef is not fit.coef_


def test_init_scale_matches_init_coef() -> None:
    """The recorded init scale must be the M-scale at the recorded init coef.

    Cross-checks the two halves of ``init_`` against each other, which is
    what exposed the aliasing above.
    """
    from pylmrob.formula import model_matrix
    from pylmrob.scale import m_scale

    formula = "Y ~ X1 + X2 + X3"
    df = _load("salinity")
    design = model_matrix(formula, df)
    X = np.ascontiguousarray(design.X, dtype=np.float64)
    y = np.ascontiguousarray(design.y, dtype=np.float64)

    for engine_c in (True, False):
        fit = lmrob(formula, df, control=Control(nResample=500, engine_c=engine_c), seed=42)
        init_coef = np.asarray(fit.init_["coef"], dtype=np.float64)
        recomputed = m_scale(
            y - X @ init_coef, "bisquare", (1.54764,), 0.5, 200, 1e-10, p=X.shape[1]
        )
        np.testing.assert_allclose(
            recomputed,
            float(fit.init_["scale"]),  # type: ignore[arg-type]
            rtol=1e-6,
            err_msg=f"engine_c={engine_c}: init scale does not match init coef",
        )
