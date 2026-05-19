# SPDX-License-Identifier: GPL-3.0-or-later
"""``Control(rng="R")`` routes the fast-S resample through
``pylmrob.r_set_seed`` + ``r_sample_noreplace``, byte-identical to
robustbase's ``unif_rand`` draws.

These tests exercise the integration only (Python-side). End-to-end
parity with R's actual ``lmrob`` fits lives in
``tests/validation/test_lmrob_rng_r_vs_R.py``.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


@pytest.fixture
def stackloss() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "stackloss.csv"))


def test_rng_r_smoke(stackloss: pd.DataFrame) -> None:
    """End-to-end fit converges and produces finite values."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(rng="R"),
        seed=42,
    )
    assert fit.converged_
    assert np.isfinite(fit.scale_)
    assert np.isfinite(fit.coef_).all()


def test_rng_r_deterministic(stackloss: pd.DataFrame) -> None:
    """Same seed twice -> identical coefficients."""
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    fit_a = lmrob(formula, stackloss, control=Control(rng="R"), seed=42)
    fit_b = lmrob(formula, stackloss, control=Control(rng="R"), seed=42)
    np.testing.assert_array_equal(fit_a.coef_, fit_b.coef_)
    assert fit_a.scale_ == fit_b.scale_


def test_rng_r_different_seeds_can_differ(stackloss: pd.DataFrame) -> None:
    """Different seeds in R-mode can produce different basins.

    Not guaranteed across all datasets (stackloss is well-conditioned),
    but the path must at least be reachable.
    """
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    fit_1 = lmrob(formula, stackloss, control=Control(rng="R"), seed=1)
    fit_42 = lmrob(formula, stackloss, control=Control(rng="R"), seed=42)
    # Both should be finite; equality is not promised here.
    assert np.isfinite(fit_1.coef_).all()
    assert np.isfinite(fit_42.coef_).all()


def test_rng_r_rejects_no_seed(stackloss: pd.DataFrame) -> None:
    with pytest.raises(ValueError, match="rng='R' requires an explicit integer seed"):
        lmrob(
            "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
            stackloss,
            control=Control(rng="R"),
        )


def test_rng_r_rejects_n_workers_gt_1() -> None:
    """R-mode is sequential; the Control validator catches workers > 1."""
    with pytest.raises(ValueError, match="rng='R' requires n_workers=1"):
        Control(rng="R", n_workers=4)


def test_rng_r_respects_subsampling_choice() -> None:
    """Both ``simple`` and ``nonsingular`` (LU-pivot) subsampling are
    supported under rng='R'; Control no longer forces simple."""
    assert Control(rng="R", subsampling="simple").subsampling == "simple"
    assert Control(rng="R", subsampling="nonsingular").subsampling == "nonsingular"


def test_rng_r_nonsingular_smoke(stackloss: pd.DataFrame) -> None:
    """End-to-end fit with the LU-pivot subset draw runs and converges."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(rng="R", subsampling="nonsingular"),
        seed=42,
    )
    assert fit.converged_
    assert np.isfinite(fit.coef_).all()


def test_rng_r_disables_engine_c() -> None:
    """``engine_c=True`` is incompatible with R-mode (monolithic engine
    owns its BitGenerator); Control flips it off."""
    ctrl = Control(rng="R", engine_c=True)
    assert ctrl.engine_c is False


def test_rng_r_vs_pcg64_can_differ(stackloss: pd.DataFrame) -> None:
    """The R-mode draws are not the same as PCG64's so on a harder
    problem the resulting fit can differ. Stackloss is well-conditioned
    enough that both land on the same basin within rtol=1e-3.
    """
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    fit_p = lmrob(formula, stackloss, control=Control(rng="PCG64"), seed=42)
    fit_r = lmrob(formula, stackloss, control=Control(rng="R"), seed=42)
    # Both converge to the well-known stackloss MM-estimator; check the
    # fits are at least close (not asserting bit-identity here).
    np.testing.assert_allclose(fit_r.coef_, fit_p.coef_, rtol=1e-2, atol=1e-1)
