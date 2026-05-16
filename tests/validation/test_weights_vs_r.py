# SPDX-License-Identifier: GPL-3.0-or-later
"""Per-case weights match R's lmrob element-wise.

R's implementation applies a sqrt(w)-transform at the R level before
calling the unweighted C kernel (robustbase/R/lmrob.R:96-98). pyrobustlm
mirrors this in ``_lmrob_impl``.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from pyrobustlm import Control, lmrob

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


@pytest.fixture
def stackloss() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "stackloss.csv"))


def test_weights_none_matches_weights_ones(stackloss: pd.DataFrame) -> None:
    """``weights=None`` and ``weights=np.ones(n)`` produce byte-identical fits."""
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    fit_a = lmrob(formula, stackloss, control=Control(nResample=500), seed=1)
    fit_b = lmrob(
        formula,
        stackloss,
        weights=np.ones(len(stackloss)),
        control=Control(nResample=500),
        seed=1,
    )
    np.testing.assert_array_equal(fit_a.coef_, fit_b.coef_)
    assert fit_a.scale_ == fit_b.scale_


def test_uniform_weight_scaling(stackloss: pd.DataFrame) -> None:
    """Uniform weight scaling w -> c*w leaves the coefficients unchanged and
    scales the M-scale by sqrt(c) (sqrt(w)-transform on the design)."""
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    n = len(stackloss)
    fit_1 = lmrob(formula, stackloss, weights=np.ones(n), control=Control(nResample=500), seed=1)
    fit_2 = lmrob(
        formula, stackloss, weights=2.0 * np.ones(n), control=Control(nResample=500), seed=1
    )
    np.testing.assert_allclose(fit_2.coef_, fit_1.coef_, rtol=1e-10)
    np.testing.assert_allclose(fit_2.scale_ / fit_1.scale_, np.sqrt(2.0), rtol=1e-10)


def test_mixed_weights_match_r(stackloss: pd.DataFrame) -> None:
    """``weights=c(rep(1,10), rep(2,11))`` on stackloss matches R element-wise.

    R reference values from ``robustbase::lmrob(stack.loss ~ .,
    data=stackloss, weights=c(rep(1,10),rep(2,11)))`` with ``set.seed(1)``.
    """
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    w = np.concatenate([np.ones(10), 2.0 * np.ones(11)])
    fit = lmrob(formula, stackloss, weights=w, control=Control(nResample=500), seed=1)
    r_coef = np.array([-43.007449, 0.878068, 0.734852, -0.088959])
    np.testing.assert_allclose(fit.coef_, r_coef, atol=1e-4)
    np.testing.assert_allclose(fit.scale_, 2.340538, atol=1e-4)


def test_zero_weights_drop_rows(stackloss: pd.DataFrame) -> None:
    """Rows with ``weight=0`` are dropped (R's behaviour)."""
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    n = len(stackloss)
    w = np.ones(n)
    w[0] = 0.0
    fit = lmrob(formula, stackloss, weights=w, control=Control(nResample=500), seed=1)
    assert fit.nobs_ == n - 1
    # Should match fitting on stackloss[1:] without weights:
    fit_ref = lmrob(
        formula, stackloss.iloc[1:].reset_index(drop=True), control=Control(nResample=500), seed=1
    )
    np.testing.assert_allclose(fit.coef_, fit_ref.coef_, rtol=1e-10)
    np.testing.assert_allclose(fit.scale_, fit_ref.scale_, rtol=1e-10)


def test_weights_validation(stackloss: pd.DataFrame) -> None:
    """Bad weights raise a clear ValueError."""
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    n = len(stackloss)
    with pytest.raises(ValueError, match="non-negative"):
        lmrob(formula, stackloss, weights=-np.ones(n))
    with pytest.raises(ValueError, match="non-negative"):
        lmrob(formula, stackloss, weights=np.full(n, np.nan))
    with pytest.raises(ValueError, match="length"):
        lmrob(formula, stackloss, weights=np.ones(n - 1))


def test_residuals_on_original_scale(stackloss: pd.DataFrame) -> None:
    """``residuals_`` and ``fitted_`` are reported on the original scale
    (y - X @ coef), not on the sqrt(w)-transformed scale, matching R."""
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    w = np.concatenate([np.ones(10), 2.0 * np.ones(11)])
    fit = lmrob(formula, stackloss, weights=w, control=Control(nResample=500), seed=1)
    # residuals should match y - fitted on the original (untransformed) y
    y = stackloss["stack.loss"].to_numpy(dtype=np.float64)
    np.testing.assert_allclose(fit.residuals_, y - fit.fitted_, rtol=1e-12)
