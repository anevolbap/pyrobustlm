# SPDX-License-Identifier: GPL-3.0-or-later
"""Bootstrap inference tests."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, bootstrap, lmrob

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


@pytest.fixture
def stackloss() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "stackloss.csv"))


def test_confint_method_bootstrap(stackloss: pd.DataFrame) -> None:
    """``fit.confint(method='bootstrap')`` returns the percentile CIs
    from the underlying bootstrap object."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(),
        seed=42,
    )
    ci = fit.confint(method="bootstrap", n_boot=200, seed=42)
    assert ci.shape == (len(fit.coef_), 2)
    # Each interval brackets the point estimate.
    assert (ci[:, 0] <= fit.coef_).all()
    assert (fit.coef_ <= ci[:, 1]).all()


def test_confint_method_bootstrap_basic_differs_from_percentile(
    stackloss: pd.DataFrame,
) -> None:
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(),
        seed=42,
    )
    pct = fit.confint(method="bootstrap", n_boot=200, seed=42, kind="percentile")
    basic = fit.confint(method="bootstrap", n_boot=200, seed=42, kind="basic")
    # Basic CI is the reflection of percentile around 2 * theta_hat;
    # they shouldn't be equal in general.
    assert not np.allclose(pct, basic)


def test_confint_invalid_method_raises(stackloss: pd.DataFrame) -> None:
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(),
        seed=42,
    )
    with pytest.raises(ValueError, match="must be 'wald' or 'bootstrap'"):
        fit.confint(method="bca")


def test_bootstrap_returns_result(stackloss: pd.DataFrame) -> None:
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(nResample=200),
        seed=1,
    )
    boot = bootstrap(fit, n_boot=200, seed=42)
    p = fit.coef_.size
    assert boot.coefs.shape == (boot.n_converged, p)
    assert boot.percentile_ci.shape == (p, 2)
    assert boot.basic_ci.shape == (p, 2)
    assert boot.se.shape == (p,)
    assert boot.bias.shape == (p,)
    assert boot.n_converged > 0
    assert boot.n_converged <= 200
    assert boot.level == 0.95
    assert boot.term_names == fit.term_names_


def test_bootstrap_ci_brackets_estimate(stackloss: pd.DataFrame) -> None:
    """The bootstrap percentile CI should contain the point estimate for
    most coefficients (basic sanity check on the resampling)."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(nResample=200),
        seed=1,
    )
    boot = bootstrap(fit, n_boot=300, seed=42)
    inside = (boot.percentile_ci[:, 0] <= fit.coef_) & (fit.coef_ <= boot.percentile_ci[:, 1])
    # On a 95% CI we expect at least 3 of 4 coefs to bracket the estimate.
    assert inside.sum() >= 3


def test_bootstrap_deterministic_under_seed(stackloss: pd.DataFrame) -> None:
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(nResample=200),
        seed=1,
    )
    b1 = bootstrap(fit, n_boot=50, seed=123)
    b2 = bootstrap(fit, n_boot=50, seed=123)
    np.testing.assert_array_equal(b1.coefs, b2.coefs)
    np.testing.assert_array_equal(b1.percentile_ci, b2.percentile_ci)


def test_bootstrap_no_design_raises() -> None:
    """Bootstrap requires the stashed design matrix."""
    from pylmrob.results import LmRobResults

    fake = LmRobResults(
        coef_=np.zeros(3),
        scale_=1.0,
        weights_=np.ones(10),
        rweights_=np.ones(10),
        residuals_=np.zeros(10),
        fitted_=np.zeros(10),
        cov_=np.eye(3),
        df_residual_=7,
        converged_=True,
        n_iter_=5,
        nobs_=10,
        term_names_=["a", "b", "c"],
        control=Control(),
        design_x_=None,
        design_y_=None,
    )
    with pytest.raises(RuntimeError, match="design matrix"):
        bootstrap(fake, n_boot=5)


def test_fit_bootstrap_method_matches_function(stackloss: pd.DataFrame) -> None:
    """``fit.bootstrap(...)`` and ``bootstrap(fit, ...)`` produce the same result."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(nResample=200),
        seed=1,
    )
    method_result = fit.bootstrap(n_boot=50, seed=99)
    fn_result = bootstrap(fit, n_boot=50, seed=99)
    np.testing.assert_array_equal(method_result.coefs, fn_result.coefs)
    np.testing.assert_array_equal(method_result.percentile_ci, fn_result.percentile_ci)


def test_bootstrap_se_consistent_with_quantiles(stackloss: pd.DataFrame) -> None:
    """The reported ``se`` matches numpy.std of the coef draws (ddof=1)."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(nResample=200),
        seed=1,
    )
    boot = bootstrap(fit, n_boot=100, seed=42)
    np.testing.assert_allclose(boot.se, np.std(boot.coefs, axis=0, ddof=1), rtol=1e-12)
