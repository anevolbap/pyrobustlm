# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 5: M-S init for designs with categorical predictors."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from pyrobustlm import Control, lmrob
from pyrobustlm.formula import model_matrix
from pyrobustlm.ms_estimator import _l1_fit, m_s_fit

REPO_ROOT = Path(__file__).resolve().parents[2]


def _education_df() -> pd.DataFrame:
    out = REPO_ROOT / "tests" / "data" / "education.csv"
    if not out.exists():
        subprocess.run(
            [
                "Rscript",
                "-e",
                f"library(robustbase); data(education); write.csv(education, '{out}', row.names=FALSE)",
            ],
            capture_output=True,
            check=True,
        )
    df = pd.read_csv(out)
    df["Region"] = pd.Categorical(df["Region"], categories=[1, 2, 3, 4])
    return df


def test_factor_columns_detected():
    """formulaic encodes categorical columns as 'Var[T.level]'; we tag those."""
    df = _education_df()
    design = model_matrix("Y ~ Region + X1 + X2 + X3", df)
    # Expect 3 factor cols (Region has 4 levels, drop reference -> 3 contrasts).
    assert design.is_factor_col.sum() == 3
    # Plus intercept (not factor) and 3 numeric -> 7 total.
    assert design.X.shape[1] == 7


def test_l1_fit_matches_lstsq_on_clean_data():
    """L1 ~ OLS on clean Gaussian data."""
    rng = np.random.default_rng(0)
    n, p = 200, 4
    X = rng.standard_normal((n, p))
    beta_true = rng.standard_normal(p)
    y = X @ beta_true + rng.standard_normal(n) * 0.1
    b_l1 = _l1_fit(X, y)
    b_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
    np.testing.assert_allclose(b_l1, b_ols, atol=0.05)


def test_l1_fit_robust_to_outliers():
    """L1 should be far less affected by outliers than OLS."""
    rng = np.random.default_rng(1)
    n, p = 100, 3
    X = rng.standard_normal((n, p))
    beta_true = np.array([1.0, -0.5, 2.0])
    y = X @ beta_true + rng.standard_normal(n) * 0.1
    # Inject 10 wild outliers
    y[:10] += 50.0
    b_l1 = _l1_fit(X, y)
    b_ols, *_ = np.linalg.lstsq(X, y, rcond=None)
    err_l1 = np.linalg.norm(b_l1 - beta_true)
    err_ols = np.linalg.norm(b_ols - beta_true)
    assert err_l1 < err_ols  # L1 is more robust


def test_m_s_fit_basic():
    """M-S should produce a sensible coefficient vector on a factor design."""
    df = _education_df()
    design = model_matrix("Y ~ Region + X1 + X2 + X3", df)
    is_cat = design.is_factor_col
    X_cat = design.X[:, is_cat]
    X_cont = design.X[:, ~is_cat]
    y = design.y

    res = m_s_fit(
        X_cat=X_cat,
        X_cont=X_cont,
        y=y,
        nResample=100,
        k_m_s=15,
        seed=0,
    )
    # All coefficients finite
    assert np.isfinite(res.coef).all()
    # Scale is positive
    assert res.scale > 0


def test_lmrob_with_init_M_S_factor_design():
    """End-to-end: lmrob(init="M-S") on a categorical+continuous design."""
    df = _education_df()
    fit = lmrob(
        "Y ~ Region + X1 + X2 + X3",
        df,
        control=Control(init="M-S", nResample=100, k_m_s=15),
        seed=0,
    )
    # MM should converge
    assert fit.converged_
    assert fit.scale_ > 0
    # init metadata
    assert fit.init_["method"] == "M-S"


def test_init_auto_picks_M_S_for_factor_designs():
    df = _education_df()
    fit = lmrob(
        "Y ~ Region + X1 + X2 + X3",
        df,
        control=Control(init="auto", nResample=100, k_m_s=10),
        seed=0,
    )
    assert fit.init_["method"] == "M-S"


def test_init_auto_picks_S_for_continuous_designs():
    """Without factors, init='auto' should pick S."""
    rng = np.random.default_rng(0)
    n, p = 50, 3
    X = rng.standard_normal((n, p))
    df = pd.DataFrame(X, columns=["x0", "x1", "x2"])
    df["y"] = X @ np.array([1.0, -1.0, 2.0]) + rng.standard_normal(n)
    fit = lmrob(
        "y ~ x0 + x1 + x2",
        df,
        control=Control(init="auto", nResample=100),
        seed=0,
    )
    assert fit.init_["method"] == "S"
