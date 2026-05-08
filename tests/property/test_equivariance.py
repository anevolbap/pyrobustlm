# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 10 property tests: regression equivariance.

Three theorems an MM-estimator must satisfy:

1. **Translation equivariance** in y: ``lmrob(y + c, X)`` shifts only the
   intercept by c, leaving slopes unchanged.
2. **Scale equivariance** in y: ``lmrob(s * y, X)`` scales both beta and
   sigma by ``s``.
3. **Linear equivariance** in y wrt X: ``lmrob(y + X @ c, X)`` shifts beta
   by ``c`` and leaves sigma unchanged.

We don't test affine equivariance in X (more delicate; covered separately).
"""

from __future__ import annotations

import hypothesis.strategies as st
import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings

from pyrobustlm import Control, lmrob


def _toy_dataset(n: int = 80, p: int = 3, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, p))
    beta = np.array([1.0, -0.5, 2.0])
    eps = rng.standard_normal(n)
    y = X @ beta + eps
    # 10% gross outliers
    bad = rng.choice(n, size=n // 10, replace=False)
    y[bad] += rng.uniform(15, 25, bad.size) * np.where(rng.random(bad.size) < 0.5, -1, 1)
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(p)])
    df["y"] = y
    return df


@pytest.fixture(scope="module")
def base_fit():
    df = _toy_dataset()
    fit = lmrob("y ~ x0 + x1 + x2", df, control=Control(nResample=200), seed=7)
    return df, fit


@settings(
    max_examples=10, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(c=st.floats(min_value=-10.0, max_value=10.0, allow_nan=False, allow_infinity=False))
def test_translation_equivariance_intercept(base_fit, c):
    """y -> y + c shifts intercept by c, slopes unchanged."""
    df, fit = base_fit
    df2 = df.copy()
    df2["y"] = df2["y"] + c
    fit2 = lmrob("y ~ x0 + x1 + x2", df2, control=Control(nResample=200), seed=7)

    # Intercept shifts by c
    np.testing.assert_allclose(
        fit2.coef_[0],
        fit.coef_[0] + c,
        rtol=1e-6,
        atol=1e-6,
    )
    # Slopes unchanged
    np.testing.assert_allclose(fit2.coef_[1:], fit.coef_[1:], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(fit2.scale_, fit.scale_, rtol=1e-6, atol=1e-6)


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(s=st.floats(min_value=0.5, max_value=3.0, allow_nan=False, allow_infinity=False))
def test_scale_equivariance(base_fit, s):
    """y -> s*y scales both beta and sigma by s."""
    df, fit = base_fit
    df2 = df.copy()
    df2["y"] = df2["y"] * s
    fit2 = lmrob("y ~ x0 + x1 + x2", df2, control=Control(nResample=200), seed=7)
    np.testing.assert_allclose(fit2.coef_, fit.coef_ * s, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(fit2.scale_, fit.scale_ * abs(s), rtol=1e-5, atol=1e-5)


@settings(
    max_examples=8, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    c=st.lists(
        st.floats(min_value=-3.0, max_value=3.0, allow_nan=False, allow_infinity=False),
        min_size=3,
        max_size=3,
    )
)
def test_regression_equivariance(base_fit, c):
    """y -> y + X @ c shifts only the slopes by c (intercept unchanged because c has no
    intercept component, X already has columns x0 x1 x2)."""
    df, fit = base_fit
    c_arr = np.array(c)
    df2 = df.copy()
    X = df[["x0", "x1", "x2"]].to_numpy()
    df2["y"] = df2["y"] + X @ c_arr
    fit2 = lmrob("y ~ x0 + x1 + x2", df2, control=Control(nResample=200), seed=7)
    # Slopes shift by c, intercept and sigma stay
    np.testing.assert_allclose(fit2.coef_[1:], fit.coef_[1:] + c_arr, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(fit2.coef_[0], fit.coef_[0], rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(fit2.scale_, fit.scale_, rtol=1e-5, atol=1e-5)
