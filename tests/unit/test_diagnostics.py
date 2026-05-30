# SPDX-License-Identifier: GPL-3.0-or-later
"""Shape, sign, and outlier-flagging tests for pylmrob.diagnostics.

Covers ``hatvalues``, ``dfbetas``, and ``cooks_distance`` on a small
synthetic design with a planted Y-outlier so the outlier-flagging
assertions are non-trivial.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob
from pylmrob.diagnostics import cooks_distance, dfbetas, hatvalues


@pytest.fixture(scope="module")
def fit_with_outlier():
    """Fit on a clean Gaussian design with one obvious Y-outlier at row 0."""
    rng = np.random.default_rng(0)
    n, p = 60, 3
    X = rng.standard_normal((n, p))
    beta = np.array([1.0, -0.5, 2.0])
    y = X @ beta + 0.5 * rng.standard_normal(n)
    y[0] += 30.0  # planted vertical outlier

    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(p)])
    df["y"] = y
    fit = lmrob("y ~ x0 + x1 + x2", df, control=Control(nResample=200), seed=42)
    # Design matrix actually used by the fit (incl. intercept column).
    design = np.column_stack([np.ones(n), X])
    return fit, design, n, p + 1  # p+1 because of intercept


def test_hatvalues_shape_and_range(fit_with_outlier):
    fit, X, n, _ = fit_with_outlier
    h = hatvalues(fit, X)
    assert h.shape == (n,)
    assert (h >= 0.0).all()
    assert (h <= 1.0).all()


def test_hatvalues_trace_approximates_p(fit_with_outlier):
    """sum(h_i) approximates p (the rank of the weighted design)."""
    fit, X, _, p = fit_with_outlier
    h = hatvalues(fit, X)
    # Robust hat is computed on sqrt(w) * X, so its trace is <= p; on a
    # mostly-clean fit (one outlier downweighted) it should be very near p.
    assert h.sum() == pytest.approx(p, abs=0.5)


def test_dfbetas_shape(fit_with_outlier):
    fit, X, n, p = fit_with_outlier
    d = dfbetas(fit, X)
    assert d.shape == (n, p)


def test_dfbetas_is_finite(fit_with_outlier):
    fit, X, _, _ = fit_with_outlier
    d = dfbetas(fit, X)
    assert np.isfinite(d).all()


def test_dfbetas_downweights_outlier_row(fit_with_outlier):
    """The outlier row carries near-zero robust weight, so its dfbetas
    influence is small (the robust fit ignores it). This is the *robust*
    flavor of dfbetas, not the classical one."""
    fit, X, _, _ = fit_with_outlier
    d = dfbetas(fit, X)
    # Median absolute dfbetas across all rows (a typical-row scale).
    typical = float(np.median(np.abs(d)))
    outlier_norm = float(np.linalg.norm(d[0]))
    # The outlier row's contribution should be at or below the typical row
    # scale because the M-estimator downweighted it (rweights_[0] ~= 0).
    assert outlier_norm <= max(typical * 5.0, 1e-6), (
        f"outlier row should not dominate dfbetas; got {outlier_norm:.3g}, typical {typical:.3g}"
    )


def test_cooks_distance_shape(fit_with_outlier):
    fit, X, n, _ = fit_with_outlier
    d = cooks_distance(fit, X)
    assert d.shape == (n,)


def test_cooks_distance_nonnegative(fit_with_outlier):
    fit, X, _, _ = fit_with_outlier
    d = cooks_distance(fit, X)
    assert (d >= 0.0).all()
    assert np.isfinite(d).all()


def test_cooks_distance_robust_false_raises(fit_with_outlier):
    fit, X, _, _ = fit_with_outlier
    with pytest.raises(NotImplementedError):
        cooks_distance(fit, X, robust=False)


def test_diagnostics_match_results_method(fit_with_outlier):
    """fit.diagnostics() should produce the same hat values as the
    standalone hatvalues() helper, since the method delegates."""
    fit, X, _, _ = fit_with_outlier
    h_direct = hatvalues(fit, X)
    diag = fit.diagnostics()
    # The Diagnostics dataclass carries leverage; for the unweighted form
    # both should agree to numerical precision.
    np.testing.assert_allclose(h_direct, diag.leverage, atol=1e-10)
