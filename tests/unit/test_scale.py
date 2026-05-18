# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 3: M-scale validated against robustbase::lmrob.mscale."""

from __future__ import annotations

import numpy as np
import pytest

from pylmrob.scale import _mad, m_scale


@pytest.mark.parametrize("seed", [0, 1, 2, 7, 42])
def test_mscale_clean_normal(r_session, seed):
    rng = np.random.default_rng(seed)
    r = rng.standard_normal(200)
    py = m_scale(r, family="bisquare")
    rv = r_session.lmrob_mscale(r)
    np.testing.assert_allclose(py, rv, rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize("seed,contam", [(0, 0.1), (1, 0.2), (5, 0.3)])
def test_mscale_with_outliers(r_session, seed, contam):
    rng = np.random.default_rng(seed)
    n = 200
    r = rng.standard_normal(n)
    n_bad = int(np.ceil(n * contam))
    r[:n_bad] += rng.uniform(15, 30, n_bad) * np.where(rng.random(n_bad) < 0.5, -1, 1)
    py = m_scale(r, family="bisquare")
    rv = r_session.lmrob_mscale(r)
    np.testing.assert_allclose(py, rv, rtol=1e-9, atol=1e-12)


def test_mscale_short_vector(r_session):
    """Small n is the regime that exercises (n - p) division most."""
    rng = np.random.default_rng(11)
    for n in (10, 21, 30):
        r = rng.standard_normal(n)
        py = m_scale(r, family="bisquare")
        rv = r_session.lmrob_mscale(r)
        np.testing.assert_allclose(py, rv, rtol=1e-9, atol=1e-12)


def test_mad_helper():
    x = np.array([-3.0, -1.0, 0.0, 1.0, 3.0])
    # Median = 0; abs-deviations: [3,1,0,1,3]; median = 1; *1.4826 = 1.4826
    assert abs(_mad(x) - 1.4826) < 1e-12


def test_mscale_zero_residuals_returns_zero():
    s = m_scale(np.zeros(20), family="bisquare")
    assert s == 0.0


def test_mscale_nonfinite_raises():
    with pytest.raises(ValueError, match="must be finite"):
        m_scale(np.array([1.0, np.nan, 3.0]))
