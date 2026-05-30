# SPDX-License-Identifier: GPL-3.0-or-later
"""Algebraic invariants of the psi/rho/chi/weight functions.

Six redescending families plus huber. For each, check the basic
algebraic properties an M-estimator loss must satisfy:

- rho(0) = 0, psi(0) = 0, wgt(0) = 1
- rho is even, psi is odd, wgt is even
- wgt(u) = psi(u) / u for u != 0 (definition)
- redescending families have rho bounded and psi(|u| -> inf) -> 0
- huber is non-redescending: psi saturates at +/- k

These are not curve-fitting tests; they're sanity checks that catch any
implementation drift or copy/paste error in the dispatch tables.
"""

from __future__ import annotations

import hypothesis.strategies as st
import numpy as np
import pytest
from hypothesis import assume, given, settings

from pylmrob import psi as P

# Tuning constants taken directly from pylmrob.control (robustbase
# defaults). M-step targets 95% Gaussian efficiency.
TUNING_M = {
    "bisquare": (4.685061,),
    "welsh":    (2.11,),
    "optimal":  (1.060158,),
    "hampel":   (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085),
    "lqq":      (1.4734061, 0.9826779, 1.5),
    "ggw":      (1.0, 1.5, 0.5, 1.694),
    "huber":    (1.345,),
}
REDESCENDING = [fam for fam in TUNING_M if fam != "huber"]
ALL_FAMILIES = list(TUNING_M)


# A grid of test points dense near the origin and the bend, sparse far
# out. Hypothesis adds random samples on top.
GRID = np.concatenate(
    [
        np.linspace(-10.0, -2.0, 9),
        np.linspace(-2.0, 2.0, 41),
        np.linspace(2.0, 10.0, 9),
    ]
)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_rho_at_zero_is_zero(family):
    """Loss is zero at the origin for every family."""
    k = TUNING_M[family]
    assert P.rho(0.0, family, k) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_psi_at_zero_is_zero(family):
    """Influence function is zero at the origin for every family."""
    k = TUNING_M[family]
    assert P.psi(0.0, family, k) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_wgt_at_zero_is_one(family):
    """IRWLS weight at the origin is one (limit of psi(u)/u as u -> 0)."""
    k = TUNING_M[family]
    assert P.wgt(0.0, family, k) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_rho_is_even(family):
    """rho(-u) = rho(u)."""
    k = TUNING_M[family]
    np.testing.assert_allclose(
        P.rho(GRID, family, k), P.rho(-GRID, family, k), atol=1e-12
    )


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_psi_is_odd(family):
    """psi(-u) = -psi(u)."""
    k = TUNING_M[family]
    np.testing.assert_allclose(
        P.psi(GRID, family, k), -P.psi(-GRID, family, k), atol=1e-12
    )


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_wgt_is_even(family):
    """wgt(-u) = wgt(u)."""
    k = TUNING_M[family]
    np.testing.assert_allclose(
        P.wgt(GRID, family, k), P.wgt(-GRID, family, k), atol=1e-12
    )


@pytest.mark.parametrize("family", ALL_FAMILIES)
@settings(max_examples=50, deadline=None)
@given(u=st.floats(min_value=-50.0, max_value=50.0, allow_nan=False, allow_infinity=False))
def test_wgt_equals_psi_over_u(family, u):
    """wgt(u) = psi(u) / u for u sufficiently away from zero.

    Both sides are computed independently in pylmrob; if they disagree,
    one of them is buggy.
    """
    # Near-zero is the limit case, covered by test_wgt_at_zero_is_one.
    assume(abs(u) >= 1e-6)
    k = TUNING_M[family]
    w = float(P.wgt(u, family, k))
    p = float(P.psi(u, family, k))
    assert w == pytest.approx(p / u, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("family", REDESCENDING)
def test_rho_is_bounded(family):
    """Redescending rho/chi are normalised: rho(inf) = 1."""
    k = TUNING_M[family]
    far = np.array([20.0, 50.0, 100.0, 1000.0])
    vals = P.rho(far, family, k)
    # Allow a small tolerance because the asymptote is approached
    # exponentially in some families (ggw).
    assert (vals <= 1.0 + 1e-9).all()
    # And at very large |u| we should be at the asymptote.
    assert P.rho(1000.0, family, k) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("family", REDESCENDING)
def test_psi_redescends(family):
    """psi(|u| -> large) -> 0 for redescending families."""
    k = TUNING_M[family]
    far = np.array([20.0, 50.0, 100.0])
    np.testing.assert_allclose(P.psi(far, family, k), 0.0, atol=1e-9)


def test_huber_psi_saturates():
    """huber is the counter-example: psi plateaus at +/- k instead of zero."""
    (k,) = TUNING_M["huber"]
    far = np.array([10.0, 100.0, 1000.0])
    np.testing.assert_allclose(P.psi(far, "huber", (k,)), k, atol=1e-12)
    np.testing.assert_allclose(P.psi(-far, "huber", (k,)), -k, atol=1e-12)


@pytest.mark.parametrize("family", ALL_FAMILIES)
def test_rho_monotone_nondecreasing_on_positive(family):
    """rho is monotonically non-decreasing on u >= 0."""
    k = TUNING_M[family]
    u_pos = np.linspace(0.0, 20.0, 201)
    rho_vals = P.rho(u_pos, family, k)
    diffs = np.diff(rho_vals)
    # Allow a tiny numerical wobble at the asymptote.
    assert (diffs >= -1e-12).all(), f"rho not monotone on positive u for {family}"
