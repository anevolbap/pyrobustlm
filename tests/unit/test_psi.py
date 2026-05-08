# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 2: psi/chi/wgt families validated against R's Mpsi/Mchi/Mwgt.

For each family we check rho, psi, psi', and wgt at a 1001-point grid
against R, at strict tolerances (rtol=1e-12).
"""

from __future__ import annotations

import numpy as np
import pytest

from pyrobustlm import psi as P

GRID = np.linspace(-6.0, 6.0, 1001)

# (family_name, R_family_name, default_tuning_constants)
CASES_SINGLE_K = [
    ("huber", "huber", 1.345),
    ("bisquare", "bisquare", 4.685061),
    ("optimal", "optimal", 1.060158),
]

CASES_HAMPEL = [
    ("hampel", "hampel", (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085)),
]

CASES_LQQ = [
    ("lqq", "lqq", (1.4734061, 0.9826779, 1.5)),
]

# GGW tests pass user-facing form to R and the corresponding internal
# (case_idx,) form to Python. The mapping is:
#   case 1: c(-0.5, 1, 0.95, NA)   b=1,   95% eff
#   case 2: c(-0.5, 1, 0.85, NA)   b=1,   85% eff
#   case 3: c(-0.5, 1, NA, 0.5)    b=1,   bp=0.5
#   case 4: c(-0.5, 1.5, 0.95, NA) b=1.5, 95% eff
#   case 5: c(-0.5, 1.5, 0.85, NA) b=1.5, 85% eff
#   case 6: c(-0.5, 1.5, NA, 0.5)  b=1.5, bp=0.5
CASES_GGW = [
    ("case1", (-0.5, 1.0, 0.95, np.nan), (1,)),
    ("case4", (-0.5, 1.5, 0.95, np.nan), (4,)),
    ("case3", (-0.5, 1.0, np.nan, 0.5), (3,)),
    ("case6", (-0.5, 1.5, np.nan, 0.5), (6,)),
]


def _flatten_k(k):
    arr = np.atleast_1d(np.asarray(k, dtype=float))
    return tuple(float(v) for v in arr)


@pytest.mark.parametrize("py_family,r_family,k", CASES_SINGLE_K + CASES_HAMPEL + CASES_LQQ)
def test_psi_matches_r(r_session, py_family, r_family, k):
    """Mpsi: psi(x; k) should match R's Mpsi to 1e-12."""
    py = P.psi(GRID, py_family, k)
    r = r_session.Mpsi(GRID, _flatten_k(k), r_family)
    np.testing.assert_allclose(py, r, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("py_family,r_family,k", CASES_SINGLE_K + CASES_HAMPEL + CASES_LQQ)
def test_chi_matches_r(r_session, py_family, r_family, k):
    """Mchi: rho(x; k) should match R's Mchi to 1e-12."""
    py = P.rho(GRID, py_family, k)
    r = r_session.Mchi(GRID, _flatten_k(k), r_family)
    np.testing.assert_allclose(py, r, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("py_family,r_family,k", CASES_SINGLE_K + CASES_HAMPEL + CASES_LQQ)
def test_wgt_matches_r(r_session, py_family, r_family, k):
    """Mwgt: wgt(x; k) should match R's Mwgt to 1e-12."""
    py = P.wgt(GRID, py_family, k)
    r = r_session.Mwgt(GRID, _flatten_k(k), r_family)
    np.testing.assert_allclose(py, r, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("py_family,r_family,k", CASES_SINGLE_K + CASES_HAMPEL + CASES_LQQ)
def test_psi_prime_matches_r(r_session, py_family, r_family, k):
    """Mpsi(deriv=1): psi'(x; k) should match R to 1e-12 except at the
    discontinuities (the breakpoints in hampel)."""
    py = P.psi_prime(GRID, py_family, k)
    r = r_session.Mpsi(GRID, _flatten_k(k), r_family, deriv=1)

    if py_family == "hampel":
        # psi' is discontinuous at +-a, +-b, +-r. Drop a small window around
        # each break before comparing.
        breaks = np.array([k[0], k[1], k[2]])
        mask = np.ones_like(GRID, dtype=bool)
        for bk in breaks:
            mask &= np.abs(np.abs(GRID) - bk) > 1e-3
        np.testing.assert_allclose(py[mask], r[mask], rtol=1e-12, atol=1e-14)
    else:
        np.testing.assert_allclose(py, r, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("name,k_user,k_py", CASES_GGW)
def test_ggw_psi_matches_r(r_session, name, k_user, k_py):
    """ggw uses tabulated polynomial coefficients. Diff against R."""
    py = P.psi(GRID, "ggw", k_py)
    r = r_session.Mpsi(GRID, k_user, "ggw")
    np.testing.assert_allclose(py, r, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("name,k_user,k_py", CASES_GGW)
def test_ggw_wgt_matches_r(r_session, name, k_user, k_py):
    py = P.wgt(GRID, "ggw", k_py)
    r = r_session.Mwgt(GRID, k_user, "ggw")
    np.testing.assert_allclose(py, r, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("name,k_user,k_py", CASES_GGW)
def test_ggw_chi_matches_r(r_session, name, k_user, k_py):
    py = P.rho(GRID, "ggw", k_py)
    r = r_session.Mchi(GRID, k_user, "ggw")
    # ggw rho uses a polynomial approximation; cross-zero region has the
    # most discretisation error.
    np.testing.assert_allclose(py, r, rtol=1e-7, atol=1e-9)


def test_default_tuning_constants():
    """Default tuning constants exposed by psi.tuning_for_efficiency match
    the values pinned in psi._PSI_TUNING_DEFAULT_PSI."""
    bisq = P.tuning_for_efficiency("bisquare", 0.95)
    assert bisq == (4.685061,)
    huber = P.tuning_for_efficiency("huber", 0.95)
    assert huber == (1.345,)


def test_dispatch_unknown_family_raises():
    with pytest.raises(ValueError, match="unknown psi family"):
        P.psi(np.array([0.0]), "doesnotexist", 1.0)
