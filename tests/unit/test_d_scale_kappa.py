# SPDX-License-Identifier: GPL-3.0-or-later
"""kappa for the D-step must reproduce R's quadrature, not the exact integral.

``robustbase:::lmrob.kappa`` solves ``E[psi(Z) Z - kappa wgt(Z)] = 0``. The
expectation comes from ``robustbase:::lmrob.E``, which applies an
``numpoints``-node Gauss-Hermite rule rather than integrating. With
``lmrob.control(numpoints = 10)`` that rule resolves the kinked families
(hampel, ggw) only roughly, so R's kappa is up to 8e-3 away from the exact
integral. Matching R means matching the *rule*.

Reference values produced with robustbase 0.99-7 / R 4.2.2::

    lk <- robustbase:::lmrob.kappa
    lk(control = lmrob.control(psi = "ggw", tuning.psi = c(-0.5, 1.5, 0.95, NA)))
"""

from __future__ import annotations

import numpy as np
import pytest

from pylmrob.d_scale import kappa

# (label, family, internal tuning, R's kappa)
_R_KAPPA = [
    ("bisquare", "bisquare", (4.685061,), 0.8280907301726164),
    (
        "hampel",
        "hampel",
        (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085),
        0.8504151399552571,
    ),
    ("optimal", "optimal", (1.060158,), 0.9361829181830277),
    ("lqq", "lqq", (1.4734061, 0.9822707, 1.5), 0.8618074039274787),
    ("ggw_b1", "ggw", (1.0,), 0.8989804359465622),
    ("ggw_b1.5", "ggw", (4.0,), 0.8590698460697594),
    # welsh had no entry in the old hardcoded Cython table, so the D-step
    # silently fell back to the NumPy path for it. Computing kappa per fit
    # means every family works, tabulated or not.
    ("welsh", "welsh", (2.11,), 0.8165683255302613),
]


@pytest.mark.parametrize("label,family,tuning,r_value", _R_KAPPA)
def test_kappa_matches_r(label: str, family: str, tuning: tuple[float, ...], r_value: float):
    """kappa agrees with R to near machine precision for every family."""
    got = kappa(family, tuning)
    assert got == pytest.approx(r_value, rel=1e-12), (
        f"{label}: kappa={got!r} R={r_value!r} rerr={abs(got - r_value) / r_value:.2e}"
    )


def test_kappa_is_not_the_exact_integral():
    """Guard the reason the Gauss-Hermite rule is used at all.

    If someone 'improves' :func:`kappa` back into an exact quadrature, the
    values move away from R by ~8e-3 on ggw. This test pins the fact that
    R's answer and the exact integral genuinely differ, so a future reader
    does not treat the difference as our bug.
    """
    quad = pytest.importorskip("scipy.integrate")
    from pylmrob import psi as _psi

    family, tuning = "ggw", (1.0,)

    def num(t: float) -> float:
        return float(_psi.psi(np.array([t]), family, tuning)[0]) * t * _phi(t)

    def den(t: float) -> float:
        return float(_psi.wgt(np.array([t]), family, tuning)[0]) * _phi(t)

    a, _ = quad.quad(num, -40.0, 40.0, limit=400)
    b, _ = quad.quad(den, -40.0, 40.0, limit=400)
    exact = a / b

    assert abs(exact - kappa(family, tuning)) / exact > 1e-3


def _phi(t: float) -> float:
    return float(np.exp(-0.5 * t * t) / np.sqrt(2.0 * np.pi))


def test_numpoints_changes_kappa():
    """``numpoints`` is part of the answer, not just its accuracy."""
    k10 = kappa("hampel", (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085), numpoints=10)
    k50 = kappa("hampel", (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085), numpoints=50)
    assert abs(k10 - k50) / k10 > 1e-4
