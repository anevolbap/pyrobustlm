# SPDX-License-Identifier: GPL-3.0-or-later
"""The compiled psi kernels and the NumPy reference must agree.

:mod:`pylmrob.psi` now prefers ``pylmrob._core._psi`` and falls back to
:mod:`pylmrob._psifuns`. That makes the NumPy code a path most users
never take -- and therefore one the suite would stop exercising if
nothing asked for it explicitly. It is still the fallback for ``welsh``
(no compiled kernel), for ``ggw`` chi (a tabulated polynomial that lives
in ``_psifuns``), and for any source build without the extension, so it
has to stay covered and has to stay equal.

This is the same shape of test that caught the D-step kappa bug: two
implementations of one thing, fed identical inputs, compared directly.
The whole-fit tests cannot see a divergence here because both paths
produce plausible numbers.
"""

from __future__ import annotations

import numpy as np
import pytest

from pylmrob import _cykernels, _psifuns
from pylmrob.control import _DEFAULT_TUNING_CHI, _DEFAULT_TUNING_PSI

# A grid dense enough to land inside every branch: the psi families are
# piecewise, so the interesting behaviour is at the knots.
_X = np.linspace(-8.0, 8.0, 2001)

_KINDS = ("rho", "psi", "psi_prime", "wgt")
_FAMILIES = sorted(set(_DEFAULT_TUNING_PSI) & set(_DEFAULT_TUNING_CHI))

_CASES = [(fam, kind, label) for fam in _FAMILIES for kind in _KINDS for label in ("psi", "chi")]


def _tuning(family: str, which: str) -> tuple[float, ...]:
    table = _DEFAULT_TUNING_PSI if which == "psi" else _DEFAULT_TUNING_CHI
    return tuple(np.atleast_1d(np.asarray(table[family], dtype=float)).ravel())


@pytest.mark.parametrize(
    "family,kind,tuning_kind", _CASES, ids=[f"{f}-{k}-{t}" for f, k, t in _CASES]
)
def test_cython_matches_numpy(family: str, kind: str, tuning_kind: str) -> None:
    """Where a compiled kernel exists, it must equal the NumPy reference."""
    k = _tuning(family, tuning_kind)
    cy = _cykernels.evaluate(kind, _X, family, k)
    if cy is None:
        pytest.skip(f"no compiled kernel for {family}/{kind}")

    ref = _psifuns._dispatch(family, kind)(_X, np.asarray(k, dtype=float))
    np.testing.assert_allclose(
        cy,
        ref,
        rtol=1e-12,
        atol=1e-14,
        err_msg=f"{family}/{kind} ({tuning_kind} tuning): Cython and NumPy disagree",
    )


@pytest.mark.parametrize("family", _FAMILIES)
@pytest.mark.parametrize("kind", _KINDS)
def test_numpy_reference_is_finite(family: str, kind: str) -> None:
    """Exercise the NumPy path directly for every family.

    Without this the fallback goes uncovered the moment the compiled
    kernels are preferred, which is how a fallback quietly rots.
    """
    for which in ("psi", "chi"):
        k = _tuning(family, which)
        out = _psifuns._dispatch(family, kind)(_X, np.asarray(k, dtype=float))
        assert out.shape == _X.shape
        assert np.all(np.isfinite(out)), f"{family}/{kind} ({which}) produced non-finite values"


def test_public_api_agrees_with_the_numpy_reference() -> None:
    """``pylmrob.psi.*`` must not change answers by preferring Cython."""
    from pylmrob import psi as psi_mod

    for family in _FAMILIES:
        k = _tuning(family, "psi")
        for kind, fn in (
            ("rho", psi_mod.rho),
            ("psi", psi_mod.psi),
            ("psi_prime", psi_mod.psi_prime),
            ("wgt", psi_mod.wgt),
        ):
            got = fn(_X, family, k)
            ref = _psifuns._dispatch(family, kind)(_X, np.asarray(k, dtype=float))
            np.testing.assert_allclose(
                got,
                ref,
                rtol=1e-12,
                atol=1e-14,
                err_msg=f"pylmrob.psi.{kind}({family!r}) diverged from the NumPy reference",
            )


def test_families_without_a_kernel_fall_back() -> None:
    """Pin the two documented gaps, so adding a kernel updates this test.

    welsh has no compiled kernel at all; ggw has no compiled ``rho``
    because its chi is a tabulated polynomial in ``_psifuns``.
    """
    assert _cykernels.evaluate("psi", _X, "welsh", _tuning("welsh", "psi")) is None
    assert _cykernels.evaluate("rho", _X, "ggw", _tuning("ggw", "chi")) is None
