# SPDX-License-Identifier: GPL-3.0-or-later
"""``Control(rng="MT19937")`` opts in to the Mersenne Twister BitGenerator
(R's RNG family) instead of NumPy's default PCG64.

Strict byte-identical agreement with R's ``lmrob`` fits is not promised:
R's ``set.seed`` scrambles an integer seed through Marsaglia's PRNG to
fill 624 MT state words, and ``pylmrob`` doesn't replicate that
scrambling. The MT19937 option still gives users a way to use the same
RNG family as R and is a starting point for tighter agreement.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob, make_generator

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


@pytest.fixture
def stackloss() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "stackloss.csv"))


def test_make_generator_pcg64_default() -> None:
    rng = make_generator(42)
    assert isinstance(rng.bit_generator, np.random.PCG64)


def test_make_generator_mt19937() -> None:
    rng = make_generator(42, kind="MT19937")
    assert isinstance(rng.bit_generator, np.random.MT19937)


def test_make_generator_unknown_kind() -> None:
    with pytest.raises(ValueError, match="Unknown rng kind"):
        make_generator(42, kind="XORSHIFT")


def test_make_generator_passes_through_existing_generator() -> None:
    """If the user passes an existing Generator, the function returns it
    unchanged when the kind matches its BitGenerator family."""
    g = np.random.default_rng(42)
    assert make_generator(g, kind="PCG64") is g
    g2 = np.random.Generator(np.random.MT19937(42))
    assert make_generator(g2, kind="MT19937") is g2


def test_lmrob_mt19937_produces_finite_fit(stackloss: pd.DataFrame) -> None:
    """``Control(rng='MT19937')`` runs end-to-end without errors and
    produces a finite, sensible fit."""
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(rng="MT19937", nResample=200),
        seed=1,
    )
    assert fit.converged_
    assert np.isfinite(fit.scale_)
    assert np.isfinite(fit.coef_).all()


def test_lmrob_mt19937_different_basin_than_pcg64(stackloss: pd.DataFrame) -> None:
    """PCG64 and MT19937 with the same integer seed produce different
    resampling sequences (different seed -> state algorithms), so the
    fits land in slightly different basins.

    Both fits should be near R's reference, but they won't be
    byte-identical to each other.
    """
    fit_pcg = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(rng="PCG64", nResample=200),
        seed=1,
    )
    fit_mt = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(rng="MT19937", nResample=200),
        seed=1,
    )
    # Same data, different RNG seed-scrambling -> different basin.
    # Numerically close (within rtol=1e-2) but not identical.
    np.testing.assert_allclose(fit_mt.coef_, fit_pcg.coef_, rtol=1e-2, atol=1e-1)


def test_control_rng_default_is_pcg64() -> None:
    """The default Control() preserves the PCG64 backend for backward
    compatibility with v0.5.x fits."""
    assert Control().rng == "PCG64"
