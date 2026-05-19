# SPDX-License-Identifier: GPL-3.0-or-later
"""Unit tests for :mod:`pylmrob.rng`.

Locks in the algorithmic structure of R's seeding and MT19937 tempering.
Ground-truth comparison against R's own output lives in
``tests/validation/test_r_rng_vs_R.py`` (requires ``rpy2``).
"""

from __future__ import annotations

import numpy as np
import pytest

from pylmrob.rng import RState, r_sample_noreplace, r_set_seed

_UINT32_MASK = 0xFFFF_FFFF


def _lcg_iter(seed: int, n: int) -> int:
    """Reference: run R's seed-scramble LCG ``n`` times."""
    s = int(seed) & _UINT32_MASK
    for _ in range(n):
        s = (69069 * s + 1) & _UINT32_MASK
    return s


def test_state_shape_and_dtype() -> None:
    rng = r_set_seed(42)
    assert rng.state.shape == (624,)
    assert rng.state.dtype == np.uint32
    assert rng.pos == 624


def test_first_state_word_matches_lcg() -> None:
    """state[0] should be the 52nd LCG iteration (51 scramble + 1 more).

    Cross-checked against R 4.2's ``.Random.seed[3]`` for seeds
    0, 1, 42, 12345, and ``2**31 - 1``.
    """
    seed = 42
    expected_first = _lcg_iter(seed, 52)
    rng = r_set_seed(seed)
    assert int(rng.state[0]) == expected_first


def test_last_state_word_matches_lcg() -> None:
    """state[623] should be the (51 + 624) = 675th LCG iteration."""
    seed = 1
    expected_last = _lcg_iter(seed, 51 + 624)
    rng = r_set_seed(seed)
    assert int(rng.state[623]) == expected_last


def test_state_matches_R_for_seed_42() -> None:
    """First, second, and last state words for ``set.seed(42)`` against R.

    Captured from ``Rscript -e 'set.seed(42); .Random.seed[c(3, 4, 626)]'``
    (R 4.2.2). Two's-complement: negative R values become large uint32.
    """
    rng = r_set_seed(42)
    assert int(rng.state[0]) == 507561766
    assert int(rng.state[1]) == 1260545903
    assert int(rng.state[623]) == 705745481


def test_reproducible() -> None:
    a = r_set_seed(42).unif_rand_n(10)
    b = r_set_seed(42).unif_rand_n(10)
    np.testing.assert_array_equal(a, b)


def test_different_seeds_differ() -> None:
    a = r_set_seed(1).unif_rand_n(5)
    b = r_set_seed(2).unif_rand_n(5)
    assert not np.array_equal(a, b)


def test_unif_rand_in_unit_interval() -> None:
    rng = r_set_seed(42)
    draws = rng.unif_rand_n(1000)
    assert (draws >= 0.0).all()
    assert (draws < 1.0).all()


def test_unif_rand_advances_pos() -> None:
    rng = r_set_seed(42)
    rng.unif_rand()
    assert rng.pos == 1
    rng.unif_rand()
    assert rng.pos == 2


def test_pos_wraps_after_624_draws() -> None:
    """After 624 calls, pos should be at 624 again (regen on the next call)."""
    rng = r_set_seed(42)
    for _ in range(624):
        rng.unif_rand()
    assert rng.pos == 624
    rng.unif_rand()
    assert rng.pos == 1


def test_negative_seed_wraps() -> None:
    """R casts seed to ``unsigned int``; ``set.seed(-1)`` == ``set.seed(2**32-1)``."""
    a = r_set_seed(-1).unif_rand_n(5)
    b = r_set_seed(2**32 - 1).unif_rand_n(5)
    np.testing.assert_array_equal(a, b)


def test_state_is_copy() -> None:
    """The ``state`` property returns a copy so callers can't mutate internals."""
    rng = r_set_seed(42)
    s1 = rng.state
    s1[0] = 0
    assert int(rng.state[0]) != 0


def test_constructor_rejects_wrong_shape() -> None:
    with pytest.raises(ValueError, match="624"):
        RState(np.zeros(100, dtype=np.uint32))


def test_constructor_rejects_wrong_dtype() -> None:
    with pytest.raises(ValueError, match="uint32"):
        RState(np.zeros(624, dtype=np.int32))


@pytest.mark.parametrize("seed", [0, 1, 42, 12345, 2**31 - 1, 2**32 - 1])
def test_unif_rand_finite(seed: int) -> None:
    """No edge-case seeds produce NaN or inf."""
    draws = r_set_seed(seed).unif_rand_n(100)
    assert np.isfinite(draws).all()


def test_sample_noreplace_is_permutation_when_k_equals_n() -> None:
    out = r_sample_noreplace(r_set_seed(42), 20, 20)
    assert out.shape == (20,)
    assert sorted(out.tolist()) == list(range(20))


def test_sample_noreplace_partial_distinct() -> None:
    """k < n: output has k distinct values from 0..n-1."""
    out = r_sample_noreplace(r_set_seed(42), 100, 10)
    assert out.shape == (10,)
    assert len(set(out.tolist())) == 10
    assert all(0 <= v < 100 for v in out.tolist())


def test_sample_noreplace_matches_robustbase_c() -> None:
    """Bit-identical to robustbase's ``sample_noreplace`` in ``lmrob.c``.

    Captured from a direct ``.C(robustbase:::R_subsample, ..., sample=TRUE)``
    call with ``set.seed(42)``, ``n = 10``: ``ind_space`` after the call.
    """
    out = r_sample_noreplace(r_set_seed(42), 10, 10)
    expected = [9, 8, 2, 5, 3, 7, 4, 0, 1, 6]
    assert out.tolist() == expected


def test_sample_noreplace_rejects_bad_k() -> None:
    rng = r_set_seed(42)
    with pytest.raises(ValueError, match="0 <= k <= n"):
        r_sample_noreplace(rng, 5, 6)
    with pytest.raises(ValueError, match="0 <= k <= n"):
        r_sample_noreplace(rng, 5, -1)


def test_sample_noreplace_consumes_exactly_k_draws() -> None:
    """``r_sample_noreplace(rng, n, k)`` calls ``unif_rand`` exactly k times."""
    rng_a = r_set_seed(42)
    r_sample_noreplace(rng_a, 50, 7)
    next_a = rng_a.unif_rand()

    rng_b = r_set_seed(42)
    for _ in range(7):
        rng_b.unif_rand()
    next_b = rng_b.unif_rand()

    assert next_a == next_b
