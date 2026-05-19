# SPDX-License-Identifier: GPL-3.0-or-later
"""R-compatible MT19937 RNG.

Replicates R's ``set.seed`` / ``unif_rand`` path from ``src/main/RNG.c`` so a
caller can produce the exact same uniform stream R produces from the same
integer seed.

This module is the first step toward bit-identical agreement with R's
``lmrob`` fits. The other pieces still required for byte-identical fits:

- The fast-S resample loop must call ``unif_rand`` in the same order R's
  ``lmrob.c`` does (subset draw, restart, ...). ``pylmrob``'s NumPy /
  Cython kernels currently don't match that order.
- Some R code paths use ``norm_rand`` or ``runif`` with extra rejection
  steps; those aren't covered here.

See ``docs/faq.md`` for the longer caveat list.
"""

from __future__ import annotations

import numpy as np

__all__ = ["RState", "r_sample_noreplace", "r_set_seed"]


_N = 624
_M = 397
_MATRIX_A = 0x9908_B0DF
_UPPER_MASK = 0x8000_0000
_LOWER_MASK = 0x7FFF_FFFF
_UINT32_MASK = 0xFFFF_FFFF


class RState:
    """MT19937 RNG with R's seeding and one-draw-per-``unif_rand`` convention.

    Construct via :func:`r_set_seed`. Then call :meth:`unif_rand` to pull
    uniforms; each call consumes one 32-bit MT word and converts it via
    ``y * 2**-32``, matching R's ``RNG.c``.
    """

    __slots__ = ("_pos", "_state")

    def __init__(self, state: np.ndarray, pos: int = _N) -> None:
        if state.shape != (_N,) or state.dtype != np.uint32:
            raise ValueError("state must be a (624,) uint32 ndarray")
        self._state = state.copy()
        self._pos = int(pos)

    @property
    def state(self) -> np.ndarray:
        """Return a copy of the 624-word MT state (R's ``.Random.seed[2:]``)."""
        return self._state.copy()

    @property
    def pos(self) -> int:
        """Current index into the 624-word state (R's ``mti``)."""
        return self._pos

    def _regenerate(self) -> None:
        s = self._state
        for i in range(_N - _M):
            y = (int(s[i]) & _UPPER_MASK) | (int(s[i + 1]) & _LOWER_MASK)
            s[i] = (int(s[i + _M]) ^ (y >> 1) ^ (_MATRIX_A if (y & 1) else 0)) & _UINT32_MASK
        for i in range(_N - _M, _N - 1):
            y = (int(s[i]) & _UPPER_MASK) | (int(s[i + 1]) & _LOWER_MASK)
            s[i] = (int(s[i + _M - _N]) ^ (y >> 1) ^ (_MATRIX_A if (y & 1) else 0)) & _UINT32_MASK
        y = (int(s[_N - 1]) & _UPPER_MASK) | (int(s[0]) & _LOWER_MASK)
        s[_N - 1] = (int(s[_M - 1]) ^ (y >> 1) ^ (_MATRIX_A if (y & 1) else 0)) & _UINT32_MASK
        self._pos = 0

    def next_uint32(self) -> int:
        """Return the next 32-bit MT output (tempered)."""
        if self._pos >= _N:
            self._regenerate()
        y = int(self._state[self._pos])
        self._pos += 1
        y ^= y >> 11
        y ^= (y << 7) & 0x9D2C_5680
        y ^= (y << 15) & 0xEFC6_0000
        y ^= y >> 18
        return y & _UINT32_MASK

    def unif_rand(self) -> float:
        """One uniform in ``[0, 1)``, matching R's ``unif_rand()``.

        Pulls one 32-bit MT word and converts via ``y * 2**-32``. The fixup
        from R (which avoids exact 0 and 1) is applied; in practice it
        almost never triggers.
        """
        y = self.next_uint32()
        u = y * 2.0**-32
        if u <= 0.0:
            return 0.5 * 2.0**-32
        if u >= 1.0:
            return 1.0 - 0.5 * 2.0**-32
        return u

    def unif_rand_n(self, n: int) -> np.ndarray:
        """Vectorised wrapper around :meth:`unif_rand`."""
        out = np.empty(int(n), dtype=np.float64)
        for i in range(int(n)):
            out[i] = self.unif_rand()
        return out


def r_set_seed(seed: int) -> RState:
    """Build an :class:`RState` from ``seed`` via R's ``set.seed`` path.

    The path (from R's ``RNG.c``):

    1. Treat ``seed`` as a 32-bit unsigned int.
    2. Run an LCG scramble: ``seed = (69069 * seed + 1) mod 2**32``.
       Empirically R does 51 of these before the first state word
       (verified against ``.Random.seed`` for several integer seeds);
       this is one more than the inline comment "Initial scrambling"
       block in ``RNG.c`` suggests, presumably because of an extra
       advance somewhere else in the path.
    3. Run 624 more LCG iterations; each becomes one state word.
    4. Position is set to 624 so the first draw triggers regeneration.

    Parameters
    ----------
    seed
        Integer seed. Negative values wrap modulo ``2**32`` to match R's
        ``unsigned int`` cast.
    """
    s = int(seed) & _UINT32_MASK
    for _ in range(51):
        s = (69069 * s + 1) & _UINT32_MASK
    state = np.empty(_N, dtype=np.uint32)
    for j in range(_N):
        s = (69069 * s + 1) & _UINT32_MASK
        state[j] = s
    return RState(state, pos=_N)


def r_sample_noreplace(rng: RState, n: int, k: int) -> np.ndarray:
    """Draw ``k`` distinct indices from ``0..n-1`` in robustbase's order.

    Ports ``sample_noreplace`` from robustbase's ``src/lmrob.c``: a
    Knuth-style swap-and-replace that pulls one :meth:`RState.unif_rand`
    draw per output index. Fast-S calls this with ``k = n`` to get a
    full permutation; smaller ``k`` is supported for completeness.

    Parameters
    ----------
    rng
        Source :class:`RState`. Advanced by ``k`` draws.
    n
        Population size; indices range over ``0..n-1``.
    k
        Number of indices to draw. Must satisfy ``0 <= k <= n``.

    Returns
    -------
    numpy.ndarray
        ``int64`` array of length ``k``.

    Notes
    -----
    The reference C is::

        for (i = 0; i < n; i++) ind_space[i] = i;
        for (i = 0; i < k; i++) {
            j = nn * unif_rand();
            x[i] = ind_space[j];
            ind_space[j] = ind_space[--nn];
        }
    """
    n_i = int(n)
    k_i = int(k)
    if k_i < 0 or k_i > n_i:
        raise ValueError(f"need 0 <= k <= n, got n={n_i}, k={k_i}")
    ind_space = np.arange(n_i, dtype=np.int64)
    out = np.empty(k_i, dtype=np.int64)
    nn = n_i
    for i in range(k_i):
        j = int(nn * rng.unif_rand())
        out[i] = ind_space[j]
        nn -= 1
        ind_space[j] = ind_space[nn]
    return out
