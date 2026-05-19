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

import importlib

import numpy as np

try:
    _cy_rng = importlib.import_module("pylmrob._core._r_rng")
except ImportError:
    _cy_rng = None

__all__ = [
    "RState",
    "r_norm_rand",
    "r_qnorm",
    "r_sample_noreplace",
    "r_set_seed",
    "r_subsample_nonsingular",
]


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
        if _cy_rng is not None:
            out, new_pos = _cy_rng.cy_r_unif_rand_n(self._state, self._pos, int(n))
            self._pos = int(new_pos)
            return out
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
    if _cy_rng is not None:
        out, new_pos = _cy_rng.cy_r_sample_noreplace(rng._state, rng._pos, n_i, k_i)
        rng._pos = int(new_pos)
        return out
    ind_space = np.arange(n_i, dtype=np.int64)
    out = np.empty(k_i, dtype=np.int64)
    nn = n_i
    for i in range(k_i):
        j = int(nn * rng.unif_rand())
        out[i] = ind_space[j]
        nn -= 1
        ind_space[j] = ind_space[nn]
    return out


def r_subsample_nonsingular(
    rng: RState,
    X: np.ndarray,
    p: int,
    mts: int = 1000,
    tol_inv: float = 1e-7,
) -> np.ndarray | None:
    """Pick a non-singular ``p``-row subset of ``X`` matching robustbase's
    ``subsample()`` with ``ss=1`` (nonsingular).

    Ports the LU-pivot-with-row-skip block from robustbase's ``lmrob.c``.
    On each attempt, draws one full permutation of the rows via
    :func:`r_sample_noreplace` and then walks through it column by column,
    extending an incremental LU factorization with partial column pivot.
    When the pivot at column ``j`` falls below ``tol_inv``, the candidate
    row is dropped and the next row from the permutation is tried; after
    ``mts`` failures the whole resample fails.

    Parameters
    ----------
    rng
        Source :class:`RState`. Advanced by ``n`` draws per attempt.
    X
        Design matrix of shape ``(n, m)``. Must be C-contiguous float64
        for the inner loops; an as-needed copy is made if not.
    p
        Number of rows to select. Must satisfy ``1 <= p <= n``.
    mts
        Maximum number of redraws before giving up.
    tol_inv
        Robustbase's ``tolInverse``; minimum |pivot| to accept a column.

    Returns
    -------
    numpy.ndarray | None
        ``int64`` array of length ``p`` with the chosen row indices, in
        the same order robustbase's ``idc`` records them. Returns
        ``None`` if no non-singular subset is found within ``mts``
        attempts.

    Notes
    -----
    The column permutation ``idr`` that robustbase tracks is purely
    internal to the LU factorization; the row selection ``idc`` is
    independent of it, and the linear system ``X[idc] @ beta = y[idc]``
    has the same solution regardless of column order.
    """
    if not isinstance(rng, RState):
        raise TypeError("rng must be an RState")
    Xa = np.ascontiguousarray(X, dtype=np.float64)
    n_i, m_i = Xa.shape
    p_i = int(p)
    if p_i < 1 or p_i > m_i:
        raise ValueError(f"need 1 <= p <= ncol(X)={m_i}, got p={p_i}")
    if p_i > n_i:
        raise ValueError(f"need p <= nrow(X)={n_i}, got p={p_i}")

    if _cy_rng is not None:
        result, new_pos = _cy_rng.cy_r_subsample_nonsingular(
            rng._state, rng._pos, Xa, p_i, int(mts), float(tol_inv)
        )
        rng._pos = int(new_pos)
        return result

    attempt = 0
    while True:
        ind_space = r_sample_noreplace(rng, n_i, n_i)
        idc = np.empty(p_i, dtype=np.int64)
        idr = np.arange(p_i, dtype=np.int64)
        lu = np.zeros((p_i, p_i), dtype=np.float64)
        v = np.empty(p_i, dtype=np.float64)
        i = 0
        failed = False
        for j in range(p_i):
            while True:
                if i + j >= n_i:
                    failed = True
                    break
                idc[j] = ind_space[i + j]
                if j == 0:
                    for k in range(p_i):
                        v[k] = Xa[idc[j], idr[k]]
                else:
                    for k in range(j):
                        lu[k, j] = Xa[idc[j], idr[k]]
                    # Forward solve L[0:j, 0:j] @ z = lu[0:j, j], result in lu[0:j, j].
                    # L is unit lower triangular (implicit diagonal of 1s).
                    for kk in range(j):
                        s = lu[kk, j]
                        for ll in range(kk):
                            s -= lu[kk, ll] * lu[ll, j]
                        lu[kk, j] = s
                    for k in range(j, p_i):
                        s = Xa[idc[j], idr[k]]
                        for ll in range(j):
                            s -= lu[k, ll] * lu[ll, j]
                        v[k] = s
                if j < p_i - 1:
                    tmpd = abs(v[j])
                    mu = j
                    for k in range(j + 1, p_i):
                        if abs(v[k]) > tmpd:
                            mu = k
                            tmpd = abs(v[k])
                    if tmpd >= tol_inv:
                        v[j], v[mu] = v[mu], v[j]
                        idr[j], idr[mu] = idr[mu], idr[j]
                        pivot = v[j]
                        for k in range(j + 1, p_i):
                            lu[k, j] = v[k] / pivot
                        if j > 0 and mu != j:
                            for k in range(j):
                                lu[j, k], lu[mu, k] = lu[mu, k], lu[j, k]
                if abs(v[j]) < tol_inv:
                    # ss = 1 path: drop this candidate row, try the next.
                    i += 1
                    continue
                lu[j, j] = v[j]
                break
            if failed:
                break
        if not failed:
            return idc
        attempt += 1
        if attempt >= int(mts):
            return None


def r_qnorm(p: float) -> float:
    """Standard normal quantile, byte-identical to R's ``qnorm()``.

    Ports ``qnorm5`` from R's ``src/nmath/qnorm.c`` (Wichura's AS 241,
    Applied Statistics 1988). The 7-coefficient minimax rational
    approximation is accurate to about 1 part in 10^16, which is the
    only way to match R bit-for-bit (``scipy.stats.norm.ppf`` uses a
    different algorithm and differs in the last ULP).
    """
    if not (0.0 < p < 1.0):
        if p == 0.0:
            return float("-inf")
        if p == 1.0:
            return float("inf")
        return float("nan")
    q = p - 0.5
    if abs(q) <= 0.425:
        r = 0.180625 - q * q
        val = (
            q
            * (
                (
                    (
                        (
                            (
                                (
                                    (r * 2509.0809287301226727 + 33430.575583588128105) * r
                                    + 67265.770927008700853
                                )
                                * r
                                + 45921.953931549871457
                            )
                            * r
                            + 13731.693765509461125
                        )
                        * r
                        + 1971.5909503065514427
                    )
                    * r
                    + 133.14166789178437745
                )
                * r
                + 3.387132872796366608
            )
            / (
                (
                    (
                        (
                            (
                                (
                                    (r * 5226.495278852854561 + 28729.085735721942674) * r
                                    + 39307.89580009271061
                                )
                                * r
                                + 21213.794301586595867
                            )
                            * r
                            + 5394.1960214247511077
                        )
                        * r
                        + 687.1870074920579083
                    )
                    * r
                    + 42.313330701600911252
                )
                * r
                + 1.0
            )
        )
        return val
    r = p if q <= 0 else 1.0 - p
    import math as _math

    r = _math.sqrt(-_math.log(r))
    if r <= 5.0:
        r -= 1.6
        val = (
            (
                (
                    (
                        (
                            (
                                (r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r
                                + 0.24178072517745061177
                            )
                            * r
                            + 1.27045825245236838258
                        )
                        * r
                        + 3.64784832476320460504
                    )
                    * r
                    + 5.7694972214606914055
                )
                * r
                + 4.6303378461565452959
            )
            * r
            + 1.42343711074968357734
        ) / (
            (
                (
                    (
                        (
                            (
                                (r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r
                                + 0.0151986665636164571966
                            )
                            * r
                            + 0.14810397642748007459
                        )
                        * r
                        + 0.68976733498510000455
                    )
                    * r
                    + 1.6763848301838038494
                )
                * r
                + 2.05319162663775882187
            )
            * r
            + 1.0
        )
    else:
        r -= 5.0
        val = (
            (
                (
                    (
                        (
                            (
                                (r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r
                                + 0.0012426609473880784386
                            )
                            * r
                            + 0.026532189526576123093
                        )
                        * r
                        + 0.29656057182850489123
                    )
                    * r
                    + 1.7848265399172913358
                )
                * r
                + 5.4637849111641143699
            )
            * r
            + 6.6579046435011037772
        ) / (
            (
                (
                    (
                        (
                            (
                                (r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r
                                + 1.8463183175100546818e-5
                            )
                            * r
                            + 7.868691311456132591e-4
                        )
                        * r
                        + 0.0148753612908506148525
                    )
                    * r
                    + 0.13692988092273580531
                )
                * r
                + 0.59983220655588793769
            )
            * r
            + 1.0
        )
    return -val if q < 0.0 else val


# Precision-extension constant from R's snorm.c (BIG = 2^27): combines two
# 32-bit unif_rand draws into a single >27-bit uniform before qnorm.
_NORM_BIG = 1 << 27


def r_norm_rand(rng: RState) -> float:
    """One standard normal draw, byte-identical to R's ``rnorm(1)``.

    Replicates the ``Inversion`` branch of ``norm_rand`` in R's
    ``src/nmath/snorm.c`` (R's default ``RNGkind()[2]``):

    1. ``u1 = unif_rand()``
    2. ``u1 = floor(2**27 * u1) + unif_rand()``
    3. ``return qnorm5(u1 / 2**27, 0, 1, 1, 0)``

    Two ``unif_rand`` draws per normal output.
    """
    u1 = rng.unif_rand()
    u1 = int(_NORM_BIG * u1) + rng.unif_rand()
    return r_qnorm(u1 / _NORM_BIG)
