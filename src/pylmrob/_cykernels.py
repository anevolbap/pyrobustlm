# SPDX-License-Identifier: GPL-3.0-or-later
"""One dispatcher from ``(kind, family)`` to the compiled psi kernels.

``pylmrob._core._psi`` has a complete ``{family}_{rho,psi,psi_prime,wgt}``
set, but nothing called most of it: :mod:`pylmrob.psi` dispatched
straight to the NumPy reference in :mod:`pylmrob._psifuns`, and its
module docstring still said the Cython version "will be mirrored ... in
a future". Line coverage on ``_psi.pyx`` sat at 55% for that reason --
roughly half the module was compiled, shipped and unreachable.

This module is the missing link. It is deliberately the *only* place
that maps a family name to a compiled kernel: ``scale.py`` grew its own
copy of that mapping, and this codebase has been bitten repeatedly by
the same table living in several files and drifting (the ggw case index
in ``_DEFAULT_K_CHI``, the lqq constant in ``inference.py``, the D-step
kappa table in ``_lmrob.pyx``).

Every entry point returns ``None`` rather than raising when the compiled
module is absent or the family is not supported, so callers keep their
NumPy fallback and a source build without the extension still works.
"""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

# GGW (case_idx 1..6) -> (a, b, c). Mirrors SET_ABC_GGW in
# robustbase/src/lmrob.c and the table in _psi_kernels.pxi.
_GGW_ABC: dict[int, tuple[float, float, float]] = {
    1: (0.648, 1.0, 1.694),
    2: (0.4760508, 1.0, 1.2442567),
    3: (0.1674046, 1.0, 0.4375470),
    4: (1.387, 1.5, 1.063),
    5: (0.8372485, 1.5, 0.7593544),
    6: (0.2036741, 1.5, 0.2959132),
}

# Families whose kernels take a single tuning constant.
_ONE_ARG = ("bisquare", "biweight", "huber", "optimal")
# Families whose kernels take three.
_THREE_ARG = ("hampel", "lqq")

_KINDS = ("rho", "psi", "psi_prime", "wgt")


def _load() -> Any | None:
    try:
        return importlib.import_module("pylmrob._core._psi")
    except ImportError:  # pragma: no cover - source build without the extension
        return None


_CPSI = _load()


def _kernel_name(family: str, kind: str) -> str:
    """``_psi.pyx`` names its chi function ``rho`` for every family."""
    fam = "bisquare" if family == "biweight" else family
    return f"{fam}_{kind}"


def evaluate(
    kind: str,
    x: np.ndarray,
    family: str,
    k: ArrayLike,
) -> np.ndarray | None:
    """Evaluate ``kind`` for ``family`` via the compiled kernel.

    ``kind`` is one of ``rho``, ``psi``, ``psi_prime``, ``wgt``. Returns
    ``None`` when there is no compiled kernel for this combination, in
    which case the caller should use the NumPy path.

    ``welsh`` has no Cython kernel (it was added to ``_psifuns`` and the
    shared ``.pxi`` but never to ``_psi.pyx``), and ``ggw`` has no
    compiled ``rho``: its chi is a tabulated polynomial that lives in
    ``_psifuns``. Both fall through.
    """
    if _CPSI is None or kind not in _KINDS:
        return None

    fam = family.lower()
    x_buf = np.ascontiguousarray(x, dtype=np.float64)
    if x_buf.ndim != 1:
        return None
    out = np.empty_like(x_buf)
    k_arr = np.atleast_1d(np.asarray(k, dtype=np.float64)).ravel()

    if fam in _ONE_ARG:
        fn = getattr(_CPSI, _kernel_name(fam, kind), None)
        if fn is None or k_arr.size < 1:
            return None
        fn(x_buf, float(k_arr[0]), out)
        return out

    if fam in _THREE_ARG:
        fn = getattr(_CPSI, _kernel_name(fam, kind), None)
        if fn is None or k_arr.size < 3:
            return None
        fn(x_buf, float(k_arr[0]), float(k_arr[1]), float(k_arr[2]), out)
        return out

    if fam == "ggw":
        # No compiled ggw_rho: the chi is a tabulated polynomial in
        # _psifuns, so rho falls back by design.
        fn = getattr(_CPSI, f"ggw_{kind}", None)
        if fn is None:
            return None
        case_idx = int(k_arr[0]) if k_arr.size else -1
        if 1 <= case_idx <= 6:
            a, b, c = _GGW_ABC[case_idx]
        elif case_idx == 0 and k_arr.size >= 4:
            a, b, c = float(k_arr[1]), float(k_arr[2]), float(k_arr[3])
        else:
            return None
        fn(x_buf, a, b, c, out)
        return out

    return None


__all__ = ["evaluate"]
