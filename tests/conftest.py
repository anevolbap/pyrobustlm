# SPDX-License-Identifier: GPL-3.0-or-later
"""Shared pytest fixtures for the pylmrob test suite.

Phase 1 fills in the rpy2 bridge and the JSON reference loader.
For Phase 0 we only set up the basics (RNG, tolerances) and a graceful
skip when rpy2/R are unavailable.
"""

from __future__ import annotations

import json
import os
import warnings
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = REPO_ROOT / "tests" / "reference"
DATA_DIR = REPO_ROOT / "tests" / "data"


# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
# Default rtol/atol per output kind, as agreed in plan.md section 5.1.
# Tightening or loosening must be justified in a comment in the test.
DEFAULT_TOLERANCES: dict[str, dict[str, float]] = {
    "psi": {"rtol": 1e-12, "atol": 1e-14},
    "chi": {"rtol": 1e-12, "atol": 1e-14},
    "wgt": {"rtol": 1e-12, "atol": 1e-14},
    "mscale": {"rtol": 1e-9, "atol": 1e-12},
    "init_S_beta": {"rtol": 1e-4, "atol": 1e-6},
    "init_S_sigma": {"rtol": 1e-5, "atol": 1e-8},
    "mm_beta": {"rtol": 1e-6, "atol": 1e-8},
    "mm_sigma": {"rtol": 1e-8, "atol": 1e-10},
    "weights": {"rtol": 1e-6, "atol": 1e-8},
    "cov": {"rtol": 1e-6, "atol": 1e-8},
}


@pytest.fixture(scope="session")
def tol() -> dict[str, dict[str, float]]:
    return DEFAULT_TOLERANCES


# ---------------------------------------------------------------------------
# Deterministic RNG
# ---------------------------------------------------------------------------
@pytest.fixture
def rng() -> np.random.Generator:
    """Per-test deterministic RNG."""
    return np.random.default_rng(seed=20260508)


# ---------------------------------------------------------------------------
# rpy2 bridge
# ---------------------------------------------------------------------------
def _rpy2_available() -> bool:
    if os.environ.get("PYROBUSTLM_SKIP_RPY2"):
        return False
    try:
        import rpy2.robjects  # noqa: F401
    except Exception:
        return False
    return True


class _RBridge:
    """Convenience wrapper around an rpy2 session.

    Exposes the R primitives that pylmrob tests need to diff against:

    - ``Mpsi``, ``Mchi``, ``Mwgt`` (+ ``deriv`` argument for higher derivatives)
    - ``lmrob_mscale``  → ``robustbase::lmrob.mscale``
    - ``lmrob_control`` → ``robustbase::lmrob.control``
    - ``rdataset(name)`` → fetch a ``robustbase`` (or base ``datasets``) data frame
      as a pandas DataFrame.
    """

    def __init__(self) -> None:
        import numpy as np
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri
        from rpy2.robjects.packages import importr

        try:
            self.robustbase = importr("robustbase")
        except Exception as exc:  # pragma: no cover - surfaced via skip
            pytest.skip(f"R package 'robustbase' not available: {exc}")

        self.ro = ro
        self.r = ro.r
        self._numpy2ri = numpy2ri
        self._np = np

    # ----- numeric primitives ----------------------------------------------
    def _to_r_vec(self, x: np.ndarray | float) -> Any:
        if hasattr(x, "__len__"):
            return self.ro.FloatVector(self._np.asarray(x, dtype=float).ravel())
        return self.ro.FloatVector([float(x)])

    def Mpsi(
        self,
        x: np.ndarray | float,
        cc: float | tuple[float, ...],
        family: str,
        deriv: int = 0,
    ) -> np.ndarray:
        cc_vec = (
            self.ro.FloatVector([float(cc)])
            if isinstance(cc, (int, float))
            else self.ro.FloatVector([float(c) for c in cc])
        )
        out = self.r["Mpsi"](self._to_r_vec(x), cc=cc_vec, psi=family, deriv=deriv)
        return self._np.asarray(out)

    def Mchi(
        self,
        x: np.ndarray | float,
        cc: float | tuple[float, ...],
        family: str,
        deriv: int = 0,
    ) -> np.ndarray:
        cc_vec = (
            self.ro.FloatVector([float(cc)])
            if isinstance(cc, (int, float))
            else self.ro.FloatVector([float(c) for c in cc])
        )
        out = self.r["Mchi"](self._to_r_vec(x), cc=cc_vec, psi=family, deriv=deriv)
        return self._np.asarray(out)

    def Mwgt(
        self,
        x: np.ndarray | float,
        cc: float | tuple[float, ...],
        family: str,
    ) -> np.ndarray:
        cc_vec = (
            self.ro.FloatVector([float(cc)])
            if isinstance(cc, (int, float))
            else self.ro.FloatVector([float(c) for c in cc])
        )
        out = self.r["Mwgt"](self._to_r_vec(x), cc=cc_vec, psi=family)
        return self._np.asarray(out)

    # ----- M-scale ---------------------------------------------------------
    def lmrob_mscale(
        self,
        r: np.ndarray,
        family: str = "bisquare",
        k: float | tuple[float, ...] = 1.54764,  # 50% bdp default
        b0: float = 0.5,
        max_iter: int = 200,
        tol: float = 1e-10,
        p: int = 0,
    ) -> float:
        """Compute the M-scale via robustbase::Mchi inside R, mirroring
        ``find_scale`` in lmrob.c. ``lmrob.mscale`` is not exported by the
        package, so we replicate the iteration here.
        """
        cc_vec = (
            self.ro.FloatVector([float(k)])
            if isinstance(k, (int, float))
            else self.ro.FloatVector([float(c) for c in k])
        )
        rscript = """
        function(r, cc, family, b0, max_it, tol, p) {
            s <- mad(r)
            if (s == 0) return(0)
            for (it in 1:max_it) {
                chi <- robustbase::Mchi(r/s, cc, family)
                s_new <- s * sqrt(sum(chi) / (length(r) - p) / b0)
                if (abs(s_new - s) <= tol * s) return(s_new)
                s <- s_new
            }
            s_new
        }
        """
        fn = self.r(rscript)
        out = fn(
            self._to_r_vec(r),
            cc_vec,
            family,
            float(b0),
            int(max_iter),
            float(tol),
            int(p),
        )
        return float(out[0])

    def lmrob_control(self, **kwargs: Any) -> Any:
        return self.robustbase.lmrob_control(**kwargs)


@pytest.fixture(scope="session")
def r_session() -> Iterator[_RBridge]:
    """A live rpy2 session with ``robustbase`` loaded.

    Skipped automatically when rpy2, R, or ``robustbase`` are unavailable.
    """
    if not _rpy2_available():
        pytest.skip("rpy2 not available")
    yield _RBridge()


# ---------------------------------------------------------------------------
# JSON reference loader
# ---------------------------------------------------------------------------
def load_reference(name: str) -> dict[str, Any]:
    """Load a JSON reference file produced by ``scripts/generate_r_reference.R``.

    Parameters
    ----------
    name :
        Reference name (no extension), e.g. ``"stackloss_default"``.
    """
    path = REFERENCE_DIR / f"{name}.json"
    if not path.exists():
        warnings.warn(
            f"Reference {name!r} missing; run scripts/generate_r_reference.R",
            stacklevel=2,
        )
        pytest.skip(f"reference {name!r} not generated yet")
    with path.open() as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def reference_loader() -> Any:
    return load_reference
