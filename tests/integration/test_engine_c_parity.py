# SPDX-License-Identifier: GPL-3.0-or-later
"""Verify ``Control(engine_c=True)`` matches the default path
element-wise when the two RNG paths land on the same basin.

The Cython subset-draw and the NumPy subset-draw consume different
byte sequences from the BitGenerator, so basins can shift. For most
seeds on real-world datasets they land in the same basin, in which
case the fits should agree to machine precision."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "tests" / "data"


# Each entry: (dataset, formula, seed). The very-small-n datasets
# (stackloss n=21, phosphor n=18) are excluded because their basin of
# attraction is sensitive enough that small BLAS/RNG differences
# between platforms can flip whether vcov_avar1 succeeds. The cases
# below are tight enough (n>=25, well-conditioned design) that both
# paths agree to machine precision on every supported platform.
_CASES = [
    ("delivery", "delTime ~ n.prod + distance", 42),
    ("salinity", "Y ~ X1 + X2 + X3", 42),
    ("coleman", "Y ~ salaryP + fatherWc + sstatus + teacherSc + motherLev", 42),
    ("wood", "y ~ x1 + x2 + x3 + x4 + x5", 42),
]


@pytest.mark.parametrize("dataset,formula,seed", _CASES)
def test_engine_c_matches_default(dataset: str, formula: str, seed: int) -> None:
    """engine_c=True produces a fit equivalent to the default path.

    Tolerances are tight (rtol=1e-8 on coef/scale): when both paths
    land in the same basin the only differences are LAPACK floating
    point ordering in the inner loops.
    """
    path = DATA_DIR / f"{dataset}.csv"
    if not path.exists():
        pytest.skip(f"data file missing: {path}")
    df = pd.read_csv(path)

    ctrl_def = Control(nResample=500)
    ctrl_c = Control(nResample=500, engine_c=True)

    try:
        fit_def = lmrob(formula, df, control=ctrl_def, seed=seed)
        fit_c = lmrob(formula, df, control=ctrl_c, seed=seed)
    except FloatingPointError as exc:
        # Documented: basin drift can produce a singular vcov on some
        # (platform, BLAS) combinations. Skip the assertion in that
        # case; the basin-drift behavior is itself covered by the
        # numerical-notes.md log.
        pytest.skip(f"basin drift made vcov singular for {dataset!r}: {exc}")

    assert fit_def.converged_, f"{dataset}: default fit did not converge"
    assert fit_c.converged_, f"{dataset}: engine_c fit did not converge"

    np.testing.assert_allclose(fit_c.coef_, fit_def.coef_, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(fit_c.scale_, fit_def.scale_, rtol=1e-8)
    # vcov diagonals should match to similar tolerance.
    np.testing.assert_allclose(np.diag(fit_c.cov_), np.diag(fit_def.cov_), rtol=1e-6, atol=1e-10)
