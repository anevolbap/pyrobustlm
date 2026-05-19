# SPDX-License-Identifier: GPL-3.0-or-later
"""Parity tests: pylmrob ``Control(rng="R")`` vs R's ``lmrob``.

With ``rng="R"``, the resample draws are byte-identical to robustbase's
``unif_rand`` stream. Final coefficients should agree with R's
``lmrob(..., subsampling="simple")`` after ``set.seed(seed)`` to a
much tighter tolerance than the PCG64 path.

The remaining gap (rtol around 1e-5 rather than truly byte-identical)
comes from the refinement step: numpy LAPACK calls and robustbase's
inline IRWLS aren't required to produce bit-identical doubles.
"""

from __future__ import annotations

import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")


def _have_r() -> bool:
    if shutil.which("Rscript") is None:
        return False
    out = subprocess.run(
        ["Rscript", "-e", 'cat(requireNamespace("robustbase", quietly=TRUE))'],
        capture_output=True, text=True, timeout=10,
    )
    return out.stdout.strip() == "TRUE"


pytestmark = pytest.mark.skipif(not _have_r(), reason="R + robustbase not available")


@pytest.fixture
def stackloss() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "stackloss.csv"))


def _r_lmrob_simple(seed: int) -> tuple[np.ndarray, float]:
    """Run R's ``lmrob`` on stackloss with simple subsampling and
    ``set.seed(seed)``. Returns (coef, scale)."""
    script = f"""
        suppressMessages(library(robustbase))
        data(stackloss)
        ctrl <- lmrob.control(subsampling="simple")
        set.seed({seed})
        fit <- lmrob(stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.,
                     data = stackloss, control = ctrl)
        cat(sprintf("%.17g\\n", coef(fit)), sep = "")
        cat(sprintf("%.17g\\n", fit$scale))
    """
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True, check=True, text=True, timeout=30,
    )
    lines = [float(x) for x in out.stdout.splitlines() if x]
    return np.array(lines[:-1], dtype=np.float64), lines[-1]


@pytest.mark.parametrize("seed", [1, 42, 12345])
def test_stackloss_coef_close_to_R(stackloss: pd.DataFrame, seed: int) -> None:
    r_coef, r_scale = _r_lmrob_simple(seed)
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss,
        control=Control(rng="R"),
        seed=seed,
    )
    # R-mode draws match exactly; remaining drift is from refinement
    # floating-point order.
    np.testing.assert_allclose(fit.coef_, r_coef, rtol=1e-4, atol=1e-4)
    np.testing.assert_allclose(fit.scale_, r_scale, rtol=1e-4)


def test_rng_r_tighter_than_pcg64(stackloss: pd.DataFrame) -> None:
    """``rng="R"`` should land closer to R than the PCG64 default does.

    Asserts a specific inequality on |coef - R_coef|; the PCG64 path
    diverges in the 3rd-4th decimal on stackloss while R-mode lands in
    the 5th-6th.
    """
    r_coef, _ = _r_lmrob_simple(42)
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    fit_r = lmrob(formula, stackloss, control=Control(rng="R"), seed=42)
    fit_p = lmrob(formula, stackloss, control=Control(rng="PCG64"), seed=42)
    err_r = np.max(np.abs(fit_r.coef_ - r_coef))
    err_p = np.max(np.abs(fit_p.coef_ - r_coef))
    # R-mode shouldn't be worse than PCG64; allow equality on
    # well-conditioned datasets where both land in the same basin.
    assert err_r <= err_p + 1e-12
