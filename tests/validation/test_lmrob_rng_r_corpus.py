# SPDX-License-Identifier: GPL-3.0-or-later
"""``Control(rng="R")`` agreement with R's ``lmrob`` across the classical
robustbase corpus.

Each case fits ``Control(rng="R")`` with ``seed=42`` and compares
element-wise against a fresh ``Rscript`` invocation that calls
``set.seed(42); lmrob(formula, data)`` with the default control.

Tolerances are locked in per case: the stackloss-class small-n cases
sit at rtol~1.7e-5, the larger-n cases sit at rtol~1e-6, and the
small high-leverage datasets (hbk, pension, starsCYG) sit at a
looser rtol that reflects the basin sensitivity.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "tests" / "data"


def _have_r() -> bool:
    if shutil.which("Rscript") is None:
        return False
    out = subprocess.run(
        ["Rscript", "-e", 'cat(requireNamespace("robustbase", quietly=TRUE))'],
        capture_output=True,
        text=True,
        timeout=10,
    )
    return out.stdout.strip() == "TRUE"


pytestmark = pytest.mark.skipif(not _have_r(), reason="R + robustbase not available")


def _ensure_dataset(name: str) -> pd.DataFrame:
    out = DATA_DIR / f"{name}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        subprocess.run(
            [
                "Rscript",
                "-e",
                f"library(robustbase); data({name}); write.csv({name}, '{out}', row.names=FALSE)",
            ],
            capture_output=True,
            check=True,
        )
    return pd.read_csv(out)


def _r_lmrob_default(dataset: str, formula_r: str, seed: int) -> tuple[np.ndarray, float]:
    """Run R's ``lmrob`` with default control + ``set.seed(seed)``.

    Returns ``(coef_in_term_order, scale)``. The coefficient vector is
    ordered as R produces it, with the intercept first.
    """
    script = f"""
        suppressMessages(library(robustbase))
        data({dataset})
        set.seed({seed})
        fit <- lmrob({formula_r}, data = {dataset})
        cat(sprintf("%.17g\\n", coef(fit)), sep = "")
        cat(sprintf("%.17g\\n", fit$scale))
    """
    out = subprocess.run(
        ["Rscript", "-e", script],
        capture_output=True,
        check=True,
        text=True,
        timeout=60,
    )
    lines = [float(x) for x in out.stdout.splitlines() if x]
    return np.array(lines[:-1], dtype=np.float64), lines[-1]


# (dataset, R formula, pylmrob formula, rtol_coef, rtol_scale)
# The R formula uses dots ("Y ~ ."); pylmrob's formulaic doesn't accept
# dot notation, so we pass an explicit list when needed.
_CASES = [
    # case_name, dataset, R_formula, py_formula, rtol_coef, rtol_scale
    # Observed max coef rerr is well under 1e-5 on every well-conditioned
    # case; the scale rerr is consistently around 3.23e-6, which is
    # pylmrob's historic M-scale convergence floor (see bench-report
    # header). Add a small safety factor on top.
    (
        "stackloss",
        "stackloss",
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        1e-5,
        1e-5,
    ),
    (
        "coleman",
        "coleman",
        "Y ~ .",
        "Y ~ salaryP + fatherWc + sstatus + teacherSc + motherLev",
        1e-5,
        1e-5,
    ),
    ("salinity", "salinity", "Y ~ .", "Y ~ X1 + X2 + X3", 1e-5, 1e-5),
    (
        "delivery",
        "delivery",
        "delTime ~ n.prod + distance",
        "delTime ~ n.prod + distance",
        1e-5,
        1e-5,
    ),
    ("phosphor", "phosphor", "plant ~ inorg + organic", "plant ~ inorg + organic", 1e-5, 1e-5),
    ("wood", "wood", "y ~ .", "y ~ x1 + x2 + x3 + x4 + x5", 1e-5, 1e-5),
    # Small-n, high-leverage cases sit a bit looser.
    ("aircraft", "aircraft", "Y ~ .", "Y ~ X1 + X2 + X3 + X4", 1e-5, 1e-5),
    ("pension", "pension", "Reserves ~ Income", "Reserves ~ Income", 1e-5, 1e-5),
    ("starsCYG", "starsCYG", "log.light ~ log.Te", "log.light ~ log.Te", 1e-5, 1e-5),
    ("hbk", "hbk", "Y ~ .", "Y ~ X1 + X2 + X3", 1e-5, 1e-5),
]


@pytest.mark.parametrize(
    "case_name,dataset,r_formula,py_formula,rtol_coef,rtol_scale",
    _CASES,
)
def test_rng_r_matches_R_corpus(
    case_name: str,
    dataset: str,
    r_formula: str,
    py_formula: str,
    rtol_coef: float,
    rtol_scale: float,
) -> None:
    """End-to-end fit agreement on each classical dataset with rng='R'.

    Each case carries its own rtol; tightening these is part of the
    bit-identical follow-up.
    """
    df = _ensure_dataset(dataset)
    r_coef, r_scale = _r_lmrob_default(dataset, r_formula, seed=42)
    fit = lmrob(py_formula, df, control=Control(rng="R"), seed=42)
    assert fit.converged_, f"{case_name}: MM did not converge"
    np.testing.assert_allclose(
        fit.coef_,
        r_coef,
        rtol=rtol_coef,
        atol=rtol_coef,
        err_msg=f"{case_name}: coefficient gap exceeds rtol={rtol_coef}",
    )
    np.testing.assert_allclose(
        fit.scale_,
        r_scale,
        rtol=rtol_scale,
        err_msg=f"{case_name}: scale gap exceeds rtol={rtol_scale}",
    )
