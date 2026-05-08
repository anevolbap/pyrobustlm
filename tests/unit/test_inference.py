# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 7: covariance estimators validated against R."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyrobustlm import Control, lmrob

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = REPO_ROOT / "tests" / "reference"


def _ensure_dataset(name: str) -> pd.DataFrame:
    out = REPO_ROOT / "tests" / "data" / f"{name}.csv"
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


_COV_CASES = [
    ("stackloss_default", "stackloss", "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", 1e-3),
    ("delivery_default", "delivery", "delTime ~ n.prod + distance", 5e-3),
    ("phosphor_default", "phosphor", "plant ~ inorg + organic", 5e-3),
]


@pytest.mark.parametrize("ref_name,dataset,formula,rtol", _COV_CASES)
def test_vcov_avar1_matches_r(ref_name, dataset, formula, rtol):
    """Covariance from .vcov.avar1 should match R element-wise."""
    ref = json.loads((REFERENCE_DIR / f"{ref_name}.json").read_text())
    df = _ensure_dataset(dataset)
    fit = lmrob(formula, df, control=Control(nResample=1000), seed=42)
    R_cov = np.asarray(ref["cov"], dtype=float)
    np.testing.assert_allclose(fit.cov_, R_cov, rtol=rtol, atol=1e-6)


def test_standard_errors_positive():
    df = _ensure_dataset("stackloss")
    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        df,
        control=Control(nResample=500),
        seed=0,
    )
    se = fit.standard_errors_
    assert (se > 0).all()
    # Confidence interval contains coef
    ci = fit.confint(0.95)
    assert (ci[:, 0] <= fit.coef_).all()
    assert (ci[:, 1] >= fit.coef_).all()
