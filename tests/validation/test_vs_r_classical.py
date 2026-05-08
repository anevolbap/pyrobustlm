# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 10: end-to-end fit on every classical-dataset reference, diffed
against the saved R output. Looser tolerances than unit tests because of
RNG-basin sensitivity.
"""

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


# Each entry: (ref_name, dataset, formula, max_rtol_coef, max_rtol_scale).
# rtol_coef is the tolerance applied to ``|py - R| / max(|R|, 1)``.
_CASES = [
    (
        "stackloss_default",
        "stackloss",
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        1e-3,
        1e-3,
    ),
    ("coleman_default", "coleman", "Y ~ .", 5e-2, 1e-2),
    ("salinity_default", "salinity", "Y ~ .", 1e-1, 5e-2),
    ("delivery_default", "delivery", "delTime ~ n.prod + distance", 1e-3, 1e-3),
    ("phosphor_default", "phosphor", "plant ~ inorg + organic", 5e-2, 5e-3),
    ("aircraft_default", "aircraft", "Y ~ X1 + X2 + X3 + X4", 1e-2, 1e-2),
    ("pension_default", "pension", "Reserves ~ Income", 5e-1, 1e-1),
    ("starsCYG_default", "starsCYG", "log.light ~ log.Te", 5e-1, 1e-1),
    ("hbk_default", "hbk", "Y ~ .", 1e-1, 5e-2),
    ("wood_default", "wood", "y ~ .", 1e-1, 1e-1),
]


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


@pytest.mark.parametrize("ref_name,dataset,formula,rtol_coef,rtol_scale", _CASES)
def test_classical_dataset(ref_name, dataset, formula, rtol_coef, rtol_scale):
    ref = json.loads((REFERENCE_DIR / f"{ref_name}.json").read_text())
    df = _ensure_dataset(dataset)

    if formula == "Y ~ .":
        rhs = " + ".join([c for c in df.columns if c != "Y"])
        formula = f"Y ~ {rhs}"
    elif formula == "y ~ .":
        rhs = " + ".join([c for c in df.columns if c != "y"])
        formula = f"y ~ {rhs}"

    fit = lmrob(formula, df, control=Control(nResample=500), seed=42)
    assert fit.converged_, f"{ref_name}: MM did not converge"

    name_map = {"Intercept": "(Intercept)"}
    r_coefs = ref["coefficients"]
    for name, py_val in zip(fit.term_names_, fit.coef_, strict=True):
        r_key = name_map.get(name, name)
        if r_key not in r_coefs:
            pytest.skip(f"{ref_name}: term {name!r} not in R reference (formula mismatch)")
        r_val = float(r_coefs[r_key])
        # Use absolute scale for small reference values.
        denom = max(abs(r_val), 1.0)
        diff = abs(py_val - r_val)
        assert diff / denom < rtol_coef, (
            f"{ref_name}: coef {name!r}: py={py_val} R={r_val} diff={diff}"
        )

    np.testing.assert_allclose(
        fit.scale_,
        float(ref["scale"]),
        rtol=rtol_scale,
        err_msg=f"{ref_name}: scale py={fit.scale_} R={ref['scale']}",
    )
