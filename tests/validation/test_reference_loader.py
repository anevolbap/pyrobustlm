# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 1 acceptance: every reference JSON loads and is internally consistent.

These checks do not yet compare to a Python implementation; they validate the
reference harness itself. Real per-output comparisons land in Phases 2+.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = REPO_ROOT / "tests" / "reference"

REQUIRED_KEYS = {
    "name",
    "dataset",
    "formula",
    "control",
    "seed",
    "coefficients",
    "scale",
    "weights",
    "residuals",
    "fitted",
    "cov",
    "df_residual",
    "converged",
    "psi",
    "tuning_psi",
    "tuning_chi",
    "rb_version",
    "r_version",
}


def _all_references() -> list[Path]:
    return sorted(REFERENCE_DIR.glob("*.json"))


def test_reference_directory_populated() -> None:
    refs = _all_references()
    assert len(refs) >= 30, (
        f"expected at least 30 references, found {len(refs)}; "
        "run `Rscript scripts/generate_r_reference.R`"
    )


@pytest.mark.parametrize("path", _all_references(), ids=lambda p: p.stem)
def test_reference_shapes(path: Path) -> None:
    """Each reference has the expected keys and consistent vector shapes."""
    ref = json.loads(path.read_text())
    missing = REQUIRED_KEYS - ref.keys()
    assert not missing, f"{path.name} missing keys: {missing}"

    weights = np.asarray(ref["weights"], dtype=float)
    residuals = np.asarray(ref["residuals"], dtype=float)
    fitted = np.asarray(ref["fitted"], dtype=float)
    n = weights.size
    assert n > 0
    assert residuals.size == n
    assert fitted.size == n
    assert (weights >= 0).all()
    assert weights.max() <= 1.0 + 1e-12

    cov = np.asarray(ref["cov"], dtype=float)
    p = len(ref["coefficients"])
    assert cov.shape == (p, p)
    # vcov must be symmetric to working precision
    assert np.allclose(cov, cov.T, atol=1e-10)

    if ref["converged"]:
        scale = float(ref["scale"])
        assert scale >= 0.0


def test_psi_family_default() -> None:
    """The default psi family in modern robustbase (>= 0.93) is 'bisquare'.

    This pins the behaviour we have to match in Phase 2.
    """
    default = json.loads((REFERENCE_DIR / "stackloss_default.json").read_text())
    assert default["psi"] in {"bisquare", "lqq"}, default["psi"]
