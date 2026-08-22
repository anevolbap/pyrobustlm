# SPDX-License-Identifier: GPL-3.0-or-later
"""``vcov_avar1`` against R, with the resampling RNG removed.

``docs/bench-report.md`` shows cov-diagonal errors of 0.4-0.5 for
``psi_hampel`` and ``psi_ggw`` on stackloss, and ``numerical-notes.md``
entry 1 used to attribute both to basin drift. Whole-fit comparisons
cannot tell a wrong formula from a different fit, which is the same
blind spot that hid the D-step kappa bug.

This file removes the search entirely: it feeds R's *own* final
residuals, initial-S residuals and scale into :func:`vcov_avar1` and
compares against R's ``vcov()``. Any disagreement here is the formula.
Both families come back at 1e-9 or better, so the formula is right and
the bench-report numbers are about which fit was found, not how the
covariance is computed.

Reference: ``tests/reference/vcov/stackloss_avar1_inputs.json``,
regenerate with the snippet in that file's ``meta`` block (robustbase
0.99-7, R 4.2.2).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control
from pylmrob.formula import model_matrix
from pylmrob.inference import vcov_avar1

REPO_ROOT = Path(__file__).resolve().parents[2]
# Lives in a subdirectory: tests/reference/*.json is globbed by
# test_reference_loader, which expects the single-fit schema, and this
# is a multi-family bundle of vcov inputs.
REFERENCE = REPO_ROOT / "tests" / "reference" / "vcov" / "stackloss_avar1_inputs.json"
DATA = REPO_ROOT / "tests" / "data" / "stackloss.csv"

_FAMILIES = ["bisquare", "hampel", "optimal", "lqq", "ggw"]

# R's lmrob.S returns an ``init$residuals`` field that does not match
# ``y - X %*% init$coefficients`` (about 2.5 on this dataset), because
# ``fast_s()`` returns the last refined survivor's residuals next to the
# best survivor's coefficients. ``.vcov.avar1`` reads that field, so R's
# own covariance for an affected fit is built from residuals
# inconsistent with the coefficients it reports. R-Forge bug #6873, see
# numerical-notes entry 11.
#
# Which families land on it is seed-dependent: hampel and welsh show the
# same gap on other seeds. This set is what the committed fixture's seed
# produced, not a statement about the families.
_R_INIT_INCONSISTENT = {"ggw", "lqq"}
_FLOAT_NOISE = 1e-4


def _load_reference() -> dict:
    if not REFERENCE.exists():  # pragma: no cover - fixture is committed
        pytest.skip(f"reference missing: {REFERENCE}")
    return json.loads(REFERENCE.read_text())


def _design() -> tuple[np.ndarray, np.ndarray]:
    if not DATA.exists():  # pragma: no cover
        pytest.skip(f"data missing: {DATA}")
    ref = _load_reference()
    design = model_matrix(ref["meta"]["formula"], pd.read_csv(DATA))
    return (
        np.ascontiguousarray(design.X, dtype=np.float64),
        np.ascontiguousarray(design.y, dtype=np.float64),
    )


def _tuning(family: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Internal-form tuning constants.

    R reports ggw/lqq tuning in the user-facing NA-padded 4-vector form;
    our kernels take the internal form, which ``Control`` already holds.
    """
    ctrl = Control(psi=family)  # type: ignore[arg-type]
    return (
        tuple(np.atleast_1d(np.asarray(ctrl.tuning_psi, dtype=float)).ravel()),
        tuple(np.atleast_1d(np.asarray(ctrl.tuning_chi, dtype=float)).ravel()),
    )


@pytest.mark.parametrize("family", _FAMILIES)
def test_vcov_avar1_matches_r_given_r_inputs(family: str) -> None:
    """Same inputs in, same covariance out. No RNG anywhere."""
    ref = _load_reference()["families"][family]
    X, _y = _design()
    psi_k, chi_k = _tuning(family)

    cov_py = vcov_avar1(
        X=X,
        residuals=np.asarray(ref["residuals"], dtype=np.float64),
        sigma=float(ref["scale"]),
        psi_family=family,
        psi_k=psi_k,
        init_residuals=np.asarray(ref["init_residuals"], dtype=np.float64),
        chi_family=family,
        chi_k=chi_k,
        bb=float(ref["bb"]),
    )
    p = X.shape[1]
    cov_r = np.asarray(ref["cov"], dtype=np.float64).reshape(p, p)

    np.testing.assert_allclose(cov_py, cov_r, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("family", _FAMILIES)
def test_r_init_residual_self_consistency(family: str) -> None:
    """Pin the upstream quirk the divergence note depends on.

    If a future robustbase makes ``lmrob.S`` self-consistent, this fails
    for the families in ``_R_INIT_INCONSISTENT`` and the note in
    ``numerical-notes.md`` (entry 11) needs revisiting, because our
    ``init_residuals`` would then agree with R's and the ggw covariance
    gap should close on its own.
    """
    ref = _load_reference()["families"][family]
    gap = float(ref["init_resid_vs_coef_maxabs"])

    if family in _R_INIT_INCONSISTENT:
        assert gap > _FLOAT_NOISE, (
            f"{family}: R's init$residuals now matches y - X b_init "
            f"(gap={gap:.2e}); revisit numerical-notes entry 11"
        )
    else:
        assert gap <= _FLOAT_NOISE, f"{family}: unexpected init inconsistency, gap={gap:.2e}"
