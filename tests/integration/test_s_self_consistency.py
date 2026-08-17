# SPDX-License-Identifier: GPL-3.0-or-later
"""The reported residuals must belong to the reported coefficients.

Guards against the analogue of robustbase R-Forge bug
[#6873](https://r-forge.r-project.org/tracker/index.php?func=detail&aid=6873&group_id=59&atid=302):
``lmrob.S()`` refines the ``best_r`` survivors through one shared ``res``
scratch array and only copies ``beta`` when a survivor improves the
scale, so the returned ``residuals`` are the last survivor's while the
returned ``coefficients`` are the best one's. On stackloss the two differ
by up to 2.99 for hampel, lqq, ggw and welsh, and ``.vcov.avar1()`` reads
that stale field. See ``docs/numerical-notes.md`` entry 11.

pylmrob never carries residuals out of the search: ``fast_s`` returns
only ``(coef, scale)`` and every residual vector is rebuilt from a
coefficient vector. These tests pin that, so an optimisation that starts
passing a residual buffer around cannot reintroduce the upstream bug
silently.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob
from pylmrob.formula import model_matrix
from pylmrob.inference import vcov_avar1
from pylmrob.scale import m_scale

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "tests" / "data" / "stackloss.csv"
FORMULA = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."

# All six families the upstream sweep covered. bisquare and optimal never
# showed the gap in R; the other four did, on 2 to 7 of 20 seeds.
_FAMILIES = ["bisquare", "optimal", "welsh", "hampel", "lqq", "ggw"]
_SEEDS = (2, 9, 12, 14, 42)

# Measured on the reference machine: the largest relative gap between the
# reported init scale and the M-scale at the reported init coefficients is
# 7.2e-10 over 6 families x 20 seeds x both engines. R's bug shows up here
# as 2.5e-3 or more, so this threshold separates the two by six orders of
# magnitude.
_SCALE_RTOL = 1e-7


def _design() -> tuple[np.ndarray, np.ndarray]:
    if not DATA.exists():  # pragma: no cover - fixture is committed
        pytest.skip(f"data missing: {DATA}")
    design = model_matrix(FORMULA, pd.read_csv(DATA))
    return (
        np.ascontiguousarray(design.X, dtype=np.float64),
        np.ascontiguousarray(design.y, dtype=np.float64),
    )


@pytest.mark.parametrize("family", _FAMILIES)
@pytest.mark.parametrize("seed", _SEEDS)
@pytest.mark.parametrize("engine_c", [True, False])
def test_residuals_match_coefficients(family: str, seed: int, engine_c: bool) -> None:
    """``residuals_`` and ``fitted_`` are the ones implied by ``coef_``."""
    X, y = _design()
    fit = lmrob(
        FORMULA, pd.read_csv(DATA), control=Control(psi=family, engine_c=engine_c), seed=seed
    )  # type: ignore[arg-type]
    coef = np.asarray(fit.coef_, dtype=np.float64)

    np.testing.assert_allclose(fit.residuals_, y - X @ coef, rtol=0, atol=1e-12)
    np.testing.assert_allclose(fit.fitted_, X @ coef, rtol=0, atol=1e-12)


@pytest.mark.parametrize("family", _FAMILIES)
@pytest.mark.parametrize("seed", _SEEDS)
@pytest.mark.parametrize("engine_c", [True, False])
def test_init_scale_matches_init_coef(family: str, seed: int, engine_c: bool) -> None:
    """The reported S scale is the M-scale at the reported S coefficients.

    This is the check that catches the upstream mix-up: when the scale
    belongs to one candidate and the residuals to another, the M-scale of
    ``y - X b_init`` no longer reproduces the reported scale.
    """
    X, y = _design()
    control = Control(psi=family, engine_c=engine_c)  # type: ignore[arg-type]
    fit = lmrob(FORMULA, pd.read_csv(DATA), control=control, seed=seed)
    init_coef = np.asarray(fit.init_["coef"], dtype=np.float64)
    chi_k = tuple(np.atleast_1d(np.asarray(control.tuning_chi, dtype=float)).ravel())

    recomputed = m_scale(y - X @ init_coef, family, chi_k, control.bb, 200, 1e-10, p=X.shape[1])
    np.testing.assert_allclose(
        recomputed,
        float(fit.init_["scale"]),  # type: ignore[arg-type]
        rtol=_SCALE_RTOL,
        err_msg=f"{family}/seed={seed}/engine_c={engine_c}: init scale does not match init coef",
    )


@pytest.mark.parametrize("family", _FAMILIES)
@pytest.mark.parametrize("seed", _SEEDS)
def test_engine_vcov_uses_consistent_init_residuals(family: str, seed: int) -> None:
    """The kernel's inline covariance rebuilds the init residuals too.

    ``cy_lmrob_fit`` computes ``vcov_avar1`` without leaving Cython, so
    the init residuals it feeds the formula are not visible from Python.
    Recomputing the covariance from ``y - X b_init`` pins that they are
    the same ones. This is where the upstream bug does its damage: R's
    ``.vcov.avar1()`` reads ``obj$init$resid``, which is the stale field.
    """
    X, y = _design()
    control = Control(psi=family, engine_c=True, cov=".vcov.avar1")  # type: ignore[arg-type]
    fit = lmrob(FORMULA, pd.read_csv(DATA), control=control, seed=seed)
    init_coef = np.asarray(fit.init_["coef"], dtype=np.float64)

    cov_py = vcov_avar1(
        X=X,
        residuals=np.asarray(fit.residuals_, dtype=np.float64),
        sigma=float(fit.scale_),
        psi_family=family,
        psi_k=tuple(np.atleast_1d(np.asarray(control.tuning_psi, dtype=float)).ravel()),
        init_residuals=y - X @ init_coef,
        chi_family=family,
        chi_k=tuple(np.atleast_1d(np.asarray(control.tuning_chi, dtype=float)).ravel()),
        bb=control.bb,
    )
    np.testing.assert_allclose(cov_py, fit.cov_, rtol=1e-9, atol=1e-12)
