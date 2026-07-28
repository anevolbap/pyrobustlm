# SPDX-License-Identifier: GPL-3.0-or-later
"""Differential tests between the Cython kernels and the NumPy reference.

``test_engine_c_parity`` compares whole fits, which means it also compares
two different resampling streams: when the engines disagree there, you
cannot tell an arithmetic bug from a different basin of attraction.

These tests remove the RNG. Every case feeds *identical* inputs to both
implementations, so any disagreement is a real difference in the code.
That is what caught the D-step kappa bug: the Cython path carried its own
hardcoded kappa table whose ggw case-4 entry was a copy of case 1, which
put the D-scale 3.8% off on every ggw fit while whole-fit parity tests
attributed it to basin drift.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import pylmrob._fast_s as fast_s_mod
from pylmrob import psi as _psi
from pylmrob._mm import mm_iterate
from pylmrob.d_scale import d_scale, kappa
from pylmrob.formula import model_matrix
from pylmrob.inference import vcov_avar1

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "tests" / "data"

_cy_lmrob = importlib.import_module("pylmrob._core._lmrob")
_cy_fast_s = importlib.import_module("pylmrob._core._fast_s")

_FAMILY_IDS = {"bisquare": 0, "hampel": 1, "optimal": 2, "lqq": 3, "ggw": 4}

_TUNING_PSI: dict[str, tuple[float, ...]] = {
    "bisquare": (4.685061,),
    "hampel": (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085),
    "optimal": (1.060158,),
    "lqq": (1.4734061, 0.9822707, 1.5),
    "ggw": (4.0,),
}
_TUNING_CHI: dict[str, tuple[float, ...]] = {
    "bisquare": (1.547645,),
    "hampel": (1.5 * 0.2119163, 3.5 * 0.2119163, 8.0 * 0.2119163),
    "optimal": (0.4047,),
    "lqq": (0.4015457, 0.2676971, 1.5),
    "ggw": (6.0,),
}

_DATASETS = [
    ("stackloss", "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."),
    ("coleman", "Y ~ salaryP + fatherWc + sstatus + teacherSc + motherLev"),
    ("salinity", "Y ~ X1 + X2 + X3"),
    ("delivery", "delTime ~ n.prod + distance"),
    ("phosphor", "plant ~ inorg + organic"),
    ("aircraft", "Y ~ X1 + X2 + X3 + X4"),
    ("pension", "Reserves ~ Income"),
    ("starsCYG", "log.light ~ log.Te"),
    ("hbk", "Y ~ X1 + X2 + X3"),
    ("wood", "y ~ x1 + x2 + x3 + x4 + x5"),
]

_CASES = [(ds, f, fam) for ds, f in _DATASETS for fam in _FAMILY_IDS]
_IDS = [f"{ds}-{fam}" for ds, _f, fam in _CASES]


def _pad3(t: tuple[float, ...]) -> np.ndarray:
    out = np.zeros(3, dtype=np.float64)
    for i, v in enumerate(t[:3]):
        out[i] = float(v)
    return out


def _fixture(dataset: str, formula: str):
    """Design plus an RNG-free starting point (OLS + MAD)."""
    path = DATA_DIR / f"{dataset}.csv"
    if not path.exists():
        pytest.skip(f"data file missing: {path}")
    design = model_matrix(formula, pd.read_csv(path))
    X = np.ascontiguousarray(design.X, dtype=np.float64)
    y = np.ascontiguousarray(design.y, dtype=np.float64)
    beta0 = np.linalg.lstsq(X, y, rcond=None)[0]
    r0 = y - X @ beta0
    sigma0 = float(1.4826 * np.median(np.abs(r0 - np.median(r0))))
    if sigma0 <= 0:
        pytest.skip(f"{dataset}: degenerate OLS residual scale")
    return X, y, beta0, r0, sigma0


@pytest.mark.parametrize("dataset,formula,family", _CASES, ids=_IDS)
def test_cy_mm_matches_numpy(dataset: str, formula: str, family: str) -> None:
    X, y, beta0, _r0, sigma0 = _fixture(dataset, formula)
    beta_cy = np.ascontiguousarray(beta0).copy()
    _cy_lmrob.cy_lmrob_mm(
        X, y, beta_cy, sigma0, _FAMILY_IDS[family], _pad3(_TUNING_PSI[family]), 50, 1e-7
    )
    beta_np = mm_iterate(
        X=X,
        y=y,
        beta_init=beta0,
        sigma=sigma0,
        psi_family=family,
        psi_k=_TUNING_PSI[family],
        max_it=50,
        rel_tol=1e-7,
    ).coef
    np.testing.assert_allclose(beta_cy, beta_np, rtol=1e-9, atol=1e-11)


@pytest.mark.parametrize("dataset,formula,family", _CASES, ids=_IDS)
def test_cy_refine_matches_numpy(dataset: str, formula: str, family: str) -> None:
    X, y, beta0, _r0, sigma0 = _fixture(dataset, formula)
    cfg = fast_s_mod.FastSConfig(psi_chi=family, k_chi=_TUNING_CHI[family])

    beta_cy = np.ascontiguousarray(beta0).copy()
    scale_cy, _conv, _it, status = _cy_fast_s.cy_refine_to_convergence(
        X,
        y,
        beta_cy,
        sigma0,
        _FAMILY_IDS[family],
        _pad3(_TUNING_CHI[family]),
        0.5,
        50,
        1e-7,
        200,
        1e-10,
    )
    if status != 0:
        pytest.skip(f"{dataset}/{family}: Cython refine reported status {status}")

    saved = fast_s_mod._CY_REFINE
    fast_s_mod._CY_REFINE = None  # force the NumPy implementation
    try:
        beta_np, scale_np, _c, _i = fast_s_mod._refine_to_convergence(X, y, beta0, sigma0, cfg)
    finally:
        fast_s_mod._CY_REFINE = saved

    np.testing.assert_allclose(beta_cy, beta_np, rtol=1e-9, atol=1e-11)
    np.testing.assert_allclose(scale_cy, scale_np, rtol=1e-9)


@pytest.mark.parametrize("dataset,formula,family", _CASES, ids=_IDS)
def test_cy_vcov_matches_numpy(dataset: str, formula: str, family: str) -> None:
    X, y, beta0, r0, sigma0 = _fixture(dataset, formula)
    p = X.shape[1]
    residuals = y - X @ beta0

    cov_cy = np.zeros((p, p), dtype=np.float64)
    status = _cy_lmrob.cy_lmrob_vcov_avar1(
        X,
        np.ascontiguousarray(residuals),
        np.ascontiguousarray(r0),
        sigma0,
        _FAMILY_IDS[family],
        _pad3(_TUNING_PSI[family]),
        _pad3(_TUNING_CHI[family]),
        0.5,
        cov_cy,
    )
    if status != 0:
        pytest.skip(f"{dataset}/{family}: Cython vcov reported status {status}")

    cov_np = vcov_avar1(
        X=X,
        residuals=residuals,
        sigma=sigma0,
        psi_family=family,
        psi_k=_TUNING_PSI[family],
        init_residuals=r0,
        chi_family=family,
        chi_k=_TUNING_CHI[family],
        bb=0.5,
    )
    np.testing.assert_allclose(cov_cy, cov_np, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("dataset,formula,family", _CASES, ids=_IDS)
def test_cy_d_scale_matches_numpy(dataset: str, formula: str, family: str) -> None:
    """The regression test for the kappa table bug.

    Before the fix this failed for every ggw case at rtol=1e-2.
    """
    X, y, beta0, _r0, sigma0 = _fixture(dataset, formula)
    n = X.shape[0]
    residuals = y - X @ beta0
    rweights = _psi.wgt(residuals / sigma0, family, _TUNING_PSI[family])

    tau_cy = np.empty(n, dtype=np.float64)
    scale_cy, _conv_cy, status = _cy_lmrob.cy_lmrob_d_scale(
        X,
        np.ascontiguousarray(residuals),
        np.ascontiguousarray(rweights),
        sigma0,
        _FAMILY_IDS[family],
        _pad3(_TUNING_PSI[family]),
        200,
        1e-7,
        tau_cy,
        kappa(family, _TUNING_PSI[family]),
    )
    if status != 0:
        pytest.skip(f"{dataset}/{family}: Cython d_scale reported status {status}")

    scale_np, _conv_np, tau_np, _h = d_scale(
        X=X,
        residuals=residuals,
        rweights=rweights,
        init_scale=sigma0,
        family=family,
        c_psi=_TUNING_PSI[family],
        max_iter=200,
        tol=1e-7,
    )
    np.testing.assert_allclose(scale_cy, scale_np, rtol=1e-9)
    np.testing.assert_allclose(tau_cy, tau_np, rtol=1e-9, atol=1e-12)
