# SPDX-License-Identifier: GPL-3.0-or-later
"""Phase 7+: KS2011 D-scale validation against R's lmrob..D..fit."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyrobustlm.d_scale import d_scale, find_d_scale, kappa, tau

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def stackloss_df() -> pd.DataFrame:
    out = REPO_ROOT / "tests" / "data" / "stackloss.csv"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "Rscript",
                "-e",
                f"library(robustbase); data(stackloss); write.csv(stackloss, '{out}', row.names=FALSE)",
            ],
            capture_output=True,
            check=True,
        )
    return pd.read_csv(out)


def test_kappa_bisquare_matches_r():
    """robustbase:::lmrob.kappa for bisquare 4.685 returns ~0.8281."""
    kp = kappa("bisquare", (4.685061,))
    np.testing.assert_allclose(kp, 0.8280907, rtol=1e-4)


def test_tau_fast_table_round_trip():
    """tau() with the fast table should reproduce sqrt(1-tfact h)(tcorr h + 1)."""
    h = np.array([0.05, 0.10, 0.20, 0.32])
    t = tau(h, "bisquare", (4.685061,))
    tfact, tcorr = 0.9473684, -0.0900833
    expected = np.sqrt(1.0 - tfact * h) * (tcorr * h + 1.0)
    np.testing.assert_allclose(t, expected, rtol=1e-12)


def test_find_d_scale_idempotent_at_fixed_point():
    """Iterating from the converged D-scale should return the same value."""
    rng = np.random.default_rng(0)
    n = 80
    r = rng.standard_normal(n)
    tau_vec = np.full(n, 0.95)
    sgma, conv = find_d_scale(
        r,
        tau_vec,
        kappa_val=0.86,
        family="bisquare",
        c_psi=(4.685061,),
        init_scale=1.0,
    )
    assert conv
    sgma2, _conv2 = find_d_scale(
        r,
        tau_vec,
        kappa_val=0.86,
        family="bisquare",
        c_psi=(4.685061,),
        init_scale=sgma,
    )
    np.testing.assert_allclose(sgma2, sgma, rtol=1e-6)


def _r_has_robustbase() -> bool:
    """True iff Rscript is on PATH and the robustbase package is loadable."""
    import shutil

    if shutil.which("Rscript") is None:
        return False
    try:
        r = subprocess.run(
            ["Rscript", "-e", "library(robustbase)"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        return r.returncode == 0
    except Exception:
        return False


@pytest.mark.skipif(
    not _r_has_robustbase(),
    reason="Rscript + robustbase not available (test drives R's R_find_D_scale C kernel)",
)
def test_d_iteration_matches_r_kernel(stackloss_df):
    """Our D-iteration matches R's ``R_find_D_scale`` C kernel to 5 decimals
    when given identical (r, w_MM, h, tau, kappa) inputs.

    Notes
    -----
    On stackloss, calling our full ``lmrob(setting='KS2011')`` and R's
    full ``lmrob(setting='KS2011')`` produces ~8% different scales; that
    gap comes from RNG basin drift in MM (different residuals and
    rweights), which propagates through tau via the hat diagonal. This
    test isolates the D-iteration itself by feeding identical inputs.
    """
    from pyrobustlm import Control, lmrob

    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        stackloss_df,
        control=Control(nResample=1000),
        seed=42,
    )
    r_path = REPO_ROOT / "tests" / "data" / "_d_test_residuals.csv"
    w_path = REPO_ROOT / "tests" / "data" / "_d_test_rweights.csv"
    np.savetxt(r_path, fit.residuals_)
    np.savetxt(w_path, fit.rweights_)

    # Drive R's R_find_D_scale C kernel directly with our (r, w, h, tau, kappa).
    rscript = f"""
suppressWarnings(library(robustbase))
data(stackloss)
r <- scan('{r_path}', quiet=TRUE)
w_mm <- scan('{w_path}', quiet=TRUE)
X <- model.matrix(~Air.Flow+Water.Temp+Acid.Conc., data=stackloss)
h <- robustbase:::.lmrob.hat(X, w_mm)
tfact <- 0.9473684; tcorr <- -0.0900833
tau_vec <- sqrt(1 - tfact*h) * (tcorr*h + 1)
ctrl <- lmrob.control(psi='bisquare'); ctrl$method <- 'SMDM'
fake <- list(control=ctrl); class(fake) <- 'lmrob'
kp <- robustbase:::lmrob.kappa(fake, ctrl)
ipsi <- robustbase:::.psi2ipsi('bisquare')
sigma_init <- sqrt(sum(w_mm * r^2) / kp / sum(tau_vec^2 * w_mm))
ret <- .C(robustbase:::R_find_D_scale, r=as.double(r), kappa=as.double(kp),
          tau=as.double(tau_vec), length=as.integer(length(r)),
          scale=as.double(sigma_init), c=as.double(4.685061),
          ipsi=as.integer(ipsi), type=3L, rel.tol=as.double(1e-7),
          k.max=as.integer(200), converged=logical(1))
cat(ret$scale, '\\n')
"""
    out = subprocess.run(["Rscript", "-e", rscript], capture_output=True, text=True, check=True)
    r_d_scale = float(out.stdout.strip().split("\n")[-1])

    # Compute py D-scale on identical inputs.
    X = np.column_stack(
        [
            np.ones(len(stackloss_df)),
            stackloss_df[["Air.Flow", "Water.Temp", "Acid.Conc."]].to_numpy(dtype=float),
        ]
    )
    py_d_scale, _conv, _tau, _h = d_scale(
        X=X,
        residuals=fit.residuals_,
        rweights=fit.rweights_,
        init_scale=fit.scale_,
        family="bisquare",
        c_psi=(4.685061,),
    )
    np.testing.assert_allclose(py_d_scale, r_d_scale, rtol=1e-4)


def test_lmrob_setting_KS2014_runs_d_step():
    """End-to-end: lmrob(setting='KS2014') runs the SMDM pipeline and
    populates init_['d_scale'] (R's setting='KS2014' uses method='SMDM').
    """
    df = pd.read_csv(REPO_ROOT / "tests" / "data" / "stackloss.csv")
    from pyrobustlm import Control, lmrob

    fit = lmrob(
        "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
        df,
        control=Control(setting="KS2014", nResample=500),
        seed=42,
    )
    assert "d_scale" in fit.init_
    assert fit.init_["d_converged"] is True
    # The D-scale typically inflates the S-scale on contaminated data.
    assert fit.scale_ > fit.init_["scale"] * 1.05
