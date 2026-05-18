# SPDX-License-Identifier: GPL-3.0-or-later
"""Validation: ``summary`` and ``anova`` against R's ``robustbase``."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, anova, lmrob

REPO_ROOT = Path(__file__).resolve().parents[2]


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


# ---------------------------------------------------------------------------
# summary()
# ---------------------------------------------------------------------------


# We compare R's summary.lmrob outputs (R-squared, t-values, p-values) for
# stackloss and delivery. Both fits must converge in the same RNG basin as R,
# which they do at nResample=500 default.
@pytest.mark.parametrize(
    "dataset,formula",
    [
        ("stackloss", "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."),
        ("delivery", "delTime ~ n.prod + distance"),
    ],
)
def test_summary_matches_r(dataset, formula, r_session):
    """summary() coefficient table and R^2 match R element-wise."""
    df = _ensure_dataset(dataset)

    # Run R's lmrob and capture summary outputs.
    ro = r_session.ro
    rcode = f"""
    suppressMessages(library(robustbase))
    df <- read.csv("{REPO_ROOT}/tests/data/{dataset}.csv")
    set.seed(42)
    fit <- lmrob({formula}, data = df)
    s <- summary(fit)
    list(
        coef = s$coefficients,
        r2 = s$r.squared,
        adj_r2 = s$adj.r.squared,
        sigma = s$sigma
    )
    """
    res = ro.r(rcode)
    r_coef = np.asarray(res.rx2("coef"))
    r_r2 = float(res.rx2("r2")[0])
    r_adj = float(res.rx2("adj_r2")[0])
    r_sigma = float(res.rx2("sigma")[0])

    fit = lmrob(formula, df, control=Control(nResample=500), seed=42)
    summ = fit.summary()

    # R-squared and sigma match closely.
    np.testing.assert_allclose(summ.r_squared, r_r2, rtol=1e-3)
    np.testing.assert_allclose(summ.adj_r_squared, r_adj, rtol=1e-3)
    np.testing.assert_allclose(summ.scale, r_sigma, rtol=1e-3)

    # Coefficient table: estimate / std error / t / p, all element-wise.
    np.testing.assert_allclose(summ.coefficients[:, 0], r_coef[:, 0], rtol=1e-3)
    np.testing.assert_allclose(summ.coefficients[:, 1], r_coef[:, 1], rtol=2e-3)
    np.testing.assert_allclose(summ.coefficients[:, 2], r_coef[:, 2], rtol=2e-3)
    # P-values can be tiny; use atol for the absolute floor.
    np.testing.assert_allclose(summ.coefficients[:, 3], r_coef[:, 3], rtol=1e-2, atol=1e-6)


def test_summary_str_contains_terms():
    """The summary printout includes the term names and 'Robust'."""
    df = _ensure_dataset("stackloss")
    fit = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, seed=0)
    s = fit.summary()
    out = str(s)
    assert "Air.Flow" in out
    assert "Robust residual standard error" in out
    assert "Multiple R-squared" in out


# ---------------------------------------------------------------------------
# anova()
# ---------------------------------------------------------------------------


def test_anova_wald_pair_matches_r() -> None:
    """anova(full, reduced) Wald chi-squared and p-value match R."""
    df = _ensure_dataset("stackloss")
    ctrl = Control(nResample=500)
    full = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, control=ctrl, seed=42)
    red = lmrob("stack.loss ~ Air.Flow + Water.Temp", df, control=ctrl, seed=42)
    red2 = lmrob("stack.loss ~ Air.Flow", df, control=ctrl, seed=42)

    # Pair full vs red (drops Acid.Conc.).
    tbl = anova(full, red).table
    # Row 0 is the full model header; row 1 is the test row.
    assert tbl[1, 2] == 1  # df = 1 dropped term
    # R values from `anova(full, red)`:
    np.testing.assert_allclose(tbl[1, 1], 2.6105, rtol=2e-3)
    np.testing.assert_allclose(tbl[1, 3], 0.1062, rtol=2e-3)

    # Pair full vs red2 (drops Water.Temp + Acid.Conc.).
    tbl2 = anova(full, red2).table
    assert tbl2[1, 2] == 2
    np.testing.assert_allclose(tbl2[1, 1], 7.0281, rtol=2e-3)
    np.testing.assert_allclose(tbl2[1, 3], 0.02978, rtol=2e-3)


def test_anova_chained_sequential() -> None:
    """anova(m1, m2, m3) compares each model against the previous one
    (matches anova.lm semantics; differs from robustbase's chained
    output which has a bug with the largest-model reference)."""
    df = _ensure_dataset("stackloss")
    ctrl = Control(nResample=500)
    full = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, control=ctrl, seed=42)
    red = lmrob("stack.loss ~ Air.Flow + Water.Temp", df, control=ctrl, seed=42)
    red2 = lmrob("stack.loss ~ Air.Flow", df, control=ctrl, seed=42)

    tbl = anova(full, red, red2).table
    # Row 1 (full vs red): df=1.
    assert tbl[1, 2] == 1
    # Row 2 (red vs red2): df=1.
    assert tbl[2, 2] == 1


def test_anova_rejects_non_nested():
    """Models that aren't strictly nested should raise."""
    df = _ensure_dataset("stackloss")
    ctrl = Control(nResample=500)
    full = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, control=ctrl, seed=42)
    other = lmrob("stack.loss ~ Air.Flow + Acid.Conc.", df, control=ctrl, seed=42)

    # full ⊃ other (drops Water.Temp). That's a valid nesting → should not raise.
    anova(full, other)

    # Two completely different models with one term each (not nested).
    a = lmrob("stack.loss ~ Air.Flow", df, control=ctrl, seed=42)
    b = lmrob("stack.loss ~ Water.Temp", df, control=ctrl, seed=42)
    with pytest.raises(ValueError):
        anova(a, b)


def test_anova_requires_two_fits():
    df = _ensure_dataset("stackloss")
    fit = lmrob("stack.loss ~ Air.Flow", df, seed=0)
    with pytest.raises(ValueError, match="at least two"):
        anova(fit)


def test_anova_unknown_test():
    df = _ensure_dataset("stackloss")
    fit = lmrob("stack.loss ~ Air.Flow + Water.Temp", df, seed=0)
    fit2 = lmrob("stack.loss ~ Air.Flow", df, seed=0)
    with pytest.raises(NotImplementedError, match="'Wald' or 'Deviance'"):
        anova(fit, fit2, test="LRT")


def test_anova_deviance_pair_matches_r() -> None:
    """anova(full, reduced, test='Deviance') matches R element-wise."""
    df = _ensure_dataset("stackloss")
    ctrl = Control(nResample=500)
    full = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, control=ctrl, seed=42)
    red = lmrob("stack.loss ~ Air.Flow + Water.Temp", df, control=ctrl, seed=42)
    red2 = lmrob("stack.loss ~ Air.Flow", df, control=ctrl, seed=42)

    tbl = anova(full, red, test="Deviance").table
    # R values from `anova(full, red, test='Deviance')`:
    np.testing.assert_allclose(tbl[1, 1], 1.5978, rtol=2e-3)
    np.testing.assert_allclose(tbl[1, 3], 0.2062, rtol=2e-3)
    assert tbl[1, 2] == 1

    tbl2 = anova(full, red2, test="Deviance").table
    np.testing.assert_allclose(tbl2[1, 1], 6.4836, rtol=2e-3)
    np.testing.assert_allclose(tbl2[1, 3], 0.03909, rtol=2e-3)
    assert tbl2[1, 2] == 2


def test_anova_factor_design_matches_r() -> None:
    """anova() on a multi-column factor drop (Wald + Deviance) matches R.

    Drops the ``RegionF`` factor (3 dummies) from the education
    regression. R reference values are from
    ``anova(lmrob(Y ~ RegionF + X1 + X2 + X3), lmrob(Y ~ X1 + X2 + X3))``.
    """
    df = _ensure_dataset("education")
    df["RegionF"] = df["Region"].astype("category")
    ctrl = Control(nResample=500)
    full = lmrob("Y ~ RegionF + X1 + X2 + X3", df, control=ctrl, seed=42)
    red = lmrob("Y ~ X1 + X2 + X3", df, control=ctrl, seed=42)

    wald = anova(full, red).table
    assert wald[1, 2] == 3  # three Region dummies dropped
    np.testing.assert_allclose(wald[1, 1], 7.7867, rtol=5e-3)
    np.testing.assert_allclose(wald[1, 3], 0.05063, rtol=5e-3)

    dev = anova(full, red, test="Deviance").table
    assert dev[1, 2] == 3
    np.testing.assert_allclose(dev[1, 1], 11.797, rtol=5e-3)
    np.testing.assert_allclose(dev[1, 3], 0.008113, rtol=5e-2)
