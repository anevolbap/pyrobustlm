# SPDX-License-Identifier: GPL-3.0-or-later
"""``outlier_stats`` against ``robustbase::outlierStats``.

The fixture is a factor design where 70% of one level is contaminated
with non-constant offsets, so the level dummy cannot absorb it and the
fit breaks down *locally*: overall it looks healthy (23% rejected, mean
robustness weight 0.74) while within level ``c`` it has rejected 14 of
20 and the mean weight is 0.29. That is the situation R warns about and
that ``setting="KS2014"`` is meant to avoid.

R's robustness weights are fed in directly, so this compares the
statistic rather than the fit; whether we find the same S optimum as R
is a separate question covered elsewhere.

Reference generated with robustbase 0.99-7 / R 4.2.2.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylmrob import Control, lmrob, outlier_stats

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE = REPO_ROOT / "tests" / "reference" / "outliers" / "local_breakdown.json"
DATA = REPO_ROOT / "tests" / "data" / "outlier_local_breakdown.csv"

_FORMULA = "y ~ x + g"


def _r_name(name: str) -> str:
    """formulaic's ``g[T.c]`` is R's ``gc``."""
    m = re.match(r"^(\w+)\[T\.(.+)\]$", name)
    return f"{m.group(1)}{m.group(2)}" if m else name


@pytest.fixture(scope="module")
def case():
    if not REFERENCE.exists() or not DATA.exists():  # pragma: no cover
        pytest.skip("outlier reference fixture missing")
    ref = json.loads(REFERENCE.read_text())
    df = pd.read_csv(DATA)
    df["g"] = pd.Categorical(df["g"])
    fit = lmrob(_FORMULA, df, control=Control(nResample=500), seed=1)
    # Compare the statistic, not the fit: use R's own robustness weights.
    fit.rweights_ = np.asarray(ref["rweights"], dtype=np.float64)
    return ref, outlier_stats(fit, shout=False)


def test_rows_match_r(case) -> None:
    ref, stats = case
    r_rows = {row["name"]: row for row in ref["rows"]}
    assert len(stats) == len(r_rows), (
        f"row count {len(stats)} vs R {len(r_rows)}: "
        f"ours={[r.name for r in stats]} R={list(r_rows)}"
    )
    for row in stats:
        r = r_rows.get(_r_name(row.name))
        assert r is not None, f"{row.name}: no R counterpart in {list(r_rows)}"
        assert row.n_nonzero == r["n_nonzero"], f"{row.name}: N.nonzero"
        assert row.n_rejected == r["n_rejected"], f"{row.name}: N.rejected"
        assert row.ratio == pytest.approx(r["ratio"], abs=1e-12), f"{row.name}: Ratio"
        assert row.mean_robweight == pytest.approx(r["mean_rw"], abs=1e-12), (
            f"{row.name}: Mean.RobWeight"
        )


def test_flags_the_broken_level(case) -> None:
    _ref, stats = case
    assert [_r_name(n) for n in stats.flagged] == ["gc"]


def test_eps_outlier_default_matches_r(case) -> None:
    ref, stats = case
    assert stats.eps_outlier == pytest.approx(float(ref["eps_outlier"]), rel=1e-12)


def test_intercept_and_continuous_columns_are_excluded(case) -> None:
    """Only indicator-like columns carry local information.

    R selects on ``colSums(xnz) < NROW(xnz)``, which drops the intercept
    and any continuous predictor.
    """
    _ref, stats = case
    names = [row.name for row in stats]
    assert "Intercept" not in names and "(Intercept)" not in names
    assert "x" not in names


def test_warns_on_local_breakdown() -> None:
    ref = json.loads(REFERENCE.read_text())
    df = pd.read_csv(DATA)
    df["g"] = pd.Categorical(df["g"])
    fit = lmrob(_FORMULA, df, control=Control(nResample=500), seed=1)
    fit.rweights_ = np.asarray(ref["rweights"], dtype=np.float64)
    with pytest.warns(RuntimeWarning, match="local breakdown"):
        outlier_stats(fit)


def test_clean_fit_does_not_warn() -> None:
    """A design with no local breakdown must stay quiet."""
    rng = np.random.default_rng(0)
    n = 90
    df = pd.DataFrame(
        {
            "x": rng.standard_normal(n),
            "g": pd.Categorical(np.repeat(["a", "b", "c"], n // 3)),
        }
    )
    df["y"] = 1.0 + 2.0 * df["x"] + rng.standard_normal(n) * 0.3
    fit = lmrob("y ~ x + g", df, control=Control(nResample=500), seed=3)
    stats = outlier_stats(fit, shout=None)
    assert stats.flagged == [], stats.flagged


def test_shout_forces_and_suppresses() -> None:
    rng = np.random.default_rng(1)
    n = 60
    df = pd.DataFrame(
        {
            "x": rng.standard_normal(n),
            "g": pd.Categorical(np.repeat(["a", "b"], n // 2)),
        }
    )
    df["y"] = 1.0 + 2.0 * df["x"] + rng.standard_normal(n) * 0.3
    fit = lmrob("y ~ x + g", df, control=Control(nResample=500), seed=4)

    with pytest.warns(RuntimeWarning, match="local breakdown"):
        outlier_stats(fit, shout=True)

    import warnings as _w

    with _w.catch_warnings():
        _w.simplefilter("error")
        outlier_stats(fit, shout=False)
