# SPDX-License-Identifier: GPL-3.0-or-later
"""Controls that are accepted but not implemented must say so.

``Control`` mirrors ``lmrob.control()``'s field list so that porting R
code is mechanical, but a few of those fields are not wired to anything
yet. Accepting a value and ignoring it silently is indistinguishable
from a bug at the call site, so setting one warns.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import pytest

from pylmrob import Control, lmrob

# Absolute, not relative: cibuildwheel runs ``pytest {project}/tests/unit``
# from the wheel build's working directory, not the repo root, so a
# relative "tests/data/..." raises FileNotFoundError there. That broke
# the v0.5.25 and v0.5.26 wheel builds, so neither release published.
REPO_ROOT = Path(__file__).resolve().parents[2]
STACKLOSS = REPO_ROOT / "tests" / "data" / "stackloss.csv"


def _stackloss() -> pd.DataFrame:
    if not STACKLOSS.exists():  # pragma: no cover - sdist without test data
        pytest.skip(f"data file missing: {STACKLOSS}")
    return pd.read_csv(STACKLOSS)


_UNIMPLEMENTED = [
    ("trace_lev", 4),
    ("eps_outlier", 1e-3),
    ("eps_x", 1e-9),
    ("solve_tol", 1e-9),
]


@pytest.mark.parametrize("field,value", _UNIMPLEMENTED)
def test_unimplemented_control_warns(field: str, value: float) -> None:
    with pytest.warns(UserWarning, match="not implemented yet"):
        Control(**{field: value})


@pytest.mark.parametrize("field,value", _UNIMPLEMENTED)
def test_warning_names_the_field(field: str, value: float) -> None:
    with pytest.warns(UserWarning) as record:
        Control(**{field: value})
    assert field in str(record[0].message)


def test_default_control_is_silent() -> None:
    """The common path must not warn, or the warning becomes noise."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        Control()
        Control(setting="KS2014")
        Control(nResample=1000, psi="lqq")


def test_nonconvergence_warns_and_does_not_raise() -> None:
    """plan.md 5.3: warn, set converged_=False, do not raise."""
    df = _stackloss()
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    # One IRWLS iteration cannot converge from the S estimate.
    with pytest.warns(RuntimeWarning, match="did not converge"):
        fit = lmrob(formula, df, control=Control(max_it=1, nResample=100), seed=42)
    assert fit.converged_ is False
    assert fit.coef_.shape == (4,)


def test_converged_fit_does_not_warn_about_convergence() -> None:
    df = _stackloss()
    formula = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        fit = lmrob(formula, df, control=Control(nResample=500), seed=42)
    assert fit.converged_
