# SPDX-License-Identifier: GPL-3.0-or-later
"""``diagnostics.plot`` actually runs.

``docs/r-source-map.md`` advertises ``pylmrob.diagnostics.plot(fit)`` as
the ``plot.lmrob(fit)`` equivalent, but nothing imported or called it:
matplotlib is an optional runtime dependency, so it sat outside every
test path. A user following the mapping was the first to execute it.

These are smoke tests, not image comparisons. They check the function
builds the four panels R's ``plot.lmrob`` provides and does not raise on
the shapes it will actually be handed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")  # headless: no display on CI

from pylmrob import Control, lmrob  # noqa: E402
from pylmrob.diagnostics import plot  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
STACKLOSS = REPO_ROOT / "tests" / "data" / "stackloss.csv"
_FORMULA = "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc."


@pytest.fixture
def fit():
    if not STACKLOSS.exists():  # pragma: no cover
        pytest.skip(f"data file missing: {STACKLOSS}")
    df = pd.read_csv(STACKLOSS)
    return lmrob(_FORMULA, df, control=Control(nResample=200), seed=42)


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


def test_plot_returns_a_figure(fit) -> None:
    import matplotlib.pyplot as plt

    out = plot(fit)
    assert isinstance(out, plt.Figure), f"expected a Figure, got {type(out)!r}"


def test_plot_draws_four_panels(fit) -> None:
    """R's plot.lmrob is a four-panel diagnostic; so is ours."""
    fig = plot(fit)
    assert len(fig.axes) == 4, f"expected 4 panels, got {len(fig.axes)}"
    for ax in fig.axes:
        assert ax.get_xlabel() or ax.get_title(), "panel has neither title nor x label"


def test_plot_panels_contain_data(fit) -> None:
    """Guard against a panel that is set up but never drawn into."""
    fig = plot(fit)
    for i, ax in enumerate(fig.axes):
        drawn = len(ax.lines) + len(ax.collections) + len(ax.patches)
        assert drawn > 0, f"panel {i} ({ax.get_title()!r}) has nothing drawn on it"


def test_plot_survives_a_perfect_fit() -> None:
    """scale == 0 makes the standardized residuals degenerate.

    An exact fit through more than half the data is a legitimate S
    solution, so plot() has to cope rather than divide by zero.
    """
    n = 20
    x = np.arange(n, dtype=float)
    df = pd.DataFrame({"x": x, "y": 1.0 + 2.0 * x})
    fit = lmrob("y ~ x", df, control=Control(nResample=100), seed=0)
    fig = plot(fit)
    assert len(fig.axes) == 4


def test_plot_raises_a_useful_error_without_matplotlib(fit, monkeypatch) -> None:
    """The ImportError should name the extra, not just fail on import."""
    import importlib

    def _blocked(name, *a, **k):
        if name.startswith("matplotlib"):
            raise ImportError("No module named 'matplotlib'")
        return importlib.import_module(name, *a, **k)

    monkeypatch.setattr(importlib, "import_module", _blocked)
    with pytest.raises(ImportError, match="matplotlib"):
        plot(fit)
