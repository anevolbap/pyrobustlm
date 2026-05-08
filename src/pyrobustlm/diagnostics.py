# SPDX-License-Identifier: GPL-3.0-or-later
"""Diagnostic plots and influence statistics for ``lmrob`` fits.

Mirrors a subset of R's ``plot.lmrob``. matplotlib is lazily imported so the
package does not require it at import time.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pyrobustlm.results import LmRobResults


def plot(results: LmRobResults) -> object:
    """Four-panel diagnostic plot.

    Panels: (a) Residuals vs Fitted, (b) Standardized residuals Q-Q,
    (c) Robust weights vs index, (d) Residuals vs Leverage proxy.
    """
    import matplotlib.pyplot as plt
    from scipy.stats import norm

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))

    # (a)
    ax = axes[0, 0]
    ax.scatter(results.fitted_, results.residuals_, s=15, alpha=0.7)
    ax.axhline(0.0, color="grey", lw=0.5)
    ax.set_xlabel("Fitted values")
    ax.set_ylabel("Residuals")
    ax.set_title("Residuals vs Fitted")

    # (b)
    ax = axes[0, 1]
    std = results.scale_
    z = results.residuals_ / std if std > 0 else results.residuals_
    z_sorted = np.sort(z)
    n = z.size
    q_theo = norm.ppf((np.arange(1, n + 1) - 0.5) / n)
    ax.scatter(q_theo, z_sorted, s=15, alpha=0.7)
    lo, hi = float(q_theo.min()), float(q_theo.max())
    ax.plot([lo, hi], [lo, hi], color="grey", lw=0.5)
    ax.set_xlabel("Theoretical quantiles")
    ax.set_ylabel("Standardized residuals")
    ax.set_title("Normal Q-Q")

    # (c)
    ax = axes[1, 0]
    ax.bar(np.arange(n), results.rweights_, color="steelblue")
    ax.set_xlabel("Observation index")
    ax.set_ylabel("Robust weight")
    ax.set_title("Robust weights")

    # (d) Residuals vs proxy leverage (sum of weight*x_i^2 for fitted points)
    ax = axes[1, 1]
    ax.scatter(results.rweights_, np.abs(z), s=15, alpha=0.7)
    ax.set_xlabel("Robust weight")
    ax.set_ylabel("|standardized residual|")
    ax.set_title("|residual| vs weight")

    fig.tight_layout()
    return fig


def cooks_distance(results: LmRobResults, robust: bool = True) -> np.ndarray:
    """Robust analogue of Cook's distance.

    For a robust fit, we compute the *weighted* analogue
    ``D_i = (w_i * r_i^2) / (p * sigma^2 * (1 - h_i))``
    where ``h_i`` is the diagonal of the hat-like matrix
    ``H = X (X^T W X)^-1 X^T W`` evaluated at the fit's robust weights.
    """
    if not robust:
        raise NotImplementedError("Only robust=True is supported in v1.")
    raise NotImplementedError(
        "diagnostics.cooks_distance — Phase 9+ (deferred); only the plot "
        "wrapper is implemented in v1. PRs welcome."
    )


def hatvalues(results: LmRobResults) -> np.ndarray:
    """Robust analogue of hat values."""
    raise NotImplementedError("diagnostics.hatvalues — Phase 9+ (deferred). PRs welcome.")
