# SPDX-License-Identifier: GPL-3.0-or-later
"""Diagnostic plots and influence statistics for ``lmrob`` fits.

Mirrors a subset of R's ``plot.lmrob``. matplotlib is lazily imported so the
package does not require it at import time.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pylmrob.results import LmRobResults


def plot(results: LmRobResults) -> object:
    """Four-panel diagnostic plot.

    Panels: (a) Residuals vs Fitted, (b) Standardized residuals Q-Q,
    (c) Robust weights vs index, (d) Residuals vs Leverage proxy.

    Requires ``matplotlib`` (an optional install dep). Install with
    ``pip install pylmrob[plot]`` or just add matplotlib to your env.
    """
    import importlib

    try:
        plt = importlib.import_module("matplotlib.pyplot")
    except ImportError as exc:
        raise ImportError(
            "diagnostics.plot requires matplotlib; install with "
            "`pip install matplotlib` or `pip install pylmrob[plot]`."
        ) from exc
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


def hatvalues(results: LmRobResults, X: np.ndarray) -> np.ndarray:
    """Diagonal of the (robust) hat matrix ``H = X (X' W X)^-1 X' W``.

    Mirrors ``robustbase:::.lmrob.hat``: compute the weighted QR of
    ``sqrt(w) * X`` and read the row sums of squares of the orthogonal
    factor's leading columns. The result is clamped to [0, 1].

    Parameters
    ----------
    results :
        Fitted ``LmRobResults`` object (uses ``rweights_``).
    X :
        Design matrix used in the fit, shape ``(n, p)``.
    """
    w = results.rweights_
    sw = np.sqrt(np.maximum(w, 0.0))
    Xw = np.asarray(X, dtype=np.float64) * sw[:, None]
    Q, _ = np.linalg.qr(Xw, mode="reduced")
    h = np.minimum(1.0, np.sum(Q * Q, axis=1))
    return h


def cooks_distance(
    results: LmRobResults,
    X: np.ndarray,
    robust: bool = True,
) -> np.ndarray:
    """Robust Cook's distance per observation.

    .. math::
        D_i = \\frac{w_i\\,r_i^2}{p\\,\\sigma^2\\,(1 - h_i)^2}\\,
              \\frac{h_i}{1 - h_i}

    where ``h_i`` are the robust hat values from :func:`hatvalues`.
    """
    if not robust:
        raise NotImplementedError("Only robust=True is supported in v1.")
    p = int(results.coef_.size)
    h = hatvalues(results, X)
    h_safe = np.minimum(h, 1.0 - 1e-12)
    sigma = max(results.scale_, 1e-300)
    r = results.residuals_
    w = results.rweights_
    return (w * r * r) / (p * sigma * sigma * (1.0 - h_safe) ** 2) * (h_safe / (1.0 - h_safe))
