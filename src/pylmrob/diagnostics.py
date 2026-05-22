# SPDX-License-Identifier: GPL-3.0-or-later
"""Diagnostic plots and influence statistics for ``lmrob`` fits.

Mirrors a subset of R's ``plot.lmrob``. matplotlib is lazily imported so the
package does not require it at import time.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pylmrob.results import LmRobResults


@dataclass
class DiagnosticsTable:
    """Per-observation diagnostic statistics returned by ``fit.diagnostics()``.

    Attributes
    ----------
    leverage :
        Diagonal of the robust hat matrix ``H = X (X' W X)^-1 X' W``.
        Values in ``[0, 1]``; the average is ``p / n``.
    cooks_distance :
        Robust Cook's distance per observation; large values mark
        influential points.
    std_residuals :
        Residuals divided by the M-scale ``sigma``. ``|z| > 2.5`` is a
        common informal outlier flag.
    rweights :
        Robustness weights ``psi(r/sigma) / (r/sigma)`` from the fit.
        Observations with weight close to zero were effectively dropped.
    outliers :
        Boolean mask, ``|std_residuals| > outlier_threshold``. Default
        threshold is ``2.5``.
    masked_outliers :
        Boolean mask flagging high-leverage observations that were
        also nearly fully downweighted by the robust fit, i.e.
        ``leverage > leverage_threshold`` AND ``rweights < weight_eps``.
        These rows are simultaneously influential in the design and
        rejected as outliers. They merit a closer look than rows
        flagged by either criterion alone.
    dfbetas :
        ``(n, p)`` matrix where ``dfbetas[i, j]`` approximates the
        change in coefficient ``j`` (in SE units) when observation ``i``
        is removed. See :func:`dfbetas` for the formula.
    """

    leverage: np.ndarray
    cooks_distance: np.ndarray
    std_residuals: np.ndarray
    rweights: np.ndarray
    outliers: np.ndarray
    masked_outliers: np.ndarray
    dfbetas: np.ndarray


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


def dfbetas(results: LmRobResults, X: np.ndarray) -> np.ndarray:
    """Per-observation influence on each regression coefficient.

    Returns an ``(n, p)`` matrix where ``dfbetas[i, j]`` approximates
    ``(beta_full[j] - beta_minus_i[j]) / SE(beta[j])``. Large absolute
    values mark observations whose deletion would visibly move
    coefficient ``j``.

    Uses Belsley/Kuh/Welsch's closed-form OLS-style formula, adapted
    with the M-estimator's robust weights so the influence reflects
    only the rows that the fit actually uses:

    .. math::
        \\text{dfbetas}_{i,j} \\approx
            \\frac{w_i \\cdot (X' W X)^{-1}_{j,\\cdot} \\, x_i \\, r_i}
                  {(1 - h_i) \\, \\sigma \\, \\sqrt{(X' W X)^{-1}_{j,j}}}

    where ``w_i`` are the robust weights, ``h_i`` the robust hat values,
    ``r_i`` the residual, and ``sigma`` the M-scale.

    This is an approximation; a leave-one-out refit is the exact
    quantity but costs ``n`` extra fits.
    """
    w = results.rweights_
    h = hatvalues(results, X)
    r = results.residuals_
    sigma = max(results.scale_, 1e-300)

    Xw = X * np.sqrt(np.maximum(w, 0.0))[:, None]
    # (X' W X)^-1
    xtwx_inv = np.linalg.inv(Xw.T @ Xw)
    # Standard errors are sqrt(sigma^2 * diag((X' W X)^-1)).
    se = sigma * np.sqrt(np.maximum(np.diag(xtwx_inv), 0.0))

    one_minus_h = np.maximum(1.0 - h, 1e-12)
    # Influence vector per i: ``w_i * (X'WX)^-1 @ x_i * r_i / ((1 - h_i) * sigma)``
    # Result shape: (n, p). Broadcasting form.
    influence = (w * r / (one_minus_h * sigma))[:, None] * (X @ xtwx_inv)
    # Scale each column by 1 / se[j] to get the dfbetas form.
    with np.errstate(divide="ignore", invalid="ignore"):
        scaled = influence / np.where(se > 0, se, 1.0)
    return scaled


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
