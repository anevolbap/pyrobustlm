# SPDX-License-Identifier: GPL-3.0-or-later
"""Per-coefficient outlier statistics, port of ``robustbase::outlierStats``.

A high-breakdown fit can look healthy overall and still have broken down
*locally*: within one level of a factor, most observations are rejected
and the level's coefficient is determined by the handful that remain.
The global scale and the overall robustness weights hide this, because
the other levels dominate the averages.

``outlier_stats`` reports, per indicator column of the design, how many
observations that column touches, how many of those the fit rejected,
and their mean robustness weight. Koller & Stahel's ``setting="KS2014"``
exists partly to avoid this failure mode, which is what R's warning
message points users toward.

Port of ``robustbase::outlierStats`` (robustbase 0.99-7).
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from pylmrob.results import LmRobResults

# R accepts either a number or a function of one argument here.
_EpsSpec = float | Callable[[float], float] | None


def default_eps_outlier(nobs: float) -> float:
    """R's ``lmrob.control()`` default: ``0.1 / nobs``.

    The threshold below which a robustness weight counts as "rejected".
    """
    return 0.1 / float(nobs)


def default_eps_x(max_abs_x: float) -> float:
    """R's default: ``.Machine$double.eps^0.75 * max(abs(x))``."""
    return float(np.finfo(float).eps ** 0.75 * max_abs_x)


@dataclass(frozen=True)
class OutlierStatsRow:
    """One row of the report."""

    name: str
    n_nonzero: int
    n_rejected: int
    ratio: float
    mean_robweight: float


@dataclass(frozen=True)
class OutlierStats:
    """``outlierStats`` report.

    ``flagged`` lists the coefficients whose rejection ratio or mean
    robustness weight crossed the warning limits.
    """

    rows: list[OutlierStatsRow]
    flagged: list[str]
    eps_outlier: float
    eps_x: float

    def __iter__(self):
        return iter(self.rows)

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, key: str) -> OutlierStatsRow:
        for row in self.rows:
            if row.name == key:
                return row
        raise KeyError(key)

    def __str__(self) -> str:
        head = f"{'':<16}{'N.nonzero':>10}{'N.rejected':>12}{'Ratio':>9}{'Mean.RobWeight':>16}"
        body = [
            f"{r.name:<16}{r.n_nonzero:>10d}{r.n_rejected:>12d}"
            f"{r.ratio:>9.4f}{r.mean_robweight:>16.4f}"
            for r in self.rows
        ]
        out = "\n".join([head, *body])
        if self.flagged:
            out += "\nPossible local breakdown of " + ", ".join(f"'{n}'" for n in self.flagged)
        return out

    def __repr__(self) -> str:  # pragma: no cover - convenience
        return self.__str__()


def _resolve(value: _EpsSpec, default: Callable[[float], float], arg: float) -> float:
    """R accepts a number or a function of one argument for these."""
    if value is None:
        return default(arg)
    # Narrow by type rather than ``callable()``: the latter widens to an
    # unknown signature, which the type checker refuses to call.
    if isinstance(value, (int, float)):
        return float(value)
    return float(value(arg))


def outlier_stats(
    results: LmRobResults,
    X: np.ndarray | None = None,
    eps_outlier: _EpsSpec = None,
    eps_x: _EpsSpec = None,
    warn_limit_reject: float | None = 0.5,
    warn_limit_meanrw: float | None = 0.5,
    shout: bool | None = None,
) -> OutlierStats:
    """Per-coefficient outlier statistics.

    Parameters
    ----------
    results :
        A fitted :class:`~pylmrob.results.LmRobResults`.
    X :
        Design matrix. Defaults to the one stashed on the fit.
    eps_outlier :
        Robustness weights below this count as rejected. A number, or a
        callable of ``nobs``. Defaults to ``0.1 / nobs`` (R's default).
    eps_x :
        Design entries with ``abs(x) <= eps_x`` are treated as zero, i.e.
        the observation is not "touched" by that column. A number, or a
        callable of ``max(abs(X))``. Defaults to
        ``eps**0.75 * max(abs(X))`` (R's default).
    warn_limit_reject, warn_limit_meanrw :
        Warn when a column's rejection ratio reaches the first, or its
        mean robustness weight falls to the second. ``None`` disables
        that half of the check.
    shout :
        ``True`` always warns, ``False`` never warns, ``None`` (default)
        warns only when a limit is crossed. Mirrors R's ``shout``.

    Returns
    -------
    OutlierStats

    Notes
    -----
    Only indicator-like columns are reported: a column is included when
    at least one observation has ``abs(x) <= eps_x`` there. That drops
    the intercept and continuous predictors, which touch every row and
    so carry no local information, and keeps factor dummies.
    """
    rw = np.asarray(results.rweights_, dtype=np.float64).ravel()

    if X is None:
        X = results.design_x_
    if X is None:  # pragma: no cover - design is stashed by lmrob()
        raise ValueError("outlier_stats needs the design matrix; pass X explicitly")
    Xa = np.asarray(X, dtype=np.float64)
    if Xa.shape[0] != rw.size:
        raise ValueError(f"X has {Xa.shape[0]} rows but there are {rw.size} robustness weights")

    n = rw.size
    epsw = _resolve(eps_outlier, default_eps_outlier, float(n))
    max_abs = float(np.max(np.abs(Xa))) if Xa.size else 0.0
    epsx = _resolve(eps_x, default_eps_x, max_abs)

    rejected = np.abs(rw) < epsw
    touches = np.abs(Xa) > epsx

    names = list(results.term_names_)
    if len(names) != Xa.shape[1]:  # pragma: no cover - defensive
        names = [f"x{i}" for i in range(Xa.shape[1])]

    def _row(name: str, idx: np.ndarray) -> OutlierStatsRow:
        nnz = int(np.count_nonzero(idx))
        n_rej = int(np.count_nonzero(rejected & idx))
        return OutlierStatsRow(
            name=name,
            n_nonzero=nnz,
            n_rejected=n_rej,
            ratio=(n_rej / nnz) if nnz else float("nan"),
            mean_robweight=float(np.mean(rw[idx])) if nnz else float("nan"),
        )

    rows = [_row("Overall", np.ones(n, dtype=bool))]
    # Columns that touch every observation carry no local information;
    # R selects on ``colSums(xnz) < NROW(xnz)``.
    for j, name in enumerate(names):
        col = touches[:, j]
        if int(np.count_nonzero(col)) < n:
            rows.append(_row(name, col))

    flagged: list[str] = []
    for row in rows:
        hit = False
        if warn_limit_reject is not None and not np.isnan(row.ratio):
            hit = hit or row.ratio >= warn_limit_reject
        if warn_limit_meanrw is not None and not np.isnan(row.mean_robweight):
            hit = hit or row.mean_robweight <= warn_limit_meanrw
        if hit:
            flagged.append(row.name)

    stats = OutlierStats(rows=rows, flagged=flagged, eps_outlier=epsw, eps_x=epsx)

    if shout is not False and (shout is True or flagged):
        method = getattr(results.control, "method", None) or "MM"
        what = (
            (f"{len(flagged)} coefficients" if len(flagged) > 1 else f"coefficient {flagged[0]!r}")
            if flagged
            else "the fit"
        )
        warnings.warn(
            f"Detected possible local breakdown of the {method}-estimate in {what}. "
            "Use setting='KS2014' to avoid this problem.",
            RuntimeWarning,
            stacklevel=2,
        )
    return stats


__all__ = ["OutlierStats", "OutlierStatsRow", "outlier_stats"]
