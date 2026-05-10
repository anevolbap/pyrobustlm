# SPDX-License-Identifier: GPL-3.0-or-later
"""Robust Wald test for nested ``LmRobResults`` models.

Mirrors ``robustbase:::anova.lmrob`` (Wald variant). Compares two or more
nested fits via

    chi2 = beta_drop' V[drop, drop]^{-1} beta_drop

where ``beta_drop`` is the subvector of the largest model's coefficients
that the reduced model omits and ``V`` is the largest model's robust
covariance. P-value is from ``chi2(df=ndrop)``.

The Deviance variant (which refits the reduced model with M-iteration on
the larger model's scale) is not yet ported.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pyrobustlm.results import LmRobResults


@dataclass
class AnovaTable:
    """Robust Wald-test table comparing nested ``lmrob`` fits.

    Columns mirror R's ``anova.lmrob`` output: ``pseudoDf``, ``Test.Stat``,
    ``Df``, ``Pr(>chisq)``.
    """

    table: np.ndarray  # shape (nmodels, 4)
    term_lists: list[list[str]]
    test: str
    method: str

    def __str__(self) -> str:
        title = f"Robust {self.test} Test Table"
        rows: list[str] = [title, ""]
        for i, terms in enumerate(self.term_lists, start=1):
            rows.append(f"Model {i}: {' + '.join(terms) if terms else '(no terms)'}")
        rows.append(f"Largest model fitted by lmrob(), i.e. {self.method}")
        rows.append("")
        header = f"{'':<3} {'pseudoDf':>9} {'Test.Stat':>11} {'Df':>4} {'Pr(>chisq)':>11}"
        rows.append(header)
        for i, row in enumerate(self.table, start=1):
            pdf, ts, df, pv = row
            if np.isnan(ts):
                rows.append(f"{i:<3} {pdf:9.0f} {'':>11} {'':>4} {'':>11}")
            else:
                rows.append(f"{i:<3} {pdf:9.0f} {ts:11.4g} {df:4.0f} {pv:11.4g}")
        return "\n".join(rows)

    def __repr__(self) -> str:
        return self.__str__()


def _h0_indices(full_terms: list[str], reduced_terms: list[str]) -> list[int]:
    """Indices in ``full_terms`` that are absent from ``reduced_terms``.

    Reduced model must be strictly nested in the full one (every term in the
    reduced fit must appear in the full fit).
    """
    full_set = set(full_terms)
    missing = [t for t in reduced_terms if t not in full_set]
    if missing:
        raise ValueError(f"models are not nested; reduced has terms not in full: {missing!r}")
    drop = [i for i, t in enumerate(full_terms) if t not in set(reduced_terms)]
    if not drop:
        raise ValueError("models are not strictly nested (no terms dropped)")
    return drop


def _wald_pair(full: LmRobResults, reduced: LmRobResults) -> tuple[int, float, int, float]:
    """Wald chi-sq for one nested pair.

    Returns ``(pseudoDf, chi2, df_drop, p_value)`` where ``pseudoDf`` is
    R's ``n - p_full + df_drop``.
    """
    from scipy.stats import chi2 as chi2_dist

    drop = _h0_indices(list(full.term_names_), list(reduced.term_names_))
    coef = np.asarray(full.coef_, dtype=np.float64)
    cov = np.asarray(full.cov_, dtype=np.float64)

    h0_coef = coef[drop]
    h0_cov = cov[np.ix_(drop, drop)]
    sol = np.linalg.solve(h0_cov, h0_coef)
    chi2 = float(h0_coef @ sol)

    df = len(drop)
    p = full.coef_.size
    n = full.nobs_
    pseudo_df = n - p + df
    pval = float(chi2_dist.sf(chi2, df=df))
    return pseudo_df, chi2, df, pval


def anova(*fits: LmRobResults, test: str = "Wald") -> AnovaTable:
    """Robust Wald test on a sequence of nested ``LmRobResults``.

    The first argument is the largest (full) model; subsequent arguments are
    progressively reduced models. Each pair must be strictly nested via term
    names.

    Parameters
    ----------
    fits :
        Two or more ``LmRobResults`` ordered from largest to smallest.
    test :
        Currently only ``"Wald"`` is supported.
    """
    if test != "Wald":
        raise NotImplementedError(
            f"anova test={test!r}: only 'Wald' is supported. "
            "The Deviance variant requires re-fitting the reduced model and is not "
            "yet ported."
        )
    if len(fits) < 2:
        raise ValueError("anova needs at least two fits")

    full = fits[0]
    nmodels = len(fits)
    table = np.full((nmodels, 4), np.nan, dtype=np.float64)
    table[0, 0] = full.nobs_ - full.coef_.size

    cur_full = full
    for k in range(1, nmodels):
        reduced = fits[k]
        if reduced.nobs_ != full.nobs_:
            raise ValueError("all fits in anova() must be on the same data (different nobs)")
        pseudo_df, chi2, df, pval = _wald_pair(cur_full, reduced)
        table[k, 0] = pseudo_df
        table[k, 1] = chi2
        table[k, 2] = df
        table[k, 3] = pval
        cur_full = reduced

    method = full.control.method or "MM"
    term_lists = [list(f.term_names_) for f in fits]
    return AnovaTable(table=table, term_lists=term_lists, test=test, method=method)


__all__ = ["AnovaTable", "anova"]
