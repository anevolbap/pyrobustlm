# SPDX-License-Identifier: GPL-3.0-or-later
"""Robust Wald and Deviance tests for nested ``LmRobResults`` models.

Mirrors ``robustbase:::anova.lmrob``.

- **Wald** (default): ``chi2 = beta_drop' V[drop, drop]^{-1} beta_drop`` where
  ``beta_drop`` is the subvector of the largest model's coefficients that the
  reduced model omits and ``V`` is the largest model's robust covariance.
- **Deviance**: refits the reduced model via M-iteration at the full model's
  scale, then computes ``T = 2 * tauStar * (sum_rho_reduced - sum_rho_full)``
  with ``tauStar = mean(psi'(r/s0)) / mean(psi(r/s0)^2)`` evaluated at the
  full residuals.

P-value is from the chi-squared distribution with ``df = ndrop``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pylmrob.results import LmRobResults


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
        title = "Robust Wald Test Table" if self.test == "Wald" else "Robust Deviance Table"
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


def _deviance_pair(full: LmRobResults, reduced: LmRobResults) -> tuple[int, float, int, float]:
    """Deviance chi-sq for one nested pair.

    Refits the reduced model via M-iteration at the full's scale, then
    computes ``T = 2 * tauStar * (RM_sRho - FM_sRho)``.
    """
    from scipy.stats import chi2 as chi2_dist

    from pylmrob import psi as _psi_mod
    from pylmrob._mm import mm_iterate

    if full.design_x_ is None or full.design_y_ is None:
        raise RuntimeError(
            "anova(test='Deviance') requires the full fit to carry its design "
            "matrix; refit with the current pylmrob to populate design_x_."
        )
    if reduced.design_x_ is None:
        raise RuntimeError("reduced fit is missing design_x_; refit and try again")

    method = full.control.method or "MM"
    if not method.endswith("M"):
        raise ValueError(
            f"anova(test='Deviance') requires the full model's method to end "
            f"with 'M'; got method={method!r}."
        )

    drop = _h0_indices(list(full.term_names_), list(reduced.term_names_))
    df = len(drop)
    df_full = full.coef_.size
    n = int(full.nobs_)

    keep_idx = [i for i, t in enumerate(full.term_names_) if t in set(reduced.term_names_)]
    x_full = np.asarray(full.design_x_, dtype=np.float64)
    y = np.asarray(full.design_y_, dtype=np.float64)
    x_red = x_full[:, keep_idx]

    s0 = float(full.scale_)
    if s0 <= 0:
        raise ValueError("full model has non-positive scale; deviance test undefined")

    psi_family = full.control.psi or "bisquare"
    psi_k = tuple(np.atleast_1d(np.asarray(full.control.tuning_psi, dtype=float)).ravel())

    beta_red_init = np.asarray(reduced.coef_, dtype=np.float64).copy()
    if beta_red_init.size != x_red.shape[1]:
        raise ValueError(
            f"reduced fit has {beta_red_init.size} coefficients but the kept "
            f"columns of the full design have {x_red.shape[1]} (term ordering mismatch)"
        )

    rm = mm_iterate(
        X=x_red,
        y=y,
        beta_init=beta_red_init,
        sigma=s0,
        psi_family=psi_family,
        psi_k=psi_k,
        max_it=full.control.max_it,
        rel_tol=full.control.rel_tol,
    )
    rm_res = y - x_red @ rm.coef
    fm_res = np.asarray(full.residuals_, dtype=np.float64)

    # Our ``psi.rho`` is the normalised chi (chi(inf)=1); R's ``Mpsi(deriv=-1)``
    # is the unnormalised rho. Scale by ``rho_inf = 1 / chi_prime_factor`` to
    # recover R's convention before differencing.
    from pylmrob._psifuns import _chi_prime_factor

    psi_k_arr = np.asarray(psi_k, dtype=float)
    rho_inf = 1.0 / _chi_prime_factor(psi_family, psi_k_arr)
    fm_sRho = float(np.sum(_psi_mod.rho(fm_res / s0, psi_family, psi_k))) * rho_inf
    rm_sRho = float(np.sum(_psi_mod.rho(rm_res / s0, psi_family, psi_k))) * rho_inf

    psi_vals = _psi_mod.psi(fm_res / s0, psi_family, psi_k)
    psi_p_vals = _psi_mod.psi_prime(fm_res / s0, psi_family, psi_k)
    denom = float(np.mean(psi_vals**2))
    if denom <= 0:
        raise ValueError("E[psi^2] is zero on full residuals; deviance test undefined")
    tau_star = float(np.mean(psi_p_vals) / denom)

    chi2 = 2.0 * tau_star * (rm_sRho - fm_sRho)
    pval = float(chi2_dist.sf(chi2, df=df))
    pseudo_df = n - df_full + df
    return pseudo_df, chi2, df, pval


def anova(*fits: LmRobResults, test: str = "Wald") -> AnovaTable:
    """Robust nested-model anova on a sequence of ``LmRobResults``.

    The first argument is the largest (full) model; subsequent arguments are
    progressively reduced models. Each adjacent pair must be strictly nested
    via term names.

    Parameters
    ----------
    fits :
        Two or more ``LmRobResults`` ordered from largest to smallest.
    test :
        ``"Wald"`` (default) or ``"Deviance"``. The Deviance test refits the
        reduced model via M-iteration at the full model's scale; it requires
        the full model's method to end with ``"M"``.
    """
    if test not in ("Wald", "Deviance"):
        raise NotImplementedError(f"anova test={test!r}: use 'Wald' or 'Deviance'")
    if len(fits) < 2:
        raise ValueError("anova needs at least two fits")

    full = fits[0]
    nmodels = len(fits)
    table = np.full((nmodels, 4), np.nan, dtype=np.float64)
    table[0, 0] = full.nobs_ - full.coef_.size

    pair_fn = _wald_pair if test == "Wald" else _deviance_pair
    cur_full = full
    for k in range(1, nmodels):
        reduced = fits[k]
        if reduced.nobs_ != full.nobs_:
            raise ValueError("all fits in anova() must be on the same data (different nobs)")
        pseudo_df, chi2, df, pval = pair_fn(cur_full, reduced)
        table[k, 0] = pseudo_df
        table[k, 1] = chi2
        table[k, 2] = df
        table[k, 3] = pval
        cur_full = reduced

    method = full.control.method or "MM"
    term_lists = [list(f.term_names_) for f in fits]
    return AnovaTable(table=table, term_lists=term_lists, test=test, method=method)


__all__ = ["AnovaTable", "anova"]
