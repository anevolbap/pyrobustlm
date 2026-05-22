# SPDX-License-Identifier: GPL-3.0-or-later
"""``summary.lmrob``-equivalent output for ``LmRobResults``.

Mirrors ``robustbase:::summary.lmrob`` and ``robustbase:::print.summary.lmrob``:

- Coefficient table with ``Estimate``, ``Std. Error``, ``t value``, ``Pr(>|t|)``
  using the t-distribution with ``df.residual`` degrees of freedom.
- Robust residual standard error = ``scale``.
- Robust multiple R-squared and adjusted R-squared computed from the rweights-
  weighted total and residual sum-of-squares with the family-specific
  correction factor ``E[wgt(r)] / E[r psi(r)]``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pylmrob.control import Control
    from pylmrob.results import LmRobResults


# R's summary.lmrob: tabulated correc factor for default tunings of each family.
# Values copied verbatim from robustbase R/lmrob.MM.R::summary.lmrob.
_R_SQUARED_CORREC_TABLE: dict[str, float] = {
    "bisquare": 1.207617,
    "biweight": 1.207617,
    "welsh": 1.224617,
    "optimal": 1.068939,
    "hampel": 1.166891,
    "lqq": 1.159232,
}


def _ggw_correc(tuning: tuple[float, ...]) -> float | None:
    """R's summary.lmrob ggw special-case: four hand-tabulated tunings, else None."""
    if len(tuning) != 4:
        return None
    table = {
        (-0.5, 1.0, 0.95): 1.121708,
        (-0.5, 1.5, 0.95): 1.163192,
        (-0.5, 1.0, 0.85): 1.33517,
        (-0.5, 1.5, 0.85): 1.395828,
    }
    key = (float(tuning[0]), float(tuning[1]), float(tuning[2]))
    return table.get(key)


def _numeric_correc(psi_family: str, tuning: tuple[float, ...]) -> float:
    """Generic ``E[wgt(r)] / E[r psi(r)]`` under the standard normal."""
    from scipy import integrate

    from pylmrob.psi import psi as _psi
    from pylmrob.psi import wgt as _wgt

    inv_sqrt2pi = 1.0 / np.sqrt(2.0 * np.pi)

    def num(t: float) -> float:
        return (
            float(_wgt(np.array([t]), psi_family, tuning)[0]) * np.exp(-0.5 * t * t) * inv_sqrt2pi
        )

    def den(t: float) -> float:
        psi_val = float(_psi(np.array([t]), psi_family, tuning)[0])
        return t * psi_val * np.exp(-0.5 * t * t) * inv_sqrt2pi

    enum, _ = integrate.quad(num, -20.0, 20.0, epsabs=1e-10)
    eden, _ = integrate.quad(den, -20.0, 20.0, epsabs=1e-10)
    if eden <= 0:
        return 1.0
    return enum / eden


def _is_default_tuning(psi_family: str, tuning: tuple[float, ...]) -> bool:
    from pylmrob.psi import _PSI_TUNING_DEFAULT_PSI

    default = _PSI_TUNING_DEFAULT_PSI.get(psi_family.lower())
    if default is None or len(default) != len(tuning):
        return False
    return all(abs(a - b) < 1e-10 for a, b in zip(default, tuning, strict=True))


def _correc_factor(psi_family: str, tuning: tuple[float, ...]) -> float:
    fam = psi_family.lower()
    if fam == "ggw":
        # R uses a literal-tuning lookup for ggw, independent of our default-table check.
        ggw = _ggw_correc(tuning)
        if ggw is not None:
            return ggw
        return _numeric_correc(psi_family, tuning)
    if fam in _R_SQUARED_CORREC_TABLE and _is_default_tuning(psi_family, tuning):
        return _R_SQUARED_CORREC_TABLE[fam]
    return _numeric_correc(psi_family, tuning)


@dataclass
class SummaryLmRob:
    """``summary.lmrob`` analogue. Stringifies to an R-style printout."""

    coefficients: np.ndarray
    term_names: list[str]
    scale: float
    r_squared: float
    adj_r_squared: float
    df_residual: int
    nobs: int
    residuals: np.ndarray
    rweights: np.ndarray
    cov: np.ndarray
    converged: bool
    n_iter: int
    control: Control
    has_intercept: bool

    style: str = "r"  # "r" (default) or "statsmodels"
    # ``detail="brief"`` (default) matches the historic output. ``"full"``
    # appends a small footer with init method, init scale, init iter count,
    # MM iter count, and the engine/RNG that drove the fit; useful for
    # debugging convergence issues without ``trace_lev=3``.
    detail: str = "brief"
    init_info: dict[str, object] | None = None
    engine_c: bool | None = None
    rng: str | None = None

    def __str__(self) -> str:
        return self.render(style=self.style, detail=self.detail)

    def render(self, *, style: str = "r", detail: str = "brief") -> str:
        """Render the summary table.

        Parameters
        ----------
        style
            ``"r"`` (default): R-style ``summary.lmrob`` output.
            ``"statsmodels"``: a fixed-width table matching the
            ``statsmodels.iolib.summary.Summary`` layout for users who
            pipe pylmrob fits into statsmodels-shaped reporting code.
        detail
            ``"brief"`` (default) emits the standard summary. ``"full"``
            appends a footer with the init method, init scale, MM iter
            count, and engine settings (engine_c, rng). Use this when
            debugging convergence.
        """
        if detail not in ("brief", "full"):
            raise ValueError(f"detail must be 'brief' or 'full', got {detail!r}")
        if style == "statsmodels":
            text = self._render_statsmodels()
        elif style == "r":
            text = self._render_r()
        else:
            raise ValueError(f"style must be 'r' or 'statsmodels', got {style!r}")
        if detail == "full":
            text = text + "\n\n" + self._render_full_footer()
        return text

    def _render_full_footer(self) -> str:
        """Diagnostic footer appended when ``detail='full'``."""
        rows: list[str] = ["--- Fit details ---"]
        if self.init_info is not None:
            method = self.init_info.get("method", "?")
            iscale = self.init_info.get("scale", float("nan"))
            iiter = self.init_info.get("n_iter", "?")
            rows.append(f"  init method      : {method}")
            try:
                rows.append(f"  init scale       : {float(iscale):.6g}")  # ty: ignore[invalid-argument-type]
            except (TypeError, ValueError):
                rows.append(f"  init scale       : {iscale}")
            rows.append(f"  init n_iter      : {iiter}")
        rows.append(f"  MM n_iter        : {self.n_iter}")
        rows.append(f"  scale (M-est)    : {self.scale:.6g}")
        if self.engine_c is not None:
            rows.append(f"  engine_c         : {self.engine_c}")
        if self.rng is not None:
            rows.append(f"  rng              : {self.rng}")
        rows.append(f"  psi family       : {self.control.psi or '?'}")
        if self.control.tuning_psi is not None:
            rows.append(f"  tuning_psi       : {self.control.tuning_psi}")
        if self.control.tuning_chi is not None:
            rows.append(f"  tuning_chi       : {self.control.tuning_chi}")
        return "\n".join(rows)

    def _render_r(self) -> str:
        rows: list[str] = []
        rows.append(
            f"lmrob(method={self.control.method!r}, psi={self.control.psi!r}, "
            f"setting={self.control.setting!r})"
        )
        rows.append("")
        rows.append("Residuals:")
        if self.df_residual > 5:
            qs = np.quantile(self.residuals, [0.0, 0.25, 0.5, 0.75, 1.0])
            qline = "    Min      1Q  Median      3Q     Max"
            vline = "  ".join(f"{q:7.4g}" for q in qs)
            rows.append(qline)
            rows.append(vline)
        else:
            rows.append("  ".join(f"{r:7.4g}" for r in self.residuals))
        rows.append("")
        rows.append("Coefficients:")
        header = f"{'':<20} {'Estimate':>11} {'Std. Error':>12} {'t value':>9} {'Pr(>|t|)':>10}"
        rows.append(header)
        for name, row in zip(self.term_names, self.coefficients, strict=True):
            est, se, tv, pv = row
            rows.append(f"{name:<20} {est:11.5g} {se:12.5g} {tv:9.4g} {pv:10.4g}")
        rows.append("")
        rows.append(f"Robust residual standard error: {self.scale:.6g}")
        if self.has_intercept:
            rows.append(
                f"Multiple R-squared: {self.r_squared:.6g},  "
                f"Adjusted R-squared: {self.adj_r_squared:.6g}"
            )
        rows.append(
            f"Convergence in {self.n_iter} IRWLS iterations"
            if self.converged
            else "Algorithm did not converge"
        )
        rows.append(f"n = {self.nobs}, df residual = {self.df_residual}")
        return "\n".join(rows)

    def _render_statsmodels(self) -> str:
        """statsmodels-style fixed-width table."""
        from scipy.stats import norm

        z = norm.ppf(0.975)
        bar = "=" * 78
        sep = "-" * 78
        rows: list[str] = []
        rows.append(bar)
        rows.append(f"{'lmrob (MM-estimator)':^78}")
        rows.append(bar)
        rows.append(
            f"{'Method:':<20}{self.control.method or '?':<19}"
            f"{'Psi:':<20}{self.control.psi or '?':<19}"
        )
        rows.append(
            f"{'No. Observations:':<20}{self.nobs:<19}{'Df Residuals:':<20}{self.df_residual:<19}"
        )
        rows.append(
            f"{'R-squared (robust):':<20}{self.r_squared:<19.4f}"
            f"{'Adj. R-squared:':<20}{self.adj_r_squared:<19.4f}"
        )
        rows.append(
            f"{'Robust scale:':<20}{self.scale:<19.6g}{'Converged:':<20}{self.converged!s:<19}"
        )
        rows.append(sep)
        rows.append(
            f"{'':<20}{'coef':>10}{'std err':>10}{'t':>8}{'P>|t|':>9}{'[0.025':>10}{'0.975]':>10}"
        )
        rows.append(sep)
        for name, row in zip(self.term_names, self.coefficients, strict=True):
            est, se, tv, pv = row
            lo = est - z * se
            hi = est + z * se
            rows.append(
                f"{name:<20}{est:>10.4f}{se:>10.4f}{tv:>8.3f}{pv:>9.4f}{lo:>10.4f}{hi:>10.4f}"
            )
        rows.append(bar)
        return "\n".join(rows)

    def __repr__(self) -> str:
        return self.__str__()

    def __contains__(self, item: object) -> bool:
        if not isinstance(item, str):
            return False
        return item in str(self)


def make_summary(result: LmRobResults) -> SummaryLmRob:
    """Build a ``SummaryLmRob`` from an ``LmRobResults``."""
    from scipy.stats import t as t_dist

    coef = np.asarray(result.coef_, dtype=np.float64)
    cov = np.asarray(result.cov_, dtype=np.float64)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    df = int(result.df_residual_)

    if result.converged_:
        with np.errstate(divide="ignore", invalid="ignore"):
            tval = np.where(se > 0, coef / np.where(se > 0, se, 1.0), np.nan)
            pval = 2.0 * t_dist.sf(np.abs(tval), df=df)
    else:
        tval = np.full_like(coef, np.nan)
        pval = np.full_like(coef, np.nan)

    coef_table = np.column_stack([coef, se, tval, pval])

    has_intercept = bool(
        result.term_names_ and result.term_names_[0] in {"Intercept", "(Intercept)"}
    )

    psi_family = result.control.psi or "bisquare"
    psi_k = tuple(np.atleast_1d(np.asarray(result.control.tuning_psi, dtype=float)).ravel())
    correc = _correc_factor(psi_family, psi_k)

    y = result.fitted_ + result.residuals_
    if has_intercept and df > 0:
        wsum = float(np.sum(result.rweights_))
        resp_mean = float(np.sum(result.rweights_ * y) / wsum) if wsum > 0 else 0.0
        yMy = float(np.sum(result.rweights_ * (y - resp_mean) ** 2))
        rMr = float(np.sum(result.rweights_ * result.residuals_**2))
        denom = yMy + rMr * (correc - 1.0)
        r2 = (yMy - rMr) / denom if denom > 0 else 0.0
        df_int = 1
        adj = 1.0 - (1.0 - r2) * ((result.nobs_ - df_int) / df) if df > 0 else 0.0
    else:
        r2 = 0.0
        adj = 0.0

    return SummaryLmRob(
        coefficients=coef_table,
        term_names=list(result.term_names_),
        scale=float(result.scale_),
        r_squared=float(r2),
        adj_r_squared=float(adj),
        df_residual=df,
        nobs=int(result.nobs_),
        residuals=np.asarray(result.residuals_, dtype=np.float64),
        rweights=np.asarray(result.rweights_, dtype=np.float64),
        cov=cov,
        converged=bool(result.converged_),
        n_iter=int(result.n_iter_),
        control=result.control,
        has_intercept=has_intercept,
        init_info=dict(result.init_) if result.init_ else None,
        engine_c=bool(result.control.engine_c),
        rng=str(result.control.rng),
    )


__all__ = ["SummaryLmRob", "make_summary"]
