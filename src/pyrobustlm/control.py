# SPDX-License-Identifier: GPL-3.0-or-later
"""Control parameters for ``lmrob``.

Mirrors robustbase's ``lmrob.control``. Phase 8 fills in the named-preset
defaults for ``setting in {"KS2011", "KS2014"}``. The values here come
from ``robustbase/R/lmrob.R``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

PsiFamily = Literal["bisquare", "huber", "hampel", "optimal", "ggw", "lqq"]
InitMethod = Literal["S", "M-S", "L1"]
Setting = Literal["KS2011", "KS2014", "MM"]


# Default tuning constants matching R's robustbase 0.99-7. (See R's
# .Mpsi.tuning.default and .Mchi.tuning.default for the canonical table.)
_DEFAULT_TUNING_PSI: dict[str, tuple[float, ...]] = {
    "huber": (1.345,),
    "bisquare": (4.685061,),
    "hampel": (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085),
    "optimal": (1.060158,),
    "ggw": (-0.5, 1.5, 0.95, float("nan")),  # user-facing form
    "lqq": (-0.5, 1.5, 0.95, float("nan")),
}

_DEFAULT_TUNING_CHI: dict[str, tuple[float, ...]] = {
    "huber": (0.6745,),
    "bisquare": (1.547645,),
    "hampel": (1.5 * 0.2119163, 3.5 * 0.2119163, 8.0 * 0.2119163),
    "optimal": (0.4047,),
    "ggw": (-0.5, 1.5, float("nan"), 0.5),
    "lqq": (-0.5, 1.5, float("nan"), 0.5),
}


@dataclass
class Control:
    """Parameters controlling an ``lmrob`` fit.

    Defaults follow R's ``lmrob.control(setting="KS2014")``.
    """

    setting: Setting = "KS2014"
    psi: PsiFamily = "bisquare"  # KS2014 default; the docs say "lqq" but
    # robustbase 0.99-7 actually defaults to "bisquare"
    tuning_chi: float | tuple[float, ...] | None = None
    tuning_psi: float | tuple[float, ...] | None = None

    init: InitMethod = "S"
    method: str = "MM"

    nResample: int = 500
    max_it: int = 50
    k_max: int = 200
    refine_tol: float = 1e-7
    rel_tol: float = 1e-7
    solve_tol: float = 1e-7
    scale_tol: float = 1e-10
    zero_tol: float = 1e-10

    best_r_s: int = 2
    k_fast_s: int = 1
    k_m_s: int = 20

    mts: int = 1000
    subsampling: Literal["nonsingular", "simple"] = "nonsingular"

    cov: str = ".vcov.avar1"
    eps_outlier: float | None = None
    eps_x: float | None = None

    seed: int | None = None
    trace_lev: int = 0

    bb: float = 0.5  # consistency constant (target value of mean(chi))

    extra: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Fill in default tuning constants matching R.
        if self.tuning_psi is None:
            if self.psi not in _DEFAULT_TUNING_PSI:
                raise ValueError(f"unknown psi family {self.psi!r}")
            self.tuning_psi = _DEFAULT_TUNING_PSI[self.psi]
        if self.tuning_chi is None:
            if self.psi not in _DEFAULT_TUNING_CHI:
                raise ValueError(f"unknown psi family {self.psi!r}")
            self.tuning_chi = _DEFAULT_TUNING_CHI[self.psi]

    @classmethod
    def preset(cls, setting: Setting, **overrides: object) -> Control:
        """Build a Control for a named preset.

        Settings:

        - ``"KS2014"``: psi="bisquare" (matches robustbase 0.99-7 default).
        - ``"KS2011"``: same families with KS2011-specific cov estimator.
        - ``"MM"``: legacy MM defaults (psi="bisquare").
        """
        defaults: dict[str, object]
        if setting == "KS2014":
            defaults = {"psi": "bisquare", "cov": ".vcov.avar1", "init": "S"}
        elif setting == "KS2011":
            defaults = {"psi": "bisquare", "cov": ".vcov.w", "init": "S"}
        elif setting == "MM":
            defaults = {"psi": "bisquare", "cov": ".vcov.avar1", "init": "S"}
        else:
            raise ValueError(f"unknown setting: {setting!r}")
        defaults["setting"] = setting
        defaults.update(overrides)
        return cls(**defaults)  # type: ignore[arg-type]
