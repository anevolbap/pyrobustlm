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
InitMethod = Literal["auto", "S", "M-S", "L1"]
Setting = Literal["KS2011", "KS2014", "MM"]


# Default tuning constants matching R's robustbase 0.99-7. (See R's
# .Mpsi.tuning.default and .Mchi.tuning.default for the canonical table.)
# Default tuning constants in **R's internal form** (post-``.psi.conv.cc``).
# Source: ``robustbase:::.psi.conv.cc`` applied to the user-facing defaults
# from ``robustbase::.Mpsi.tuning.default()``.
_DEFAULT_TUNING_PSI: dict[str, tuple[float, ...]] = {
    "huber": (1.345,),
    "bisquare": (4.685061,),
    "hampel": (1.5 * 0.9016085, 3.5 * 0.9016085, 8.0 * 0.9016085),
    "optimal": (1.060158,),
    # ggw: case index 4 = (b=1.5, 95% efficiency).
    "ggw": (4,),
    # lqq: (b, c, s) - converted from user-facing (-0.5, 1.5, 0.95, NA).
    "lqq": (1.4734061, 0.9822707, 1.5),
}

_DEFAULT_TUNING_CHI: dict[str, tuple[float, ...]] = {
    "huber": (0.6745,),
    "bisquare": (1.547645,),
    "hampel": (1.5 * 0.2119163, 3.5 * 0.2119163, 8.0 * 0.2119163),
    "optimal": (0.4047,),
    # ggw: case index 6 = (b=1.5, breakdown=0.5).
    "ggw": (6,),
    # lqq: (b, c, s) - converted from user-facing (-0.5, 1.5, NA, 0.5).
    "lqq": (0.4015457, 0.2676971, 1.5),
}


@dataclass
class Control:
    """Parameters controlling an ``lmrob`` fit.

    Defaults follow R's ``lmrob.control(setting="KS2014")``.
    """

    # ``setting=None`` (and ``"MM"``) follow R's plain default: psi="bisquare",
    # method="MM", cov=".vcov.avar1". ``setting="KS2014"`` switches to
    # psi="lqq", method="SMDM" (with D-step), cov=".vcov.w".
    # ``setting="KS2011"`` is psi="lqq", method="MM", cov=".vcov.w".
    # Pass ``psi`` / ``method`` / ``cov`` explicitly to override.
    setting: Setting | None = None
    psi: PsiFamily | None = None
    tuning_chi: float | tuple[float, ...] | None = None
    tuning_psi: float | tuple[float, ...] | None = None

    init: InitMethod = "S"
    method: str | None = None

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

    cov: str | None = None
    eps_outlier: float | None = None
    eps_x: float | None = None

    seed: int | None = None
    trace_lev: int = 0

    # Worker threads for the fast-S resampling loop.
    # 1 = serial (default; deterministic). 0 = auto (only enables threading
    # for problems large enough that BLAS dominates Python/GIL overhead).
    # Any positive integer = explicit worker count.
    n_workers: int = 1
    # When True, draw resampling subsets inside the Cython kernel using
    # numpy's BitGenerator C API (Floyd's combination algorithm). About
    # 1.4-2x faster at n<=500 because the rank-check SVD goes away.
    # Off by default: the draw sequence is not byte-identical with
    # ``np.random.Generator.choice`` so the basin of attraction can shift
    # slightly, which changes which fits beyond rtol=1e-3 vs R on some
    # small-n datasets. Opt in when you want raw speed and can tolerate
    # the drift; not recommended for reproducibility-sensitive workloads.
    fast_rng: bool = False

    bb: float = 0.5  # consistency constant (target value of mean(chi))

    extra: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Apply named-setting defaults, mirroring R's lmrob.control(setting=...)
        # for fields the user passed as None.
        if self.setting in ("KS2014", "KS2011"):
            if self.psi is None:
                self.psi = "lqq"
            if self.method is None:
                # Both KS2011 and KS2014 use the SMDM pipeline in
                # robustbase. They differ in tuning constants and
                # ``cov.corrfact`` defaults.
                self.method = "SMDM"
            if self.cov is None:
                self.cov = ".vcov.w"
        else:
            if self.psi is None:
                self.psi = "bisquare"
            if self.method is None:
                self.method = "MM"
            if self.cov is None:
                self.cov = ".vcov.avar1"

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
        if setting == "KS2014":
            ctrl = cls(setting="KS2014", psi="bisquare", cov=".vcov.avar1", init="S")
        elif setting == "KS2011":
            ctrl = cls(setting="KS2011", psi="bisquare", cov=".vcov.w", init="S")
        elif setting == "MM":
            ctrl = cls(setting="MM", psi="bisquare", cov=".vcov.avar1", init="S")
        else:
            raise ValueError(f"unknown setting: {setting!r}")
        for key, value in overrides.items():
            if not hasattr(ctrl, key):
                raise TypeError(f"unknown Control field: {key!r}")
            setattr(ctrl, key, value)
        return ctrl
