# SPDX-License-Identifier: GPL-3.0-or-later
"""Control parameters for ``lmrob``.

Mirrors R's ``robustbase::lmrob.control``. The full set of parameters and the
preset values for ``setting in {"KS2011", "KS2014"}`` will be filled in during
Phase 8. This file currently provides the dataclass shape and a stub
constructor so the rest of the package can import it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

PsiFamily = Literal["bisquare", "huber", "hampel", "optimal", "ggw", "lqq"]
InitMethod = Literal["S", "M-S", "L1"]
Setting = Literal["KS2011", "KS2014", "MM"]


@dataclass
class Control:
    """Parameters controlling an ``lmrob`` fit.

    Defaults follow R's ``lmrob.control()`` for ``setting="KS2014"`` once
    Phase 8 fills them in. Until then the dataclass exists only to support
    the import surface; constructing it does not yet validate its contents.
    """

    setting: Setting = "KS2014"
    psi: PsiFamily = "lqq"
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

    extra: dict[str, object] = field(default_factory=dict)

    @classmethod
    def preset(cls, setting: Setting, **overrides: object) -> Control:
        """Build a Control for one of the named presets in R.

        Will be implemented in Phase 8 to match R's defaults exactly. For now
        this returns the dataclass defaults with ``setting`` set.
        """

        ctrl = cls(setting=setting)
        for key, value in overrides.items():
            if not hasattr(ctrl, key):
                raise TypeError(f"unknown Control field: {key!r}")
            setattr(ctrl, key, value)
        return ctrl
