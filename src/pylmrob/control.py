# SPDX-License-Identifier: GPL-3.0-or-later
"""Control parameters for ``lmrob``.

Mirrors robustbase's ``lmrob.control``. Phase 8 fills in the named-preset
defaults for ``setting in {"KS2011", "KS2014"}``. The values here come
from ``robustbase/R/lmrob.R``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any, Literal

PsiFamily = Literal["bisquare", "huber", "hampel", "optimal", "ggw", "lqq", "welsh"]
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
    "welsh": (2.11,),
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
    "welsh": (0.5773502691896258,),  # = 1/sqrt(3); gives E[chi(Z)]=0.5 under Z~N(0,1).
    # ggw: case index 6 = (b=1.5, breakdown=0.5).
    "ggw": (6,),
    # lqq: (b, c, s) - converted from user-facing (-0.5, 1.5, NA, 0.5).
    "lqq": (0.4015457, 0.2676971, 1.5),
}


@dataclass
class Control:
    """Parameters controlling an ``lmrob`` fit.

    Defaults follow R's plain ``lmrob.control()``: psi="bisquare",
    method="MM", cov=".vcov.avar1". Pass ``setting="KS2011"`` or
    ``setting="KS2014"`` for the Koller-Stahel presets (both are
    psi="lqq", method="SMDM", cov=".vcov.w" in robustbase 0.99-7).
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
    # Gauss-Hermite nodes used for the E[...] expectations in the D-step
    # (kappa). Matches R's ``lmrob.control(numpoints = 10)``; R's
    # ``lmrob.E`` uses this rule rather than exact integration, so the
    # node count is part of the answer, not just its accuracy.
    numpoints: int = 10
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
    # BitGenerator backing the resample RNG.
    #
    # - ``"PCG64"`` (default): the modern NumPy default. Fast and
    #   statistically better than MT19937, but produces a different
    #   subset-draw sequence than R.
    # - ``"MT19937"``: NumPy's Mersenne Twister. Closer to R's family
    #   but the seed-to-state path and the 64-bit-vs-32-bit raw output
    #   still differ from R's ``RNG.c``.
    # - ``"R"``: drive the resample loop from ``pylmrob.r_set_seed`` and
    #   ``r_sample_noreplace``, which is byte-identical to robustbase's
    #   draw sequence. Implies ``n_workers=1``, ``engine_c=False``, and
    #   ``subsampling="simple"`` (set automatically); R's
    #   ``subsampling="nonsingular"`` LU-pivot path isn't yet ported.
    rng: Literal["PCG64", "MT19937", "R"] = "PCG64"
    # Use the monolithic Cython engine (pylmrob._core._lmrob). When True,
    # fast-S + MM + vcov_avar1 run in one nogil C block with one workspace
    # allocation. On small n this is 5-10x the default Python path; on
    # large n ``lmrob()`` auto-falls-back to the threaded default path
    # (since the monolithic kernel is a single C call that does not
    # parallelise). The Cython subset-draw is not byte-identical to
    # ``np.random.choice``, so the two engines can land in different
    # basins of attraction; they agree element-wise on 98 of 100 fits
    # over the classical corpus. ``lmrob()`` still catches a singular
    # ``X'WX`` from this path and refits with ``engine_c=False``, but
    # that is now a warned last resort rather than a routine event.
    engine_c: bool = True

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

        # rng="R" implies a serial, R-call-order resample path. Force
        # the matching control values now so the rest of the pipeline
        # doesn't have to special-case combinations that don't make sense.
        # Both "simple" and "nonsingular" subsampling are now supported
        # via r_set_seed + r_sample_noreplace / r_subsample_nonsingular,
        # so we honour whichever the user requested.
        if self.rng == "R":
            if self.n_workers != 1:
                raise ValueError(
                    "rng='R' requires n_workers=1 (R's unif_rand stream is"
                    " sequential); got n_workers=" + str(self.n_workers)
                )
            self.engine_c = False

        # Fill in default tuning constants matching R.
        if self.tuning_psi is None:
            if self.psi not in _DEFAULT_TUNING_PSI:
                raise ValueError(f"unknown psi family {self.psi!r}")
            self.tuning_psi = _DEFAULT_TUNING_PSI[self.psi]
        if self.tuning_chi is None:
            if self.psi not in _DEFAULT_TUNING_CHI:
                raise ValueError(f"unknown psi family {self.psi!r}")
            self.tuning_chi = _DEFAULT_TUNING_CHI[self.psi]

        self._warn_unimplemented()

    def _warn_unimplemented(self) -> None:
        """Tell the caller when a control they set is not wired up yet.

        These fields exist so the surface matches ``lmrob.control()``, but
        nothing reads them. Accepting a value and silently ignoring it is
        worse than not offering it: a user who passes ``trace_lev=4`` and
        sees no trace has no way to tell the difference between "no
        output" and "not implemented".
        """
        import warnings

        unimplemented = {
            "trace_lev": (self.trace_lev, 0),
            "eps_outlier": (self.eps_outlier, None),
            "eps_x": (self.eps_x, None),
            "solve_tol": (self.solve_tol, 1e-7),
        }
        set_anyway = [name for name, (got, default) in unimplemented.items() if got != default]
        if set_anyway:
            warnings.warn(
                f"Control({', '.join(sorted(set_anyway))}) is accepted for "
                "lmrob.control() compatibility but not implemented yet; the "
                "value is ignored.",
                UserWarning,
                stacklevel=3,
            )

    # sklearn-compatibility shim. Lets ``GridSearchCV`` reach Control fields
    # via the standard ``estimator__nested__field`` parameter syntax, e.g.
    # ``param_grid={"control__nResample": [200, 500, 1000]}``.
    def get_params(self, deep: bool = True) -> dict[str, object]:
        """Return ``Control``'s public fields as a sklearn-style dict."""
        import dataclasses

        return {f.name: getattr(self, f.name) for f in dataclasses.fields(self)}

    def set_params(self, **params: object) -> Control:
        """Set ``Control``'s public fields in place, sklearn convention."""
        import dataclasses

        valid = {f.name for f in dataclasses.fields(self)}
        for key, value in params.items():
            if key not in valid:
                raise ValueError(f"Invalid parameter {key!r} for Control")
            setattr(self, key, value)
        # Re-run post-init so derived defaults stay consistent. We invoke
        # __post_init__ directly because Control is a regular @dataclass.
        self.__post_init__()
        return self

    @classmethod
    def preset(cls, setting: Setting, **overrides: Any) -> Control:
        """Build a Control for a named preset.

        Equivalent to ``Control(setting=..., **overrides)``; the presets
        follow ``robustbase::lmrob.control(setting=)``:

        - ``"KS2014"``: psi="lqq", method="SMDM", cov=".vcov.w".
        - ``"KS2011"``: same as KS2014 in robustbase 0.99-7.
        - ``"MM"``: psi="bisquare", method="MM", cov=".vcov.avar1".

        This used to build the object and then ``setattr`` the overrides
        on top, which ran after ``__post_init__``. Two bugs came out of
        that: ``preset("KS2014")`` returned psi="bisquare" with
        cov=".vcov.avar1" (so it disagreed with ``Control(setting=
        "KS2014")`` and with R), and an override like ``psi="lqq"``
        left the tuning constants at the previous family's values.
        Going through the constructor keeps one code path.
        """
        if setting not in ("KS2014", "KS2011", "MM"):
            raise ValueError(f"unknown setting: {setting!r}")
        valid = {f.name for f in fields(cls)}
        unknown = sorted(set(overrides) - valid)
        if unknown:
            raise TypeError(f"unknown Control field: {unknown[0]!r}")
        return cls(setting=setting, **overrides)
