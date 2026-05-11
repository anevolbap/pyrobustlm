# Changelog

All notable changes to `pyrobustlm` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **``Control(engine_c=True)``: monolithic Cython lmrob engine.** New
  module ``pyrobustlm._core._lmrob`` runs the whole fast-S + survivor
  refinement + MM iteration in one nogil C block with a single
  workspace allocation. Mirrors the structure of
  ``robustbase/src/lmrob.c::R_lmrob_S`` and ``rwls``. Subset draws use
  numpy's ``bitgen_t`` capsule (Floyd's combination algorithm);
  LAPACK via ``scipy.linalg.cython_lapack`` (``dgesv`` for the
  p-subset solve, ``dgels`` for IRWLS). All five lmrob-supported psi
  families dispatch via a family enum.

  Stackloss seed=0 timing (Default MM, n=21, p=4):

  | family   | default   | engine_c   | speedup |
  |----------|-----------|------------|---------|
  | bisquare | 22.4 ms   | 6.8 ms     | 3.3x    |
  | hampel   | 22.4 ms   | 5.3 ms     | 4.2x    |
  | optimal  | 23.4 ms   | 6.5 ms     | 3.6x    |
  | lqq      | 26.3 ms   | 8.2 ms     | 3.2x    |
  | ggw      | 25.5 ms   | 7.2 ms     | 3.5x    |

  R baseline on the same fit is 3-5 ms, so engine_c lands at 1.3-2.0x
  R end-to-end for the default MM pipeline.

  Off by default because the bitgen draw sequence is not byte-identical
  with ``np.random.Generator.choice``; basin-of-attraction drift can
  occasionally produce a degenerate vcov on tiny-n problems. The
  remaining stages of the monolithic port (D-scale, vcov_avar1,
  vcov_w) stay on the NumPy path for now, so ``setting="KS2014"`` /
  ``"KS2011"`` still pays the Python cost there.

## [0.4.1] - 2026-05-11

### Added

- **``Control(fast_rng=True)``** (opt-in). Routes resampling subset
  draws through a Cython BitGenerator path (Floyd's combination
  algorithm via numpy's ``bitgen_t`` C API) and replaces the SVD-based
  singularity check with ``dgesv`` info inspection. End-to-end speedup
  on serial fits: 2.3x at n=100/p=5/nR=500, 1.3-1.4x at n=500-2000.
  Off by default because the draw sequence is not byte-identical with
  ``np.random.Generator.choice``, so the basin of attraction can shift
  slightly. See the field docstring for the trade-off.
- **Cython survivor-refinement.** ``_refine_to_convergence``, which runs
  ``best_r`` candidates from the resampling pool to convergence, now uses
  the same nogil Cython kernel as the resampling loop
  (``cy_refine_to_convergence``). Removes the Python overhead from the
  IRWLS + m_scale iteration. Roughly 20-30% serial speedup end-to-end on
  problems where refinement was a significant fraction of wall-clock
  (n=100..2000).
- **Read the Docs URL** wired into ``[project.urls]`` and a badge added
  to the README.

### Fixed

- ``lqq`` IRWLS weight formula in the new Cython kernel matches
  ``_psi.pyx::lqq_wgt`` exactly. Stackloss lqq coefficient drift seen
  during the refinement port (rerr 1.7e-3) is gone.
- ``optimal`` IRWLS weight no longer divides by the (cancelling)
  constant ``3.25``; cosmetic match with ``_psi.pyx::optimal_wgt``.

## [0.4.0] - 2026-05-11

### Added

- **Cython fast path extended to hampel, optimal, lqq, ggw.** The
  Cython resampling kernel from v0.3.0 now dispatches on a family enum;
  the bisquare-only ``cy_resample_iter_bisquare`` wrapper is kept for
  compatibility. ggw brings over the six polynomial chi coefficient
  tables from ``_psi.pyx`` plus the ``_GGW_ABC`` (a, b, c) lookup, so
  the resampling iteration stays nogil. End-to-end runtime on
  n=500/p=10/nR=500 is now uniform across all five families
  (~110-180 ms vs. 200-280 ms on the v0.3.0 NumPy path). Stackloss
  element-wise parity vs R is preserved for all families (rtol=1e-4
  on coef and scale; rtol=5e-3 for ggw, limited by R's polynomial
  chi approximation).
- **Read the Docs config (``.readthedocs.yaml``).** Builds the Sphinx
  site on every push and PR. Uses ubuntu-22.04 + Python 3.11 with
  libopenblas/liblapack apt packages so meson-python can link against
  system BLAS/LAPACK. ``fail_on_warning: true`` keeps doc breakage
  visible in CI.

## [0.3.0] - 2026-05-10

### Added

- **Sphinx documentation scaffolding.** ``docs/`` now builds with
  ``sphinx-build -b html -W docs docs/_build/html``. Pages: index,
  quickstart, API reference (autodoc on
  ``lmrob`` / ``Control`` / ``LmRobResults`` / ``SummaryLmRob`` /
  ``anova``), R-to-Python porting cheatsheet, and the existing
  numerical-notes log. Theme: ``sphinx-rtd-theme`` with
  ``myst-parser``, ``napoleon`` (numpy-style), ``intersphinx`` to
  numpy/scipy/pandas, and ``sphinx-autodoc-typehints``.
- **Cython resampling kernel for the bisquare default.** New module
  ``pyrobustlm._core._fast_s`` runs the per-iteration body
  (initial p-subset solve, k-step refinement, m-scale, IRWLS) in
  ``nogil`` C with LAPACK calls via ``scipy.linalg.cython_lapack``
  (``dgesv`` for the subset solve, ``dgels`` for IRWLS). This is the
  hot path for ``Control(psi="bisquare")``, which is the default.
  Subset draws still happen in Python (RNG path is Python-side); the
  iteration body is fully nogil so ``n_workers > 1`` now scales on
  small problems too.
  Measured serial speedup vs the NumPy implementation: 2.4x at
  n=100/p=5/nR=500, 2.0x at n=500/p=10. With 8 worker threads:
  3-5x at n>=500 problems and ~5x faster than R wall-clock at
  n=2000/p=20.
- **Parallel fast-S resampling.** ``Control(n_workers=...)`` distributes
  the resampling loop across a ``ThreadPoolExecutor``. ``n_workers=0`` is
  auto (only kicks in for problems where ``n*p^2 >= 1e6`` and
  ``nResample >= 250``). Per-thread RNG via ``SeedSequence.spawn`` keeps
  results deterministic for a given ``(seed, n_workers, nResample)``.
  Measured ~1.7x speedup at n=5000, p=30, nResample=2000 on a 16-core
  OpenBLAS box. Small-n problems are GIL-bound and not affected; see
  ``docs/numerical-notes.md`` entry 5.
- **``anova(test="Deviance")``** ports the Deviance variant from
  ``robustbase::anovaLmrobPair``. Refits the reduced model via
  ``mm_iterate`` at the full's scale, then computes
  ``T = 2 * tauStar * (sum_rho_reduced - sum_rho_full)`` with
  ``tauStar = mean(psi'(r/s0)) / mean(psi(r/s0)^2)``. Matches R element-wise
  on stackloss (chi-sq and p-value within rtol=2e-3).
- ``LmRobResults`` now stashes the design matrix (``design_x_``) and
  response (``design_y_``) for downstream operations (currently used by
  the Deviance anova refit).

## [0.2.0] - 2026-05-10

### Added

- **``summary()``** mirrors ``robustbase:::summary.lmrob``. Returns a
  ``SummaryLmRob`` with the coefficient table (``Estimate``, ``Std. Error``,
  ``t value``, ``Pr(>|t|)`` from the t-distribution with ``df_residual``
  degrees of freedom), robust multiple R-squared and adjusted R-squared
  (using the family-specific ``E[wgt(r)] / E[r psi(r)]`` correction
  factor), and an R-style ``__str__``.
- **``anova(full, *reduced, test="Wald")``** runs the robust Wald test
  for nested ``LmRobResults``. Pair-wise output matches R element-wise
  on stackloss (chi-sq and p-value within rtol=2e-3). Chained mode does
  sequential ``Model k`` vs ``Model k-1`` comparisons (matches
  ``anova.lm`` semantics).
- **Full ``vcov_w`` port.** All five ``cov.corrfact`` branches
  (``asympt``, ``empirical``, ``tau``, ``hybrid``, ``tauold``) plus the
  Huber finite-sample correction (``cov.hubercorr``), the five
  ``cov.dfcorr`` branches (``mean``, ``mn.vc``, ``none``, ``varc``,
  ``mn.df``), and the three ``cov.resid`` modes (``final``, ``initial``,
  ``trick``). Setting-driven defaults match R's
  ``lmrob.control(setting="KS2014" / "KS2011")`` exactly.
  ``setting="KS2014"`` and ``setting="KS2011"`` covariance now matches R
  **element-wise to rtol=1e-3** on stackloss (was 0.10 / 0.18 before).
- **Full M-S port.** Replaces the simplified L1/S alternation with a
  direct port of robustbase's ``R_lmrob_M_S`` (4 phases: orthogonalize
  via L1, subsample many candidates in orth space, transform back via
  ``b1 = ot1 + b1_orth - oT2 @ b2_orth``, descent alternating L1 + WLS).
  On the ``education`` reference, M-S now matches R within ~2% on
  coefficients and ~0.6% on scale at ``nResample=2000`` (was 30%+ off).
- **Design-adaptive D-scale (Koller & Stahel 2014).** New module
  ``pyrobustlm.d_scale`` ports ``robustbase::lmrob..D..fit`` and the
  ``R_find_D_scale`` C iteration. ``setting="KS2014"`` and
  ``setting="KS2011"`` now run the full SMDM pipeline (S-init, MM,
  D-scale, MM with new scale). Matches R's C kernel to ``rtol=1e-4``
  on identical inputs (verified in ``tests/unit/test_d_scale.py``).
- **Setting-driven defaults.** ``Control(setting=None)`` (or omitted)
  follows R's plain ``lmrob.control()`` (psi=bisquare, method=MM,
  cov=.vcov.avar1). ``setting="KS2014"`` and ``setting="KS2011"``
  default to psi=lqq, method=SMDM, cov=.vcov.w (matching R 0.99-7).

### Changed

- ``Control`` now treats ``psi``, ``method``, and ``cov`` as
  ``None``-defaulted "auto from setting" fields. Pass them explicitly
  to override.

- **Cython kernels for all six psi families.** `huber`, `hampel`, `optimal`,
  `lqq`, and `ggw` now have inlined Cython implementations of
  `psi`/`wgt`/`rho`/`psi_prime` and a fully-inlined `m_scale_<family>` that
  keeps the iteration in C. End-to-end runtime is now uniform across
  families (~26x R for small problems, ~1.6x R for n=2000).
- **`vcov_avar1` parity for all five lmrob-supported families.** Fixed
  per-family `chi_prime_factor` constants verified against R's
  `Mchi(deriv=1) / Mpsi`. Hampel uses `1/nc`, optimal uses
  `1/(3.25 c²)`, lqq uses `6(s-1)/denom`, ggw uses tabulated factors per
  case (1..6).
- **Per-family parity tests** in `tests/validation/test_vs_r_psi_families.py`:
  bisquare/optimal/hampel/lqq match R within `rtol=1e-3`; ggw within `5e-3`
  (limited by R's polynomial chi approximation).
- **Benchmark harness.** `scripts/benchmark.R` and `scripts/benchmark.py`
  run identical fits in R and Python; `scripts/build_bench_report.py`
  merges them into `docs/bench-report.md` with side-by-side accuracy and
  runtime tables. 20 cases: 10 classical datasets, 5 per-psi-family,
  5 synthetic timing grids.

### Fixed

- Default `Control.tuning_psi` for `lqq` (was `0.9826779`, now R's
  `0.9822707`) and for `ggw` (was case-1 parameters, now case 4 = b=1.5,
  95% efficiency, matching R's default).

## [0.1.0] - 2026-05-09

First public release. End-to-end MM regression that matches R's
`robustbase::lmrob` element-wise on the classical datasets within
`rtol=1e-3` for both coefficients and covariance.

### Added

- **Build system.** `pyproject.toml` (meson-python backend), GPL-3 license,
  `NOTICE` crediting upstream `robustbase` authors. CI matrix
  (Linux/macOS/Windows × Python 3.10-3.13) and `cibuildwheel` workflow.
- **R reference harness.** `scripts/generate_r_reference.R` produces 35
  JSON references covering 10 classical datasets and 18 synthetic
  configurations. Pytest fixtures with rpy2 bridge for live R comparisons
  during development.
- **Psi/chi/wgt.** All six robustbase families (`bisquare`, `huber`,
  `hampel`, `optimal`, `lqq`, `ggw`). Match R's `Mpsi`/`Mchi`/`Mwgt` to
  `rtol=1e-12` (1e-7 for ggw rho's polynomial approximation).
- **M-scale.** `pyrobustlm.scale.m_scale` matches R's `find_scale` to
  `rtol=1e-9`.
- **Fast-S** resampling estimator and **MM** IRWLS iteration.
- **M-S** initial estimator for designs with categorical predictors
  (Maronna-Yohai 2000 alternating L1/S). Auto-detection via
  `init="auto"`.
- **Inference.** `vcov_avar1` (sandwich) ported from
  `robustbase/R/lmrob.MM.R:510-577`, matching R element-wise to `rtol=1e-3`.
  `vcov_w` implements the asymptotic-correction-factor branch.
- **Public API.** `pyrobustlm.lmrob(formula, data, ...)`, `LmRob`
  scikit-learn-style estimator, `Control` with KS2014/KS2011/MM presets.
- **Predict on new data.** `LmRobResults.predict(DataFrame)` re-applies
  the formula's design transformation (factor encoding, `I(x**2)`, etc.)
  via formulaic's `ModelSpec`. Also accepts a raw NumPy design matrix.
- **Diagnostics.** Four-panel `diagnostics.plot`, robust `hatvalues`,
  robust `cooks_distance`.
- **Cython acceleration.** Bisquare psi/wgt/rho/psi_prime and M-scale
  inner loop. End-to-end `lmrob()` is 3-5x faster than the pure-NumPy
  baseline; on n=2000, p=20 it lands at 1.1x R's wall time.
- **Validation.** 10 classical datasets pass at `rtol=1e-3` (one at
  `rtol=5e-3`). Three Hypothesis property tests for affine/scale/regression
  equivariance.

### Known limitations

- `init="M-S"` is a simplified port; not bit-identical with the C
  `m_s_descent` (which does sophisticated multi-restart subsampling).
- `vcov_w` Huber finite-sample correction and `empirical`/`hybrid`/
  `tau`/`tauold` branches are deferred.
- Cython kernels currently exist only for the bisquare family; other
  families use the slower NumPy path.
- RNG is not bit-identical with R (we use PCG64; R uses Mersenne Twister).
- `setting="KS2011"` runs the SM path but not the additional D-step
  scale refinement that R's SMDM method does.

See [`docs/numerical-notes.md`](docs/numerical-notes.md) for the full list
of documented divergences from R.

[Unreleased]: https://github.com/anevolbap/pyrobustlm/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/anevolbap/pyrobustlm/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/anevolbap/pyrobustlm/releases/tag/v0.1.0
