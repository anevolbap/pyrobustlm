# Changelog

All notable changes to `pyrobustlm` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.5.4] - 2026-05-11

### Added

- **vcov_avar1 folded into ``cy_lmrob_fit``.** The fast-S kernel now
  optionally takes ``cov_out`` and computes the sandwich covariance
  inline using the same workspace it just used for fast-S + MM. The
  per-element ``psi_prime`` / ``psi`` / ``chi`` / ``chi_prime``
  kernels and the ``X' diag(w_pp) X`` build + ``dgesv`` inversion are
  factored into a shared ``cdef`` helper used by both
  ``cy_lmrob_vcov_avar1`` and the new ``cov_out`` parameter.
  Removes one Python/C boundary cross and one extra
  ``np.zeros((p, p))`` per fit. Small absolute gain (~25-30 µs) but
  the architecture is cleaner: a single Cython call now produces
  beta, scale, residuals, rweights, init residuals (beta_init_out),
  and (optionally) cov.

### Build / dev experience

- Pin ``ruff>=0.15,<0.16`` in ``[dev]``. Each ruff minor can subtly
  change format rules; pinning keeps local dev and CI agreeing on
  ``ruff format --check``.

### Perf misc

- ``model_matrix`` simple-formula fast path: batch dtype lookup via
  ``data.dtypes`` (one call) instead of ``data[c].dtype`` per column
  (N calls), and use ``Series.values`` (~25% faster than
  ``.to_numpy()``).
- Tighter engine_c block in ``lmrob.py``: drop unneeded
  ``np.atleast_1d`` round-trips, skip ``np.ascontiguousarray`` when
  the arrays are already C-contig float64, drop unnecessary copies.

## [0.5.3] - 2026-05-11

### Fixed

- **CI green across Linux/macOS/Windows × Python 3.10-3.13.** Three CI-only
  failures fixed:
  - ``_fast_s.pyx`` declared ``ndarray[long, ndim=1]`` for the subset-
    indices buffer, which broke on Windows LLP64 (``long`` is 32-bit on
    Windows but our ``np.int64`` array is 64-bit). Switched to
    ``cnp.int64_t``.
  - ``test_d_iteration_matches_r_kernel`` subprocesses ``Rscript`` and
    needs ``robustbase`` installed. Skips now when either is missing
    (CI only installs R on Linux).
  - ``test_engine_c_speedup_at_small_n`` is a wall-clock assertion that
    was flaky on shared CI runners. Skips when ``CI`` env var is set.

### Build

- ``scipy`` added to ``[build-system].requires`` and its parent
  directory wired into ``cython_args`` as ``-I``. The new
  ``_fast_s.pyx`` and ``_lmrob.pyx`` modules ``cimport`` from
  ``scipy.linalg.cython_lapack``, which needs scipy at build time when
  using PEP-517 isolation (as CI does).

## [0.5.2] - 2026-05-11

### Added

- **Formula fast path.** ``model_matrix`` now detects the common
  ``y ~ x1 + x2 + ...`` pattern (numeric columns only, no factors or
  transforms) and parses it by hand, skipping formulaic's ~1.5 ms
  per-fit overhead. Falls back to formulaic for any complex formula
  (factors, interactions, ``I(x**2)``, etc.).
  ``LmRobResults.predict(DataFrame)`` learns to rebuild the design
  matrix the same way when ``rhs_spec_`` is ``None``.

  Combined with the v0.5.1 Cython vcov port, end-to-end stackloss MM
  drops from 5.6 ms to 4.4 ms with ``engine_c=True``; phosphor MM is
  3.8 ms (close to R's ~3 ms on the same fit). KS settings drop to
  6-10 ms.

  | dataset / setting | default | engine_c | speedup |
  |---|---|---|---|
  | stackloss MM     | 22.4 ms |  4.4 ms |  5.1x |
  | stackloss KS2014 | 62.8 ms |  7.3 ms |  8.7x |
  | phosphor MM      | 20.9 ms |  3.8 ms |  5.5x |
  | phosphor KS2014  | 59.1 ms |  6.3 ms |  9.5x |
  | phosphor KS2011  | 58.0 ms |  6.3 ms |  9.3x |
  | salinity KS2014  | 62.0 ms | 10.0 ms |  6.2x |

## [0.5.1] - 2026-05-11

### Added

- **``vcov_avar1`` in the monolithic Cython kernel.** New function
  ``cy_lmrob_vcov_avar1`` ports the sandwich-covariance computation
  end-to-end: per-family ``psi_prime`` / ``psi`` / ``chi`` / ``chi_prime``
  inlined nogil, ``X' diag(w_pp) X`` build + ``dgesv`` inversion, then
  the u1/u2/u3/u4 assembly in C. Matches ``inference.vcov_avar1``
  element-wise (rerr ~1e-16). Posdefify (eigendecomposition) stays in
  NumPy because p is small. ``cy_lmrob_fit`` also returns the post-S
  beta via the new optional ``beta_init_out`` parameter so the Python
  side can compute the right initial residuals.

  Side-by-side timing (single-threaded BLAS, 7-rep median):

  | dataset / setting | default | engine_c | speedup |
  |---|---|---|---|
  | stackloss MM     | 24.6 ms |  5.6 ms |  4.4x |
  | stackloss KS2014 | 66.0 ms | 10.0 ms |  6.6x |
  | phosphor MM      | 23.8 ms |  5.4 ms |  4.4x |
  | phosphor KS2014  | 63.3 ms |  7.4 ms |  8.5x |
  | phosphor KS2011  | 60.4 ms |  7.6 ms |  7.9x |
  | salinity KS2014  | 71.7 ms | 12.8 ms |  5.6x |

  R wall-clock on these fits is 3-7 ms, so ``engine_c=True`` lands at
  1.4-1.9x R across all four datasets and three settings.

## [0.5.0] - 2026-05-11

### Added

- **``Control(engine_c=True)``: monolithic Cython lmrob engine.** New
  module ``pyrobustlm._core._lmrob`` runs the whole fast-S + survivor
  refinement + MM iteration + D-scale (KS2014/KS2011) in one nogil C
  block with a single workspace allocation. Mirrors the structure of
  ``robustbase/src/lmrob.c::R_lmrob_S``, ``rwls``, and
  ``R_find_D_scale``. Subset draws use numpy's ``bitgen_t`` capsule
  (Floyd's combination algorithm); LAPACK via
  ``scipy.linalg.cython_lapack`` (``dgesv`` for the p-subset solve and
  the (X'WX) inverse, ``dgels`` for IRWLS). All five lmrob-supported
  psi families dispatch via a family enum. D-step uses tabulated
  ``kappa`` and ``(tfact, tcorr)`` per family at default tuning;
  non-default tuning falls back to the NumPy path.

  Side-by-side timing on classical datasets
  (single-threaded BLAS, 7-rep median):

  | dataset / setting | default | engine_c | speedup |
  |---|---|---|---|
  | stackloss MM       | 27.7 ms |  6.9 ms |  4.0x |
  | stackloss KS2014   | 69.2 ms |  9.5 ms |  7.3x |
  | stackloss KS2011   | 70.4 ms |  9.3 ms |  7.6x |
  | delivery MM        | 24.7 ms |  7.1 ms |  3.5x |
  | delivery KS2014    | 68.8 ms | 10.9 ms |  6.3x |
  | phosphor MM        | 23.5 ms |  5.8 ms |  4.0x |
  | phosphor KS2014    | 69.5 ms |  8.1 ms |  8.6x |
  | phosphor KS2011    | 74.2 ms |  8.1 ms |  9.2x |
  | salinity MM        | 26.8 ms |  8.8 ms |  3.0x |
  | salinity KS2014    | 76.0 ms | 13.2 ms |  5.7x |

  R baseline on these fits is 3-7 ms, so engine_c lands at 1.3-2x R
  wall-clock for both the default MM and KS pipelines.

  Off by default because the bitgen draw sequence is not byte-identical
  with ``np.random.Generator.choice``; basin-of-attraction drift can
  occasionally produce a degenerate vcov on tiny-n problems.
  ``vcov_avar1`` and ``vcov_w`` stay on the NumPy path; they're ~1 ms
  even at n=5000, so the remaining win is small.

- **``scripts/bench_engine_c.py``** for side-by-side timing of the two
  engines across the classical datasets and KS settings.

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

[Unreleased]: https://github.com/anevolbap/pyrobustlm/compare/v0.5.4...HEAD
[0.5.4]: https://github.com/anevolbap/pyrobustlm/compare/v0.5.3...v0.5.4
[0.5.3]: https://github.com/anevolbap/pyrobustlm/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/anevolbap/pyrobustlm/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/anevolbap/pyrobustlm/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/anevolbap/pyrobustlm/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/anevolbap/pyrobustlm/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/anevolbap/pyrobustlm/releases/tag/v0.1.0
