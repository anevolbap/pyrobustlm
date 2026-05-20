# Changelog

All notable changes to `pylmrob` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- **``Control(rng="R")`` survivor refinement now matches robustbase's
  ``refine_fast_s`` more closely.** Two algorithmic fixes:
  - The per-iter scale update is now ONE Newton step
    (``s = s * sqrt(sum_rho(r/s) / ((n-p) * b0))``) instead of full
    M-scale convergence. This matches the fast-S algorithm in
    Salibian-Barrera & Yohai (2006) and robustbase's C code.
  - The L2-norm convergence test uses the OLD beta's norm, matching
    R's ``refine_fast_s`` line-for-line.
  - The MM step routes through the Cython ``cy_lmrob_mm`` kernel
    (LAPACK ``dgels`` QR-based solver) instead of the NumPy
    ``np.linalg.lstsq`` (gelsd, SVD-based) fallback, even when
    ``engine_c=False``. Closer match to R's ``rwls()``.

  Net effect on stackloss vs R's ``lmrob`` with ``set.seed(seed)``:
  intercept gap improved from ~2.3e-5 to ~1.7e-5; scale gap from
  ~8.4e-6 to ~6.1e-6. The residual gap is mostly in the resample
  loop's "associated scale" computation and survivor selection, both
  of which still differ from robustbase. A full closing to bit-
  identical fits would require porting ``refine_fast_s`` and
  ``find_scale`` end-to-end.

- **Cython kernel for the R RNG** (``pylmrob._core._r_rng``). The
  inner loops of ``unif_rand_n``, ``r_sample_noreplace``, and
  ``r_subsample_nonsingular`` now run nogil in C. Pure-Python paths
  remain as fallback when the Cython module isn't built. Microbench
  on ``r_sample_noreplace(n=200, n=200)`` × 500: 133 ms → 1.0 ms
  (~130x); ``r_subsample_nonsingular(n=200, p=5)`` × 500: 162 ms →
  1.5 ms (~108x); end-to-end stackloss fit with ``rng="R"``: 33 ms
  → 9 ms (~3.6x). Result is byte-identical to the pure-Python path.

### Added

- **``pylmrob.r_qnorm(p)``**: standard normal quantile, byte-identical to
  R's ``qnorm()``. Direct port of Wichura's AS 241 from R's
  ``src/nmath/qnorm.c``.
- **``pylmrob.r_norm_rand(rng)``**: one standard normal draw, byte-identical
  to R's ``rnorm()``. Ports the ``Inversion`` branch of R's ``norm_rand``
  in ``src/nmath/snorm.c`` (two ``unif_rand`` draws combined into a
  >27-bit uniform, then ``r_qnorm``). Verified against R's actual output
  for seeds {0, 1, 42, 12345, ``2**31 - 1``} and across the regenerate
  boundary.

## [0.5.16] - 2026-05-19

### Added

- **``pylmrob.r_set_seed(seed)`` and ``pylmrob.RState``**: a pure-Python
  MT19937 RNG that replicates R's ``set.seed`` / ``unif_rand`` path
  from ``src/main/RNG.c``. The uniform stream is byte-identical to R's
  ``runif()`` after ``set.seed()`` for the same integer seed, including
  across the 624-word regenerate boundary. Verified against R 4.2 in
  ``tests/validation/test_r_rng_vs_R.py``. First step toward
  bit-identical lmrob fits; the resample kernel is not wired to it yet.
- **``pylmrob.r_sample_noreplace(rng, n, k)``**: port of
  ``sample_noreplace`` from robustbase's ``lmrob.c``. Knuth-style
  swap-and-replace; draws ``k`` distinct indices from ``0..n-1`` using
  one ``unif_rand`` per output. Verified bit-identical to robustbase's
  C implementation via ``.C(R_subsample, ..., sample = TRUE)`` for
  ``n`` in {10, 25, 100} and seeds {1, 42, 12345}.
- **``Control(rng="R")``**: route the fast-S resample loop through
  ``r_set_seed`` + (``r_sample_noreplace`` or
  ``r_subsample_nonsingular``) so subset draws are byte-identical to
  robustbase. Implies ``n_workers=1`` and ``engine_c=False`` (Control
  forces these). Supports both ``subsampling="simple"`` and
  ``subsampling="nonsingular"`` (R's default). Final fits on stackloss
  agree with R's ``lmrob`` after ``set.seed(seed)`` to ``rtol=1e-4``
  across seeds {1, 42, 12345}; the residual drift comes from
  refinement-step LAPACK ordering, not the draws.
- **``pylmrob.r_subsample_nonsingular(rng, X, p)``**: pure-Python port
  of the ``ss=1`` LU-pivot block of robustbase's ``subsample()``.
  Walks a permutation column by column with partial column pivot,
  skipping rows whose pivot falls below ``tolInverse``. Verified
  bit-identical to robustbase's C ``subsample()`` (via
  ``.C(R_subsample, ..., sample=TRUE, ss=1)``) for varying ``p``,
  collinear-row stress tests, and seeds {1, 7, 42, 12345}.

### Docs

- README decluttered: feature list and Performance section moved to
  ``docs/numerical-notes.md`` and ``docs/bench-report.md``; CI badge
  fixed; Read the Docs badge dropped (not live); inline links to the
  core papers added.
- ``CONTRIBUTING.md`` added.
- FAQ entry on bit-identical R fits expanded to point at
  ``r_set_seed``.

## [0.5.15] - 2026-05-18

### Added

- **``pylmrob.bootstrap(fit, n_boot=1000, level=0.95, ...)``** and
  the method-style spelling **``fit.bootstrap(...)``** on
  ``LmRobResults``. Case-resampling bootstrap returning a
  ``BootstrapResult`` with the full coefficient distribution,
  percentile and basic CIs, standard errors, and bias estimate.
  Useful when the asymptotic Wald CI is unreliable (small n,
  near-singular X, heavy contamination); doesn't replace
  ``fit.confint()``. Deterministic given ``(seed, n_workers)``;
  ThreadPool-parallel via the ``n_workers`` argument.

### Docs

- New **theory page** (``docs/theory.md``) with a short conceptual
  tour: M-estimator, M-scale equation, the two-stage MM-estimator,
  the Koller-Stahel D-scale refinement, guidance on picking a psi
  family.

- New **FAQ / troubleshooting page** (``docs/faq.md``) covering
  convergence failures, singular vcov, Wald vs bootstrap CIs,
  bit-identical R reproducibility, ``nResample`` tuning, performance
  levers, ``predict()`` pitfalls.

- **Bootstrap sections in quickstart and stackloss tour** showing
  how to call ``fit.bootstrap()`` and compare percentile vs Wald
  intervals.

## [0.5.14] - 2026-05-18

### Added

- **Prediction intervals.** ``LmRobResults.predict`` gains
  ``interval='confidence'|'prediction'`` and ``level=`` keyword args.
  Returns ``(n, 3)`` with ``(fit, lwr, upr)`` columns using the
  t-distribution at ``df_residual``. Matches R's ``predict.lmrob``.

- **statsmodels-style attribute aliases on ``LmRobResults``**:
  ``params``, ``bse``, ``tvalues``, ``pvalues``, ``conf_int(alpha)``.
  Property forwarders over the existing fields; pylmrob fits drop into
  ``statsmodels.regression``-shaped code without an adapter.

- **sklearn-style ``LmRob.score`` / ``get_params`` / ``set_params``.**
  Makes the estimator compatible with ``cross_val_score``,
  ``GridSearchCV``, and other sklearn pipeline tooling. ``score``
  returns OLS R² on the test set (sklearn convention).

- **``LmRobResults.diagnostics()``.** Returns a ``DiagnosticsTable``
  with per-observation leverage, robust Cook's distance, standardized
  residuals, robust weights, and an outlier flag. Wires the existing
  ``hatvalues`` / ``cooks_distance`` helpers in ``pylmrob.diagnostics``.

- **``LmRobResults.anova(*others, test=)``** as a method-style spelling
  of the existing free function. Matches R's idiom.

- **R-style ``__repr__``.** Prints a call line and an aligned
  coefficient table; ``print(fit.summary())`` still has the verbose
  per-coef stats and R-squared.

- **Stackloss walkthrough.** ``docs/examples/stackloss_tour.md`` is a
  complete robust-regression tutorial: OLS vs lmrob, outlier
  identification, CI/PI bands, nested-model test, statsmodels-style
  access, sklearn cross-validation. Hooked into the docs toctree.

## [0.5.13] - 2026-05-18

### Changed

- **Renamed the package from ``pyrobustlm`` to ``pylmrob``** to mirror
  R's ``robustbase::lmrob`` directly (the function this package ports).
  ``pyrobustlm`` was never published to PyPI, so there is no user
  migration path; if you were using the editable install or TestPyPI
  build, switch ``from pyrobustlm import ...`` to ``from pylmrob
  import ...``. The GitHub repository name is unchanged.

## [0.5.12] - 2026-05-16

### Added

- **Per-case ``weights`` argument to ``lmrob()``.** Mirrors R's lmrob
  implementation: a sqrt(w)-transform at the top level
  (robustbase/R/lmrob.R:96-98) routes the transformed design
  ``(sqrt(w)*X, sqrt(w)*y)`` through the unweighted fit. Zero-weight
  rows are dropped (R's behaviour). ``residuals_`` / ``fitted_`` are
  reported on the original scale, while ``rweights_`` / ``scale_`` /
  ``cov_`` come from the transformed fit. Matches R element-wise on
  stackloss with non-trivial weights (coef rerr 1.1e-6, scale rerr
  4.1e-6). Non-trivial weights currently force the NumPy path
  (``engine_c=False``); full Cython-engine support is a follow-up.

## [0.5.11] - 2026-05-16

### Changed

- **``Control.engine_c`` defaults to True.** The monolithic Cython
  kernel is consistently faster (median 0.93x R wall-clock across the
  34-case bench corpus, down from default-path 2.74x R) and its cov
  now agrees with the numpy path to ~1e-7. ``lmrob()`` catches the
  rare ``FloatingPointError`` from a singular ``X' W X`` (basin drift
  on small classical datasets like stackloss/hbk) and retries with
  ``engine_c=False`` so the new default never raises on cases the old
  default handled. Pass ``Control(engine_c=False)`` to force the
  legacy numpy path.

### Fixed

- **``engine_c=True`` cov for ggw on small n.** The inline
  ``vcov_avar1`` in the Cython kernel used placeholder values for
  the ggw ``chi'(x)/psi(x)`` factor (``1/1.6047`` for cases 1 and 4,
  ``6/case^2`` fallback for the rest), which blew up the cov diagonal
  by 20-100x on small-n datasets. Ported the correct case-dependent
  values from ``pylmrob._psifuns.chi_prime_over_psi``. psi_ggw on
  n=21 stackloss with engine_c=True went from cov diag
  ``[534, 0.44, 11.6, 0.013]`` to ``[29, 0.021, 0.17, 0.0046]``, in
  agreement with the numpy path.

- **``lmrob()`` auto-falls-back on engine_c singular vcov.** The
  Cython subset-draw is not byte-identical to ``np.random.choice``,
  so on a few small classical datasets (stackloss, hbk) it lands in
  a basin where ``X' W X`` is singular. ``lmrob()`` now catches the
  resulting ``FloatingPointError`` once and retries with
  ``engine_c=False``.

## [0.5.10] - 2026-05-15

### Perf

- **Hoist per-call constants out of psi/chi inner loops.** The Cython
  ``_chi_sum`` and ``_wgt_zinv`` kernels recomputed denominators,
  reciprocals, and constant subexpressions on every element of the
  per-residual loop (lqq especially: ``(3*s_l-3)/denom``, ``s_l/b_l``,
  ``s5*s5/(3*s6)``). Hoisted them above the loop and replaced the
  per-element ``r[i]/s`` with ``inv_s = 1/s`` plus a multiply.

  Bench at ``synth_n500_p10`` with ``Control(engine_c=True)``
  (single-thread OpenBLAS, A/B on the same machine):

  | psi | before | after | delta |
  |---|---|---|---|
  | lqq | 200 ms | 130 ms | -35% |
  | optimal | 152 ms | 97 ms | -36% |
  | bisquare | 107 ms | 77 ms | -29% |
  | hampel | 130 ms | 125 ms | -5% |

  Median runtime ratio ``py engine_c / R`` across the full bench corpus
  drops from 1.02x to 0.80x: pylmrob now beats R on more than half
  of cases. The worst remaining engine_c case
  (``synth_ggw_n500_p10``) is 1.91x R.

## [0.5.9] - 2026-05-15

### Perf

- **Cython MM in the default-path branch when ``engine_c=True``.** The
  engine_c block auto-falls-back at large n, but it used to drop into
  the NumPy ``mm_iterate`` after the fallback. The IRLS loop now uses
  ``cy_lmrob_mm`` whenever ``Control(engine_c=True)`` is set and the
  ``psi`` family has a Cython kernel. ``Control(engine_c=True)`` and
  ``Control(n_workers=0)`` now finish in the same ballpark on
  large-n fits.

  Bench at n=5000, p=30, nResample=500 (single-thread OpenBLAS):

  | config | wall-clock |
  |---|---|
  | ``Control()`` (default, serial) | 1882 ms |
  | ``Control(n_workers=0)`` | 730 ms |
  | ``Control(engine_c=True)`` | 719 ms |

## [0.5.8] - 2026-05-14

### Perf

- **BLAS dgemm in the Cython vcov body.** The two big matrix products
  inside ``_compute_vcov_avar1_body`` (the ``X^T diag(w) X``
  computation and the inner ``u_mat = X^T diag(psi^2) X``) were hand-
  coded triple loops. Switched to ``dgemm`` via
  ``scipy.linalg.cython_blas``: build a column-major
  ``X_w[i + j*n] = X[i,j]*w[i]`` and let dgemm do
  ``X^T @ X_w``. The sqrt(w) trick is avoided because ``psi'`` can be
  negative for redescending psi.

  Per-call vcov_avar1 wall-clock (single-thread OpenBLAS):

  | n / p | Python | Cython before | Cython after |
  |---|---|---|---|
  | 100 / 5  | 0.18 ms | 0.16 ms | 0.02 ms |
  | 500 / 10 | 0.18 ms | 0.30 ms | 0.05 ms |
  | 5000 / 30| 1.29 ms | 22 ms   | 1.20 ms |

- **Cython vcov fast path is now used at all sizes** (the
  ``_engine_c_too_big`` gate is removed from the standalone vcov
  branch). The engine_c block itself still falls back at large n
  because its embedded fast-S kernel does not parallelise.

## [0.5.7] - 2026-05-14

### Changed

- **``Control(engine_c=True)`` auto-falls-back to threaded default at
  large n.** The monolithic Cython kernel is a single non-parallel
  call; once ``n*p^2 >= 100,000`` it loses to the threaded default
  path. ``lmrob()`` now detects this regime, skips the engine_c block,
  also skips the Cython ``vcov_avar1`` fast path (hand-coded matrix
  products are ~10x slower than NumPy/BLAS at p>=20), clears
  ``cfg.engine_c`` when calling ``fast_s`` (so the embedded
  ``cy_lmrob_fast_s`` branch does not fire), and enables ``n_workers=0``
  auto-threading on the fallback. A single ``Control(engine_c=True)``
  setting now gives the fastest fit on small-n and large-n alike.

  Bench at n=5000, p=30, nResample=500 (single-thread OpenBLAS):

  | config | wall-clock |
  |---|---|
  | ``Control()`` (default, serial) | ~2230 ms |
  | ``Control(n_workers=0)`` | ~860 ms |
  | ``Control(engine_c=True)`` | ~870 ms (auto fallback) |

- **Auto worker count bumped** from ``n_iter // 250`` to ``n_iter // 64``.
  At the default ``nResample=500`` this lifts the auto picker from
  2 workers to 7, which is past the diminishing-returns knee on the
  representative configurations we measure.

### Added

- ``CITATION.cff`` so GitHub shows a "Cite this repository" button.
  Lists the four key papers: Yohai (1987), Salibian-Barrera & Yohai
  (2006), Koller & Stahel (2011), and the robustbase package itself.

### Docs

- New ``engine_c`` section in ``docs/quickstart.md``.
- ``docs/engine_c.md`` rewritten to document the auto-fallback
  behaviour (no manual path-picking for users).
- ``docs/numerical-notes.md`` entry 2 (Performance vs R) refreshed
  with engine_c numbers.

## [0.5.6] - 2026-05-14

### Added

- **Engine_c parity test corpus.** New
  ``tests/integration/test_engine_c_parity.py`` fits four classical
  datasets twice (default + ``engine_c=True``) and asserts the
  coefficients, scale, and vcov diagonals agree to rtol=1e-8 (on
  coef/scale) and rtol=1e-6 (on cov diag).

### Fixed

- ``_compute_vcov_avar1_body`` had a stray ``except 3`` clause that
  Cython interpreted as "return value 3 == Python exception
  pending"; status=3 was the legitimate non-exception status for a
  LAPACK error. Dropping the clause lets the status code flow
  through cleanly; the previously-confusing "returned NULL without
  setting an exception" is now a clean ``FloatingPointError``.

### Perf misc

- Bundle all output buffers into one ``np.empty`` allocation
  (``beta_out``, ``residuals_out``, ``rweights_out``,
  ``beta_init_out``, ``cov_out`` are now slices of one big array).
  Saves about 40 us per fit.
- ``LmRob.fit`` builds its working DataFrame in one shot from a
  dict-of-arrays instead of ``pd.DataFrame(X, columns=...)`` +
  ``df["y"] = y``. Saves about 150 us per fit.

### Docs

- README refreshed: status v0.5.6, new Performance section, "What
  works" reflects what actually ships, "What does not work yet"
  trimmed to actual gaps.
- ``docs/numerical-notes.md`` entry 2 (Performance vs R) refreshed
  with engine_c numbers.
- ``docs/bench-report.md`` refreshed against v0.5.6 default path.

## [0.5.5] - 2026-05-11

### Added

- **Posdefify (vcov projection to PSD) in Cython.** The `dsyev`
  eigendecomposition + clipping + reconstruction lives inside the
  Cython kernel now; the Python-side `numpy.linalg.eigh` round-trip is
  gone. Marginal but tidier.
- **Docs: new `engine_c` page** on the Sphinx site, explaining the
  trade-off and listing the bench numbers.
- **Anova factor-design test.** `tests/validation/test_summary_anova.py`
  now covers a multi-column factor drop (`RegionF` from the
  ``education`` dataset). Wald and Deviance match R element-wise.

### Changed

- **OpenMP auto-threshold loosened.** ``n_workers=0`` was gated on
  ``n*p^2 >= 1_000_000`` which only kicked in on quite large problems;
  measurement shows threading goes positive at ``n*p^2 >= 5_000`` and
  reaches 2x by 50,000. Threshold is now ``n*p^2 >= 10_000``, so
  ``n_workers=0`` actually parallelises on most non-trivial problems.

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
  module ``pylmrob._core._lmrob`` runs the whole fast-S + survivor
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
  ``pylmrob._core._fast_s`` runs the per-iteration body
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
  ``pylmrob.d_scale`` ports ``robustbase::lmrob..D..fit`` and the
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
- **M-scale.** `pylmrob.scale.m_scale` matches R's `find_scale` to
  `rtol=1e-9`.
- **Fast-S** resampling estimator and **MM** IRWLS iteration.
- **M-S** initial estimator for designs with categorical predictors
  (Maronna-Yohai 2000 alternating L1/S). Auto-detection via
  `init="auto"`.
- **Inference.** `vcov_avar1` (sandwich) ported from
  `robustbase/R/lmrob.MM.R:510-577`, matching R element-wise to `rtol=1e-3`.
  `vcov_w` implements the asymptotic-correction-factor branch.
- **Public API.** `pylmrob.lmrob(formula, data, ...)`, `LmRob`
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

[Unreleased]: https://github.com/anevolbap/pylmrob/compare/v0.5.8...HEAD
[0.5.8]: https://github.com/anevolbap/pylmrob/compare/v0.5.7...v0.5.8
[0.5.7]: https://github.com/anevolbap/pylmrob/compare/v0.5.6...v0.5.7
[0.5.6]: https://github.com/anevolbap/pylmrob/compare/v0.5.5...v0.5.6
[0.5.5]: https://github.com/anevolbap/pylmrob/compare/v0.5.4...v0.5.5
[0.5.4]: https://github.com/anevolbap/pylmrob/compare/v0.5.3...v0.5.4
[0.5.3]: https://github.com/anevolbap/pylmrob/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/anevolbap/pylmrob/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/anevolbap/pylmrob/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/anevolbap/pylmrob/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/anevolbap/pylmrob/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/anevolbap/pylmrob/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/anevolbap/pylmrob/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/anevolbap/pylmrob/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/anevolbap/pylmrob/releases/tag/v0.1.0
