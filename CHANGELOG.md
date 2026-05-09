# Changelog

All notable changes to `pyrobustlm` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

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

[Unreleased]: https://github.com/anevolbap/pyrobustlm/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/anevolbap/pyrobustlm/releases/tag/v0.1.0
