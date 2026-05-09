# Changelog

All notable changes to `pyrobustlm` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

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
