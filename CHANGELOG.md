# Changelog

All notable changes to `pyrobustlm` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added (initial pre-alpha)

- **Build system.** `pyproject.toml` (meson-python backend), GPL-3 license,
  `NOTICE` crediting upstream `robustbase` authors. CI matrix
  (Linux/macOS/Windows × Python 3.10-3.13) and `cibuildwheel` workflow.
- **Phase 1.** R reference harness (`scripts/generate_r_reference.R`) that
  produces 35 JSON references covering 10 classical datasets and 18
  synthetic configurations. Pytest fixtures with rpy2 bridge for live R
  comparisons.
- **Phase 2.** `pyrobustlm.psi` exposes the six robustbase families
  (`bisquare`, `huber`, `hampel`, `optimal`, `lqq`, `ggw`). All match
  R's `Mpsi`/`Mchi`/`Mwgt` to rtol=1e-12 (1e-7 for ggw rho's polynomial
  approximation).
- **Phase 3.** `pyrobustlm.scale.m_scale` matches R's `find_scale` to
  rtol=1e-9.
- **Phase 4.** Fast-S resampling estimator in pure NumPy
  (`_fast_s.fast_s`).
- **Phase 6.** MM iteration as plain IRWLS with L1-norm convergence test,
  matching robustbase's `rwls` (no step-halving).
- **Phase 7.** Asymptotic-variance sandwich estimator (`vcov_avar1`).
- **Phase 8.** End-to-end public API: `pyrobustlm.lmrob(formula, data, ...)`
  and `LmRob` scikit-learn-style class. Formula handling via
  `formulaic`. Control with KS2014/KS2011/MM presets.
- **Phase 9.** Minimal four-panel diagnostic plot via `pyrobustlm.diagnostics.plot`.
- **Phase 10.** Validation sweep against all 10 classical datasets and
  three Hypothesis property tests for affine/scale/regression equivariance.

### Known limitations

- `init="M-S"` for factor designs is not yet implemented (Phase 5).
- `vcov_w` falls back to `vcov_avar1` with a warning.
- Performance is 2x-50x slower than R; Phase 11 will Cython-accelerate
  the resampling loop.
- RNG is not bit-identical with R (we use PCG64; R uses MT).

See [`docs/numerical-notes.md`](docs/numerical-notes.md) for the full list
of documented divergences from R.
