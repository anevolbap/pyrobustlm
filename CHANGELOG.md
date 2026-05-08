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

### Added (round 2)

- **Phase 5.** ``init="M-S"`` (and ``init="auto"``) now use a real
  Maronna-Yohai 2000 alternating L1/S estimator. ``formula.model_matrix``
  exposes per-column ``is_factor_col`` metadata so ``init="auto"`` can
  pick M-S when factors are present.
- **Phase 7.** ``vcov_avar1`` now ports the full robustbase formula from
  ``R/lmrob.MM.R:510-577``, including the correct ``Mchi(deriv=1) = chi'``
  normalisation and ``posdefify`` PSD projection. Matches R element-wise
  to ``rtol=1e-3`` on the validation corpus. ``vcov_w`` now implements
  the asymptotic-correction-factor branch.
- **Phase 9.** ``hatvalues`` and ``cooks_distance`` are now real, not
  stubs.
- **Phase 11.** Cython kernels for the bisquare psi family and the
  M-scale iteration. End-to-end ``lmrob()`` is 3-5x faster:
  100×5 fits in 94 ms (was 390 ms), 2000×20 in 277 ms (was 574 ms,
  vs R's 249 ms).
- **Phase 10.** Tightened classical-dataset validation tolerances:
  9 of 10 datasets now pass at ``rtol=1e-3``; ``hbk`` at ``rtol=5e-3``.
- New tests: ``tests/integration/test_factor_designs.py`` (M-S),
  ``tests/unit/test_inference.py`` (vcov_avar1 element-wise).

### Known limitations

- ``init="M-S"`` is a simplified port; not bit-identical with the C
  ``m_s_descent``.
- ``vcov_w`` Huber finite-sample correction and ``empirical``/``hybrid``/
  ``tau``/``tauold`` branches are deferred.
- Cython kernels currently exist only for the bisquare family; other
  families use the slower NumPy path.
- RNG is not bit-identical with R (we use PCG64; R uses MT).

See [`docs/numerical-notes.md`](docs/numerical-notes.md) for the full list
of documented divergences from R.
