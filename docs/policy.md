# API stability and deprecation policy

`pylmrob` is beta as of v0.6.0. The public surface listed below is
covered by a compatibility contract; anything outside that list is
internal and may change in any release.

## What is public

The supported public API is everything in `pylmrob.__all__`:

| Symbol            | Kind        | Purpose                                           |
|-------------------|-------------|---------------------------------------------------|
| `lmrob`           | function    | Formula-based MM-estimator fit                    |
| `LmRob`           | class       | sklearn-style estimator (`fit`/`predict`/`score`) |
| `LmRobResults`    | class       | Return type of `lmrob()`                          |
| `Control`         | dataclass   | Tuning parameters (psi, nResample, setting, ...)  |
| `anova`           | function    | Nested-model Wald / Deviance test                 |
| `bootstrap`       | function    | Bootstrap inference helper                        |
| `__version__`     | string      | Installed pylmrob version                         |

Methods on `LmRobResults` that are part of the contract:
`coef_`, `scale_`, `cov_`, `residuals_`, `fitted_`, `rweights_`,
`df_residual_`, `converged_`, `summary()`, `predict()`,
`predict_std()`, `confint()`, `bootstrap()`, `diagnostics()`,
plus the statsmodels-style aliases (`params`, `bse`, `tvalues`,
`pvalues`, `conf_int`).

## What is private

Anything under a leading-underscore module
(`pylmrob._fast_s`, `pylmrob._psifuns`, `pylmrob._mm`,
`pylmrob._l1`, `pylmrob._utils`, `pylmrob._core`) is internal.
Importing from these is not supported and may break in any release,
including patch releases.

The R-compatibility helpers (`r_set_seed`, `r_norm_rand`, `r_qnorm`,
`r_sample_noreplace`, `r_subsample_nonsingular`, `RState`,
`make_generator`) stay importable from `pylmrob` and `pylmrob.rng`
for the validation workflow but are not part of the day-to-day public
API. Treat them as advanced.

## Deprecation policy

Before any breaking change to the public surface:

1. The new behavior ships in release `X.Y.Z` alongside a
   `DeprecationWarning` on the old behavior.
2. The warning is emitted for at least one full minor release cycle
   (so users have at least one `X.(Y+1).0` to see the warning).
3. The old behavior is removed no earlier than `X.(Y+2).0`.

A "breaking change" here means: removing a public symbol, removing a
keyword argument, narrowing the accepted types of a public argument,
changing default behavior in a way that alters numerical output by
more than the documented rtol vs R, or moving a symbol such that
existing imports break.

Pure additions (new keyword arguments with default values, new helper
functions, new optional `Control` fields) are not breaking and may
land in any release.

## Numerical agreement

The numerical contract with `robustbase::lmrob` is documented in
[numerical notes](numerical-notes). The headline numbers:

- Coefficients agree to `rtol=1e-3` element-wise.
- Scale agrees to `rtol=1e-3`.
- `summary()` t/p values agree to `rtol=2e-3`.

These tolerances are deliberately loose; sub-`1e-5` agreement is
available with `Control(rng="R")` but is best-effort and may shift
slightly across `pylmrob` releases as the BLAS, scipy, or robustbase
itself evolves. A measurable widening of the gap on a fixed corpus
is a regression and will be treated as a `fix:` candidate.

## Versioning

Semantic versioning (`MAJOR.MINOR.PATCH`). While the package is on
the `0.x` series:

- `PATCH` (`0.y.z` -> `0.y.(z+1)`): bug fixes, doc updates, perf
  improvements, no public-surface changes.
- `MINOR` (`0.y.z` -> `0.(y+1).0`): new features, deprecation
  announcements, deprecated-symbol removals from a prior cycle.
- `MAJOR` (`0.y.z` -> `1.0.0`): the API freezes for a real `1.0`
  release once the wider scientific Python ecosystem (conda-forge
  feedstock, downstream packages) has had a few minor cycles to
  catch up.

## Reporting an issue

If you spot a deprecation that wasn't announced, an undocumented
breaking change, or numbers drifting outside the contracted
tolerance, please open an issue. Include the previous and current
`pylmrob` version, the platform, and a minimal reproducer.
