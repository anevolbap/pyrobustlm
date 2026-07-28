# Porting from R's robustbase::lmrob

A side-by-side cheat sheet for users coming from R's `robustbase` package.

## Fit

| R | Python |
|---|---|
| `lmrob(y ~ x1 + x2, data = df)` | `lmrob("y ~ x1 + x2", df)` |
| `lmrob(y ~ ., data = df, control = lmrob.control())` | `lmrob("y ~ x1 + x2", df, control=Control())` |
| `lmrob.control(psi = "lqq", setting = "KS2014")` | `Control(psi="lqq", setting="KS2014")` |
| `set.seed(42); lmrob(...)` | `lmrob(..., seed=42)` |

## Result accessors

| R | Python |
|---|---|
| `fit$coefficients` | `fit.coef_` |
| `fit$scale` | `fit.scale_` |
| `fit$residuals` | `fit.residuals_` |
| `fit$fitted.values` | `fit.fitted_` |
| `fit$rweights` | `fit.rweights_` |
| `fit$cov` | `fit.cov_` |
| `fit$df.residual` | `fit.df_residual_` |
| `fit$converged` | `fit.converged_` |

## Inference

| R | Python |
|---|---|
| `summary(fit)` | `print(fit.summary())` |
| `summary(fit)$coefficients` | `fit.summary().coefficients` |
| `summary(fit)$r.squared` | `fit.summary().r_squared` |
| `confint(fit, level = 0.95)` | `fit.confint(level=0.95)` |
| `anova(full, red)` | `anova(full, red)` |
| `anova(full, red, test = "Deviance")` | `anova(full, red, test="Deviance")` |
| `predict(fit, newdata = nd)` | `fit.predict(nd)` |

## Control fields

The most-used R control fields and their Python equivalents:

| R control field | Python `Control` field |
|---|---|
| `nResample` | `nResample` |
| `tuning.psi` | `tuning_psi` |
| `tuning.chi` | `tuning_chi` |
| `psi` | `psi` |
| `method` | `method` |
| `cov` | `cov` |
| `setting` | `setting` |
| `max.it` | `max_it` |
| `k.max` | `k_max` |
| `refine.tol` | `refine_tol` |
| `seed` | `seed` |
| `subsampling` | `subsampling` (``"simple"`` / ``"nonsingular"``) |
| (no equivalent) | `rng` (``"PCG64"`` / ``"MT19937"`` / ``"R"``) |
| (no equivalent) | `n_workers` (parallel resampling) |

## Default behaviour

`Control()` with no arguments matches R's `lmrob.control()` (default
`setting=NULL`): `psi="bisquare"`, `method="MM"`, `cov=".vcov.avar1"`,
`tuning_psi=4.685061`, `tuning_chi=1.54764`.

`Control(setting="KS2014")` matches R's `lmrob.control(setting="KS2014")`:
`psi="lqq"`, `method="SMDM"`, `cov=".vcov.w"`, `tuning_psi=(1.4734061,
0.9822707, 1.5)`, `tuning_chi=(0.4015457, 0.2676971, 1.5)`.

`Control(setting="KS2011")` is the same as `KS2014` but with
the older D-scale tuning.

## Things that don't match exactly

The default fast-S resampling RNG uses NumPy's PCG64; R uses Mersenne
Twister. With the same seed, the two diverge at the first random draw,
so basin-of-attraction sensitivity can give small coefficient
differences on tiny-n datasets. For tighter agreement, opt in to
``Control(rng="R")`` (v0.5.16+): the resample loop runs against R's
exact ``unif_rand`` stream and stackloss fits agree with ``lmrob`` to
rtol~1.7e-5. See the [numerical notes](numerical-notes) for the
full list.
