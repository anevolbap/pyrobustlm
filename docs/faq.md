# FAQ / troubleshooting

## "Algorithm did not converge"

The MM IRWLS hit `max_it` without converging. The three usual causes:

- **Near-singular design.** Check `np.linalg.cond(fit.design_x_)`; if
  it's much above 1e10 you have near-collinear columns. Drop or merge.
- **Bad initial scale.** Increase `Control(nResample=1000)` or switch
  to `Control(setting="KS2014")` (which uses a more stable D-scale).
- **Contamination above 50%.** Past the breakdown point of the
  M-estimator. Robust regression can't fix it; you need a different
  model.

## "X' W X is singular; consider cov='.vcov.w'"

Almost every observation has been downweighted to zero, so the
asymptotic covariance is rank-deficient. Either:

- Switch to `Control(cov=".vcov.w")` (more stable when many weights
  hit zero), or
- Bootstrap: `fit.confint(method="bootstrap", n_boot=2000)`.

The default `engine_c=True` path already retries with the NumPy basin
on this error. If it bubbles up, both basins hit the same wall.

## "Wald CI looks suspiciously narrow"

For small `n` or heavy contamination the asymptotic Wald CI
underestimates. Use the bootstrap; on `n < 50` it's typically 2-3×
wider:

```python
fit.confint(method="bootstrap", n_boot=2000)
```

## Confidence vs prediction vs tolerance interval

- **Confidence** (`predict(interval="confidence")`): range for the
  *mean response* at the new row. Shrinks as `n` grows.
- **Prediction** (`predict(interval="prediction")`): range for a
  *single new observation*. Strictly wider; doesn't shrink with `n`
  because per-observation noise stays.
- **Tolerance**: covers a fraction of the population. Not built in;
  compose your own using `fit.predict_std(kind="prediction")` + a
  tolerance factor.

`fit.predict_std(new, kind=...)` returns the SE directly if you'd
rather build bands at a non-t distribution.

## How do I get bit-identical fits to R?

`pylmrob`'s default uses NumPy's PCG64 BitGenerator. R uses Mersenne
Twister, so the fast-S subsets differ. Agreement on the validation
corpus is `rtol=1e-3` on default settings.

For tighter agreement, drive the resample loop through R's exact
`unif_rand` stream:

```python
fit = lmrob(formula, df, control=Control(rng="R"), seed=42)
```

Stackloss-class agreement: `rtol=1e-4` across seeds. The residual
gap is BLAS-precision (different LAPACK implementations accumulate
different ULPs in IRWLS); see [`numerical-notes`](numerical-notes)
for the worked investigation. Truly bit-identical *fits* require
sharing R's BLAS — use rpy2 if you need that.

Lower-level helpers: `r_set_seed`, `r_sample_noreplace`,
`r_subsample_nonsingular` give you the byte-identical uniform stream
and subset draws standalone (see [`rng-r-perf`](rng-r-perf)).

## How do I pick `nResample`?

Default `nResample=500` matches R and covers the standard cases.
Bump to 1000-2000 if:

- The data is small with a difficult contamination pattern.
- You're on `init="M-S"` and want tighter R agreement.

Cost is linear; with `n_workers=4` even `2000` is sub-second.

## How do I make it faster?

`engine_c=True` (default) already runs the whole fit in one Cython
block. Beyond that:

- `Control(n_workers=4)` for OpenMP parallel resampling. 2-3× speedup
  on `n > 5000`.
- For many small fits in a loop, the formula parser adds ~3 ms each.
  Use `LmRob(...).fit(X, y)` (raw arrays) to skip it.

See [`engine_c`](engine_c) for the long version.

## `LmRob` vs `lmrob`?

Same machinery, two faces:

- `lmrob(formula, data, ...)`: R-style. Most expressive (formulas,
  factors, transforms).
- `LmRob(control).fit(X, y)`: sklearn-style. Drops into `Pipeline`,
  `cross_val_score`, `GridSearchCV`. The underlying result is at
  `est.result_`.

## `predict()` is returning weird values

Two pitfalls:

- **Array path but factor design.** A raw NumPy `X_new` must match
  the *encoded* design (intercept included). Pass a DataFrame and let
  the stored formula spec do the encoding.
- **Missing columns.** The DataFrame needs the same column names the
  formula referenced. Check `new.columns`.

## Why does `setting="KS2014"` use `psi="lqq"` by default?

It's Koller & Stahel's recommended pairing. `lqq` has slightly heavier
downweighting at the bend than bisquare, which complements the D-scale
refinement. Override with `Control(setting="KS2014", psi="bisquare")`
if you want.

## See what the fitter is doing

```python
Control(trace_lev=2)   # per-iteration output, matches R's trace.lev=2
```

`trace_lev=3` and `=4` are progressively noisier.

Or for a post-hoc diagnostic:

```python
fit.summary(detail="full")   # appends init method, scale, engine info
```

## Where do I report bugs?

[GitHub issues](https://github.com/anevolbap/pyrobustlm/issues).
Include `pylmrob.__version__`, a minimal reproducer, and (if you have
R handy) what `robustbase::lmrob` says on the same data.
