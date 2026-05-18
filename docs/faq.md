# FAQ / Troubleshooting

Short answers to recurring questions.

## "Algorithm did not converge"

The MM IRWLS hit `max_it` without converging. Likely causes:

- **Near-singular design.** `X` has nearly collinear columns. Check
  `np.linalg.cond(fit.design_x_)`; values much bigger than 1e10 usually
  mean trouble. Drop or combine columns.
- **Bad initial scale.** The MM step holds `σ` fixed from the S-step.
  If the S-step landed in a basin where `σ_S` is too small (perfect
  fit on a subset of clean points, ignoring real signal), MM never
  recovers. Increase `Control(nResample=)` to get more S candidates,
  or switch to `Control(setting="KS2014")` which uses a more stable
  D-scale.
- **Heavy contamination beyond 50%.** The M-estimator has a 50%
  breakdown point. Past that, nothing in `lmrob`'s family can save
  you; you need a different model.

## "X' W X is singular; consider cov='.vcov.w'"

The asymptotic vcov `(X' diag(psi'(r/σ)) X)^{-1}` is singular at the
fitted point, usually because nearly every observation has been
downweighted to zero. Two fixes:

- Use `Control(cov=".vcov.w")` instead of the default `.vcov.avar1`.
  `.vcov.w` uses a different correction factor that's stable when
  many weights are at zero.
- Bootstrap the inference instead: `bootstrap(fit, n_boot=2000)`.
  Robust to the singular-vcov issue.

The `engine_c=True` path catches this `FloatingPointError` and
automatically retries with `engine_c=False`, which has a different
RNG-driven basin and usually lands on a non-singular fit. If you see
this raised it means both basins hit the same wall.

## "Wald CI looks suspiciously narrow"

For small `n`, the asymptotic Wald CI underestimates real variability.
Use `bootstrap(fit, n_boot=2000)`; the percentile CIs are usually
wider and more honest. Typical pattern on `n < 50` with outliers:
Wald CI is 30-50% the width of the bootstrap CI.

## How do I get bit-identical fits to R?

You can't, because `pylmrob` uses NumPy's PCG64 BitGenerator while R
uses the Mersenne Twister. The set of resampled p-subsets in fast-S
differs between the two RNGs, so they can land on different basins.
Coefficient agreement on the validation corpus is within `rtol=1e-3`
for the well-conditioned cases. See
{doc}`numerical-notes` for the per-case tolerances.

If you really need bit-identical R output, drive R via `rpy2` and skip
`pylmrob` entirely.

## How do I pick `nResample`?

The default `nResample=500` matches R's `lmrob.control()` default and
is fine for all standard datasets. Bumping it (1000, 2000) helps when:

- You're on a small dataset with a difficult contamination pattern and
  the fit is sensitive to the initial S basin.
- You're using `init="M-S"` for a factor design and want tighter
  agreement with R.

Cost is linear in `nResample`. With `engine_c=True` and `n_workers=4`
even `nResample=2000` runs in well under a second on the standard
datasets.

## How do I make `lmrob` faster?

It's already pretty fast:

- `engine_c=True` is the default since v0.5.11.
- `Control(n_workers=4)` (or some other `>1` value) parallelises the
  resample loop via OpenMP. Get 2-3x speedup on n > 5000.
- For batch fits in a loop, the per-call Python overhead is ~3 ms; if
  you're fitting many small models, that dominates. Use the array API
  (`LmRob(...).fit(X, y)`) to skip the formula parser.

## Should I use `LmRob` or `lmrob`?

Different APIs for different audiences:

- `lmrob(formula, data, control, weights, seed)`: R-style, takes a
  formula and a DataFrame. Most expressive.
- `LmRob(control).fit(X, y)`: sklearn-style, takes raw arrays. Fits
  into `Pipeline`, `cross_val_score`, `GridSearchCV`. The underlying
  result is at `est.result_`.

Both end up calling the same fit machinery.

## Why is `predict()` returning weird values?

Two common pitfalls:

- **You passed a NumPy array but the fit was on a formula with
  factors.** The array path expects a design matrix shaped exactly
  like the original (intercept column included). Use the DataFrame
  path: the stored formula spec re-applies the encoding.
- **The DataFrame is missing columns.** `predict(new_data)` looks for
  the same column names the formula referenced; missing columns raise
  `ValueError`. Check `new_data.columns`.

## Why does `setting="KS2014"` use `psi="lqq"`?

That's Koller & Stahel's recommended pairing. `lqq` (the
linear-quadratic-quadratic family) has slightly heavier downweighting
at the bend than bisquare, which complements the D-scale refinement
used by the KS2014 setting. You can override it: `Control(setting="KS2014", psi="bisquare")`.

## Where can I see what's actually happening?

Set `Control(trace_lev=2)` to get per-iteration output (matches R's
`lmrob.control(trace.lev = 2)`). Useful for debugging convergence
problems. `trace_lev=3` is more verbose; `trace_lev=4` is debug-level
and very noisy.

## Where do I report bugs?

GitHub issues: https://github.com/anevolbap/pyrobustlm/issues.

Include:
- `pylmrob.__version__`, your Python version, OS.
- A minimal reproducer (the data + the call).
- What R says on the same data (if applicable), so we can tell whether
  it's a `pylmrob` bug or a basin-drift footnote.
