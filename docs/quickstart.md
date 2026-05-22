# Quickstart

Fit a robust regression on `stackloss`, look at the summary, run a Wald
test, predict on new data. Five minutes.

## Install

```bash
pip install pylmrob
```

## Fit

```python
import pandas as pd
from pylmrob import Control, lmrob

df = pd.read_csv("stackloss.csv")
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    seed=42,
)

print(fit.summary())
```

`fit` is a [`LmRobResults`](api) object. The interesting fields are
`fit.coef_`, `fit.scale_`, `fit.rweights_` (the per-row robust weights),
and `fit.cov_`. Everything else is computed from those.

## Compare nested models

```python
from pylmrob import anova

full = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, seed=42)
red  = lmrob("stack.loss ~ Air.Flow + Water.Temp", df, seed=42)

print(anova(full, red))                       # Wald (default)
print(anova(full, red, test="Deviance"))      # Deviance variant
```

## Predict on new data

```python
new = pd.DataFrame({
    "Air.Flow":   [60.0, 70.0],
    "Water.Temp": [20.0, 25.0],
    "Acid.Conc.": [85.0, 88.0],
})

y_hat = fit.predict(new)

# With a confidence band (mean response) or prediction band (single obs):
band = fit.predict(new, interval="confidence", level=0.95)
# band[:, 0] = fitted, band[:, 1] = lower, band[:, 2] = upper
```

Both bands use the t-distribution at `fit.df_residual_`, matching R's
`predict.lmrob`. Skip the band scaling with `fit.predict_std(new,
kind="confidence")` if you want the raw standard errors.

## Inference

```python
fit.confint(level=0.95)                       # Wald (asymptotic normal)
fit.confint(method="bootstrap", n_boot=1000)  # bootstrap percentile

# statsmodels-style aliases also work:
fit.params, fit.bse, fit.tvalues, fit.pvalues, fit.conf_int(0.05)
```

When `n` is small or contamination is heavy, prefer the bootstrap CI:
the asymptotic Wald can be 30–50% too narrow. The bootstrap costs more
(each replicate is a full fit), so use `n_workers > 1` to parallelise:

```python
boot = fit.bootstrap(n_boot=2000, n_workers=4, seed=42)
boot.percentile_ci      # (p, 2) lower/upper
boot.se                 # bootstrap standard error per coef
boot.n_converged
```

## sklearn integration

```python
from pylmrob import LmRob
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

est = Pipeline([("scale", StandardScaler()), ("lmrob", LmRob())])
est.fit(X, y, lmrob__seed=42)
est.score(X_test, y_test)
```

`LmRob` is a `BaseEstimator + RegressorMixin`, so `GridSearchCV`,
`cross_val_score`, and friends work without adapters. The underlying
`LmRobResults` is at `est.named_steps["lmrob"].result_`. Full example
at [`examples/sklearn_pipeline`](examples/sklearn_pipeline).

## Pick a different psi family

```python
fit = lmrob("stack.loss ~ ...", df, control=Control(psi="lqq"), seed=42)
```

Available: `bisquare` (default), `hampel`, `optimal`, `lqq`, `ggw`,
`welsh`. `huber` is also available for separate M-scale work but isn't
permitted for the S-step (it's not redescending).

## R compatibility presets

```python
fit = lmrob("stack.loss ~ ...", df, control=Control(setting="KS2014"), seed=42)
```

`"KS2014"` and `"KS2011"` configure the SMDM pipeline (S init → MM →
D-scale → MM-on-new-scale), matching R's
`lmrob.control(setting="KS20XX")`.

## Speed levers

The default `Control()` already runs through the monolithic Cython
kernel; you don't need to opt in. Knobs worth knowing:

- `n_workers=N`: parallel fast-S resampling. `0` = auto, `1` = serial
  (deterministic), `>1` = explicit thread count.
- `nResample=N`: more resamples = more robust fit, linear runtime cost.
  Default 500 is enough for most problems.
- `engine_c=False`: forces the legacy NumPy path. Slower but useful
  for debugging.

See [Performance](engine_c) for the long version.

## Where to next

- Worked example: [stackloss tour](examples/stackloss_tour)
- Why robust at all: [robustness vs OLS](examples/robustness_vs_ols)
- Coming from R: [porting from R](porting-from-r)
- Got an error: [FAQ](faq)
