# Walkthrough: the `stackloss` dataset

A complete robust-regression analysis on Brownlee's `stackloss` data, the
canonical small example used by every robust-regression paper since
Hampel et al. (1986). The dataset has 21 observations and one outlier
cluster (observations 1, 3, 4, 21) that OLS fits very badly. We will:

1. Compare OLS and `lmrob` coefficient estimates.
2. Identify outliers from the robust weights and from `fit.diagnostics()`.
3. Predict with confidence and prediction bands.
4. Run a nested-model Wald test.

The full script runs in under a second on a laptop.

## Setup

```python
import numpy as np
import pandas as pd
from pylmrob import lmrob, Control

# robustbase ships stackloss; for portability we use R's copy via
# scripts/benchmark.py's CSV cache, or just type it in.
stackloss = pd.DataFrame({
    "Air.Flow":   [80, 80, 75, 62, 62, 62, 62, 62, 58, 58, 58, 58, 58, 58, 50, 50, 50, 50, 50, 56, 70],
    "Water.Temp": [27, 27, 25, 24, 22, 23, 24, 24, 23, 18, 18, 17, 18, 19, 18, 18, 19, 19, 20, 20, 20],
    "Acid.Conc.": [89, 88, 90, 87, 87, 87, 93, 93, 87, 80, 89, 88, 82, 93, 89, 86, 72, 79, 80, 82, 91],
    "stack.loss": [42, 37, 37, 28, 18, 18, 19, 20, 15, 14, 14, 13, 11, 12,  8,  7,  8,  8,  9, 15, 15],
})
```

## OLS baseline

```python
import statsmodels.api as sm

X = sm.add_constant(stackloss[["Air.Flow", "Water.Temp", "Acid.Conc."]])
ols = sm.OLS(stackloss["stack.loss"], X).fit()
print(ols.params)
# const         -39.92
# Air.Flow        0.72
# Water.Temp      1.30
# Acid.Conc.     -0.15
```

OLS pulls the slopes toward the outlier cluster. `Water.Temp` is
inflated and `Air.Flow` is suppressed compared to the robust fit.

## Robust fit

```python
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    stackloss,
    control=Control(setting="KS2014"),
    seed=42,
)
print(fit)
```

```
lmrob(method="MM", psi="lqq")

Coefficients:
   Intercept     Air.Flow   Water.Temp   Acid.Conc.
      -42.30        0.939        0.580       -0.113
```

The `Air.Flow` slope is about 30% larger, `Water.Temp` 55% smaller, and
`Acid.Conc.` is essentially unchanged. The robust fit landed on the
"true" linear relationship; OLS landed on a compromise contaminated by
the 4-observation outlier cluster.

For the full coefficient table with std errors, t and p values, and the
robust R², use `summary()`:

```python
print(fit.summary())
```

## Identify outliers

The robust weights `fit.rweights_` are `psi(r / sigma) / (r / sigma)`;
values near zero mean the observation was effectively dropped.
`diagnostics()` packages this together with leverage and Cook's distance:

```python
diag = fit.diagnostics(outlier_threshold=2.5)
flagged = np.where(diag.outliers)[0] + 1  # 1-indexed for R compatibility
print("flagged observations:", flagged.tolist())
# flagged observations: [1, 3, 4, 21]

# Inspect why each was flagged
for i in flagged - 1:
    print(
        f"obs {i+1:>2}: rweight={diag.rweights[i]:.3f}  "
        f"std_resid={diag.std_residuals[i]:+.2f}  "
        f"leverage={diag.leverage[i]:.3f}  "
        f"cook's D={diag.cooks_distance[i]:.3f}"
    )
```

These four observations are the documented outliers in stackloss; the
robust fit identified them automatically without any user input. OLS
would not have flagged them (its residuals get smeared by the outliers
themselves).

## Confidence and prediction intervals

`predict(interval=...)` gives the band for the conditional mean
(`"confidence"`) or for a single new observation (`"prediction"`):

```python
new = pd.DataFrame({
    "Air.Flow":   [60, 75],
    "Water.Temp": [20, 24],
    "Acid.Conc.": [85, 88],
})
band = fit.predict(new, interval="confidence", level=0.95)
for x, (yhat, lo, hi) in zip(new.itertuples(index=False), band):
    print(f"x={tuple(x)} -> y_hat={yhat:.2f}  CI95=[{lo:.2f}, {hi:.2f}]")
```

The prediction interval (`interval="prediction"`) is wider because it
includes residual noise from the M-scale `fit.scale_`:

```python
pi = fit.predict(new, interval="prediction", level=0.95)
print("PI width:", pi[:, 2] - pi[:, 1])
print("CI width:", band[:, 2] - band[:, 1])
```

## Nested-model test

`anova()` runs a Wald or Deviance test between two nested fits.
Method-style is `full.anova(reduced)`:

```python
red = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp",
    stackloss,
    control=Control(setting="KS2014"),
    seed=42,
)
print(fit.anova(red))
# Wald: large test statistic on Acid.Conc. with p ~ 0.7, suggesting
# the coefficient is not significantly different from zero in this
# small dataset.
```

## Statsmodels-style access

If your downstream code expects a statsmodels-shaped result, the same
attributes work without an adapter:

```python
fit.params         # coef_
fit.bse            # standard_errors_
fit.tvalues        # coef / se
fit.pvalues
fit.conf_int(alpha=0.05)
```

## Sklearn pipeline

For grid search or cross-validation:

```python
from pylmrob import LmRob
from sklearn.model_selection import cross_val_score

est = LmRob()
X = stackloss[["Air.Flow", "Water.Temp", "Acid.Conc."]].to_numpy()
y = stackloss["stack.loss"].to_numpy()
scores = cross_val_score(est, X, y, cv=5)
print(scores.mean())
```

## What to read next

- {doc}`../numerical-notes` for the details on RNG basin drift and how
  `pylmrob` matches `robustbase` numerically.
- {doc}`../engine_c` for the fast Cython kernel and `n_workers`.
- {doc}`../porting-from-r` if you're coming from `robustbase`.
