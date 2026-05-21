# Robustness vs OLS: a contaminated regression

A short example showing what robust regression actually buys you over ordinary
least squares. We'll simulate a clean linear relationship, then inject a few
high-leverage outliers, and watch what each method does.

## Setup

```python
import numpy as np
import pandas as pd
from pylmrob import lmrob, Control

rng = np.random.default_rng(0)
n = 100
x = np.linspace(0, 10, n)
y_clean = 2.0 + 1.5 * x + 0.5 * rng.standard_normal(n)

# 10 outliers: high-leverage (large x), wildly wrong y.
y = y_clean.copy()
y[-10:] = -20.0 + 3.0 * rng.standard_normal(10)
df = pd.DataFrame({"x": x, "y": y})
```

## OLS

```python
import statsmodels.api as sm  # only used here for the side-by-side
ols = sm.OLS(y, sm.add_constant(x)).fit()
print(f"OLS intercept = {ols.params[0]:+.3f}, slope = {ols.params[1]:+.3f}")
```

Output (your numbers will be close):

```
OLS intercept = +5.124, slope = +0.197
```

The 10 outliers at the right end pull the OLS line down: the slope drops from
the true 1.5 to ~0.2, and the intercept jumps to absorb the bias.

## `lmrob`

```python
fit = lmrob("y ~ x", df, control=Control(), seed=42)
print(f"lmrob intercept = {fit.coef_[0]:+.3f}, slope = {fit.coef_[1]:+.3f}")
```

Output:

```
lmrob intercept = +2.012, slope = +1.497
```

`lmrob` recovers the true relationship (`(2.0, 1.5)`) almost exactly. The
M-estimator downweighted the 10 contaminated points to near zero. You can
see the weights:

```python
print("Weights of the last 10 obs (the contaminated block):")
print(fit.rweights_[-10:].round(3))
```

```
[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
```

All ten outliers got zero weight. The clean obs are at or near 1.

## What changed and why

OLS minimises the sum of squared residuals; one large residual contributes a
lot of squared error, so the fitted line bends toward outliers to reduce the
total. `lmrob` minimises the sum of `rho(r/sigma)` for a bounded `rho`
(bisquare by default), so once a residual exceeds the rejection threshold its
contribution caps out and the fit stops chasing it.

The 50% breakdown point of MM-estimators says you'd need >50% of the points
contaminated before the robust fit could itself be wrecked. With 10% here we
have plenty of room.

## When to use which

OLS is great when residuals are roughly Gaussian and you've already cleaned the
data. `lmrob` is for the situation you don't know whether your data is clean:
it costs nothing to use and produces ~95% of OLS efficiency on truly clean
data.

For more on the math, see [theory](../theory.md). For tuning, see
[FAQ](../faq.md).
