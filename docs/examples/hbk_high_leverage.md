# High-leverage outliers: the `hbk` dataset

`hbk` (Hawkins-Bradu-Kass) is a synthetic dataset designed to stress robust
estimators. It has 75 observations, 3 predictors, and 14 observations that
are simultaneously outliers *and* high-leverage. OLS-style methods that down-
weight by residual size alone can be fooled because the high-leverage points
have small residuals at the OLS fit.

## Loading the data

`hbk` ships with `robustbase`. The repo cached a copy under
`tests/data/hbk.csv`:

```python
import pandas as pd
hbk = pd.read_csv("tests/data/hbk.csv")
print(hbk.head(3))
print(f"shape: {hbk.shape}")
```

Columns: `X1`, `X2`, `X3`, `Y`. The first 14 rows are the contamination.

## OLS lies

```python
import numpy as np
import statsmodels.api as sm
X = sm.add_constant(hbk[["X1", "X2", "X3"]])
ols = sm.OLS(hbk["Y"], X).fit()
print("OLS coef:", dict(zip(ols.params.index, ols.params.round(3))))
```

```
OLS coef: {'const': -0.388, 'X1': 0.239, 'X2': -0.335, 'X3': 0.383}
```

These coefficients describe the contamination, not the underlying structure.
The first 14 obs are large and aligned in a way that gives them low OLS
residuals (they "look" like a valid trend).

## `lmrob` identifies the outliers

```python
from pylmrob import lmrob, Control
fit = lmrob("Y ~ X1 + X2 + X3", hbk, control=Control(nResample=1000), seed=42)
print("lmrob coef:", dict(zip(fit.term_names_, fit.coef_.round(3))))
print("Obs with weight < 0.1:", np.where(fit.rweights_ < 0.1)[0])
```

```
lmrob coef: {'Intercept': -0.181, 'X1': 0.082, 'X2': 0.039, 'X3': -0.052}
Obs with weight < 0.1: [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13]
```

The robust fit's coefficients are near zero (correct: by construction the
remaining 61 observations are noise around zero); the 14 contaminated points
are identified by their zero weights.

## Why this case is hard

The contamination is high-leverage (extreme in X-space) and aligned to look
like a fit. Robust methods that downweight purely by residual size, with no
account for leverage, fail on hbk. MM-estimators downweight by `psi(r/sigma)`
where `sigma` itself is robust (from the S-step); the S-step's high-breakdown
scale catches the contamination before the M-step locks in.

Setting `nResample=1000` (above the default 500) helps the fast-S resampling
hit a clean subset; on `hbk` the default 500 is sometimes enough but the
extra robustness is cheap.

## Diagnostics

```python
diag = fit.diagnostics()
print(f"Robust R-squared: {diag.r_squared_robust:.3f}")
print(f"Robust scale: {fit.scale_:.3f}")
```

A high robust R² and a small `scale_` (close to the clean-data noise level)
confirms the fit captured the underlying structure.

## R-side reference

```r
library(robustbase)
data(hbk)
fit <- lmrob(Y ~ ., data = hbk)
summary(fit)
```

pylmrob's output agrees with R's to `rtol=5e-3` on the coefficients (`hbk` is
the most challenging case in our
[bench corpus](https://github.com/anevolbap/pyrobustlm/blob/main/docs/bench-report.md)).
