# pyrobustlm

`pyrobustlm` is a Python port of R's [robustbase::lmrob](https://cran.r-project.org/package=robustbase),
the MM-estimator for robust linear regression. It targets feature parity
with R's reference implementation and ships Cython-accelerated kernels for
all six psi families.

## Why use this

Robust regression downweights outliers automatically. Where ordinary least
squares can be wrecked by a few bad rows, MM-estimators give bounded influence
and 95% Gaussian efficiency at the same time. R's `lmrob` is the standard
implementation; `pyrobustlm` brings that algorithm to Python with element-wise
agreement on coefficients, scale, and covariance for the classical datasets.

## At a glance

```python
import pandas as pd
from pyrobustlm import lmrob

df = pd.read_csv("stackloss.csv")
fit = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, seed=42)

print(fit.summary())
# lmrob(method='MM', psi='bisquare', setting=None)
# Coefficients:
#   (Intercept) -41.5246  ...
#   Air.Flow      0.9388  ...
#   ...
# Multiple R-squared: 0.9593,  Adjusted R-squared: 0.9521
```

## Status

`pyrobustlm` is at v0.2.0 and not yet on PyPI. Installation is from source:

```bash
pip install git+https://github.com/anevolbap/pyrobustlm
```

R parity (vs `robustbase` 0.99-7) on the classical datasets:

| Output | Tolerance |
|---|---|
| Coefficients | rtol=1e-3 (rtol=1e-6 on most) |
| Scale | rtol=1e-3 |
| Covariance (`vcov.avar1`) | rtol=1e-3 element-wise |
| Covariance (`vcov.w`) | rtol=1e-3 element-wise |
| `summary()` t-values, p-values | rtol=2e-3 |
| `anova()` chi-squared, p-value | rtol=2e-3 |

See [Numerical notes](numerical-notes) for the full list of documented
divergences from R.

## Contents

```{toctree}
:maxdepth: 2

quickstart
api
porting-from-r
numerical-notes
```

## Indices

- {ref}`genindex`
- {ref}`modindex`
