# pylmrob

`pylmrob` is a Python port of R's [robustbase::lmrob](https://cran.r-project.org/package=robustbase),
the MM-estimator for robust linear regression. It targets feature parity
with R's reference implementation and ships Cython-accelerated kernels for
all six psi families.

## Why use this

Robust regression downweights outliers automatically. Where ordinary least
squares can be wrecked by a few bad rows, MM-estimators give bounded influence
and 95% Gaussian efficiency at the same time. R's `lmrob` is the standard
implementation; `pylmrob` brings that algorithm to Python with element-wise
agreement on coefficients, scale, and covariance for the classical datasets.

## At a glance

```python
import pandas as pd
from pylmrob import lmrob

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

Alpha. Available on PyPI:

```bash
pip install pylmrob
```

Feature parity with R's `lmrob` for the common path (all five psi
families, S/M-S init, MM, KS2014/KS2011 settings, per-case weights,
both vcov flavours, anova). Median wall-clock 0.93x R across the bench
corpus.

R parity (vs `robustbase` 0.99-7) on the classical datasets, default
`Control(rng="PCG64")`:

| Output | Tolerance |
|---|---|
| Coefficients | rtol=1e-3 (rtol=1e-6 on most) |
| Scale | rtol=1e-3 |
| Covariance (`vcov.avar1`) | rtol=1e-3 element-wise |
| Covariance (`vcov.w`) | rtol=1e-3 element-wise |
| `summary()` t-values, p-values | rtol=2e-3 |
| `anova()` chi-squared, p-value | rtol=2e-3 |

`Control(rng="R")` (v0.5.16+) drives the resample loop through R's exact
`unif_rand` stream and tightens coefficient agreement to rtol~1.7e-5 on
stackloss. See [Numerical notes](numerical-notes) and
[`rng-r-perf`](rng-r-perf) for details.

## Contents

```{toctree}
:maxdepth: 2

quickstart
theory
examples/stackloss_tour
examples/robustness_vs_ols
examples/m_s_factors
examples/hbk_high_leverage
examples/sklearn_pipeline
api
porting-from-r
engine_c
rng-r-perf
numerical-notes
faq
```

## Indices

- {ref}`genindex`
- {ref}`modindex`
