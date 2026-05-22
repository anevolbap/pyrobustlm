# pylmrob

Python port of R's
[robustbase::lmrob](https://cran.r-project.org/package=robustbase) MM-estimator
for robust linear regression. Same algorithm, same numbers (to `rtol=1e-3` on
the classical datasets), Cython-accelerated, and roughly as fast as R itself.

```bash
pip install pylmrob
```

```python
from pylmrob import lmrob

fit = lmrob("stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.", df, seed=42)
print(fit.summary())
```

## Why robust regression

OLS gets pulled toward outliers. A few bad rows in `y` (or in `X`) can move
your slope by an order of magnitude. MM-estimators downweight those rows
automatically and still recover ~95% of OLS efficiency on clean data. The
canonical reference is `robustbase::lmrob`; pylmrob is that algorithm,
unchanged where it matters, in Python.

## Quick numbers

| Output | Agreement with R |
|---|---|
| Coefficients | `rtol=1e-3` (`rtol=1e-6` on most rows) |
| Scale | `rtol=1e-3` |
| `vcov.avar1`, `vcov.w` | `rtol=1e-3` element-wise |
| `summary()` t/p values | `rtol=2e-3` |
| `anova()` chi-sq, p | `rtol=2e-3` |

For sub-`1e-5` reproducibility against R's exact RNG stream, opt in with
`Control(rng="R")`. See [`rng-r-perf`](rng-r-perf) and
[Numerical notes](numerical-notes) for the details.

Median wall-clock on the bench corpus: **0.93× R**. Status: alpha; the
public API is stable enough to depend on but may change.

## Contents

```{toctree}
:maxdepth: 2
:caption: Getting started

quickstart
theory
```

```{toctree}
:maxdepth: 2
:caption: Theory in pictures

notebooks/01_ols_vs_robust
notebooks/02_efficiency
notebooks/03_breakdown
notebooks/04_s_estimator
```

```{toctree}
:maxdepth: 2
:caption: How-to examples

examples/stackloss_tour
examples/m_s_factors
examples/sklearn_pipeline
```

```{toctree}
:maxdepth: 2
:caption: Reference

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
