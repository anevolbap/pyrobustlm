# M-S init for factor designs

Robust regression with a mix of continuous and categorical predictors uses a
different initial estimator than the pure-continuous case. `lmrob` switches
automatically when it sees factor columns in the design; `Control(init="auto")`
does it for you, but understanding the mechanics helps when reading the trace
output or debugging convergence.

## The problem

A full S-estimator on a design with `k` factor levels needs subsamples that
include at least one observation from each level. With small per-level counts
this becomes unlikely; the resampling either fails or lands on degenerate
subsets.
[Maronna and Yohai (2000)](https://doi.org/10.1016/S0378-3758(99)00208-6)
solve this with a two-step estimator: an L1 estimate for the factor effects,
then an M estimate over the continuous part conditional on the factor
effects.

## A small worked example

```python
import numpy as np
import pandas as pd
from pylmrob import Control, lmrob

rng = np.random.default_rng(0)
n = 80
group = np.repeat(["A", "B", "C", "D"], n // 4)
x = rng.standard_normal(n)
# True model: y = 1.5 * x + group effects (-2, 0, +1, +3)
group_effect = {"A": -2.0, "B": 0.0, "C": 1.0, "D": 3.0}
y = 1.5 * x + np.array([group_effect[g] for g in group]) + 0.4 * rng.standard_normal(n)
# Inject 6 outliers in group D.
y[60:66] += 25.0
df = pd.DataFrame({"x": x, "group": pd.Categorical(group), "y": y})
```

## Default `lmrob` (auto-picks M-S)

```python
fit = lmrob("y ~ x + group", df, control=Control(), seed=42)
print(fit.summary())
```

The output includes `method='MM'` and `init='M-S'` in the header. `lmrob`
saw the factor column and chose the M-S path automatically.

## Forcing pure-S init

```python
fit_s = lmrob("y ~ x + group", df, control=Control(init="S"), seed=42)
```

This may or may not work, depending on the seed: fast-S needs subsets that
span the factor levels. On a small design (`n=80`, 4 levels), the failure
rate is non-trivial. The M-S init avoids that.

## Inspecting which observations got downweighted

```python
print("Observations with weight < 0.1:")
print(np.where(fit.rweights_ < 0.1)[0])
```

The 6 outliers in group D show up here.

## When this matters

- Whenever your design has categorical predictors with few levels per group.
- Whenever fast-S is reporting "no non-singular subsample found" warnings.
- ANCOVA-style models with a continuous covariate and a treatment factor.

R-side reference: `lmrob.control(method="MM", init="M-S")`. pylmrob's
default `Control(init="auto")` makes the same choice.
