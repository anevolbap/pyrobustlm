# Quickstart

A 5-minute tour. We will fit the classic `stackloss` dataset with
`lmrob`, look at the summary, run a Wald test, and predict on new data.

## Install

From source (PyPI publication is pending):

```bash
pip install git+https://github.com/anevolbap/pylmrob
```

## Fit

```python
import pandas as pd
from pylmrob import Control, lmrob

df = pd.read_csv("stackloss.csv")
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    control=Control(nResample=500),
    seed=42,
)

print(f"converged: {fit.converged_}")
print(f"scale:     {fit.scale_:.4f}")
print(f"coef:      {dict(zip(fit.term_names_, fit.coef_))}")
```

## Inspect

`summary()` returns a `SummaryLmRob` whose `__str__` prints the R-style
table:

```python
print(fit.summary())
```

You can also pull the table out programmatically:

```python
summ = fit.summary()
summ.coefficients      # ndarray of shape (p, 4): est, se, t, p
summ.r_squared
summ.adj_r_squared
```

## Compare nested models

`anova` runs a robust Wald or Deviance test on nested fits:

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
```

The fit stores the formula's design transformation, so factor encoding
and `I(x**2)` style transforms are re-applied to `new`.

## Use a different psi family

```python
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    control=Control(psi="lqq"),
    seed=42,
)
```

Supported families: `bisquare` (default), `huber`, `hampel`, `optimal`,
`lqq`, `ggw`. Note that `lmrob` requires a redescending psi for the
S-step; `huber` is rejected.

## R compatibility presets

`Control(setting="KS2014")` and `Control(setting="KS2011")` configure the
SMDM pipeline (S init, MM, design-adaptive D-scale, MM with new scale)
matching R's `lmrob.control(setting="KS2014" / "KS2011")`:

```python
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    control=Control(setting="KS2014"),
    seed=42,
)
```

## Threading

The fast-S resampling loop can run in parallel via a thread pool:

```python
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    control=Control(n_workers=0),  # 0 = auto
    seed=42,
)
```

`n_workers=1` (default) is serial and bit-identical with single-threaded
runs. `n_workers=0` is auto and only enables threading when the problem
is large enough that BLAS dominates Python overhead.

## The fast Cython engine

The default `Control()` runs the fit through a monolithic Cython
kernel that does fast-S, MM, D-scale, and vcov in one nogil C block.
Median wall-clock is 0.93x R across the bench corpus; small-n
classical cases drop from 5-10x R to about 1x R. To force the legacy
numpy-based path, pass `Control(engine_c=False)`:

```python
fit = lmrob(
    "stack.loss ~ Air.Flow + Water.Temp + Acid.Conc.",
    df,
    control=Control(engine_c=False),
    seed=42,
)
```

On rare small datasets the Cython subset-draw lands in a basin
where the vcov is singular; `lmrob()` catches that and falls back
automatically. See [the engine_c page](engine_c) for the full
trade-off, which mainly concerns RNG byte-level reproducibility.
