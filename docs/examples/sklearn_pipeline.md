# `LmRob` in sklearn pipelines

`pylmrob.LmRob` is an sklearn-shaped estimator. It implements `fit`,
`predict`, `score`, `get_params`, `set_params`, and follows the
sklearn regressor contract for `score` (returns OLS R² on the test
set, so `cross_val_score` and `GridSearchCV` work as written).

This page shows three patterns that come up in practice.

## 1. Drop-in inside a `Pipeline`

`LmRob` doesn't centre or standardise the design (just like
`LinearRegression`), so pairing it with a `StandardScaler` or
`PolynomialFeatures` is a natural fit:

```python
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from pylmrob import LmRob, Control

rng = np.random.default_rng(0)
n_train, n_test = 300, 200
X_train = rng.standard_normal((n_train, 5))
true_beta = np.array([1.0, -0.5, 2.0, 0.0, 0.3])
y_train = X_train @ true_beta + 0.2 * rng.standard_normal(n_train)
y_train[:30] += 30  # contaminate 10% of the training rows

# Clean held-out test set — no contamination.
X_test = rng.standard_normal((n_test, 5))
y_test = X_test @ true_beta + 0.2 * rng.standard_normal(n_test)

pipe = Pipeline([
    ("scale", StandardScaler()),
    ("lmrob", LmRob(control=Control(nResample=500))),
])
pipe.fit(X_train, y_train, lmrob__seed=42)
print(f"OLS R² on clean test: {pipe.score(X_test, y_test):.3f}")
```

The `lmrob__seed=42` kwarg propagates via sklearn's
double-underscore step-parameter convention. `LmRob.fit` accepts
`seed` and forwards it to `lmrob()`.

Comparing against `LinearRegression` on the same data:

```python
from sklearn.linear_model import LinearRegression
ols = LinearRegression().fit(X_train, y_train)
print(f"OLS test R²: {ols.score(X_test, y_test):.3f}")
```

Output (deterministic with the seeds above):

```
LmRob test R²: 0.991
OLS test R²:   -1.115
```

The OLS fit chases the 30 contaminated training rows; its
coefficients are nowhere near the true `[1, -0.5, 2, 0, 0.3]`, so
predictions on clean test data are worse than the
constant-mean predictor. `lmrob` downweights the outliers and
recovers the underlying linear model.

Inspecting the fitted estimator:

```python
fit = pipe.named_steps["lmrob"].result_
print(fit.coef_)        # coefficients on the standardised design
print(fit.scale_)
```

`coef_` is on the standardised scale (because `StandardScaler` ran
first). To recover the original-scale coefficients, undo the scaling
with `pipe["scale"].scale_`.

## 2. `cross_val_score` for honest performance estimates

Robust regression's main selling point is "the fit doesn't blow up
on outliers." A held-out CV score quantifies how the fit
generalises:

```python
from sklearn.model_selection import cross_val_score

scores = cross_val_score(
    LmRob(control=Control(nResample=500)),
    X, y,
    cv=5,
    fit_params={"seed": 42},
)
print(f"5-fold R² (mean ± std): {scores.mean():.3f} ± {scores.std():.3f}")
```

For comparison, the same call with `sklearn.linear_model.LinearRegression`
would score much lower on the contaminated data, since the OLS fit
chases the outliers and predicts poorly on clean held-out rows.

## 3. `GridSearchCV` over `Control` settings

`Control` exposes `get_params` / `set_params` so sklearn's nested
parameter syntax (`estimator__nested__field`) works directly:

```python
from sklearn.model_selection import GridSearchCV

grid = GridSearchCV(
    LmRob(),
    param_grid={
        "control__nResample": [200, 500, 1000],
        "control__psi": ["bisquare", "optimal", "lqq"],
    },
    cv=5,
    scoring="r2",
)
grid.fit(X_train, y_train, seed=42)
print("best:", grid.best_params_)
print("best score:", grid.best_score_)
```

`GridSearchCV` does the cross-validation, leaderboard, and refit on
the full training set for you. The fitted best estimator is at
`grid.best_estimator_`.

If you'd rather pass entire `Control` objects (e.g., to test KS2014
vs MM presets), use the older form:

```python
candidates = [
    Control(setting="KS2014"),
    Control(setting="KS2011"),
    Control(),  # default MM
]
grid = GridSearchCV(LmRob(), param_grid={"control": candidates}, cv=5)
```

## When sklearn integration doesn't make sense

- **Reading coefficients on the original scale.** If you need
  named, untransformed coefficients (`fit.summary()` style),
  `LmRob` inside a `Pipeline` complicates things because each
  preprocessor changes the feature space. Skip the pipeline and
  call `pylmrob.lmrob("y ~ a + b + c", df, ...)` directly.
- **`predict` with confidence intervals.** `LmRob.predict`
  returns point predictions only. For
  `predict(interval='confidence')` use `LmRob().fit(...).result_.predict(...)`.
  See the [stackloss tour](stackloss_tour.md) for the API.
- **rng="R" reproducibility.** Inside `Pipeline`,
  `Control(rng="R")` still works, but cross-validation folds shuffle
  rows, so byte-identical agreement with R becomes per-fold not
  per-fit.

## See also

- [`stackloss_tour`](stackloss_tour.md): the canonical worked
  example for `fit.predict(interval=...)`, `fit.diagnostics()`,
  `fit.anova()`.
- [`01_ols_vs_robust`](../notebooks/01_ols_vs_robust.md): the same
  kind of contamination this page uses, with plots and the M /
  Theil-Sen comparison.
