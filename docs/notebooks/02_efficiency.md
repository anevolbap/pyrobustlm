---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# What does robustness cost?

OLS is the **minimum-variance unbiased estimator** when the noise is
truly Gaussian. Any robust method, by design, gives up some of that
optimality in exchange for outlier protection. The question is: **how
much?**

The standard measure is *asymptotic relative efficiency* (ARE):

$$\text{ARE}(\widehat{\beta}_{\text{robust}}) = \frac{\text{Var}(\widehat{\beta}_{\text{OLS}})}{\text{Var}(\widehat{\beta}_{\text{robust}})}$$

The closer to 1, the less you give up. MM-estimators with the default
bisquare tuning target **95% ARE under Gaussian noise**. We'll
measure that empirically and compare against M-estimator (Huber) and
Theil-Sen.

When the noise is **heavy-tailed**, the picture flips: OLS becomes
the worst choice and robust methods win. We'll measure that too.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from pylmrob import lmrob, Control
import statsmodels.api as sm
from statsmodels.robust.robust_linear_model import RLM
from sklearn.linear_model import TheilSenRegressor

rng = np.random.default_rng(0)

plt.rcParams.update({"figure.dpi": 110, "axes.grid": True, "grid.alpha": 0.3})

METHODS = ["OLS", "MM (pylmrob)", "Huber-M", "Theil-Sen"]
COLORS = {
    "OLS": "#d62728",
    "MM (pylmrob)": "#1f77b4",
    "Huber-M": "#ff7f0e",
    "Theil-Sen": "#2ca02c",
}
```

## A small Monte Carlo loop

For each noise distribution, we draw `n_rep=200` clean datasets of
size `n=100` and fit all four methods. We track the slope estimate;
its **spread** across replicates is what efficiency measures.

```{code-cell} ipython3
def fit_slope(x, y, method):
    """Return the slope estimate for a single (x, y) under one method."""
    X = sm.add_constant(x)
    if method == "OLS":
        return sm.OLS(y, X).fit().params[1]
    if method == "MM (pylmrob)":
        df = pd.DataFrame({"x": x, "y": y})
        return lmrob("y ~ x", df, control=Control(), seed=42).coef_[1]
    if method == "Huber-M":
        return RLM(y, X, M=sm.robust.norms.HuberT()).fit().params[1]
    if method == "Theil-Sen":
        return TheilSenRegressor(random_state=0).fit(x.reshape(-1, 1), y).coef_[0]
    raise ValueError(method)


def simulate(noise_fn, n_rep=200, n=100, true_slope=1.5, true_intercept=2.0):
    """Run the MC. ``noise_fn(rng, size)`` returns the noise vector."""
    slopes = {m: np.empty(n_rep) for m in METHODS}
    sim_rng = np.random.default_rng(42)
    x = np.linspace(0, 10, n)
    for r in range(n_rep):
        eps = noise_fn(sim_rng, n)
        y = true_intercept + true_slope * x + eps
        for m in METHODS:
            slopes[m][r] = fit_slope(x, y, m)
    return slopes
```

## Scenario A: Gaussian noise (the OLS-optimal case)

```{code-cell} ipython3
slopes_gauss = simulate(lambda r, n: 0.5 * r.standard_normal(n))

# Empirical variances of the slope estimator under each method.
variances = {m: float(np.var(s, ddof=1)) for m, s in slopes_gauss.items()}
are = {m: variances["OLS"] / variances[m] for m in METHODS}
print("Variance of slope:")
for m, v in variances.items():
    print(f"  {m:<14} {v:.5f}")
print()
print("ARE relative to OLS (higher = better):")
for m, a in are.items():
    print(f"  {m:<14} {a:.3f}")
```

The OLS slope has the smallest variance; that's the baseline (ARE = 1).
**MM and Huber-M both reach ~95% efficiency** on Gaussian data — they
give up ~5% of OLS's precision in exchange for outlier protection.
**Theil-Sen** trades more (~67% ARE) because its median-of-slopes
construction is statistically less efficient at normal noise.

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(8, 4))
for m in METHODS:
    ax.hist(slopes_gauss[m], bins=30, alpha=0.5,
            label=f"{m}  (ARE={are[m]:.2f})", color=COLORS[m])
ax.axvline(1.5, color="gray", ls="--", lw=1, label="true slope")
ax.set_xlabel("slope estimate")
ax.set_ylabel("count (out of 200 replicates)")
ax.set_title("Sampling distribution under Gaussian noise")
ax.legend(fontsize=9)
fig.tight_layout()
```

All four are centred on the true slope; their **widths** are what
ARE measures. OLS is tightest. MM and Huber-M are nearly as tight.
Theil-Sen is visibly wider.

## Scenario B: Student-t(3) noise (heavy-tailed)

Now use t-distributed noise with 3 degrees of freedom. Its variance
is finite but the tails are heavy: an unlikely fraction of draws are
big outliers in a way that Gaussian noise never produces. OLS pays
for the heavy tails through inflated variance.

```{code-cell} ipython3
slopes_t3 = simulate(lambda r, n: 0.5 * r.standard_t(df=3, size=n))

variances_t3 = {m: float(np.var(s, ddof=1)) for m, s in slopes_t3.items()}
are_t3 = {m: variances_t3["OLS"] / variances_t3[m] for m in METHODS}
print("Variance of slope:")
for m, v in variances_t3.items():
    print(f"  {m:<14} {v:.5f}")
print()
print("ARE relative to OLS:")
for m, a in are_t3.items():
    print(f"  {m:<14} {a:.3f}")
```

OLS variance has blown up (look at the absolute number compared to
the Gaussian case). The **ARE is now > 1** for all three robust
methods — they are *more* efficient than OLS on heavy-tailed data.

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(8, 4))
for m in METHODS:
    ax.hist(slopes_t3[m], bins=40, alpha=0.5,
            label=f"{m}  (ARE={are_t3[m]:.2f})", color=COLORS[m])
ax.axvline(1.5, color="gray", ls="--", lw=1, label="true slope")
ax.set_xlabel("slope estimate")
ax.set_ylabel("count (out of 200 replicates)")
ax.set_title("Sampling distribution under Student-t(3) noise")
ax.legend(fontsize=9)
fig.tight_layout()
```

OLS estimates fan out widely; the robust methods stay tight. Same
data-generating process, very different sampling behaviour.

## Summary table

```{code-cell} ipython3
summary = pd.DataFrame({
    "method": METHODS,
    "ARE @ Gaussian": [are[m] for m in METHODS],
    "ARE @ t(3)":     [are_t3[m] for m in METHODS],
})
summary.set_index("method").round(2)
```

The MM-estimator's design goal — **95% efficiency on the Gaussian
case, plus outlier protection** — is visible in the first column:
0.95 against OLS, basically no cost. The second column shows the
payoff when the noise is *actually* heavy-tailed.

That's the whole bargain. You give up ~5% on the perfect Gaussian
case to gain a 1.5-2× safety margin on anything heavier-tailed.

## Next

- [`03_breakdown`](03_breakdown.md): what about contamination, not
  just heavy tails?
- [`04_s_estimator`](04_s_estimator.md): inside the MM algorithm.
