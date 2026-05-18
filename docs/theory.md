# What is robust regression?

OLS fits a linear model by minimizing the sum of squared residuals.
That works well when the noise is roughly Gaussian, but it has no
defense against outliers: a single observation far from the trend can
swing the fitted slope arbitrarily far. Robust regression replaces the
squared loss with a function that grows slowly in `|r|` so outliers
contribute less.

This page sketches what `lmrob` actually does, so you can pick sensible
defaults without reading three papers.

## The M-estimator

An **M-estimator** minimises

$$\widehat{\beta} = \arg\min_{\beta} \sum_{i=1}^{n} \rho\!\left(\frac{y_i - x_i^\top \beta}{\sigma}\right)$$

for some loss function `ρ` (rho) and a scale `σ`. The first-order
condition gives a weighted least-squares problem with weights
`ψ(r/σ) / (r/σ)`, where `ψ = ρ'`. Iterating WLS to convergence is the
IRWLS algorithm. `lmrob` and `pylmrob` use IRWLS at every step.

For a sensible robust fit you want `ρ` to be:

- **bounded** (so far-away points have bounded influence), and
- **smooth and zero near the origin** (so well-behaved points keep
  unit weight).

The classical `bisquare` family does both:
`ρ(u) = 1 - (1 - (u/k)^2)^3` for `|u| < k`, and `ρ(u) = 1` outside.
The tuning constant `k` controls the trade-off between Gaussian
efficiency and breakdown.

## S-estimator and M-scale

The M-estimator above needs a value of `σ` (the scale of the
residuals). You can't just plug in `sd(r)`, because the sd is itself
ruined by outliers. Instead `lmrob` uses an **M-scale**: find `σ` such
that

$$\frac{1}{n - p}\sum_{i=1}^{n} \chi\!\left(\frac{r_i}{\sigma}\right) = b_0$$

for a tuning function `χ` and a target `b₀ = 0.5`. The S-estimator
*jointly* finds `(β, σ)` by minimising `σ` over candidate `β`'s, each
evaluated by the M-scale equation above.

S-estimators have a **breakdown point** of 50%: you can replace up to
half the data with arbitrary outliers and the estimator still tracks
the bulk of the data. The down-side is low Gaussian efficiency
(around 28% at the 50% breakdown tuning).

## MM-estimator

The **MM-estimator** is the workhorse of `lmrob`. It runs two stages:

1. **S-step.** Find an S-estimate `β_S` and the corresponding `σ_S` at
   50% breakdown. This is the high-breakdown but low-efficiency
   starting point.
2. **MM-step.** Starting from `β_S`, holding `σ_S` fixed, iterate
   IRWLS with a *different* tuning of `ρ`: one that gives 95% Gaussian
   efficiency. The final `β_MM` retains the S-step's 50% breakdown
   point and reaches 95% efficiency under Gaussian noise.

That's the best of both worlds in one fit. The Cython kernel runs the
S-step (with the chi tuning) and the MM IRWLS (with the psi tuning)
in a single `nogil` C block. The default psi family is `bisquare`
with the 95%-efficient tuning constant `k = 4.685`.

## Why `setting="KS2014"` is the recommended default

The basic MM-estimator is Yohai (1987). Koller & Stahel (2014)
showed that the M-scale estimate of `σ_MM` can be biased when the
contamination is heavy, and proposed a **D-scale refinement**
(robustbase calls this `setting="KS2014"`) that re-estimates `σ` using
the MM weights with a design correction. Empirically it gives more
honest standard errors when the data really is contaminated.

For pure Gaussian data the MM-estimator and the SMDM pipeline give
essentially the same answer; for contaminated data the SMDM pipeline
gives slightly better calibrated inference. `setting="KS2014"` is what
the original authors recommend; `pylmrob` follows R in keeping plain
MM as the default for backwards compatibility, but exposes
`Control(setting="KS2014")` as a one-keyword switch.

## Which psi family to pick

All five families are exposed: `bisquare` (the default), `optimal`,
`hampel`, `lqq`, `ggw`. Practical guidance:

- **Stick with `bisquare`** unless you have a specific reason to
  change. It's the original Tukey psi, smooth, easy to reason about,
  and what most papers benchmark with.
- **`lqq`** (linear-quadratic-quadratic) is the default with
  `setting="KS2014"`. Slightly heavier downweighting at the bend than
  bisquare; the Koller & Stahel paper argues this is preferable for
  the D-scale refinement.
- **`hampel`** is piecewise-linear with three breakpoints; very
  classical, occasionally numerically nicer on heavy-tail
  contamination.
- **`optimal`** is V. Yohai's optimal redescending psi for MM. Mostly
  of theoretical interest; performance is close to bisquare.
- **`ggw`** is the Generalised Gauss-Welsh family; a smooth psi with
  six tabulated cases (`tuning=(case_idx,)`). Heaviest downweighting
  of the redescending families.

If you don't know which to pick: use the default. If your data has
heavy contamination, try `setting="KS2014"`.

## Computational cost

`lmrob` is dominated by the S-step: 500 candidate p-subsets, each
refined for a few IRWLS steps, then survivor refinement on the best
few. With `engine_c=True` (the default since v0.5.11) the whole
pipeline runs in one nogil C block. With `n_workers > 1` the
candidate loop runs in parallel via OpenMP.

The wall-clock floor is the underlying BLAS work, which `pylmrob`
shares with R's `lmrob` (both call OpenBLAS under the hood). Median
wall-clock across the bench corpus is **0.93x R**.

## Further reading

- Maronna, Martin, Yohai, Salibian-Barrera, *Robust Statistics: Theory
  and Methods (with R)*, 2nd ed., Wiley 2019. The standard reference.
- Koller & Stahel, "Sharpening Wald-type inference in robust regression
  for small samples", *Computational Statistics & Data Analysis* 55
  (2011): 2504-2515. The KS2011 setting.
- Koller & Stahel, "Nonsingular subsampling for regression S
  estimators with categorical covariates", *Computational Statistics*
  32 (2017): 631-646. The M-S estimator for designs with factor
  variables.
- Yohai (1987), "High breakdown-point and high efficiency robust
  estimates for regression", *Annals of Statistics* 15: 642-656. The
  original MM-estimator paper.
- Salibian-Barrera & Yohai (2006), "A fast algorithm for S-regression
  estimates", *Journal of Computational and Graphical Statistics* 15:
  414-427. The fast-S resampling algorithm `lmrob` is built on.
