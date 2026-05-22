# What is robust regression?

OLS fits a linear model by minimizing the sum of squared residuals.
That works well when the noise is roughly Gaussian, but it has no
defense against outliers: a single observation far from the trend can
swing the fitted slope arbitrarily far. Robust regression replaces the
squared loss with a function that grows slowly in `|r|` so outliers
contribute less.

This page sketches what `lmrob` actually does, so you can pick sensible
defaults without reading three papers. The four
[notebooks](notebooks/01_ols_vs_robust.md) make each idea concrete with
plots and simulations; this page is the prose companion to them.

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
efficiency and breakdown. Notebook
[01](notebooks/01_ols_vs_robust.md) shows what an M-estimator does to
vertical vs horizontal outliers, side by side with OLS and a
Theil-Sen reference.

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
(around 28% at the 50% breakdown tuning). Notebook
[03](notebooks/03_breakdown.md) measures this empirically by sweeping
the contamination fraction; notebook
[04](notebooks/04_s_estimator.md) opens up the S step and shows why
the resampling delivers a high-breakdown starting point.

## MM-estimator

The **MM-estimator** is the workhorse of `lmrob`. It runs two stages:

1. **S-step.** Find an S-estimate `β_S` and the corresponding `σ_S` at
   50% breakdown. This is the high-breakdown but low-efficiency
   starting point.
2. **MM-step.** Starting from `β_S`, holding `σ_S` fixed, iterate
   IRWLS with a *different* tuning of `ρ`: one that gives 95% Gaussian
   efficiency. The final `β_MM` retains the S-step's 50% breakdown
   point and reaches 95% efficiency under Gaussian noise.

That's the best of both worlds in one fit. The default psi family is
`bisquare` with the 95%-efficient tuning constant `k = 4.685`.
Notebook [02](notebooks/02_efficiency.md) measures the efficiency
claim with a Monte Carlo simulation and contrasts it against an
M-estimator and Theil-Sen.

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

## Where MM fails

A 50% breakdown is the worst-case guarantee, not free insurance.
If the contamination is unusually structured — for example, more
than half the data shifted along a specific direction in design
space, mimicking a real linear trend — MM can lock onto the wrong
trend. Notebook [03](notebooks/03_breakdown.md) shows the cliff:
the estimator stays accurate until ε ≈ 0.5, then collapses
abruptly. Past that point, no single fit can recover the truth;
the right move is to identify and split the populations.

## Which psi family to pick

Six families are exposed: `bisquare` (the default), `welsh`,
`optimal`, `hampel`, `lqq`, `ggw`. **Stick with `bisquare`** unless
you have a specific reason to change. `lqq` is the default with
`setting="KS2014"`. The others are mostly of theoretical interest;
the [FAQ](faq.md) has a quick description if you need to pick.

For implementation details — the Cython kernel, parallel
resampling, BLAS interactions, and bench numbers — see
[`engine_c`](engine_c.md).

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
