# Research notes

This file is the running log for the literature reading required by
plan.md §1.1. Add notes as you read; do not mark a paper "read" without
a paragraph here.

## Reading status

| # | Paper | Read | Notes section |
|---|-------|------|---------------|
| 1 | [Yohai (1987) "High breakdown-point and high efficiency robust estimates for regression"](https://doi.org/10.1214/aos/1176350366). *Annals of Statistics* 15(2). | TODO | [§1](#1-yohai-1987-mm-estimators) |
| 2 | [Salibian-Barrera & Yohai (2006) "A fast algorithm for S-regression estimates"](https://doi.org/10.1198/106186006X113629). *JCGS* 15(2). | TODO | [§2](#2-salibian-barrera--yohai-2006-fast-s) |
| 3 | [Maronna & Yohai (2000) "Robust regression with both continuous and categorical predictors"](https://doi.org/10.1016/S0378-3758(99)00208-6). *J. Stat. Plan. Inf.* | TODO | [§3](#3-maronna--yohai-2000-m-s-for-factor-designs) |
| 4 | [Koller & Stahel (2011) "Sharpening Wald-type inference in robust regression for small samples"](https://doi.org/10.1016/j.csda.2011.02.014). *CSDA*. | TODO | [§4](#4-koller--stahel-2011-ks2011) |
| 5 | [Koller & Stahel (2017) "Nonsingular subsampling for regression S estimators with categorical predictors"](https://doi.org/10.1007/s00180-016-0679-x). *Comp. Stat.* | TODO | [§5](#5-koller--stahel-2017-ks2014--lqq) |
| 6 | [Koller (2016) `lmrob_simulation.pdf` vignette](https://cran.r-project.org/web/packages/robustbase/vignettes/lmrob_simulation.pdf). | TODO | [§6](#6-koller-2016-lmrob-simulation-vignette) |

For each paper, the notes should answer:

- What problem does it solve?
- What are the algorithm steps in pseudocode?
- What are the tuning constants and how are they derived?
- Where in `robustbase`'s R / C source is it implemented?
- Which Python module(s) in this project port that algorithm?

---

## 1. Yohai (1987): MM-estimators

TODO.

## 2. Salibian-Barrera & Yohai (2006): fast-S

TODO.

## 3. Maronna & Yohai (2000): M-S for factor designs

TODO.

## 4. Koller & Stahel (2011): KS2011

TODO.

## 5. Koller & Stahel (2017): KS2014 + lqq

TODO.

## 6. Koller (2016): `lmrob` simulation vignette

TODO.

---

## Reference benchmark numbers

Once §6 is read, copy the simulation vignette's headline numbers
(efficiency vs breakdown trade-off plots, sample-size vs bias plots)
here so we can sanity-check our own validation runs.

## Glossary

- **MM-estimator**: combines a high-breakdown initial S-estimate of scale
  with a higher-efficiency M-estimate of regression coefficients holding
  that scale fixed.
- **psi function**: derivative of the rho function; controls how outlying
  residuals are downweighted.
- **chi function**: a rho-like function used inside the M-scale equation;
  in robustbase it is the same as rho up to scaling.
- **breakdown point**: max fraction of contamination an estimator can
  tolerate before producing an arbitrarily bad result.
- **efficiency**: ratio of the variance of the OLS estimator to the
  variance of the robust estimator under the Gaussian model. Tuning
  constants are usually chosen to give ~95% efficiency.
