# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 3.49x, min 1.49x, max 6.35x
- Runtime ratio (py engine_c/R): median 1.41x, min 0.77x, max 4.97x

## Environment

- pylmrob: 0.5.33
- Python: 3.12.16
- Platform: Linux-7.0.0-1019-azure-x86_64-with-glibc2.43
- robustbase: 0.99.7
- R: R version 4.6.2 (2026-10-29)

## Numerical accuracy: max relative error vs R

| case | psi | n_x_p | max coef rerr | scale rerr | cov diag max rerr |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 1.04e-07 | 1.15e-07 | 2.78e-07 |
| classical_coleman | bisquare | 20x6 | 9.88e-08 | 4.57e-07 | 2.98e-07 |
| classical_delivery | bisquare | 25x3 | 5.02e-07 | 1.20e-05 | 2.12e-05 |
| classical_hbk | bisquare | 75x4 | 5.03e-11 | 2.51e-09 | 1.29e-08 |
| classical_pension | bisquare | 18x2 | 1.50e-06 | 2.53e-06 | 1.00e-05 |
| classical_phosphor | bisquare | 18x3 | 6.27e-08 | 4.29e-07 | 3.86e-07 |
| classical_salinity | bisquare | 28x4 | 4.50e-08 | 8.36e-07 | 6.46e-07 |
| classical_stackloss | bisquare | 21x4 | 5.68e-07 | 1.43e-06 | 4.16e-06 |
| classical_starsCYG | bisquare | 47x2 | 1.69e-12 | 4.21e-12 | 3.84e-07 |
| classical_wood | bisquare | 20x6 | 1.91e-09 | 2.12e-07 | 6.92e-08 |
| psi_bisquare | bisquare | 21x4 | 5.68e-07 | 1.43e-06 | 4.16e-06 |
| psi_ggw | ggw | 21x4 | 1.33e-03 | 2.51e-03 | 5.20e-01 |
| psi_hampel | hampel | 21x4 | 1.98e-03 | 3.45e-03 | 4.01e-01 |
| psi_lqq | lqq | 21x4 | 1.86e-06 | 3.69e-06 | 7.91e-06 |
| psi_optimal | optimal | 21x4 | 5.66e-16 | 5.80e-07 | 8.87e-14 |
| setting_KS2011_stackloss | lqq | 21x4 | 7.42e-08 | 8.73e-08 | 3.04e-07 |
| setting_KS2014_stackloss | lqq | 21x4 | 7.31e-08 | 8.61e-08 | 3.00e-07 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 3.48e-13 | 1.78e-11 | 3.74e-08 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.49e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.66e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.64e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 5.85e-11 | 5.81e-09 | 1.98e-08 |
| synth_n10000_p50 | bisquare | 10000x51 | 1.11e-11 | 1.39e-09 | 2.22e-08 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.94e-14 | 3.92e-12 | 3.39e-09 |
| synth_n100_p5 | bisquare | 100x6 | 2.14e-12 | 5.35e-11 | 9.68e-09 |
| synth_n2000_p20 | bisquare | 2000x21 | 3.48e-13 | 1.78e-11 | 3.74e-08 |
| synth_n5000_p20 | bisquare | 5000x21 | 1.13e-10 | 8.35e-09 | 9.11e-09 |
| synth_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.35e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.45e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 11 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.0 | 13.7 | 4.51x | 3.9 | 1.29x |
| classical_coleman | bisquare | 20x6 | 3.2 | 13.9 | 4.29x | 3.8 | 1.16x |
| classical_delivery | bisquare | 25x3 | 2.4 | 13.0 | 5.33x | 3.8 | 1.58x |
| classical_hbk | bisquare | 75x4 | 4.3 | 16.7 | 3.85x | 7.2 | 1.66x |
| classical_pension | bisquare | 18x2 | 1.8 | 11.5 | 6.35x | 3.2 | 1.74x |
| classical_phosphor | bisquare | 18x3 | 2.2 | 12.1 | 5.51x | 3.1 | 1.40x |
| classical_salinity | bisquare | 28x4 | 2.9 | 13.7 | 4.67x | 4.1 | 1.41x |
| classical_stackloss | bisquare | 21x4 | 2.6 | 13.3 | 5.16x | 3.3 | 1.28x |
| classical_starsCYG | bisquare | 47x2 | 2.7 | 13.1 | 4.84x | 4.3 | 1.60x |
| classical_wood | bisquare | 20x6 | 3.3 | 13.9 | 4.17x | 3.7 | 1.12x |
| psi_bisquare | bisquare | 21x4 | 2.6 | 13.1 | 5.02x | 3.3 | 1.26x |
| psi_ggw | ggw | 21x4 | 4.6 | 15.8 | 3.43x | 5.4 | 1.18x |
| psi_hampel | hampel | 21x4 | 2.9 | 14.1 | 4.80x | 4.1 | 1.39x |
| psi_lqq | lqq | 21x4 | 3.3 | 14.6 | 4.48x | 4.6 | 1.40x |
| psi_optimal | optimal | 21x4 | 2.7 | 12.9 | 4.69x | 3.4 | 1.22x |
| setting_KS2011_stackloss | lqq | 21x4 | 4.4 | 15.5 | 3.56x | 5.4 | 1.23x |
| setting_KS2014_stackloss | lqq | 21x4 | 7.0 | 15.5 | 2.23x | 5.4 | 0.77x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 342.1 | 542.0 | 1.58x | 442.8 | 1.29x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 35.2 | 71.2 | 2.03x | 56.4 | 1.60x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 483.7 | 817.7 | 1.69x | 668.2 | 1.38x |
| synth_ggw_n500_p10 | ggw | 500x11 | 66.2 | 129.9 | 1.96x | 113.1 | 1.71x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 378.5 | 751.9 | 1.99x | 595.7 | 1.57x |
| synth_hampel_n500_p10 | hampel | 500x11 | 42.0 | 105.6 | 2.51x | 89.7 | 2.13x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 410.2 | 788.5 | 1.92x | 590.8 | 1.44x |
| synth_lqq_n500_p10 | lqq | 500x11 | 46.7 | 108.9 | 2.33x | 93.1 | 1.99x |
| synth_n10000_p20 | bisquare | 10000x21 | 443.2 | 2425.0 | 5.47x | 2201.0 | 4.97x |
| synth_n10000_p50 | bisquare | 10000x51 | 1738.0 | 5612.8 | 3.23x | 5431.1 | 3.12x |
| synth_n1000_p10 | bisquare | 1000x11 | 88.2 | 215.0 | 2.44x | 113.6 | 1.29x |
| synth_n100_p5 | bisquare | 100x6 | 7.3 | 21.0 | 2.88x | 10.8 | 1.48x |
| synth_n2000_p20 | bisquare | 2000x21 | 343.1 | 544.0 | 1.59x | 442.1 | 1.29x |
| synth_n5000_p20 | bisquare | 5000x21 | 350.5 | 1363.4 | 3.89x | 1226.0 | 3.50x |
| synth_n500_p10 | bisquare | 500x11 | 35.7 | 71.0 | 1.99x | 56.5 | 1.58x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 380.7 | 566.0 | 1.49x | 475.7 | 1.25x |
| synth_optimal_n500_p10 | optimal | 500x11 | 40.6 | 79.3 | 1.96x | 64.0 | 1.58x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
