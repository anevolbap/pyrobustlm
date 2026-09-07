# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.38x, min 1.42x, max 9.56x
- Runtime ratio (py engine_c/R): median 1.52x, min 0.82x, max 5.53x

## Environment

- pylmrob: 0.5.31
- Python: 3.12.14
- Platform: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39
- robustbase: 0.99.7
- R: R version 4.6.1 (2026-06-24)

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
| psi_optimal | optimal | 21x4 | 1.44e-15 | 5.80e-07 | 9.59e-14 |
| setting_KS2011_stackloss | lqq | 21x4 | 7.42e-08 | 8.73e-08 | 3.04e-07 |
| setting_KS2014_stackloss | lqq | 21x4 | 7.31e-08 | 8.61e-08 | 3.00e-07 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 3.50e-13 | 1.78e-11 | 3.74e-08 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.48e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.67e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.64e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 5.85e-11 | 5.81e-09 | 1.98e-08 |
| synth_n10000_p50 | bisquare | 10000x51 | 1.12e-11 | 1.39e-09 | 2.22e-08 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.72e-14 | 3.92e-12 | 3.39e-09 |
| synth_n100_p5 | bisquare | 100x6 | 2.14e-12 | 5.35e-11 | 9.68e-09 |
| synth_n2000_p20 | bisquare | 2000x21 | 3.50e-13 | 1.78e-11 | 3.74e-08 |
| synth_n5000_p20 | bisquare | 5000x21 | 1.13e-10 | 8.35e-09 | 9.11e-09 |
| synth_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.34e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.44e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 11 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.6 | 26.7 | 7.33x | 5.0 | 1.37x |
| classical_coleman | bisquare | 20x6 | 3.9 | 26.9 | 6.90x | 4.9 | 1.25x |
| classical_delivery | bisquare | 25x3 | 3.0 | 25.0 | 8.31x | 4.8 | 1.59x |
| classical_hbk | bisquare | 75x4 | 7.9 | 30.0 | 3.82x | 8.9 | 1.14x |
| classical_pension | bisquare | 18x2 | 2.3 | 21.9 | 9.56x | 3.8 | 1.67x |
| classical_phosphor | bisquare | 18x3 | 2.7 | 23.8 | 8.68x | 3.9 | 1.41x |
| classical_salinity | bisquare | 28x4 | 3.5 | 26.2 | 7.41x | 5.3 | 1.50x |
| classical_stackloss | bisquare | 21x4 | 3.2 | 26.1 | 8.23x | 4.3 | 1.37x |
| classical_starsCYG | bisquare | 47x2 | 3.3 | 24.4 | 7.39x | 5.3 | 1.60x |
| classical_wood | bisquare | 20x6 | 3.9 | 26.7 | 6.81x | 4.8 | 1.23x |
| psi_bisquare | bisquare | 21x4 | 3.2 | 26.0 | 8.25x | 4.2 | 1.33x |
| psi_ggw | ggw | 21x4 | 5.7 | 29.1 | 5.15x | 6.7 | 1.19x |
| psi_hampel | hampel | 21x4 | 3.4 | 27.7 | 8.21x | 5.1 | 1.52x |
| psi_lqq | lqq | 21x4 | 3.9 | 27.9 | 7.21x | 5.7 | 1.48x |
| psi_optimal | optimal | 21x4 | 3.2 | 25.5 | 7.85x | 4.4 | 1.34x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.1 | 29.3 | 5.74x | 7.0 | 1.36x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.5 | 29.6 | 3.48x | 7.0 | 0.82x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 393.5 | 695.0 | 1.77x | 577.6 | 1.47x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 39.6 | 103.4 | 2.61x | 72.1 | 1.82x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 615.3 | 873.3 | 1.42x | 905.5 | 1.47x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.0 | 166.6 | 2.25x | 138.0 | 1.86x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 440.2 | 993.7 | 2.26x | 812.7 | 1.85x |
| synth_hampel_n500_p10 | hampel | 500x11 | 46.9 | 151.5 | 3.23x | 121.1 | 2.58x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 473.5 | 925.7 | 1.96x | 832.2 | 1.76x |
| synth_lqq_n500_p10 | lqq | 500x11 | 52.4 | 150.4 | 2.87x | 124.2 | 2.37x |
| synth_n10000_p20 | bisquare | 10000x21 | 521.1 | 3234.8 | 6.21x | 2883.3 | 5.53x |
| synth_n10000_p50 | bisquare | 10000x51 | 1905.8 | 7495.8 | 3.93x | 7243.0 | 3.80x |
| synth_n1000_p10 | bisquare | 1000x11 | 121.0 | 270.0 | 2.23x | 202.9 | 1.68x |
| synth_n100_p5 | bisquare | 100x6 | 7.9 | 36.4 | 4.61x | 13.4 | 1.70x |
| synth_n2000_p20 | bisquare | 2000x21 | 394.5 | 704.7 | 1.79x | 579.0 | 1.47x |
| synth_n5000_p20 | bisquare | 5000x21 | 403.2 | 1677.4 | 4.16x | 1427.1 | 3.54x |
| synth_n500_p10 | bisquare | 500x11 | 40.4 | 103.2 | 2.55x | 72.2 | 1.79x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 445.1 | 766.0 | 1.72x | 674.7 | 1.52x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.1 | 121.8 | 2.64x | 93.7 | 2.03x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
