# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 5.47x, min 1.34x, max 16.13x
- Runtime ratio (py engine_c/R): median 1.51x, min 0.72x, max 8.36x

## Environment

- pylmrob: 0.5.24
- Python: 3.12.13
- Platform: Linux-6.17.0-1018-azure-x86_64-with-glibc2.39
- robustbase: 0.99.7
- R: R version 4.6.0 (2026-04-24)

## Numerical accuracy: max relative error vs R

| case | psi | n_x_p | max coef rerr | scale rerr | cov diag max rerr |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.03e-06 | 3.35e-06 | 7.76e-06 |
| classical_coleman | bisquare | 20x6 | 7.97e-07 | 3.69e-06 | 2.31e-06 |
| classical_delivery | bisquare | 25x3 | 6.36e-07 | 1.53e-05 | 2.76e-05 |
| classical_hbk | bisquare | 75x4 | 6.48e-08 | 3.23e-06 | 1.08e-06 |
| classical_pension | bisquare | 18x2 | 3.43e-06 | 5.76e-06 | 2.27e-05 |
| classical_phosphor | bisquare | 18x3 | 5.36e-07 | 3.66e-06 | 2.87e-06 |
| classical_salinity | bisquare | 28x4 | 2.19e-07 | 4.07e-06 | 3.11e-06 |
| classical_stackloss | bisquare | 21x4 | 1.85e-06 | 4.67e-06 | 1.33e-05 |
| classical_starsCYG | bisquare | 47x2 | 1.78e-06 | 3.23e-06 | 5.97e-06 |
| classical_wood | bisquare | 20x6 | 3.10e-08 | 3.44e-06 | 7.53e-07 |
| psi_bisquare | bisquare | 21x4 | 1.85e-06 | 4.67e-06 | 1.33e-05 |
| psi_ggw | ggw | 21x4 | 1.33e-03 | 2.51e-03 | 5.20e-01 |
| psi_hampel | hampel | 21x4 | 1.98e-03 | 3.45e-03 | 4.01e-01 |
| psi_lqq | lqq | 21x4 | 1.86e-06 | 3.69e-06 | 7.91e-06 |
| psi_optimal | optimal | 21x4 | 1.44e-15 | 5.80e-07 | 1.23e-13 |
| setting_KS2011_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| setting_KS2014_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.48e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.66e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.48e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.66e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 3.24e-08 | 3.24e-06 | 4.66e-07 |
| synth_n10000_p50 | bisquare | 10000x51 | 2.52e-08 | 3.23e-06 | 5.10e-07 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.10e-08 | 3.23e-06 | 7.74e-07 |
| synth_n100_p5 | bisquare | 100x6 | 1.29e-07 | 3.23e-06 | 1.39e-06 |
| synth_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_n5000_p20 | bisquare | 5000x21 | 4.39e-08 | 3.24e-06 | 4.73e-07 |
| synth_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.35e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.44e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 5 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 4.0 | 29.5 | 7.39x | 5.2 | 1.30x |
| classical_coleman | bisquare | 20x6 | 5.0 | 29.5 | 5.89x | 5.1 | 1.02x |
| classical_delivery | bisquare | 25x3 | 3.0 | 28.3 | 9.44x | 4.9 | 1.63x |
| classical_hbk | bisquare | 75x4 | 6.0 | 37.1 | 6.18x | 46.5 | 7.75x |
| classical_pension | bisquare | 18x2 | 3.0 | 25.3 | 8.44x | 3.9 | 1.31x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 26.2 | 8.74x | 3.9 | 1.30x |
| classical_salinity | bisquare | 28x4 | 4.0 | 29.7 | 7.42x | 5.4 | 1.34x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 28.6 | 7.15x | 33.5 | 8.36x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 29.3 | 7.31x | 5.6 | 1.39x |
| classical_wood | bisquare | 20x6 | 5.0 | 29.3 | 5.85x | 5.1 | 1.02x |
| psi_bisquare | bisquare | 21x4 | 4.0 | 28.6 | 7.16x | 33.2 | 8.30x |
| psi_ggw | ggw | 21x4 | 6.0 | 30.2 | 5.03x | 7.0 | 1.17x |
| psi_hampel | hampel | 21x4 | 4.0 | 29.0 | 7.24x | 5.3 | 1.33x |
| psi_lqq | lqq | 21x4 | 4.0 | 32.1 | 8.03x | 6.0 | 1.49x |
| psi_optimal | optimal | 21x4 | 4.0 | 29.2 | 7.29x | 4.5 | 1.13x |
| setting_KS2011_stackloss | lqq | 21x4 | 6.0 | 96.8 | 16.13x | 7.2 | 1.20x |
| setting_KS2014_stackloss | lqq | 21x4 | 10.0 | 96.9 | 9.69x | 7.2 | 0.72x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 396.0 | 846.3 | 2.14x | 589.1 | 1.49x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.0 | 143.8 | 3.60x | 74.8 | 1.87x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 634.0 | 851.8 | 1.34x | 912.0 | 1.44x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.0 | 170.4 | 2.30x | 138.3 | 1.87x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 447.0 | 841.0 | 1.88x | 822.1 | 1.84x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.0 | 162.5 | 3.46x | 120.4 | 2.56x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 482.0 | 1144.0 | 2.37x | 803.3 | 1.67x |
| synth_lqq_n500_p10 | lqq | 500x11 | 53.0 | 247.7 | 4.67x | 119.0 | 2.25x |
| synth_n10000_p20 | bisquare | 10000x21 | 499.0 | 3898.9 | 7.81x | 2961.8 | 5.94x |
| synth_n10000_p50 | bisquare | 10000x51 | 1880.0 | 8322.1 | 4.43x | 7397.3 | 3.93x |
| synth_n1000_p10 | bisquare | 1000x11 | 122.0 | 333.9 | 2.74x | 206.6 | 1.69x |
| synth_n100_p5 | bisquare | 100x6 | 9.0 | 45.7 | 5.08x | 13.8 | 1.53x |
| synth_n2000_p20 | bisquare | 2000x21 | 400.0 | 829.2 | 2.07x | 591.4 | 1.48x |
| synth_n5000_p20 | bisquare | 5000x21 | 390.0 | 1975.5 | 5.07x | 1539.1 | 3.95x |
| synth_n500_p10 | bisquare | 500x11 | 40.0 | 143.6 | 3.59x | 75.1 | 1.88x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 452.0 | 923.9 | 2.04x | 641.3 | 1.42x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.0 | 193.2 | 4.20x | 90.6 | 1.97x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
