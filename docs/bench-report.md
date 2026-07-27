# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 5.04x, min 1.64x, max 16.18x
- Runtime ratio (py engine_c/R): median 1.65x, min 0.85x, max 7.75x

## Environment

- pylmrob: 0.5.24
- Python: 3.12.13
- Platform: Linux-6.17.0-1020-azure-x86_64-with-glibc2.39
- robustbase: 0.99.7
- R: R version 4.6.1 (2026-06-24)

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
| classical_aircraft | bisquare | 23x5 | 4.0 | 21.1 | 5.28x | 5.1 | 1.28x |
| classical_coleman | bisquare | 20x6 | 4.0 | 21.1 | 5.28x | 5.0 | 1.25x |
| classical_delivery | bisquare | 25x3 | 3.0 | 20.6 | 6.85x | 4.9 | 1.65x |
| classical_hbk | bisquare | 75x4 | 5.0 | 29.9 | 5.98x | 38.7 | 7.75x |
| classical_pension | bisquare | 18x2 | 2.0 | 17.9 | 8.96x | 3.9 | 1.95x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 18.4 | 6.12x | 3.9 | 1.31x |
| classical_salinity | bisquare | 28x4 | 4.0 | 21.8 | 5.45x | 5.4 | 1.36x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 20.1 | 5.01x | 24.2 | 6.06x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 21.9 | 7.29x | 5.4 | 1.80x |
| classical_wood | bisquare | 20x6 | 4.0 | 20.8 | 5.21x | 4.9 | 1.22x |
| psi_bisquare | bisquare | 21x4 | 4.0 | 20.2 | 5.06x | 24.2 | 6.05x |
| psi_ggw | ggw | 21x4 | 6.0 | 21.4 | 3.56x | 6.9 | 1.15x |
| psi_hampel | hampel | 21x4 | 4.0 | 20.3 | 5.08x | 5.3 | 1.31x |
| psi_lqq | lqq | 21x4 | 4.0 | 24.2 | 6.05x | 5.9 | 1.47x |
| psi_optimal | optimal | 21x4 | 3.0 | 21.3 | 7.11x | 4.5 | 1.49x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 80.9 | 16.18x | 6.8 | 1.36x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.0 | 81.4 | 10.18x | 6.8 | 0.85x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 393.0 | 912.1 | 2.32x | 610.7 | 1.55x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 43.0 | 145.1 | 3.37x | 81.2 | 1.89x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 593.0 | 972.0 | 1.64x | 933.3 | 1.57x |
| synth_ggw_n500_p10 | ggw | 500x11 | 78.0 | 165.3 | 2.12x | 143.8 | 1.84x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 432.0 | 1099.0 | 2.54x | 811.9 | 1.88x |
| synth_hampel_n500_p10 | hampel | 500x11 | 51.0 | 166.2 | 3.26x | 129.2 | 2.53x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 472.0 | 1335.3 | 2.83x | 835.6 | 1.77x |
| synth_lqq_n500_p10 | lqq | 500x11 | 57.0 | 263.3 | 4.62x | 127.0 | 2.23x |
| synth_n10000_p20 | bisquare | 10000x21 | 492.0 | 3893.0 | 7.91x | 2944.1 | 5.98x |
| synth_n10000_p50 | bisquare | 10000x51 | 1904.0 | 8694.6 | 4.57x | 7776.1 | 4.08x |
| synth_n1000_p10 | bisquare | 1000x11 | 120.0 | 346.0 | 2.88x | 199.0 | 1.66x |
| synth_n100_p5 | bisquare | 100x6 | 9.0 | 38.6 | 4.29x | 14.3 | 1.59x |
| synth_n2000_p20 | bisquare | 2000x21 | 390.0 | 920.9 | 2.36x | 611.7 | 1.57x |
| synth_n5000_p20 | bisquare | 5000x21 | 394.0 | 2132.0 | 5.41x | 1541.5 | 3.91x |
| synth_n500_p10 | bisquare | 500x11 | 43.0 | 145.2 | 3.38x | 81.0 | 1.88x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 447.0 | 1103.0 | 2.47x | 669.7 | 1.50x |
| synth_optimal_n500_p10 | optimal | 500x11 | 49.0 | 203.9 | 4.16x | 97.6 | 1.99x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
