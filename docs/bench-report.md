# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 5.74x, min 1.77x, max 16.16x
- Runtime ratio (py engine_c/R): median 1.42x, min 0.78x, max 8.15x

## Environment

- pylmrob: 0.5.24
- Python: 3.12.13
- Platform: Linux-6.17.0-1018-azure-x86_64-with-glibc2.39
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
| classical_aircraft | bisquare | 23x5 | 4.0 | 29.2 | 7.30x | 5.1 | 1.27x |
| classical_coleman | bisquare | 20x6 | 5.0 | 29.1 | 5.81x | 5.0 | 0.99x |
| classical_delivery | bisquare | 25x3 | 4.0 | 28.0 | 6.99x | 4.8 | 1.21x |
| classical_hbk | bisquare | 75x4 | 6.0 | 36.2 | 6.03x | 45.9 | 7.66x |
| classical_pension | bisquare | 18x2 | 3.0 | 24.2 | 8.08x | 3.9 | 1.28x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 26.0 | 8.66x | 3.8 | 1.28x |
| classical_salinity | bisquare | 28x4 | 4.0 | 29.2 | 7.31x | 5.4 | 1.34x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 28.3 | 7.08x | 32.4 | 8.11x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 28.5 | 7.13x | 5.5 | 1.38x |
| classical_wood | bisquare | 20x6 | 4.0 | 28.7 | 7.17x | 4.9 | 1.24x |
| psi_bisquare | bisquare | 21x4 | 4.0 | 27.7 | 6.91x | 32.6 | 8.15x |
| psi_ggw | ggw | 21x4 | 6.0 | 29.5 | 4.91x | 6.9 | 1.14x |
| psi_hampel | hampel | 21x4 | 4.0 | 28.6 | 7.14x | 5.2 | 1.30x |
| psi_lqq | lqq | 21x4 | 4.0 | 31.6 | 7.89x | 5.8 | 1.46x |
| psi_optimal | optimal | 21x4 | 4.0 | 28.5 | 7.13x | 4.4 | 1.10x |
| setting_KS2011_stackloss | lqq | 21x4 | 6.0 | 96.9 | 16.16x | 7.1 | 1.18x |
| setting_KS2014_stackloss | lqq | 21x4 | 9.0 | 97.3 | 10.81x | 7.0 | 0.78x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 394.0 | 842.6 | 2.14x | 487.8 | 1.24x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.0 | 142.4 | 3.56x | 74.7 | 1.87x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 633.0 | 1122.2 | 1.77x | 729.5 | 1.15x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.0 | 168.1 | 2.27x | 138.3 | 1.87x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 452.0 | 1017.4 | 2.25x | 682.3 | 1.51x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.0 | 163.6 | 3.48x | 121.4 | 2.58x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 489.0 | 1304.1 | 2.67x | 669.6 | 1.37x |
| synth_lqq_n500_p10 | lqq | 500x11 | 53.0 | 245.8 | 4.64x | 119.6 | 2.26x |
| synth_n10000_p20 | bisquare | 10000x21 | 503.0 | 3665.6 | 7.29x | 2937.0 | 5.84x |
| synth_n10000_p50 | bisquare | 10000x51 | 1896.0 | 8064.2 | 4.25x | 7395.9 | 3.90x |
| synth_n1000_p10 | bisquare | 1000x11 | 120.0 | 394.1 | 3.28x | 205.6 | 1.71x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 45.3 | 5.66x | 13.9 | 1.74x |
| synth_n2000_p20 | bisquare | 2000x21 | 396.0 | 757.3 | 1.91x | 580.7 | 1.47x |
| synth_n5000_p20 | bisquare | 5000x21 | 388.0 | 1991.3 | 5.13x | 1457.7 | 3.76x |
| synth_n500_p10 | bisquare | 500x11 | 41.0 | 143.1 | 3.49x | 74.9 | 1.83x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 451.0 | 957.7 | 2.12x | 553.1 | 1.23x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.0 | 196.9 | 4.28x | 90.8 | 1.97x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
