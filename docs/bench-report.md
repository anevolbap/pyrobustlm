# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 3.81x, min 1.81x, max 9.99x
- Runtime ratio (py engine_c/R): median 1.84x, min 0.89x, max 9.54x

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
| classical_starsCYG | bisquare | 47x2 | 1.68e-12 | 4.21e-12 | 3.84e-07 |
| classical_wood | bisquare | 20x6 | 1.91e-09 | 2.12e-07 | 6.92e-08 |
| psi_bisquare | bisquare | 21x4 | 5.68e-07 | 1.43e-06 | 4.16e-06 |
| psi_ggw | ggw | 21x4 | 1.33e-03 | 2.51e-03 | 5.20e-01 |
| psi_hampel | hampel | 21x4 | 1.98e-03 | 3.45e-03 | 4.01e-01 |
| psi_lqq | lqq | 21x4 | 1.86e-06 | 3.69e-06 | 7.91e-06 |
| psi_optimal | optimal | 21x4 | 4.58e-16 | 5.80e-07 | 8.87e-14 |
| setting_KS2011_stackloss | lqq | 21x4 | 7.42e-08 | 8.73e-08 | 3.04e-07 |
| setting_KS2014_stackloss | lqq | 21x4 | 7.31e-08 | 8.61e-08 | 3.00e-07 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 3.48e-13 | 1.78e-11 | 3.74e-08 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.50e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.67e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.66e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 5.85e-11 | 5.81e-09 | 1.98e-08 |
| synth_n10000_p50 | bisquare | 10000x51 | 1.12e-11 | 1.39e-09 | 2.22e-08 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.90e-14 | 3.92e-12 | 3.39e-09 |
| synth_n100_p5 | bisquare | 100x6 | 2.14e-12 | 5.35e-11 | 9.68e-09 |
| synth_n2000_p20 | bisquare | 2000x21 | 3.48e-13 | 1.78e-11 | 3.74e-08 |
| synth_n5000_p20 | bisquare | 5000x21 | 1.13e-10 | 8.35e-09 | 9.11e-09 |
| synth_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.35e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.43e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 11 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 2.8 | 13.1 | 4.65x | 4.3 | 1.52x |
| classical_coleman | bisquare | 20x6 | 3.1 | 13.3 | 4.32x | 4.3 | 1.38x |
| classical_delivery | bisquare | 25x3 | 2.3 | 12.3 | 5.25x | 4.5 | 1.92x |
| classical_hbk | bisquare | 75x4 | 4.6 | 16.6 | 3.63x | 8.0 | 1.74x |
| classical_pension | bisquare | 18x2 | 1.8 | 10.9 | 5.95x | 3.4 | 1.84x |
| classical_phosphor | bisquare | 18x3 | 2.1 | 11.4 | 5.34x | 3.4 | 1.58x |
| classical_salinity | bisquare | 28x4 | 2.8 | 13.2 | 4.76x | 4.8 | 1.73x |
| classical_stackloss | bisquare | 21x4 | 2.5 | 12.5 | 5.03x | 4.2 | 1.70x |
| classical_starsCYG | bisquare | 47x2 | 2.6 | 12.8 | 4.85x | 4.7 | 1.79x |
| classical_wood | bisquare | 20x6 | 3.1 | 13.2 | 4.32x | 4.1 | 1.34x |
| psi_bisquare | bisquare | 21x4 | 2.5 | 12.7 | 5.06x | 3.6 | 1.43x |
| psi_ggw | ggw | 21x4 | 4.3 | 14.9 | 3.45x | 5.7 | 1.32x |
| psi_hampel | hampel | 21x4 | 2.7 | 13.5 | 4.99x | 4.4 | 1.64x |
| psi_lqq | lqq | 21x4 | 3.0 | 14.5 | 4.82x | 4.9 | 1.61x |
| psi_optimal | optimal | 21x4 | 2.6 | 12.2 | 4.72x | 3.8 | 1.45x |
| setting_KS2011_stackloss | lqq | 21x4 | 3.8 | 15.2 | 4.00x | 6.1 | 1.60x |
| setting_KS2014_stackloss | lqq | 21x4 | 6.3 | 14.6 | 2.33x | 5.6 | 0.89x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 285.6 | 612.7 | 2.15x | 571.6 | 2.00x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 30.9 | 76.8 | 2.48x | 60.9 | 1.97x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 459.5 | 833.8 | 1.81x | 813.8 | 1.77x |
| synth_ggw_n500_p10 | ggw | 500x11 | 62.2 | 128.2 | 2.06x | 114.1 | 1.83x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 318.4 | 845.6 | 2.66x | 774.0 | 2.43x |
| synth_hampel_n500_p10 | hampel | 500x11 | 37.3 | 119.1 | 3.19x | 102.6 | 2.75x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 345.9 | 934.5 | 2.70x | 740.4 | 2.14x |
| synth_lqq_n500_p10 | lqq | 500x11 | 40.0 | 117.5 | 2.94x | 102.9 | 2.57x |
| synth_n10000_p20 | bisquare | 10000x21 | 388.2 | 3878.4 | 9.99x | 3702.8 | 9.54x |
| synth_n10000_p50 | bisquare | 10000x51 | 1557.3 | 7876.5 | 5.06x | 7628.5 | 4.90x |
| synth_n1000_p10 | bisquare | 1000x11 | 91.9 | 263.8 | 2.87x | 173.6 | 1.89x |
| synth_n100_p5 | bisquare | 100x6 | 6.2 | 21.3 | 3.45x | 12.0 | 1.94x |
| synth_n2000_p20 | bisquare | 2000x21 | 292.6 | 722.1 | 2.47x | 570.5 | 1.95x |
| synth_n5000_p20 | bisquare | 5000x21 | 307.4 | 1667.6 | 5.42x | 1433.1 | 4.66x |
| synth_n500_p10 | bisquare | 500x11 | 31.3 | 75.8 | 2.42x | 62.5 | 2.00x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 318.8 | 695.0 | 2.18x | 606.5 | 1.90x |
| synth_optimal_n500_p10 | optimal | 500x11 | 34.9 | 91.2 | 2.62x | 76.6 | 2.20x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
