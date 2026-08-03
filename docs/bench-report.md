# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.19x, min 1.73x, max 7.48x
- Runtime ratio (py engine_c/R): median 1.64x, min 0.85x, max 5.85x

## Environment

- pylmrob: 0.5.29
- Python: 3.12.13
- Platform: Linux-6.17.0-1020-azure-x86_64-with-glibc2.39
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
| synth_n10000_p50 | bisquare | 10000x51 | 1.11e-11 | 1.39e-09 | 2.22e-08 |
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
| classical_aircraft | bisquare | 23x5 | 3.6 | 19.2 | 5.37x | 5.2 | 1.46x |
| classical_coleman | bisquare | 20x6 | 3.8 | 19.1 | 5.02x | 5.0 | 1.33x |
| classical_delivery | bisquare | 25x3 | 2.9 | 17.8 | 6.11x | 5.2 | 1.77x |
| classical_hbk | bisquare | 75x4 | 5.2 | 22.6 | 4.36x | 9.5 | 1.83x |
| classical_pension | bisquare | 18x2 | 2.1 | 16.0 | 7.48x | 4.0 | 1.87x |
| classical_phosphor | bisquare | 18x3 | 2.6 | 16.8 | 6.45x | 4.0 | 1.53x |
| classical_salinity | bisquare | 28x4 | 3.5 | 18.9 | 5.46x | 5.5 | 1.60x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 18.1 | 5.97x | 4.4 | 1.44x |
| classical_starsCYG | bisquare | 47x2 | 3.3 | 17.7 | 5.40x | 5.5 | 1.68x |
| classical_wood | bisquare | 20x6 | 4.0 | 19.1 | 4.79x | 5.0 | 1.25x |
| psi_bisquare | bisquare | 21x4 | 3.1 | 18.2 | 5.89x | 4.3 | 1.40x |
| psi_ggw | ggw | 21x4 | 5.5 | 21.5 | 3.91x | 7.0 | 1.27x |
| psi_hampel | hampel | 21x4 | 4.7 | 19.5 | 4.18x | 5.4 | 1.15x |
| psi_lqq | lqq | 21x4 | 3.8 | 20.0 | 5.26x | 6.0 | 1.57x |
| psi_optimal | optimal | 21x4 | 3.2 | 17.7 | 5.53x | 4.6 | 1.42x |
| setting_KS2011_stackloss | lqq | 21x4 | 4.9 | 21.2 | 4.36x | 7.0 | 1.44x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.3 | 21.3 | 2.58x | 7.0 | 0.85x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 393.9 | 732.3 | 1.86x | 634.0 | 1.61x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 42.7 | 100.9 | 2.36x | 79.0 | 1.85x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 590.6 | 1024.1 | 1.73x | 956.4 | 1.62x |
| synth_ggw_n500_p10 | ggw | 500x11 | 76.8 | 165.2 | 2.15x | 144.4 | 1.88x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 438.3 | 935.0 | 2.13x | 831.9 | 1.90x |
| synth_hampel_n500_p10 | hampel | 500x11 | 50.9 | 150.6 | 2.96x | 129.3 | 2.54x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 470.9 | 1042.6 | 2.21x | 882.5 | 1.87x |
| synth_lqq_n500_p10 | lqq | 500x11 | 56.6 | 148.9 | 2.63x | 131.4 | 2.32x |
| synth_n10000_p20 | bisquare | 10000x21 | 507.3 | 3353.2 | 6.61x | 2966.9 | 5.85x |
| synth_n10000_p50 | bisquare | 10000x51 | 1939.3 | 8134.1 | 4.19x | 7951.1 | 4.10x |
| synth_n1000_p10 | bisquare | 1000x11 | 120.8 | 317.5 | 2.63x | 200.1 | 1.66x |
| synth_n100_p5 | bisquare | 100x6 | 8.2 | 28.9 | 3.51x | 14.3 | 1.74x |
| synth_n2000_p20 | bisquare | 2000x21 | 395.9 | 795.0 | 2.01x | 615.7 | 1.56x |
| synth_n5000_p20 | bisquare | 5000x21 | 406.7 | 1902.1 | 4.68x | 1554.0 | 3.82x |
| synth_n500_p10 | bisquare | 500x11 | 43.0 | 100.7 | 2.34x | 79.0 | 1.84x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 442.1 | 830.1 | 1.88x | 721.3 | 1.63x |
| synth_optimal_n500_p10 | optimal | 500x11 | 49.7 | 119.4 | 2.40x | 98.7 | 1.99x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
