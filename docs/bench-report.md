# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 3.71x, min 1.52x, max 6.24x
- Runtime ratio (py engine_c/R): median 1.73x, min 0.86x, max 5.95x

## Environment

- pylmrob: 0.5.30
- Python: 3.12.13
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
| classical_aircraft | bisquare | 23x5 | 2.1 | 9.1 | 4.36x | 3.1 | 1.47x |
| classical_coleman | bisquare | 20x6 | 2.2 | 9.0 | 4.03x | 3.0 | 1.34x |
| classical_delivery | bisquare | 25x3 | 1.7 | 8.4 | 5.07x | 2.9 | 1.77x |
| classical_hbk | bisquare | 75x4 | 3.0 | 11.2 | 3.70x | 6.2 | 2.04x |
| classical_pension | bisquare | 18x2 | 1.2 | 7.6 | 6.24x | 2.3 | 1.92x |
| classical_phosphor | bisquare | 18x3 | 1.5 | 7.7 | 5.25x | 2.3 | 1.58x |
| classical_salinity | bisquare | 28x4 | 2.0 | 9.0 | 4.48x | 3.3 | 1.64x |
| classical_stackloss | bisquare | 21x4 | 1.8 | 8.5 | 4.65x | 2.6 | 1.43x |
| classical_starsCYG | bisquare | 47x2 | 1.9 | 8.6 | 4.62x | 3.5 | 1.87x |
| classical_wood | bisquare | 20x6 | 2.4 | 8.9 | 3.73x | 2.9 | 1.23x |
| psi_bisquare | bisquare | 21x4 | 1.7 | 8.5 | 4.84x | 2.5 | 1.45x |
| psi_ggw | ggw | 21x4 | 3.2 | 10.3 | 3.25x | 4.1 | 1.30x |
| psi_hampel | hampel | 21x4 | 1.9 | 9.2 | 4.90x | 3.2 | 1.70x |
| psi_lqq | lqq | 21x4 | 2.2 | 9.5 | 4.34x | 3.5 | 1.61x |
| psi_optimal | optimal | 21x4 | 1.8 | 8.4 | 4.62x | 2.7 | 1.50x |
| setting_KS2011_stackloss | lqq | 21x4 | 2.7 | 10.1 | 3.71x | 4.1 | 1.49x |
| setting_KS2014_stackloss | lqq | 21x4 | 4.7 | 10.1 | 2.13x | 4.1 | 0.86x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 219.8 | 427.1 | 1.94x | 340.4 | 1.55x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 24.7 | 54.0 | 2.19x | 49.8 | 2.02x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 349.9 | 627.5 | 1.79x | 533.3 | 1.52x |
| synth_ggw_n500_p10 | ggw | 500x11 | 46.4 | 95.1 | 2.05x | 85.6 | 1.85x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 261.2 | 586.9 | 2.25x | 498.4 | 1.91x |
| synth_hampel_n500_p10 | hampel | 500x11 | 29.4 | 81.8 | 2.78x | 73.1 | 2.48x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 278.9 | 495.4 | 1.78x | 491.8 | 1.76x |
| synth_lqq_n500_p10 | lqq | 500x11 | 32.1 | 82.0 | 2.56x | 73.6 | 2.29x |
| synth_n10000_p20 | bisquare | 10000x21 | 305.1 | 1875.0 | 6.15x | 1814.3 | 5.95x |
| synth_n10000_p50 | bisquare | 10000x51 | 1069.0 | 4140.9 | 3.87x | 4130.5 | 3.86x |
| synth_n1000_p10 | bisquare | 1000x11 | 72.3 | 161.8 | 2.24x | 127.2 | 1.76x |
| synth_n100_p5 | bisquare | 100x6 | 4.9 | 14.6 | 2.99x | 8.9 | 1.83x |
| synth_n2000_p20 | bisquare | 2000x21 | 232.3 | 354.1 | 1.52x | 367.1 | 1.58x |
| synth_n5000_p20 | bisquare | 5000x21 | 252.3 | 1023.1 | 4.06x | 896.5 | 3.55x |
| synth_n500_p10 | bisquare | 500x11 | 25.7 | 53.8 | 2.10x | 48.1 | 1.87x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 254.9 | 558.7 | 2.19x | 424.9 | 1.67x |
| synth_optimal_n500_p10 | optimal | 500x11 | 29.7 | 65.3 | 2.20x | 58.8 | 1.98x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
