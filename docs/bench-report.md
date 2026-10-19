# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.96x, min 1.51x, max 9.70x
- Runtime ratio (py engine_c/R): median 1.53x, min 0.83x, max 5.61x

## Environment

- pylmrob: 0.5.33
- Python: 3.12.15
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
| classical_aircraft | bisquare | 23x5 | 3.7 | 26.9 | 7.36x | 5.1 | 1.38x |
| classical_coleman | bisquare | 20x6 | 3.9 | 27.0 | 6.91x | 4.9 | 1.25x |
| classical_delivery | bisquare | 25x3 | 3.8 | 25.4 | 6.73x | 4.8 | 1.27x |
| classical_hbk | bisquare | 75x4 | 5.1 | 30.6 | 5.98x | 8.9 | 1.75x |
| classical_pension | bisquare | 18x2 | 2.3 | 22.6 | 9.70x | 3.9 | 1.67x |
| classical_phosphor | bisquare | 18x3 | 2.8 | 24.1 | 8.63x | 3.9 | 1.39x |
| classical_salinity | bisquare | 28x4 | 3.7 | 26.4 | 7.13x | 5.4 | 1.46x |
| classical_stackloss | bisquare | 21x4 | 3.2 | 26.3 | 8.28x | 4.2 | 1.33x |
| classical_starsCYG | bisquare | 47x2 | 3.3 | 25.1 | 7.52x | 5.3 | 1.59x |
| classical_wood | bisquare | 20x6 | 4.0 | 26.9 | 6.77x | 4.8 | 1.21x |
| psi_bisquare | bisquare | 21x4 | 3.1 | 26.3 | 8.42x | 4.2 | 1.34x |
| psi_ggw | ggw | 21x4 | 5.6 | 29.4 | 5.22x | 6.7 | 1.19x |
| psi_hampel | hampel | 21x4 | 3.3 | 27.5 | 8.29x | 5.1 | 1.55x |
| psi_lqq | lqq | 21x4 | 3.9 | 28.0 | 7.19x | 5.7 | 1.47x |
| psi_optimal | optimal | 21x4 | 3.2 | 25.6 | 7.94x | 4.4 | 1.35x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.1 | 29.4 | 5.77x | 7.0 | 1.37x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.5 | 29.4 | 3.47x | 7.0 | 0.83x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 394.3 | 783.3 | 1.99x | 574.0 | 1.46x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.3 | 103.5 | 2.57x | 72.7 | 1.80x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 620.5 | 1081.8 | 1.74x | 901.1 | 1.45x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.1 | 168.6 | 2.27x | 138.0 | 1.86x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 443.6 | 1004.8 | 2.27x | 808.9 | 1.82x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.4 | 148.9 | 3.14x | 121.3 | 2.56x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 471.0 | 792.4 | 1.68x | 828.7 | 1.76x |
| synth_lqq_n500_p10 | lqq | 500x11 | 52.5 | 148.3 | 2.82x | 124.0 | 2.36x |
| synth_n10000_p20 | bisquare | 10000x21 | 520.4 | 3220.4 | 6.19x | 2917.1 | 5.61x |
| synth_n10000_p50 | bisquare | 10000x51 | 1928.7 | 7657.8 | 3.97x | 7408.9 | 3.84x |
| synth_n1000_p10 | bisquare | 1000x11 | 121.0 | 347.0 | 2.87x | 203.1 | 1.68x |
| synth_n100_p5 | bisquare | 100x6 | 7.8 | 36.9 | 4.70x | 13.4 | 1.71x |
| synth_n2000_p20 | bisquare | 2000x21 | 393.6 | 764.0 | 1.94x | 583.2 | 1.48x |
| synth_n5000_p20 | bisquare | 5000x21 | 403.1 | 1726.9 | 4.28x | 1436.5 | 3.56x |
| synth_n500_p10 | bisquare | 500x11 | 40.7 | 103.2 | 2.54x | 72.5 | 1.78x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 444.4 | 672.2 | 1.51x | 670.9 | 1.51x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.5 | 121.9 | 2.62x | 94.1 | 2.03x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
