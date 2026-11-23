# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 3.89x, min 1.57x, max 6.70x
- Runtime ratio (py engine_c/R): median 1.50x, min 0.79x, max 4.91x

## Environment

- pylmrob: 0.5.33
- Python: 3.12.16
- Platform: Linux-7.0.0-1020-azure-x86_64-with-glibc2.43
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
| classical_aircraft | bisquare | 23x5 | 3.0 | 14.2 | 4.72x | 3.9 | 1.29x |
| classical_coleman | bisquare | 20x6 | 3.2 | 14.3 | 4.39x | 3.8 | 1.16x |
| classical_delivery | bisquare | 25x3 | 2.4 | 13.4 | 5.54x | 3.8 | 1.58x |
| classical_hbk | bisquare | 75x4 | 4.3 | 17.1 | 3.98x | 7.2 | 1.67x |
| classical_pension | bisquare | 18x2 | 1.8 | 11.9 | 6.70x | 3.1 | 1.77x |
| classical_phosphor | bisquare | 18x3 | 2.2 | 12.5 | 5.71x | 3.1 | 1.41x |
| classical_salinity | bisquare | 28x4 | 2.9 | 14.0 | 4.91x | 4.2 | 1.45x |
| classical_stackloss | bisquare | 21x4 | 2.6 | 13.6 | 5.24x | 3.3 | 1.27x |
| classical_starsCYG | bisquare | 47x2 | 2.7 | 13.4 | 5.01x | 4.3 | 1.59x |
| classical_wood | bisquare | 20x6 | 3.3 | 14.2 | 4.29x | 3.7 | 1.12x |
| psi_bisquare | bisquare | 21x4 | 2.6 | 13.7 | 5.31x | 3.3 | 1.26x |
| psi_ggw | ggw | 21x4 | 4.3 | 16.4 | 3.85x | 5.4 | 1.27x |
| psi_hampel | hampel | 21x4 | 2.7 | 14.9 | 5.43x | 4.1 | 1.49x |
| psi_lqq | lqq | 21x4 | 3.1 | 15.2 | 4.89x | 4.6 | 1.46x |
| psi_optimal | optimal | 21x4 | 2.7 | 13.4 | 4.94x | 3.3 | 1.23x |
| setting_KS2011_stackloss | lqq | 21x4 | 3.9 | 16.0 | 4.08x | 5.3 | 1.36x |
| setting_KS2014_stackloss | lqq | 21x4 | 6.8 | 16.1 | 2.37x | 5.3 | 0.79x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 330.1 | 545.7 | 1.65x | 445.1 | 1.35x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 35.0 | 71.5 | 2.04x | 56.6 | 1.62x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 472.5 | 846.0 | 1.79x | 672.9 | 1.42x |
| synth_ggw_n500_p10 | ggw | 500x11 | 58.6 | 130.5 | 2.23x | 113.2 | 1.93x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 365.3 | 678.0 | 1.86x | 597.9 | 1.64x |
| synth_hampel_n500_p10 | hampel | 500x11 | 42.1 | 106.3 | 2.53x | 90.1 | 2.14x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 392.9 | 717.6 | 1.83x | 594.7 | 1.51x |
| synth_lqq_n500_p10 | lqq | 500x11 | 45.9 | 109.3 | 2.38x | 93.1 | 2.03x |
| synth_n10000_p20 | bisquare | 10000x21 | 445.4 | 2409.0 | 5.41x | 2188.3 | 4.91x |
| synth_n10000_p50 | bisquare | 10000x51 | 1909.0 | 5539.2 | 2.90x | 5454.9 | 2.86x |
| synth_n1000_p10 | bisquare | 1000x11 | 88.3 | 237.2 | 2.69x | 143.3 | 1.62x |
| synth_n100_p5 | bisquare | 100x6 | 6.9 | 21.4 | 3.12x | 10.6 | 1.55x |
| synth_n2000_p20 | bisquare | 2000x21 | 345.7 | 544.0 | 1.57x | 459.9 | 1.33x |
| synth_n5000_p20 | bisquare | 5000x21 | 351.5 | 1383.5 | 3.94x | 1132.9 | 3.22x |
| synth_n500_p10 | bisquare | 500x11 | 35.4 | 71.4 | 2.01x | 56.4 | 1.59x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 367.2 | 639.0 | 1.74x | 478.4 | 1.30x |
| synth_optimal_n500_p10 | optimal | 500x11 | 40.3 | 79.7 | 1.97x | 63.8 | 1.58x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
