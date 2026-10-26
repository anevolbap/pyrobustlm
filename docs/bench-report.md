# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.69x, min 1.66x, max 9.51x
- Runtime ratio (py engine_c/R): median 1.50x, min 0.78x, max 5.12x

## Environment

- pylmrob: 0.5.33
- Python: 3.12.15
- Platform: Linux-7.0.0-1017-azure-x86_64-with-glibc2.43
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
| psi_optimal | optimal | 21x4 | 1.44e-15 | 5.80e-07 | 1.72e-13 |
| setting_KS2011_stackloss | lqq | 21x4 | 7.42e-08 | 8.73e-08 | 3.04e-07 |
| setting_KS2014_stackloss | lqq | 21x4 | 7.31e-08 | 8.61e-08 | 3.00e-07 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 3.51e-13 | 1.78e-11 | 3.74e-08 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.50e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.66e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.67e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 5.85e-11 | 5.81e-09 | 1.98e-08 |
| synth_n10000_p50 | bisquare | 10000x51 | 1.12e-11 | 1.39e-09 | 2.22e-08 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.86e-14 | 3.92e-12 | 3.39e-09 |
| synth_n100_p5 | bisquare | 100x6 | 2.14e-12 | 5.35e-11 | 9.68e-09 |
| synth_n2000_p20 | bisquare | 2000x21 | 3.51e-13 | 1.78e-11 | 3.74e-08 |
| synth_n5000_p20 | bisquare | 5000x21 | 1.13e-10 | 8.35e-09 | 9.11e-09 |
| synth_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.38e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.44e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 11 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.8 | 26.2 | 6.92x | 4.8 | 1.27x |
| classical_coleman | bisquare | 20x6 | 4.1 | 26.2 | 6.41x | 4.7 | 1.15x |
| classical_delivery | bisquare | 25x3 | 3.1 | 24.8 | 8.11x | 4.7 | 1.54x |
| classical_hbk | bisquare | 75x4 | 5.2 | 29.2 | 5.58x | 9.0 | 1.71x |
| classical_pension | bisquare | 18x2 | 2.3 | 21.9 | 9.51x | 3.8 | 1.67x |
| classical_phosphor | bisquare | 18x3 | 2.8 | 23.6 | 8.51x | 3.8 | 1.37x |
| classical_salinity | bisquare | 28x4 | 3.6 | 25.6 | 7.11x | 5.2 | 1.44x |
| classical_stackloss | bisquare | 21x4 | 3.3 | 25.5 | 7.76x | 4.1 | 1.25x |
| classical_starsCYG | bisquare | 47x2 | 3.4 | 24.1 | 7.17x | 5.3 | 1.58x |
| classical_wood | bisquare | 20x6 | 4.1 | 26.1 | 6.35x | 4.7 | 1.15x |
| psi_bisquare | bisquare | 21x4 | 3.2 | 25.6 | 7.97x | 4.2 | 1.29x |
| psi_ggw | ggw | 21x4 | 5.6 | 29.1 | 5.16x | 6.8 | 1.20x |
| psi_hampel | hampel | 21x4 | 3.4 | 26.9 | 7.88x | 5.0 | 1.46x |
| psi_lqq | lqq | 21x4 | 3.9 | 27.5 | 6.95x | 5.6 | 1.42x |
| psi_optimal | optimal | 21x4 | 3.3 | 24.9 | 7.48x | 4.2 | 1.26x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.2 | 28.6 | 5.54x | 6.8 | 1.33x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.9 | 28.8 | 3.25x | 6.9 | 0.78x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 381.7 | 688.3 | 1.80x | 524.5 | 1.37x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 41.1 | 93.9 | 2.29x | 67.2 | 1.64x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 611.7 | 1012.8 | 1.66x | 857.8 | 1.40x |
| synth_ggw_n500_p10 | ggw | 500x11 | 73.6 | 163.0 | 2.22x | 133.2 | 1.81x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 439.0 | 858.0 | 1.95x | 729.6 | 1.66x |
| synth_hampel_n500_p10 | hampel | 500x11 | 48.2 | 131.4 | 2.73x | 104.8 | 2.18x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 473.4 | 829.9 | 1.75x | 737.3 | 1.56x |
| synth_lqq_n500_p10 | lqq | 500x11 | 53.4 | 137.2 | 2.57x | 111.4 | 2.09x |
| synth_n10000_p20 | bisquare | 10000x21 | 517.7 | 2906.8 | 5.61x | 2651.9 | 5.12x |
| synth_n10000_p50 | bisquare | 10000x51 | 1799.7 | 7344.4 | 4.08x | 7203.0 | 4.00x |
| synth_n1000_p10 | bisquare | 1000x11 | 110.9 | 270.0 | 2.44x | 180.7 | 1.63x |
| synth_n100_p5 | bisquare | 100x6 | 8.2 | 34.7 | 4.22x | 12.8 | 1.56x |
| synth_n2000_p20 | bisquare | 2000x21 | 381.2 | 734.4 | 1.93x | 529.8 | 1.39x |
| synth_n5000_p20 | bisquare | 5000x21 | 411.9 | 1557.4 | 3.78x | 1314.3 | 3.19x |
| synth_n500_p10 | bisquare | 500x11 | 41.6 | 93.8 | 2.26x | 67.0 | 1.61x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 440.2 | 767.4 | 1.74x | 559.9 | 1.27x |
| synth_optimal_n500_p10 | optimal | 500x11 | 47.1 | 103.8 | 2.20x | 75.5 | 1.60x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
