# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.63x, min 1.57x, max 9.38x
- Runtime ratio (py engine_c/R): median 1.48x, min 0.79x, max 5.15x

## Environment

- pylmrob: 0.5.33
- Python: 3.12.15
- Platform: Linux-7.0.0-1017-azure-x86_64-with-glibc2.43
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
| classical_aircraft | bisquare | 23x5 | 3.8 | 25.7 | 6.83x | 4.9 | 1.29x |
| classical_coleman | bisquare | 20x6 | 4.2 | 26.0 | 6.21x | 4.8 | 1.14x |
| classical_delivery | bisquare | 25x3 | 3.1 | 24.6 | 7.94x | 4.7 | 1.52x |
| classical_hbk | bisquare | 75x4 | 5.3 | 29.0 | 5.51x | 8.7 | 1.66x |
| classical_pension | bisquare | 18x2 | 2.3 | 21.7 | 9.38x | 3.9 | 1.68x |
| classical_phosphor | bisquare | 18x3 | 2.8 | 23.3 | 8.19x | 3.8 | 1.34x |
| classical_salinity | bisquare | 28x4 | 3.7 | 25.4 | 6.90x | 5.1 | 1.39x |
| classical_stackloss | bisquare | 21x4 | 3.2 | 25.3 | 7.84x | 4.1 | 1.28x |
| classical_starsCYG | bisquare | 47x2 | 3.4 | 24.0 | 7.10x | 5.3 | 1.57x |
| classical_wood | bisquare | 20x6 | 4.2 | 25.9 | 6.23x | 4.7 | 1.12x |
| psi_bisquare | bisquare | 21x4 | 3.2 | 25.5 | 7.86x | 4.1 | 1.26x |
| psi_ggw | ggw | 21x4 | 5.7 | 28.6 | 5.04x | 6.6 | 1.17x |
| psi_hampel | hampel | 21x4 | 3.5 | 26.7 | 7.69x | 5.0 | 1.44x |
| psi_lqq | lqq | 21x4 | 3.9 | 27.2 | 6.92x | 5.6 | 1.43x |
| psi_optimal | optimal | 21x4 | 3.4 | 24.7 | 7.35x | 4.2 | 1.24x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.1 | 28.6 | 5.62x | 6.8 | 1.33x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.7 | 28.5 | 3.29x | 6.8 | 0.79x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 380.5 | 690.6 | 1.81x | 522.7 | 1.37x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 41.2 | 93.9 | 2.28x | 67.0 | 1.63x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 623.4 | 980.1 | 1.57x | 820.0 | 1.32x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.7 | 163.5 | 2.19x | 132.9 | 1.78x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 434.0 | 852.7 | 1.96x | 713.0 | 1.64x |
| synth_hampel_n500_p10 | hampel | 500x11 | 48.2 | 132.0 | 2.74x | 102.9 | 2.13x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 472.8 | 844.8 | 1.79x | 736.7 | 1.56x |
| synth_lqq_n500_p10 | lqq | 500x11 | 53.2 | 137.9 | 2.59x | 108.6 | 2.04x |
| synth_n10000_p20 | bisquare | 10000x21 | 518.0 | 2914.8 | 5.63x | 2669.3 | 5.15x |
| synth_n10000_p50 | bisquare | 10000x51 | 1876.1 | 7165.0 | 3.82x | 7029.0 | 3.75x |
| synth_n1000_p10 | bisquare | 1000x11 | 110.8 | 264.0 | 2.38x | 182.2 | 1.64x |
| synth_n100_p5 | bisquare | 100x6 | 8.2 | 34.6 | 4.23x | 12.8 | 1.56x |
| synth_n2000_p20 | bisquare | 2000x21 | 380.9 | 651.9 | 1.71x | 530.0 | 1.39x |
| synth_n5000_p20 | bisquare | 5000x21 | 410.7 | 1540.5 | 3.75x | 1352.3 | 3.29x |
| synth_n500_p10 | bisquare | 500x11 | 41.3 | 93.7 | 2.27x | 67.2 | 1.63x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 434.2 | 736.5 | 1.70x | 556.3 | 1.28x |
| synth_optimal_n500_p10 | optimal | 500x11 | 47.3 | 103.5 | 2.19x | 75.5 | 1.59x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
