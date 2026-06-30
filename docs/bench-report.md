# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 6.13x, min 1.66x, max 18.59x
- Runtime ratio (py engine_c/R): median 1.61x, min 0.74x, max 10.36x

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
| classical_aircraft | bisquare | 23x5 | 4.0 | 27.8 | 6.94x | 4.9 | 1.23x |
| classical_coleman | bisquare | 20x6 | 4.0 | 27.7 | 6.92x | 4.8 | 1.20x |
| classical_delivery | bisquare | 25x3 | 3.0 | 26.8 | 8.95x | 4.7 | 1.57x |
| classical_hbk | bisquare | 75x4 | 5.0 | 34.8 | 6.97x | 44.4 | 8.89x |
| classical_pension | bisquare | 18x2 | 2.0 | 23.7 | 11.84x | 3.8 | 1.89x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 24.7 | 8.24x | 3.8 | 1.26x |
| classical_salinity | bisquare | 28x4 | 4.0 | 27.9 | 6.97x | 5.2 | 1.31x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 26.8 | 8.95x | 31.1 | 10.36x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 27.5 | 6.88x | 5.4 | 1.34x |
| classical_wood | bisquare | 20x6 | 4.0 | 27.5 | 6.88x | 4.7 | 1.18x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 26.7 | 8.89x | 30.8 | 10.28x |
| psi_ggw | ggw | 21x4 | 6.0 | 28.4 | 4.73x | 6.6 | 1.10x |
| psi_hampel | hampel | 21x4 | 4.0 | 27.3 | 6.83x | 5.1 | 1.27x |
| psi_lqq | lqq | 21x4 | 4.0 | 30.3 | 7.58x | 5.6 | 1.40x |
| psi_optimal | optimal | 21x4 | 4.0 | 27.5 | 6.87x | 4.3 | 1.06x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 92.9 | 18.59x | 6.7 | 1.34x |
| setting_KS2014_stackloss | lqq | 21x4 | 9.0 | 91.7 | 10.19x | 6.7 | 0.74x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 389.0 | 768.0 | 1.97x | 597.5 | 1.54x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.0 | 140.4 | 3.51x | 74.8 | 1.87x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 630.0 | 1043.3 | 1.66x | 728.9 | 1.16x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.0 | 167.7 | 2.27x | 137.7 | 1.86x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 448.0 | 927.7 | 2.07x | 797.1 | 1.78x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.0 | 164.4 | 3.50x | 118.6 | 2.52x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 483.0 | 1270.6 | 2.63x | 793.7 | 1.64x |
| synth_lqq_n500_p10 | lqq | 500x11 | 52.0 | 244.6 | 4.70x | 117.4 | 2.26x |
| synth_n10000_p20 | bisquare | 10000x21 | 505.0 | 3763.6 | 7.45x | 2939.4 | 5.82x |
| synth_n10000_p50 | bisquare | 10000x51 | 1919.0 | 8091.4 | 4.22x | 7243.0 | 3.77x |
| synth_n1000_p10 | bisquare | 1000x11 | 119.0 | 342.8 | 2.88x | 210.5 | 1.77x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 43.4 | 5.42x | 13.9 | 1.73x |
| synth_n2000_p20 | bisquare | 2000x21 | 389.0 | 765.2 | 1.97x | 595.5 | 1.53x |
| synth_n5000_p20 | bisquare | 5000x21 | 389.0 | 1968.0 | 5.06x | 1469.9 | 3.78x |
| synth_n500_p10 | bisquare | 500x11 | 40.0 | 139.7 | 3.49x | 74.6 | 1.86x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 442.0 | 1054.4 | 2.39x | 666.5 | 1.51x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.0 | 192.6 | 4.19x | 93.1 | 2.02x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
