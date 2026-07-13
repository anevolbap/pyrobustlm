# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 5.71x, min 1.58x, max 18.76x
- Runtime ratio (py engine_c/R): median 1.66x, min 0.74x, max 10.45x

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
| classical_aircraft | bisquare | 23x5 | 4.0 | 27.5 | 6.86x | 4.9 | 1.22x |
| classical_coleman | bisquare | 20x6 | 4.0 | 27.4 | 6.85x | 4.8 | 1.21x |
| classical_delivery | bisquare | 25x3 | 3.0 | 26.5 | 8.84x | 4.7 | 1.56x |
| classical_hbk | bisquare | 75x4 | 6.0 | 35.2 | 5.86x | 44.6 | 7.44x |
| classical_pension | bisquare | 18x2 | 2.0 | 23.2 | 11.61x | 3.8 | 1.88x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 24.5 | 8.15x | 3.8 | 1.25x |
| classical_salinity | bisquare | 28x4 | 4.0 | 27.8 | 6.95x | 5.2 | 1.31x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 26.7 | 6.68x | 31.8 | 7.94x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 27.2 | 6.79x | 5.3 | 1.32x |
| classical_wood | bisquare | 20x6 | 4.0 | 27.2 | 6.79x | 4.8 | 1.20x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 26.6 | 8.87x | 31.3 | 10.45x |
| psi_ggw | ggw | 21x4 | 6.0 | 28.3 | 4.71x | 6.6 | 1.11x |
| psi_hampel | hampel | 21x4 | 3.0 | 26.9 | 8.96x | 5.0 | 1.68x |
| psi_lqq | lqq | 21x4 | 4.0 | 30.0 | 7.50x | 5.6 | 1.40x |
| psi_optimal | optimal | 21x4 | 3.0 | 27.2 | 9.05x | 4.3 | 1.44x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 93.8 | 18.76x | 6.7 | 1.34x |
| setting_KS2014_stackloss | lqq | 21x4 | 9.0 | 93.4 | 10.38x | 6.7 | 0.74x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 395.0 | 725.5 | 1.84x | 587.5 | 1.49x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.0 | 142.0 | 3.55x | 74.6 | 1.87x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 634.0 | 1000.8 | 1.58x | 906.8 | 1.43x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.0 | 167.3 | 2.26x | 138.0 | 1.86x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 447.0 | 1055.7 | 2.36x | 816.5 | 1.83x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.0 | 161.6 | 3.44x | 120.5 | 2.56x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 483.0 | 1236.2 | 2.56x | 797.5 | 1.65x |
| synth_lqq_n500_p10 | lqq | 500x11 | 54.0 | 245.8 | 4.55x | 119.5 | 2.21x |
| synth_n10000_p20 | bisquare | 10000x21 | 504.0 | 3671.3 | 7.28x | 2947.1 | 5.85x |
| synth_n10000_p50 | bisquare | 10000x51 | 1884.0 | 7986.0 | 4.24x | 7200.6 | 3.82x |
| synth_n1000_p10 | bisquare | 1000x11 | 120.0 | 386.4 | 3.22x | 205.8 | 1.72x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 44.5 | 5.56x | 13.4 | 1.67x |
| synth_n2000_p20 | bisquare | 2000x21 | 395.0 | 874.2 | 2.21x | 582.2 | 1.47x |
| synth_n5000_p20 | bisquare | 5000x21 | 393.0 | 1941.6 | 4.94x | 1527.8 | 3.89x |
| synth_n500_p10 | bisquare | 500x11 | 40.0 | 144.1 | 3.60x | 74.6 | 1.86x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 450.0 | 1034.1 | 2.30x | 640.1 | 1.42x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.0 | 193.5 | 4.21x | 90.3 | 1.96x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
