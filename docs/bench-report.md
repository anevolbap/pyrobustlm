# Benchmark report

Element-wise comparison between `pyrobustlm` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 2.75x, min 0.89x, max 14.02x
- Runtime ratio (py engine_c/R): median 1.01x, min 0.35x, max 11.32x

## Environment

- pyrobustlm: 0.5.10
- Python: 3.11.2
- Platform: Linux-6.1.0-47-amd64-x86_64-with-glibc2.36
- robustbase: 0.99.7
- R: R version 4.2.2 Patched (2022-11-10 r83330)

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
| psi_optimal | optimal | 21x4 | 1.32e-15 | 5.80e-07 | 1.83e-13 |
| setting_KS2011_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| setting_KS2014_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.94e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.49e-13 | 8.35e-12 | 1.46e-08 |
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
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.32e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.43e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 5 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 5.0 | 25.9 | 5.18x | 4.3 | 0.86x |
| classical_coleman | bisquare | 20x6 | 5.0 | 24.5 | 4.91x | 4.3 | 0.86x |
| classical_delivery | bisquare | 25x3 | 4.0 | 24.7 | 6.17x | 4.3 | 1.08x |
| classical_hbk | bisquare | 75x4 | 6.0 | 31.2 | 5.20x | 48.8 | 8.13x |
| classical_pension | bisquare | 18x2 | 3.0 | 23.5 | 7.83x | 3.3 | 1.09x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 24.6 | 8.19x | 3.6 | 1.19x |
| classical_salinity | bisquare | 28x4 | 5.0 | 23.9 | 4.78x | 5.6 | 1.13x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 24.8 | 6.19x | 30.6 | 7.65x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 24.4 | 8.12x | 5.1 | 1.71x |
| classical_wood | bisquare | 20x6 | 5.0 | 25.8 | 5.15x | 5.3 | 1.06x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 24.9 | 8.31x | 34.0 | 11.32x |
| psi_ggw | ggw | 21x4 | 11.0 | 27.2 | 2.48x | 8.4 | 0.76x |
| psi_hampel | hampel | 21x4 | 6.0 | 25.0 | 4.17x | 5.8 | 0.97x |
| psi_lqq | lqq | 21x4 | 5.0 | 28.7 | 5.73x | 5.8 | 1.16x |
| psi_optimal | optimal | 21x4 | 5.0 | 26.9 | 5.37x | 4.3 | 0.87x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 70.1 | 14.02x | 7.9 | 1.57x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.0 | 71.6 | 8.95x | 9.3 | 1.16x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 559.0 | 533.7 | 0.95x | 211.9 | 0.38x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 70.0 | 123.5 | 1.76x | 68.9 | 0.98x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 696.0 | 788.3 | 1.13x | 334.4 | 0.48x |
| synth_ggw_n500_p10 | ggw | 500x11 | 84.0 | 179.8 | 2.14x | 152.7 | 1.82x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 588.0 | 630.8 | 1.07x | 227.7 | 0.39x |
| synth_hampel_n500_p10 | hampel | 500x11 | 64.0 | 143.6 | 2.24x | 109.5 | 1.71x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 634.0 | 862.4 | 1.36x | 260.6 | 0.41x |
| synth_lqq_n500_p10 | lqq | 500x11 | 70.0 | 202.2 | 2.89x | 126.9 | 1.81x |
| synth_n10000_p20 | bisquare | 10000x21 | 1005.0 | 2592.0 | 2.58x | 806.8 | 0.80x |
| synth_n10000_p50 | bisquare | 10000x51 | 4750.0 | 7399.9 | 1.56x | 3295.2 | 0.69x |
| synth_n1000_p10 | bisquare | 1000x11 | 128.0 | 213.5 | 1.67x | 100.9 | 0.79x |
| synth_n100_p5 | bisquare | 100x6 | 15.0 | 37.0 | 2.47x | 13.0 | 0.87x |
| synth_n2000_p20 | bisquare | 2000x21 | 576.0 | 512.3 | 0.89x | 199.2 | 0.35x |
| synth_n5000_p20 | bisquare | 5000x21 | 1005.0 | 1254.2 | 1.25x | 439.1 | 0.44x |
| synth_n500_p10 | bisquare | 500x11 | 66.0 | 117.1 | 1.77x | 68.5 | 1.04x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 575.0 | 764.8 | 1.33x | 257.0 | 0.45x |
| synth_optimal_n500_p10 | optimal | 500x11 | 69.0 | 180.8 | 2.62x | 79.5 | 1.15x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
