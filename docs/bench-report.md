# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.63e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.32x, min 1.18x, max 9.21x
- Runtime ratio (py engine_c/R): median 1.30x, min 0.82x, max 2.91x

## Environment

- pylmrob: 0.5.24
- Python: 3.12.13
- Platform: Linux-6.12.95+deb12-amd64-x86_64-with-glibc2.36
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
| psi_optimal | optimal | 21x4 | 1.55e-15 | 5.80e-07 | 9.35e-14 |
| setting_KS2011_stackloss | lqq | 21x4 | 7.42e-08 | 8.73e-08 | 3.04e-07 |
| setting_KS2014_stackloss | lqq | 21x4 | 7.31e-08 | 8.61e-08 | 3.00e-07 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.94e-12 | 8.82e-11 | 2.43e-07 |
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
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.31e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.44e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 11 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.1 | 22.4 | 7.14x | 3.7 | 1.17x |
| classical_coleman | bisquare | 20x6 | 3.1 | 21.6 | 7.08x | 4.5 | 1.46x |
| classical_delivery | bisquare | 25x3 | 2.3 | 21.3 | 9.21x | 3.7 | 1.58x |
| classical_hbk | bisquare | 75x4 | 5.1 | 27.2 | 5.32x | 7.5 | 1.46x |
| classical_pension | bisquare | 18x2 | 1.9 | 17.6 | 9.05x | 2.9 | 1.48x |
| classical_phosphor | bisquare | 18x3 | 2.2 | 19.9 | 8.89x | 2.9 | 1.30x |
| classical_salinity | bisquare | 28x4 | 3.1 | 23.0 | 7.47x | 4.8 | 1.56x |
| classical_stackloss | bisquare | 21x4 | 2.4 | 19.7 | 8.11x | 4.0 | 1.65x |
| classical_starsCYG | bisquare | 47x2 | 2.9 | 18.5 | 6.35x | 4.2 | 1.45x |
| classical_wood | bisquare | 20x6 | 3.4 | 20.5 | 6.11x | 3.8 | 1.12x |
| psi_bisquare | bisquare | 21x4 | 2.7 | 20.7 | 7.58x | 3.3 | 1.19x |
| psi_ggw | ggw | 21x4 | 4.5 | 22.0 | 4.93x | 5.1 | 1.15x |
| psi_hampel | hampel | 21x4 | 2.9 | 18.7 | 6.56x | 4.0 | 1.39x |
| psi_lqq | lqq | 21x4 | 3.0 | 21.6 | 7.15x | 4.3 | 1.43x |
| psi_optimal | optimal | 21x4 | 2.9 | 19.9 | 6.82x | 3.2 | 1.11x |
| setting_KS2011_stackloss | lqq | 21x4 | 4.6 | 23.0 | 5.03x | 5.4 | 1.19x |
| setting_KS2014_stackloss | lqq | 21x4 | 6.9 | 24.3 | 3.53x | 5.7 | 0.82x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 470.4 | 608.8 | 1.29x | 422.1 | 0.90x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 50.4 | 108.2 | 2.15x | 56.4 | 1.12x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 602.0 | 709.2 | 1.18x | 711.8 | 1.18x |
| synth_ggw_n500_p10 | ggw | 500x11 | 80.4 | 140.4 | 1.75x | 115.8 | 1.44x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 544.9 | 709.2 | 1.30x | 585.2 | 1.07x |
| synth_hampel_n500_p10 | hampel | 500x11 | 59.0 | 135.7 | 2.30x | 97.2 | 1.65x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 576.8 | 975.4 | 1.69x | 603.9 | 1.05x |
| synth_lqq_n500_p10 | lqq | 500x11 | 66.4 | 197.6 | 2.98x | 102.1 | 1.54x |
| synth_n10000_p20 | bisquare | 10000x21 | 820.3 | 3254.6 | 3.97x | 2388.8 | 2.91x |
| synth_n10000_p50 | bisquare | 10000x51 | 3662.3 | 10572.0 | 2.89x | 9478.1 | 2.59x |
| synth_n1000_p10 | bisquare | 1000x11 | 99.5 | 238.0 | 2.39x | 131.5 | 1.32x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 37.5 | 4.67x | 10.4 | 1.30x |
| synth_n2000_p20 | bisquare | 2000x21 | 494.5 | 592.6 | 1.20x | 440.8 | 0.89x |
| synth_n5000_p20 | bisquare | 5000x21 | 595.8 | 1539.1 | 2.58x | 1178.0 | 1.98x |
| synth_n500_p10 | bisquare | 500x11 | 50.6 | 115.2 | 2.28x | 59.2 | 1.17x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 537.3 | 782.7 | 1.46x | 476.5 | 0.89x |
| synth_optimal_n500_p10 | optimal | 500x11 | 62.5 | 149.2 | 2.39x | 78.6 | 1.26x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
