# Benchmark report

Element-wise comparison between `pyrobustlm` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 14.15x, min 1.25x, max 36.92x

## Environment

- pyrobustlm: 0.1.0
- Python: 3.11.2
- Platform: Linux-6.1.0-45-amd64-x86_64-with-glibc2.36
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
| psi_optimal | optimal | 21x4 | 1.22e-15 | 5.80e-07 | 2.43e-13 |
| setting_KS2011_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| setting_KS2014_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.50e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.67e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.48e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.67e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 3.24e-08 | 3.24e-06 | 4.66e-07 |
| synth_n10000_p50 | bisquare | 10000x51 | 2.52e-08 | 3.23e-06 | 5.10e-07 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.10e-08 | 3.23e-06 | 7.74e-07 |
| synth_n100_p5 | bisquare | 100x6 | 1.29e-07 | 3.23e-06 | 1.39e-06 |
| synth_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_n5000_p20 | bisquare | 5000x21 | 4.39e-08 | 3.24e-06 | 4.73e-07 |
| synth_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.31e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.45e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 5 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.0 | 78.2 | 26.08x |
| classical_coleman | bisquare | 20x6 | 4.0 | 80.5 | 20.12x |
| classical_delivery | bisquare | 25x3 | 3.0 | 75.6 | 25.19x |
| classical_hbk | bisquare | 75x4 | 5.0 | 84.5 | 16.89x |
| classical_pension | bisquare | 18x2 | 2.0 | 72.1 | 36.05x |
| classical_phosphor | bisquare | 18x3 | 2.0 | 73.8 | 36.92x |
| classical_salinity | bisquare | 28x4 | 4.0 | 79.1 | 19.77x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 75.4 | 25.13x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 79.2 | 26.41x |
| classical_wood | bisquare | 20x6 | 4.0 | 82.2 | 20.54x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 75.3 | 25.11x |
| psi_ggw | ggw | 21x4 | 5.0 | 83.4 | 16.68x |
| psi_hampel | hampel | 21x4 | 3.0 | 75.6 | 25.19x |
| psi_lqq | lqq | 21x4 | 4.0 | 77.8 | 19.45x |
| psi_optimal | optimal | 21x4 | 3.0 | 74.6 | 24.87x |
| setting_KS2011_stackloss | lqq | 21x4 | 4.0 | 116.0 | 29.00x |
| setting_KS2014_stackloss | lqq | 21x4 | 7.0 | 115.3 | 16.48x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 476.0 | 621.8 | 1.31x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 63.0 | 173.9 | 2.76x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 594.0 | 841.7 | 1.42x |
| synth_ggw_n500_p10 | ggw | 500x11 | 80.0 | 231.9 | 2.90x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 559.0 | 696.8 | 1.25x |
| synth_hampel_n500_p10 | hampel | 500x11 | 62.0 | 189.6 | 3.06x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 544.0 | 906.0 | 1.67x |
| synth_lqq_n500_p10 | lqq | 500x11 | 62.0 | 241.7 | 3.90x |
| synth_n10000_p20 | bisquare | 10000x21 | 676.0 | 3891.9 | 5.76x |
| synth_n10000_p50 | bisquare | 10000x51 | 3219.0 | 19811.9 | 6.15x |
| synth_n1000_p10 | bisquare | 1000x11 | 91.0 | 304.8 | 3.35x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 94.6 | 11.82x |
| synth_n2000_p20 | bisquare | 2000x21 | 440.0 | 682.8 | 1.55x |
| synth_n5000_p20 | bisquare | 5000x21 | 510.0 | 1577.5 | 3.09x |
| synth_n500_p10 | bisquare | 500x11 | 48.0 | 183.2 | 3.82x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 526.0 | 755.0 | 1.44x |
| synth_optimal_n500_p10 | optimal | 500x11 | 56.0 | 204.4 | 3.65x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
