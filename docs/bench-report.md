# Benchmark report

Element-wise comparison between `pyrobustlm` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 5.06x, min 1.32x, max 12.41x

## Environment

- pyrobustlm: 0.5.1
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

| case | psi | n_x_p | R (ms) | py (ms) | py/R |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.0 | 23.8 | 7.92x |
| classical_coleman | bisquare | 20x6 | 3.0 | 22.7 | 7.57x |
| classical_delivery | bisquare | 25x3 | 3.0 | 22.1 | 7.37x |
| classical_hbk | bisquare | 75x4 | 5.0 | 30.4 | 6.07x |
| classical_pension | bisquare | 18x2 | 3.0 | 19.8 | 6.59x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 21.7 | 7.24x |
| classical_salinity | bisquare | 28x4 | 3.0 | 22.8 | 7.59x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 22.3 | 7.45x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 23.3 | 7.77x |
| classical_wood | bisquare | 20x6 | 3.0 | 24.7 | 8.23x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 23.5 | 7.85x |
| psi_ggw | ggw | 21x4 | 5.0 | 28.1 | 5.61x |
| psi_hampel | hampel | 21x4 | 3.0 | 24.0 | 8.01x |
| psi_lqq | lqq | 21x4 | 3.0 | 26.9 | 8.97x |
| psi_optimal | optimal | 21x4 | 3.0 | 23.2 | 7.72x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 62.0 | 12.41x |
| setting_KS2014_stackloss | lqq | 21x4 | 7.0 | 64.6 | 9.24x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 444.0 | 585.5 | 1.32x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 45.0 | 104.8 | 2.33x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 548.0 | 903.0 | 1.65x |
| synth_ggw_n500_p10 | ggw | 500x11 | 75.0 | 167.6 | 2.23x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 513.0 | 726.4 | 1.42x |
| synth_hampel_n500_p10 | hampel | 500x11 | 57.0 | 136.7 | 2.40x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 477.0 | 988.5 | 2.07x |
| synth_lqq_n500_p10 | lqq | 500x11 | 60.0 | 196.2 | 3.27x |
| synth_n10000_p20 | bisquare | 10000x21 | 662.0 | 2846.2 | 4.30x |
| synth_n10000_p50 | bisquare | 10000x51 | 3114.0 | 9616.0 | 3.09x |
| synth_n1000_p10 | bisquare | 1000x11 | 88.0 | 210.2 | 2.39x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 36.0 | 4.50x |
| synth_n2000_p20 | bisquare | 2000x21 | 449.0 | 627.5 | 1.40x |
| synth_n5000_p20 | bisquare | 5000x21 | 562.0 | 1332.9 | 2.37x |
| synth_n500_p10 | bisquare | 500x11 | 45.0 | 106.3 | 2.36x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 499.0 | 832.1 | 1.67x |
| synth_optimal_n500_p10 | optimal | 500x11 | 52.0 | 148.9 | 2.86x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
