# Benchmark report

Element-wise comparison between `pyrobustlm` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 6.28x, min 1.37x, max 21.26x

## Environment

- pyrobustlm: 0.4.0
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
| psi_optimal | optimal | 21x4 | 1.32e-15 | 5.80e-07 | 1.24e-13 |
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
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.30e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.43e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 5 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 4.0 | 30.4 | 7.59x |
| classical_coleman | bisquare | 20x6 | 4.0 | 29.5 | 7.37x |
| classical_delivery | bisquare | 25x3 | 3.0 | 28.9 | 9.64x |
| classical_hbk | bisquare | 75x4 | 5.0 | 37.7 | 7.54x |
| classical_pension | bisquare | 18x2 | 3.0 | 25.6 | 8.52x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 30.9 | 10.30x |
| classical_salinity | bisquare | 28x4 | 3.0 | 29.3 | 9.75x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 28.0 | 9.32x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 31.4 | 10.46x |
| classical_wood | bisquare | 20x6 | 4.0 | 32.2 | 8.04x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 29.8 | 9.92x |
| psi_ggw | ggw | 21x4 | 5.0 | 35.3 | 7.05x |
| psi_hampel | hampel | 21x4 | 3.0 | 30.6 | 10.21x |
| psi_lqq | lqq | 21x4 | 4.0 | 33.9 | 8.47x |
| psi_optimal | optimal | 21x4 | 3.0 | 30.5 | 10.16x |
| setting_KS2011_stackloss | lqq | 21x4 | 4.0 | 85.1 | 21.26x |
| setting_KS2014_stackloss | lqq | 21x4 | 7.0 | 100.5 | 14.36x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 431.0 | 633.8 | 1.47x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 47.0 | 141.7 | 3.01x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 633.0 | 947.9 | 1.50x |
| synth_ggw_n500_p10 | ggw | 500x11 | 79.0 | 190.9 | 2.42x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 521.0 | 714.7 | 1.37x |
| synth_hampel_n500_p10 | hampel | 500x11 | 54.0 | 130.5 | 2.42x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 560.0 | 1008.1 | 1.80x |
| synth_lqq_n500_p10 | lqq | 500x11 | 58.0 | 197.5 | 3.40x |
| synth_n10000_p20 | bisquare | 10000x21 | 743.0 | 4087.6 | 5.50x |
| synth_n10000_p50 | bisquare | 10000x51 | 3195.0 | 10727.1 | 3.36x |
| synth_n1000_p10 | bisquare | 1000x11 | 93.0 | 280.8 | 3.02x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 42.0 | 5.25x |
| synth_n2000_p20 | bisquare | 2000x21 | 454.0 | 917.5 | 2.02x |
| synth_n5000_p20 | bisquare | 5000x21 | 516.0 | 1999.3 | 3.87x |
| synth_n500_p10 | bisquare | 500x11 | 50.0 | 127.0 | 2.54x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 466.0 | 790.4 | 1.70x |
| synth_optimal_n500_p10 | optimal | 500x11 | 51.0 | 147.4 | 2.89x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
