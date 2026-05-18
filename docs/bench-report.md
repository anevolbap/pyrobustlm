# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 2.74x, min 0.87x, max 13.79x
- Runtime ratio (py engine_c/R): median 0.93x, min 0.32x, max 8.91x

## Environment

- pylmrob: 0.5.10
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
| classical_aircraft | bisquare | 23x5 | 5.0 | 25.2 | 5.03x | 4.8 | 0.95x |
| classical_coleman | bisquare | 20x6 | 5.0 | 26.7 | 5.33x | 5.3 | 1.06x |
| classical_delivery | bisquare | 25x3 | 4.0 | 23.9 | 5.98x | 3.8 | 0.95x |
| classical_hbk | bisquare | 75x4 | 6.0 | 32.6 | 5.43x | 38.2 | 6.36x |
| classical_pension | bisquare | 18x2 | 3.0 | 20.8 | 6.93x | 3.4 | 1.13x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 21.2 | 7.06x | 3.2 | 1.06x |
| classical_salinity | bisquare | 28x4 | 5.0 | 24.5 | 4.90x | 5.4 | 1.07x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 25.1 | 6.28x | 29.6 | 7.39x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 25.3 | 8.42x | 4.4 | 1.45x |
| classical_wood | bisquare | 20x6 | 5.0 | 25.3 | 5.07x | 3.7 | 0.73x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 24.7 | 8.23x | 26.7 | 8.91x |
| psi_ggw | ggw | 21x4 | 11.0 | 30.9 | 2.81x | 6.3 | 0.57x |
| psi_hampel | hampel | 21x4 | 6.0 | 25.8 | 4.30x | 4.6 | 0.77x |
| psi_lqq | lqq | 21x4 | 5.0 | 28.6 | 5.71x | 4.6 | 0.92x |
| psi_optimal | optimal | 21x4 | 5.0 | 26.0 | 5.19x | 3.6 | 0.72x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 69.0 | 13.79x | 5.3 | 1.06x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.0 | 68.4 | 8.55x | 5.3 | 0.66x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 559.0 | 484.5 | 0.87x | 210.8 | 0.38x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 70.0 | 109.7 | 1.57x | 70.6 | 1.01x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 696.0 | 741.3 | 1.07x | 267.5 | 0.38x |
| synth_ggw_n500_p10 | ggw | 500x11 | 84.0 | 176.2 | 2.10x | 152.3 | 1.81x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 588.0 | 584.5 | 0.99x | 259.1 | 0.44x |
| synth_hampel_n500_p10 | hampel | 500x11 | 64.0 | 129.3 | 2.02x | 110.6 | 1.73x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 634.0 | 812.6 | 1.28x | 296.3 | 0.47x |
| synth_lqq_n500_p10 | lqq | 500x11 | 70.0 | 187.4 | 2.68x | 119.6 | 1.71x |
| synth_n10000_p20 | bisquare | 10000x21 | 1005.0 | 2488.2 | 2.48x | 866.4 | 0.86x |
| synth_n10000_p50 | bisquare | 10000x51 | 4750.0 | 6381.6 | 1.34x | 3165.3 | 0.67x |
| synth_n1000_p10 | bisquare | 1000x11 | 128.0 | 206.7 | 1.62x | 79.7 | 0.62x |
| synth_n100_p5 | bisquare | 100x6 | 15.0 | 38.1 | 2.54x | 11.5 | 0.76x |
| synth_n2000_p20 | bisquare | 2000x21 | 576.0 | 499.5 | 0.87x | 185.5 | 0.32x |
| synth_n5000_p20 | bisquare | 5000x21 | 1005.0 | 1285.4 | 1.28x | 494.8 | 0.49x |
| synth_n500_p10 | bisquare | 500x11 | 66.0 | 116.1 | 1.76x | 61.8 | 0.94x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 575.0 | 633.7 | 1.10x | 243.3 | 0.42x |
| synth_optimal_n500_p10 | optimal | 500x11 | 69.0 | 146.2 | 2.12x | 83.8 | 1.21x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
