# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 5.48x, min 1.94x, max 16.81x
- Runtime ratio (py engine_c/R): median 1.77x, min 0.87x, max 8.79x

## Environment

- pylmrob: 0.5.22
- Python: 3.12.13
- Platform: Linux-6.17.0-1015-azure-x86_64-with-glibc2.39
- robustbase: 0.99.7
- R: R version 4.6.0 (2026-04-24)

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
| psi_optimal | optimal | 21x4 | 7.55e-16 | 5.80e-07 | 2.89e-13 |
| setting_KS2011_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| setting_KS2014_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.50e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.66e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.64e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 3.24e-08 | 3.24e-06 | 4.66e-07 |
| synth_n10000_p50 | bisquare | 10000x51 | 2.52e-08 | 3.23e-06 | 5.10e-07 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.10e-08 | 3.23e-06 | 7.74e-07 |
| synth_n100_p5 | bisquare | 100x6 | 1.29e-07 | 3.23e-06 | 1.39e-06 |
| synth_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_n5000_p20 | bisquare | 5000x21 | 4.39e-08 | 3.24e-06 | 4.73e-07 |
| synth_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.32e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.45e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 5 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 4.0 | 22.3 | 5.57x | 5.5 | 1.37x |
| classical_coleman | bisquare | 20x6 | 4.0 | 22.2 | 5.56x | 5.4 | 1.35x |
| classical_delivery | bisquare | 25x3 | 3.0 | 21.4 | 7.14x | 5.3 | 1.76x |
| classical_hbk | bisquare | 75x4 | 5.0 | 28.9 | 5.78x | 39.1 | 7.82x |
| classical_pension | bisquare | 18x2 | 2.0 | 18.9 | 9.46x | 4.3 | 2.13x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 19.9 | 6.62x | 4.2 | 1.42x |
| classical_salinity | bisquare | 28x4 | 4.0 | 22.6 | 5.65x | 5.8 | 1.44x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 21.6 | 7.20x | 26.4 | 8.79x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 22.6 | 5.65x | 5.9 | 1.49x |
| classical_wood | bisquare | 20x6 | 4.0 | 22.1 | 5.51x | 5.3 | 1.32x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 21.4 | 7.15x | 26.1 | 8.70x |
| psi_ggw | ggw | 21x4 | 6.0 | 23.5 | 3.92x | 7.2 | 1.19x |
| psi_hampel | hampel | 21x4 | 4.0 | 22.2 | 5.56x | 5.5 | 1.37x |
| psi_lqq | lqq | 21x4 | 4.0 | 24.5 | 6.13x | 6.0 | 1.51x |
| psi_optimal | optimal | 21x4 | 4.0 | 21.8 | 5.45x | 4.7 | 1.17x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 84.1 | 16.81x | 6.9 | 1.38x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.0 | 84.5 | 10.56x | 6.9 | 0.87x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 410.0 | 864.8 | 2.11x | 680.3 | 1.66x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.0 | 125.1 | 3.13x | 76.0 | 1.90x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 607.0 | 1177.3 | 1.94x | 991.5 | 1.63x |
| synth_ggw_n500_p10 | ggw | 500x11 | 71.0 | 168.5 | 2.37x | 145.2 | 2.05x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 449.0 | 1225.4 | 2.73x | 860.3 | 1.92x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.0 | 154.4 | 3.28x | 113.9 | 2.42x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 485.0 | 1427.7 | 2.94x | 944.7 | 1.95x |
| synth_lqq_n500_p10 | lqq | 500x11 | 51.0 | 222.9 | 4.37x | 121.1 | 2.38x |
| synth_n10000_p20 | bisquare | 10000x21 | 511.0 | 5366.2 | 10.50x | 4302.5 | 8.42x |
| synth_n10000_p50 | bisquare | 10000x51 | 2088.0 | 10479.9 | 5.02x | 9230.9 | 4.42x |
| synth_n1000_p10 | bisquare | 1000x11 | 127.0 | 406.5 | 3.20x | 226.9 | 1.79x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 36.3 | 4.54x | 14.6 | 1.82x |
| synth_n2000_p20 | bisquare | 2000x21 | 407.0 | 979.8 | 2.41x | 657.8 | 1.62x |
| synth_n5000_p20 | bisquare | 5000x21 | 398.0 | 2260.3 | 5.68x | 1721.2 | 4.32x |
| synth_n500_p10 | bisquare | 500x11 | 39.0 | 126.6 | 3.25x | 76.1 | 1.95x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 448.0 | 1023.4 | 2.28x | 705.1 | 1.57x |
| synth_optimal_n500_p10 | optimal | 500x11 | 44.0 | 166.3 | 3.78x | 87.4 | 1.99x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
