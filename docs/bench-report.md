# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.87x, min 2.14x, max 15.54x
- Runtime ratio (py engine_c/R): median 1.66x, min 0.86x, max 8.26x

## Environment

- pylmrob: 0.5.24
- Python: 3.12.13
- Platform: Linux-6.17.0-1018-azure-x86_64-with-glibc2.39
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
| psi_optimal | optimal | 21x4 | 1.13e-15 | 5.80e-07 | 9.27e-14 |
| setting_KS2011_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| setting_KS2014_stackloss | lqq | 21x4 | 8.14e-04 | 9.58e-04 | 9.84e-04 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.49e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.66e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
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
| classical_aircraft | bisquare | 23x5 | 4.0 | 21.1 | 5.28x | 5.1 | 1.29x |
| classical_coleman | bisquare | 20x6 | 5.0 | 20.9 | 4.17x | 5.0 | 1.00x |
| classical_delivery | bisquare | 25x3 | 3.0 | 20.5 | 6.82x | 5.0 | 1.66x |
| classical_hbk | bisquare | 75x4 | 6.0 | 29.9 | 4.99x | 39.7 | 6.61x |
| classical_pension | bisquare | 18x2 | 3.0 | 17.8 | 5.92x | 3.9 | 1.31x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 18.5 | 6.18x | 4.0 | 1.32x |
| classical_salinity | bisquare | 28x4 | 4.0 | 21.7 | 5.44x | 5.5 | 1.37x |
| classical_stackloss | bisquare | 21x4 | 4.0 | 20.3 | 5.08x | 25.1 | 6.27x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 21.8 | 5.46x | 5.5 | 1.37x |
| classical_wood | bisquare | 20x6 | 4.0 | 20.7 | 5.18x | 5.0 | 1.24x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 20.0 | 6.67x | 24.8 | 8.26x |
| psi_ggw | ggw | 21x4 | 6.0 | 21.7 | 3.61x | 6.9 | 1.15x |
| psi_hampel | hampel | 21x4 | 3.0 | 20.5 | 6.83x | 5.3 | 1.77x |
| psi_lqq | lqq | 21x4 | 4.0 | 24.0 | 6.01x | 5.9 | 1.48x |
| psi_optimal | optimal | 21x4 | 4.0 | 21.2 | 5.30x | 4.5 | 1.13x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.0 | 77.7 | 15.54x | 6.9 | 1.38x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.0 | 77.8 | 9.72x | 6.9 | 0.86x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 395.0 | 978.7 | 2.48x | 653.2 | 1.65x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 43.0 | 144.7 | 3.36x | 77.1 | 1.79x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 592.0 | 1267.2 | 2.14x | 974.9 | 1.65x |
| synth_ggw_n500_p10 | ggw | 500x11 | 77.0 | 165.4 | 2.15x | 144.0 | 1.87x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 442.0 | 1067.6 | 2.42x | 866.1 | 1.96x |
| synth_hampel_n500_p10 | hampel | 500x11 | 51.0 | 163.3 | 3.20x | 130.8 | 2.56x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 479.0 | 1457.2 | 3.04x | 882.4 | 1.84x |
| synth_lqq_n500_p10 | lqq | 500x11 | 57.0 | 262.7 | 4.61x | 126.6 | 2.22x |
| synth_n10000_p20 | bisquare | 10000x21 | 494.0 | 3829.1 | 7.75x | 3047.1 | 6.17x |
| synth_n10000_p50 | bisquare | 10000x51 | 1892.0 | 8977.2 | 4.74x | 8148.0 | 4.31x |
| synth_n1000_p10 | bisquare | 1000x11 | 119.0 | 384.3 | 3.23x | 200.5 | 1.68x |
| synth_n100_p5 | bisquare | 100x6 | 9.0 | 38.6 | 4.29x | 14.2 | 1.58x |
| synth_n2000_p20 | bisquare | 2000x21 | 394.0 | 919.1 | 2.33x | 653.9 | 1.66x |
| synth_n5000_p20 | bisquare | 5000x21 | 393.0 | 2108.0 | 5.36x | 1646.5 | 4.19x |
| synth_n500_p10 | bisquare | 500x11 | 43.0 | 145.0 | 3.37x | 77.2 | 1.80x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 451.0 | 1172.4 | 2.60x | 721.9 | 1.60x |
| synth_optimal_n500_p10 | optimal | 500x11 | 50.0 | 202.5 | 4.05x | 100.5 | 2.01x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
