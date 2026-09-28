# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 3.07e-11, max 1.98e-03
- Scale relative error: median 4.16e-09, max 3.45e-03
- Cov diagonal max-rerr: median 2.02e-07, max 5.20e-01
- Runtime ratio (py/R): median 4.28x, min 1.43x, max 8.21x
- Runtime ratio (py engine_c/R): median 1.62x, min 0.82x, max 7.64x

## Environment

- pylmrob: 0.5.32
- Python: 3.12.14
- Platform: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39
- robustbase: 0.99.7
- R: R version 4.6.1 (2026-06-24)

## Numerical accuracy: max relative error vs R

| case | psi | n_x_p | max coef rerr | scale rerr | cov diag max rerr |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 1.04e-07 | 1.15e-07 | 2.78e-07 |
| classical_coleman | bisquare | 20x6 | 9.88e-08 | 4.57e-07 | 2.98e-07 |
| classical_delivery | bisquare | 25x3 | 5.02e-07 | 1.20e-05 | 2.12e-05 |
| classical_hbk | bisquare | 75x4 | 5.03e-11 | 2.51e-09 | 1.29e-08 |
| classical_pension | bisquare | 18x2 | 1.50e-06 | 2.53e-06 | 1.00e-05 |
| classical_phosphor | bisquare | 18x3 | 6.27e-08 | 4.29e-07 | 3.86e-07 |
| classical_salinity | bisquare | 28x4 | 4.50e-08 | 8.36e-07 | 6.46e-07 |
| classical_stackloss | bisquare | 21x4 | 5.68e-07 | 1.43e-06 | 4.16e-06 |
| classical_starsCYG | bisquare | 47x2 | 1.68e-12 | 4.21e-12 | 3.84e-07 |
| classical_wood | bisquare | 20x6 | 1.91e-09 | 2.12e-07 | 6.92e-08 |
| psi_bisquare | bisquare | 21x4 | 5.68e-07 | 1.43e-06 | 4.16e-06 |
| psi_ggw | ggw | 21x4 | 1.33e-03 | 2.51e-03 | 5.20e-01 |
| psi_hampel | hampel | 21x4 | 1.98e-03 | 3.45e-03 | 4.01e-01 |
| psi_lqq | lqq | 21x4 | 1.86e-06 | 3.69e-06 | 7.91e-06 |
| psi_optimal | optimal | 21x4 | 4.58e-16 | 5.80e-07 | 8.87e-14 |
| setting_KS2011_stackloss | lqq | 21x4 | 7.42e-08 | 8.73e-08 | 3.04e-07 |
| setting_KS2014_stackloss | lqq | 21x4 | 7.31e-08 | 8.61e-08 | 3.00e-07 |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 3.48e-13 | 1.78e-11 | 3.74e-08 |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 1.95e-12 | 8.82e-11 | 2.43e-07 |
| synth_ggw_n500_p10 | ggw | 500x11 | 1.50e-13 | 8.35e-12 | 1.46e-08 |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 2.75e-12 | 7.83e-11 | 1.62e-07 |
| synth_hampel_n500_p10 | hampel | 500x11 | 1.67e-12 | 2.71e-11 | 4.13e-07 |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 7.47e-12 | 3.30e-10 | 3.43e-07 |
| synth_lqq_n500_p10 | lqq | 500x11 | 2.66e-13 | 3.65e-12 | 1.01e-08 |
| synth_n10000_p20 | bisquare | 10000x21 | 5.85e-11 | 5.81e-09 | 1.98e-08 |
| synth_n10000_p50 | bisquare | 10000x51 | 1.12e-11 | 1.39e-09 | 2.22e-08 |
| synth_n1000_p10 | bisquare | 1000x11 | 3.90e-14 | 3.92e-12 | 3.39e-09 |
| synth_n100_p5 | bisquare | 100x6 | 2.14e-12 | 5.35e-11 | 9.68e-09 |
| synth_n2000_p20 | bisquare | 2000x21 | 3.48e-13 | 1.78e-11 | 3.74e-08 |
| synth_n5000_p20 | bisquare | 5000x21 | 1.13e-10 | 8.35e-09 | 9.11e-09 |
| synth_n500_p10 | bisquare | 500x11 | 1.04e-12 | 3.36e-11 | 7.76e-08 |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 1.35e-13 | 1.53e-12 | 1.11e-07 |
| synth_optimal_n500_p10 | optimal | 500x11 | 1.43e-13 | 3.95e-12 | 1.52e-07 |

## Runtime: median over 11 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R | py engine_c (ms) | py engine_c/R |
|---|---|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.8 | 20.8 | 5.48x | 5.5 | 1.45x |
| classical_coleman | bisquare | 20x6 | 4.0 | 20.9 | 5.25x | 5.4 | 1.34x |
| classical_delivery | bisquare | 25x3 | 3.0 | 19.7 | 6.53x | 5.3 | 1.75x |
| classical_hbk | bisquare | 75x4 | 5.2 | 24.7 | 4.80x | 10.1 | 1.96x |
| classical_pension | bisquare | 18x2 | 2.3 | 17.6 | 7.55x | 4.3 | 1.85x |
| classical_phosphor | bisquare | 18x3 | 2.8 | 18.5 | 6.57x | 4.3 | 1.53x |
| classical_salinity | bisquare | 28x4 | 3.7 | 20.7 | 5.63x | 5.8 | 1.56x |
| classical_stackloss | bisquare | 21x4 | 3.3 | 20.4 | 6.19x | 4.6 | 1.40x |
| classical_starsCYG | bisquare | 47x2 | 3.3 | 20.0 | 6.06x | 6.0 | 1.83x |
| classical_wood | bisquare | 20x6 | 4.1 | 20.7 | 5.00x | 5.2 | 1.26x |
| psi_bisquare | bisquare | 21x4 | 3.3 | 20.1 | 6.04x | 4.6 | 1.39x |
| psi_ggw | ggw | 21x4 | 5.6 | 23.2 | 4.15x | 7.2 | 1.29x |
| psi_hampel | hampel | 21x4 | 3.4 | 21.4 | 6.23x | 5.5 | 1.60x |
| psi_lqq | lqq | 21x4 | 4.0 | 21.7 | 5.41x | 6.0 | 1.50x |
| psi_optimal | optimal | 21x4 | 3.4 | 19.6 | 5.81x | 4.7 | 1.39x |
| setting_KS2011_stackloss | lqq | 21x4 | 5.2 | 22.9 | 4.40x | 7.0 | 1.35x |
| setting_KS2014_stackloss | lqq | 21x4 | 8.6 | 23.0 | 2.66x | 7.1 | 0.82x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 456.7 | 690.2 | 1.51x | 742.2 | 1.63x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 43.6 | 114.1 | 2.62x | 89.1 | 2.04x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 676.0 | 980.0 | 1.45x | 913.0 | 1.35x |
| synth_ggw_n500_p10 | ggw | 500x11 | 83.1 | 188.9 | 2.27x | 171.4 | 2.06x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 507.3 | 1119.4 | 2.21x | 950.7 | 1.87x |
| synth_hampel_n500_p10 | hampel | 500x11 | 51.7 | 162.0 | 3.13x | 139.6 | 2.70x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 547.6 | 961.8 | 1.76x | 776.0 | 1.42x |
| synth_lqq_n500_p10 | lqq | 500x11 | 59.3 | 170.2 | 2.87x | 140.2 | 2.36x |
| synth_n10000_p20 | bisquare | 10000x21 | 587.6 | 4823.2 | 8.21x | 4490.4 | 7.64x |
| synth_n10000_p50 | bisquare | 10000x51 | 2373.1 | 10118.8 | 4.26x | 9770.0 | 4.12x |
| synth_n1000_p10 | bisquare | 1000x11 | 144.2 | 269.2 | 1.87x | 268.5 | 1.86x |
| synth_n100_p5 | bisquare | 100x6 | 7.7 | 30.5 | 3.95x | 14.7 | 1.91x |
| synth_n2000_p20 | bisquare | 2000x21 | 474.5 | 677.0 | 1.43x | 751.1 | 1.58x |
| synth_n5000_p20 | bisquare | 5000x21 | 472.6 | 2031.0 | 4.30x | 1843.9 | 3.90x |
| synth_n500_p10 | bisquare | 500x11 | 46.0 | 114.0 | 2.48x | 88.9 | 1.93x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 494.2 | 880.6 | 1.78x | 796.2 | 1.61x |
| synth_optimal_n500_p10 | optimal | 500x11 | 48.9 | 129.8 | 2.66x | 105.5 | 2.16x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
