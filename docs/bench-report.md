# Benchmark report

Element-wise comparison between `pylmrob` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

## Headline (across 34 cases)

- Coefficient max-relative-error: median 5.45e-08, max 1.98e-03
- Scale relative error: median 3.23e-06, max 3.45e-03
- Cov diagonal max-rerr: median 7.95e-07, max 5.20e-01
- Runtime ratio (py/R): median 6.35x, min 1.35x, max 15.89x
- Runtime ratio (py engine_c/R): median 1.55x, min 0.75x, max 10.49x

## Environment

- pylmrob: 0.5.24
- Python: 3.12.13
- Platform: Linux-6.17.0-1020-azure-x86_64-with-glibc2.39
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
| classical_aircraft | bisquare | 23x5 | 4.0 | 27.8 | 6.95x | 4.9 | 1.22x |
| classical_coleman | bisquare | 20x6 | 4.0 | 27.9 | 6.98x | 4.8 | 1.20x |
| classical_delivery | bisquare | 25x3 | 3.0 | 26.9 | 8.96x | 4.7 | 1.56x |
| classical_hbk | bisquare | 75x4 | 5.0 | 35.7 | 7.14x | 44.6 | 8.92x |
| classical_pension | bisquare | 18x2 | 2.0 | 23.6 | 11.82x | 3.7 | 1.87x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 24.7 | 8.25x | 3.8 | 1.26x |
| classical_salinity | bisquare | 28x4 | 4.0 | 28.2 | 7.04x | 5.2 | 1.30x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 27.7 | 9.22x | 31.3 | 10.44x |
| classical_starsCYG | bisquare | 47x2 | 4.0 | 27.6 | 6.91x | 5.2 | 1.30x |
| classical_wood | bisquare | 20x6 | 4.0 | 27.6 | 6.90x | 4.7 | 1.18x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 26.9 | 8.98x | 31.5 | 10.49x |
| psi_ggw | ggw | 21x4 | 5.0 | 29.0 | 5.80x | 6.7 | 1.33x |
| psi_hampel | hampel | 21x4 | 3.0 | 27.7 | 9.24x | 5.0 | 1.67x |
| psi_lqq | lqq | 21x4 | 4.0 | 30.7 | 7.66x | 5.6 | 1.40x |
| psi_optimal | optimal | 21x4 | 4.0 | 27.7 | 6.92x | 4.2 | 1.06x |
| setting_KS2011_stackloss | lqq | 21x4 | 6.0 | 95.4 | 15.89x | 6.7 | 1.11x |
| setting_KS2014_stackloss | lqq | 21x4 | 9.0 | 95.7 | 10.63x | 6.7 | 0.75x |
| synth_bisquare_n2000_p20 | bisquare | 2000x21 | 391.0 | 856.0 | 2.19x | 490.0 | 1.25x |
| synth_bisquare_n500_p10 | bisquare | 500x11 | 40.0 | 140.4 | 3.51x | 74.5 | 1.86x |
| synth_ggw_n2000_p20 | ggw | 2000x21 | 624.0 | 842.7 | 1.35x | 729.2 | 1.17x |
| synth_ggw_n500_p10 | ggw | 500x11 | 74.0 | 166.7 | 2.25x | 138.3 | 1.87x |
| synth_hampel_n2000_p20 | hampel | 2000x21 | 444.0 | 978.0 | 2.20x | 680.5 | 1.53x |
| synth_hampel_n500_p10 | hampel | 500x11 | 47.0 | 162.3 | 3.45x | 120.5 | 2.56x |
| synth_lqq_n2000_p20 | lqq | 2000x21 | 481.0 | 1206.0 | 2.51x | 670.6 | 1.39x |
| synth_lqq_n500_p10 | lqq | 500x11 | 52.0 | 243.8 | 4.69x | 118.9 | 2.29x |
| synth_n10000_p20 | bisquare | 10000x21 | 502.0 | 3712.4 | 7.40x | 3015.3 | 6.01x |
| synth_n10000_p50 | bisquare | 10000x51 | 1871.0 | 8083.9 | 4.32x | 7254.7 | 3.88x |
| synth_n1000_p10 | bisquare | 1000x11 | 119.0 | 313.3 | 2.63x | 205.9 | 1.73x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 43.8 | 5.48x | 13.4 | 1.68x |
| synth_n2000_p20 | bisquare | 2000x21 | 393.0 | 771.0 | 1.96x | 586.9 | 1.49x |
| synth_n5000_p20 | bisquare | 5000x21 | 389.0 | 1930.1 | 4.96x | 1488.2 | 3.83x |
| synth_n500_p10 | bisquare | 500x11 | 40.0 | 140.8 | 3.52x | 74.5 | 1.86x |
| synth_optimal_n2000_p20 | optimal | 2000x21 | 444.0 | 975.9 | 2.20x | 551.8 | 1.24x |
| synth_optimal_n500_p10 | optimal | 500x11 | 46.0 | 193.1 | 4.20x | 90.4 | 1.97x |

## Coverage

- Cases in both: 34
- Only R: (none)
- Only py: (none)
