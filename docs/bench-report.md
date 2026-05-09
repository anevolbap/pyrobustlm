# Benchmark report

Element-wise comparison between `pyrobustlm` and `robustbase::lmrob` on a fixed corpus of fits. Re-generate with::

    Rscript scripts/benchmark.R
    python  scripts/benchmark.py
    python  scripts/build_bench_report.py

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
| synth_n1000_p10 | bisquare | 1000x11 | 3.10e-08 | 3.23e-06 | 7.74e-07 |
| synth_n100_p5 | bisquare | 100x6 | 1.29e-07 | 3.23e-06 | 1.39e-06 |
| synth_n2000_p20 | bisquare | 2000x21 | 5.45e-08 | 3.23e-06 | 6.16e-07 |
| synth_n5000_p20 | bisquare | 5000x21 | 4.39e-08 | 3.24e-06 | 4.73e-07 |
| synth_n500_p10 | bisquare | 500x11 | 4.91e-08 | 3.23e-06 | 7.95e-07 |

## Runtime: median over 5 reps (lower is better)

| case | psi | n_x_p | R (ms) | py (ms) | py/R |
|---|---|---|---|---|---|
| classical_aircraft | bisquare | 23x5 | 3.0 | 81.2 | 27.07x |
| classical_coleman | bisquare | 20x6 | 4.0 | 81.3 | 20.32x |
| classical_delivery | bisquare | 25x3 | 3.0 | 77.6 | 25.86x |
| classical_hbk | bisquare | 75x4 | 5.0 | 89.6 | 17.92x |
| classical_pension | bisquare | 18x2 | 2.0 | 74.4 | 37.18x |
| classical_phosphor | bisquare | 18x3 | 3.0 | 77.8 | 25.93x |
| classical_salinity | bisquare | 28x4 | 3.0 | 83.0 | 27.65x |
| classical_stackloss | bisquare | 21x4 | 3.0 | 77.4 | 25.79x |
| classical_starsCYG | bisquare | 47x2 | 3.0 | 83.7 | 27.89x |
| classical_wood | bisquare | 20x6 | 3.0 | 86.1 | 28.69x |
| psi_bisquare | bisquare | 21x4 | 3.0 | 78.8 | 26.26x |
| psi_ggw | ggw | 21x4 | 5.0 | 89.2 | 17.85x |
| psi_hampel | hampel | 21x4 | 3.0 | 82.6 | 27.53x |
| psi_lqq | lqq | 21x4 | 3.0 | 82.4 | 27.47x |
| psi_optimal | optimal | 21x4 | 3.0 | 80.3 | 26.75x |
| synth_n1000_p10 | bisquare | 1000x11 | 95.0 | 327.1 | 3.44x |
| synth_n100_p5 | bisquare | 100x6 | 8.0 | 109.8 | 13.72x |
| synth_n2000_p20 | bisquare | 2000x21 | 430.0 | 704.4 | 1.64x |
| synth_n5000_p20 | bisquare | 5000x21 | 511.0 | 1599.1 | 3.13x |
| synth_n500_p10 | bisquare | 500x11 | 47.0 | 198.0 | 4.21x |

## Coverage

- Cases in both: 20
- Only R: (none)
- Only py: (none)
