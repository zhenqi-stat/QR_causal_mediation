# Recorded runtimes

The original simulation runs used 30 worker processes. Approximate runtimes in hours are summarized from the original logs. Method and replication coverage are recorded in `runtime_batches.csv`.

## Core estimation smoke test

After installing the required packages, run:

```bash
Rscript smoke_test.R
```

A local single-worker test on macOS arm64 with R 4.5.3 took **7.6 minutes**.

This runs Setting (iv) at `N = 600`, replication 1, using one worker. It generates the simulation data and fits MQR1, MQR2, QR1, QR2, TR, O-TR, O-TR', and Oracle through the original `simu_4.R` entrypoint. The tuning grids, tree counts, sample splitting, and finite-pool target are unchanged.

At completion, the script prints `PASS` after checking that all eight methods returned finite estimates and nonnegative variances and that the finite-pool target and Oracle calculation match the archived reference. It also prints the elapsed time. The run directory and absolute log-file path are printed at startup; method-by-method progress is recorded in `logs/progress.log`.

Results are saved under `output/smoke_test/simu_4/n_600/run_<timestamp>_pid<PID>/`. This directory contains the replication result, `results/results.RData`, `results/summary.csv`, and `results/smoke_test_status.txt`.

## Recorded 100-replication nonparametric batches

These batches compute QR1, QR2, TR, and O-TR. Times are hours for one configuration, using 30 workers. Full estimation also computes the parametric methods.

| Setting | N = 600 | N = 1000 | N = 1400 | N = 1800 |
| --- | --- | --- | --- | --- |
| (i) | 1.6 | 3.1 | 6.6 | 12.3 |
| (ii) | 1.1 | 2.1 | 3.2 | 7.4 |
| (iii) | 1.5 | 2.5 | 4.1 | 6.5 |
| (iv) | 0.6 | 1.1 | 1.8 | 2.8 |

| Setting (v), signal | N | Hours |
| --- | --- | --- |
| 0.25 | 1800 | 7.8 |
| 0.5 | 1800 | 7.7 |
| 2 | 1800 | 7.7 |

| Setting (vi), sparsity | N | Hours |
| --- | --- | --- |
| (4,12) | 1800 | 6.4 |
| (7,18) | 1800 | 7.0 |

## Configurations run in separate batches

| Configuration | Methods | Saved replications | Hours |
| --- | --- | --- | --- |
| Setting (v), signal = 1 | QR1, QR2, TR, O-TR | 1-60 | 6.8 |
| Setting (v), signal = 1 | QR1, QR2, TR, O-TR | 60-89 | 2.3 |
| Setting (v), signal = 1 | QR1, QR2, TR, O-TR | 90-100 | 1.3 |
| Setting (vi), (13,30) | MQR1, MQR2, O-TR', TR, O-TR | 1-30 | 1.3 |
| Setting (vi), (13,30) | MQR1, MQR2, O-TR', TR, O-TR | 31-100 | 3.8 |
| Setting (vi), (13,30) | QR1, QR2 | 31-100 | 3.9 |
| Setting (vi), (17,39) | MQR1, MQR2, O-TR', TR, O-TR | 1-30 | 1.5 |
| Setting (vi), (17,39) | MQR1, MQR2, O-TR', TR, O-TR | 31-100 | 4.3 |
| Setting (vi), (17,39) | QR1, QR2 | 1-30 | 1.8 |
| Setting (vi), (17,39) | QR1, QR2 | 31-100 | 4.7 |

For signal = 1, replication 60 is taken once from the first batch; the next contributes 61-89. The runs overlap in time, so their durations are not added as a single wall-clock runtime. For (13,30), the separate QR log for replications 1-30 was not retained.

## ACTG175 analysis

The QR2, TR, and O-TR batch took about 2 hours for 100 sample splits. The real-data implementation specifies 31 worker processes.

## Progress during a new run

Each analysis prints its run directory and absolute log-file path before loading packages or preparing data. Method-by-method progress is written to `output/<analysis>/<configuration>/run_<timestamp>_pid<PID>/logs/progress.log`. Replication files are saved under `replications/`; the final `results.RData` and `summary.csv` are saved under `results/`. The same layout is used for simulations and ACTG175.

Regenerating all eight numerical tables from the saved results takes less than one minute. Full estimation requires more than eight hours; the recorded batch timings above provide a guide for planning individual runs.
