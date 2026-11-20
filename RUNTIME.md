# Recorded runtimes

The simulation entrypoints default to 30 worker processes. `runtime_batches.csv` lists all 135 source batches contributing to the final results across seven simulation settings and the ACTG175 analysis. Each row records the source filename and MD5, retained methods and replications, and available timing and worker information.

`elapsed_hours` is measured from matched run logs. `snapshot_elapsed_hours` is the interval from a recorded run-start timestamp to the saved result snapshot. Times are rounded to tenths of an hour; `NA` denotes an unrecorded measurement.

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

These source batches compute QR1, QR2, TR, and O-TR using 30 workers. Times are hours for one configuration. Some methods were subsequently rerun in separate batches; `runtime_batches.csv` identifies the exact contributions to the final tables.

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

For signal = 1, replication 60 is taken once from the first batch; the next contributes 61-89. The runs overlap in time, so their durations are not added as a single wall-clock runtime.

## Setting (vii)

The final results combine the following saved batches. The times below are snapshot ages, as defined above.

| N | Methods | Replications | Snapshot age (hours) |
| --- | --- | --- | --- |
| 1000 | All eight methods | 1-100 | 2.4 |
| 1800 | MQR1, MQR2, O-TR', Oracle | 1-30 | 2.5 |
| 1800 | QR1, QR2, TR, O-TR, Oracle | 1-30 | 4.1 |
| 1800 | All eight methods | 31-100 | 15.0 |
| 2500 | All eight methods | 1-50 | 4.3 |
| 2500 | All eight methods | 51-100 | 4.1 |
| 5000 | MQR1, MQR2, O-TR', Oracle | 1-30 | 6.7 |
| 5000 | MQR1, MQR2, O-TR', Oracle | 31-100 | 16.1 |
| 5000 | QR1, QR2, TR, Oracle | 1-30 | 27.7 |
| 5000 | O-TR, Oracle | 1-30 | 0.2 |
| 5000 | QR1, QR2, TR, O-TR, Oracle | 31-60 | 27.1 |
| 5000 | QR1, QR2, TR, O-TR, Oracle | 61-100 | 41.3 |

The current script uses 30 workers. The source snapshots do not record their worker counts.

## ACTG175 analysis

The final MQR2, O-TR', and QR2 effect batch records 30 workers. The corrected TR/O-TR KDE batch records 50 workers. The `real.R` entrypoint defaults to 31 workers. Earlier ACTG175 logs provide a rough runtime reference of 2-3 hours for 100 sample splits; timings and worker records for the final source files are listed in `runtime_batches.csv`.

## Progress during a new run

Each analysis prints its run directory and absolute log-file path before loading packages or preparing data. Method-by-method progress is written to `output/<analysis>/<configuration>/run_<timestamp>_pid<PID>/logs/progress.log`. Replication files are saved under `replications/`; the final `results.RData` and `summary.csv` are saved under `results/`. The same layout is used for simulations and ACTG175.

Regenerating all eight numerical tables from the saved results takes less than one minute. Full estimation requires more than eight hours; the recorded batch timings above provide a guide for planning individual runs.
