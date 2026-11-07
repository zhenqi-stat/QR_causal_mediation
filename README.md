# Reproducibility code

Code and final computational results for *Quadruply robust methods for causal mediation analysis*.

## Experiments and numerical tables

| Script | Experiment | Table | Configurations |
| --- | --- | --- | --- |
| `simu_1.R` | Setting (i) | Main Table 4 | `N = 600, 1000, 1400, 1800` |
| `simu_2.R` | Setting (ii) | Main Table 5 | `N = 600, 1000, 1400, 1800` |
| `simu_3.R` | Setting (iii), correctly specified parametric models | Supplement Table 7 | `N = 600, 1000, 1400, 1800` |
| `real.R` | ACTG175 analysis | Main Table 6 | 100 sample splits |
| `simu_4.R` | Setting (iv), binary mediator | Supplement Table 8 | `N = 600, 1000, 1400, 1800` |
| `simu_5.R` | Setting (v), varying nonlinear signal | Supplement Table 9 | `signal = 0.25, 0.5, 1, 2`; `N = 1800` |
| `simu_6.R` | Setting (vi), varying sparsity | Supplement Table 10 | `s_level = 1, 3, 5, 6`; `N = 1800` |

The four sparsity levels correspond to `(s_pi, s_q) = (4,12), (7,18), (13,30), (17,39)`. `table_map.csv` gives the same mapping in machine-readable form. The other manuscript tables summarize theoretical properties, and the figures are schematic diagrams.

## Installation and data

Set the R working directory to this folder and run:

```r
source("install_packages.R")
```

The required packages are listed in `install_packages.R`; recorded versions are in `package_versions.csv`. ACTG175 is distributed with `speff2trial` 1.0.5 at https://CRAN.R-project.org/package=speff2trial. It is loaded automatically by the real-data script. `DATA_DICTIONARY.md` documents the analysis variables, participant selection, transformations, and result-file objects.

## Regenerate the numerical tables

```r
source("reproduce_tables.R")
```

or, from a terminal:

```bash
Rscript run_all.R
```

This command rebuilds the seven numerical tables from the 25 final result files in `results/`. It takes less than one minute and requires base R. Unrounded CSVs, formatted CSVs, and LaTeX rows are written to `tables/`; filenames follow the manuscript table numbers, such as `main_table_4_display.csv`.

## Rerun estimation

Select the configuration, replication IDs, and worker count near the top of an entrypoint, then run it, for example:

```r
source("simu_1.R")
```

Each experiment script is self-contained and can also be run directly with `Rscript`. To run replications 1--30, set `replication_ids <- 1:30`; keep `tmax = 100`.

Run the full collection of 24 simulation configurations and the real-data analysis with:

```bash
Rscript run_all.R --fit
```

The default is `1:100`. Simulations report Oracle, TR, O-TR, QR1, QR2, O-TR', MQR1, and MQR2. Main Table 6 displays TR, O-TR, QR2, O-TR', and MQR2. The real-data analysis also computes and saves QR1.

Simulations use 30 worker processes. The real-data implementation specifies 31 workers. Its QR2, TR, and O-TR batch took about 2 hours for 100 sample splits. Approximate runtimes and batch coverage are summarized in `RUNTIME.md`.

### Progress and output files

At startup, the console prints the run directory and log-file path:

```text
Simulation begun with 30 parallel workers
Log file: /path/to/output/<run>/logs/Log_simu_1_N600_<timestamp>_pid<pid>.txt
```

Follow estimation progress in the log file. It records package messages, warnings, errors, and method completion.

Each completed run writes `Final_*.RData` and a summary CSV to its run directory. Logs are stored in its `logs/` subfolder.

## Summary definitions

Simulation summaries use:

```r
theta <- mean(tau_all)
bias <- apply(pred_theta - theta, 2, median)
RMSE <- sqrt(apply((pred_theta - theta)^2, 2, mean))
```

AL is the median 95% confidence-interval length; AC is empirical coverage. Real-data results aggregate estimates over 100 sample splits. Full definitions are in `DATA_DICTIONARY.md` and `tools/verify_results.R`.

## Computing environment

| Component | Configuration |
| --- | --- |
| OS | Ubuntu 20.04.4 LTS |
| Processor | AMD EPYC 9654 |
| R | 4.5.2 |
| ranger | 0.18.0 |
| randomForestSRC | 3.9.0 |
| caret | 7.0-1 |
| glmnet | 5.1 |
| BLAS | OpenBLAS 0.3.8 |

Additional package versions, compiler information, and library-linkage details are provided in `environment.txt` and `package_versions.csv`. `SOFTWARE_REFERENCES.bib` contains references for the principal packages.

## Structure and verification

`configurations/` contains the configuration-specific implementations, and `configuration_index.rds` records their parameters. `source_manifest.csv` and result metadata document the sources and replications. `tools/` contains execution and verification utilities. `results/` holds 24 simulation files with 100 replications each and `real.RData` with 100 sample splits.

The seven analysis entrypoints contain their complete implementations. The configuration-specific files provide the corresponding DGP definitions used by the verification utilities.

```bash
Rscript verify.R
Rscript verify.R --dgp
```

`verify.R` checks source files, result completeness, and summary calculations. The `--dgp` option also regenerates the 24 DGP configurations and checks their finite-pool targets and 2,400 Oracle estimates and variances. Use `run_all.R --fit` for full estimation.

## License

The supplied code is released under the MIT License; see `LICENSE`. Third-party R packages and the ACTG175 data retain their respective terms.
