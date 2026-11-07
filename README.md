# Reproducibility code

Code and final computational results for *Quadruply robust methods for causal mediation analysis*.

## Experiments and numerical tables

| Script | Experiment | Table | Configurations |
| --- | --- | --- | --- |
| `simu_1.R` | Setting (i) | Main Table 4 | `N = 600, 1000, 1400, 1800` |
| `simu_2.R` | Setting (ii) | Main Table 5 | `N = 600, 1000, 1400, 1800` |
| `simu_3.R` | Setting (iii), correctly specified parametric models | Supplement Table 7 | `N = 600, 1000, 1400, 1800` |
| `simu_4.R` | Setting (iv), binary mediator | Supplement Table 8 | `N = 600, 1000, 1400, 1800` |
| `simu_5.R` | Setting (v), varying nonlinear signal | Supplement Table 9 | `signal = 0.25, 0.5, 1, 2`; `N = 1800` |
| `simu_6.R` | Setting (vi), varying sparsity | Supplement Table 10 | `ell = 1, 2, 3, 4`; `N = 1800` |
| `real.R` | ACTG175 analysis | Main Table 6 | 100 sample splits |

`table_map.csv` lists the scripts, configurations, and corresponding tables.

## Computing environment

QR1, QR2, O-TR, and TR use the ranger package for nonparametric modeling. Repeated runs may produce small numerical differences due to differences in hardware, compilers and mathematical libraries, parallel scheduling, and floating-point arithmetic.

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

## Installation and data

Set the R working directory to this folder and run:

```r
source("install_packages.R")
```

The required R packages are listed in `install_packages.R`. ACTG175 is distributed with `speff2trial` 1.0.5 at https://CRAN.R-project.org/package=speff2trial. The real-data script loads it automatically. `DATA_DICTIONARY.md` documents the analysis variables, participant selection, transformations, and result-file objects.

## Regenerate the numerical tables

```r
source("reproduce_tables.R")
```

or, from a terminal:

```bash
Rscript run_all.R
```

This recalculates the seven numerical tables from the 25 final result files in `results/` and writes unrounded CSVs, formatted CSVs, and LaTeX rows to `tables/`. It takes less than one minute and requires only base R.

## Rerun estimation

Set the configuration, replication IDs, and worker count at the top of the corresponding script, then run it, for example:

```r
source("simu_1.R")
```

For a short test, set `N <- 600L`, `replication_ids <- 1:2`, and `workers <- 2L` in `simu_1.R`.

Run all 24 simulation configurations and the real-data analysis with:

```bash
Rscript run_all.R --fit
```

The default is 100 replications or sample splits. New fits are saved under `output/`. Simulation results report Oracle, TR, O-TR, QR1, QR2, O-TR', MQR1, and MQR2. The real-data table uses TR, O-TR, QR2, O-TR', and MQR2. The real-data script also computes QR1.

The default worker counts are 30 for simulations and 31 for real data. `RUNTIME.md` gives approximate runtimes and the methods included in each recorded batch.

### Progress and output files

Every analysis prints its run directory and absolute log-file path before loading packages or preparing data. Method-by-method progress is written to the log file. All seven analysis scripts use the same directory structure:

```text
output/<analysis>/<configuration>/run_<YYYYMMDD_HHMMSS>_pid<PID>/
  run_manifest.rds
  logs/
    progress.log
  replications/
    rep_001.rds
    rep_002.rds
    ...
  results/
    results.RData
    summary.csv
```

The analysis directory is `simu_1` through `simu_6`, or `real`. Configuration directories use `n_600`, `signal_0p25`, `ell_1`, or `dataset_actg175`, as appropriate. A numeric suffix distinguishes runs started in the same second by the same process. Replication files use the requested replication IDs. `run_manifest.rds` records the configuration, replication IDs, worker count, and output paths.

## Summary definitions

Simulation summaries use:

```r
theta <- mean(tau_all)
bias <- apply(pred_theta - theta, 2, median)
RMSE <- sqrt(apply((pred_theta - theta)^2, 2, mean))
```

AL is the median 95% confidence-interval length; AC is empirical coverage. Real-data results aggregate estimates over 100 sample splits. Full definitions are in `DATA_DICTIONARY.md` and `tools/verify_results.R`.

## Structure and verification

The seven analysis scripts are self-contained. `configurations/` supplies the corresponding DGP definitions for verification. `configuration_index.rds` lists configuration parameters; `source_manifest.csv` and result metadata identify the sources and replications. `tools/` contains execution and verification utilities. The bundled `results/` directory contains 24 simulation files with 100 replications and eight methods each, plus `real.RData` with 100 splits for six methods.

```bash
Rscript verify.R
Rscript verify.R --dgp
```

The first command checks the source files, result completeness, and summary calculations. The second also regenerates the 24 DGP configurations and checks each finite-pool target and its 100 Oracle estimates and variances.

## License

The supplied code is released under the MIT License; see `LICENSE`. Third-party R packages and the ACTG175 data retain their respective terms.
