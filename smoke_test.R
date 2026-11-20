# Fit one replication of Setting (iv), N=600, with the complete estimator code.
# Run from a terminal: Rscript smoke_test.R
.smoke_root <- function() {
  paths <- vapply(sys.frames(), function(frame) {
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) paths <- sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(paths)) stop("Run with Rscript or source().")
  dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
}

run_smoke_test <- function(archive_root) {
  source(file.path(archive_root, "tools", "archive_runtime.R"), local = TRUE)
  started <- proc.time()[["elapsed"]]
  cat("Smoke test: Setting (iv), N=600, replication 1, one worker.\n",
      "Fits MQR1, MQR2, QR1, QR2, TR, O-TR, O-TR', and Oracle.\n",
      "The analysis will print its log path; follow method progress there.\n", sep = "")
  flush.console()
  result_path <- run_archive(group = 4L, config = 600L, replication_ids = 1L,
    archive_root = archive_root, workers = 1L,
    output_dir = file.path(archive_root, "output", "smoke_test"))
  elapsed <- proc.time()[["elapsed"]] - started
  fitted <- new.env(parent = baseenv())
  load(result_path, fitted)
  reference <- new.env(parent = baseenv())
  load(file.path(archive_root, "results", "simu_4_N600.RData"), reference)
  methods <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")
  stopifnot(identical(fitted$source_id, "simu_4_N600"), fitted$N == 600,
    identical(fitted$replication_ids, 1L), identical(fitted$method_names, methods),
    identical(dim(fitted$pred_theta), c(1L, 8L)),
    identical(dim(fitted$pred_sig2), c(1L, 8L)),
    all(is.finite(fitted$pred_theta)), all(is.finite(fitted$pred_sig2)),
    all(fitted$pred_sig2 >= 0), abs(fitted$theta - reference$theta) < 1e-10,
    abs(fitted$pred_theta[1L, "Oracle"] - reference$pred_theta[1L, "Oracle"]) < 1e-10,
    abs(fitted$pred_sig2[1L, "Oracle"] - reference$pred_sig2[1L, "Oracle"]) < 1e-8)
  estimates <- data.frame(method = methods, estimate = as.numeric(fitted$pred_theta[1L, ]),
    SE = sqrt(as.numeric(fitted$pred_sig2[1L, ]) / fitted$N), row.names = NULL)
  print(estimates, row.names = FALSE)
  status <- c("PASS: all eight methods fitted and returned finite estimates and variances.",
    "PASS: the finite-pool target and Oracle calculation match the archived reference.",
    sprintf("Elapsed time: %.1f minutes (%.0f seconds).", elapsed / 60, elapsed),
    paste0("Result file: ", result_path))
  writeLines(status, file.path(dirname(result_path), "smoke_test_status.txt"))
  cat(paste(status, collapse = "\n"), "\n")
  invisible(result_path)
}

run_smoke_test(.smoke_root())
