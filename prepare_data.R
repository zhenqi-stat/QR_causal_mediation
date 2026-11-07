# Rebuild the ACTG175 analysis sample and covariate bases from the supplied CSV.
.data_root <- function() {
  paths <- vapply(sys.frames(), function(frame) {
    if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!length(paths)) stop("Run with source() or Rscript.")
  dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
}

# Retain full numeric precision in the analysis CSV.
.write_data_csv <- function(data, file) {
  encoded <- lapply(data, function(x) {
    if (is.numeric(x)) {
      y <- if (is.integer(x)) as.character(x) else sprintf("%.17g", x)
      y[is.na(x)] <- "NA"
      y
    } else as.character(x)
  })
  write.table(as.data.frame(encoded, check.names = FALSE), file,
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
}

prepare_actg175 <- function(archive_root,
                           output_directory = file.path(archive_root, "output", "data_preparation")) {
  data_directory <- file.path(archive_root, "data")
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)
  columns <- c("pidnum", "age", "wtkg", "hemo", "homo", "drugs", "karnof",
    "oprior", "z30", "zprior", "preanti", "race", "gender", "str2", "strat",
    "symptom", "treat", "offtrt", "cd40", "cd420", "cd496", "r", "cd80",
    "cd820", "cens", "days", "arms")
  column_classes <- setNames(rep("integer", length(columns)), columns)
  column_classes["wtkg"] <- "numeric"
  raw <- read.csv(file.path(data_directory, "ACTG175_raw.csv"),
    colClasses = column_classes, check.names = FALSE, na.strings = "NA")
  stopifnot(identical(names(raw), columns), identical(dim(raw), c(2139L, 27L)))
  source_rows <- which(raw$offtrt == 0 & raw$r == 1 & raw$arms != 3)
  df <- raw[source_rows, , drop = FALSE]
  X_conti <- c("age", "wtkg", "karnof", "cd40", "cd80")
  X_binar <- c("hemo", "homo", "drugs", "race", "gender", "str2", "symptom")
  M_vars <- c("cd420", "cd820")
  X_all <- df[, c(X_conti, X_binar)]
  M_all <- log(df[, M_vars])
  Y_all <- log(df$cd496)
  A_all <- df$treat
  stopifnot(nrow(df) == 780L, sum(A_all == 0L) == 253L, sum(A_all == 1L) == 527L,
    !anyNA(cbind(X_all, M_all, Y_all, A_all)),
    all(is.finite(as.matrix(cbind(X_all, M_all, Y_all, A_all)))))
  functions <- new.env(parent = asNamespace("stats"))
  sys.source(file.path(archive_root, "configurations", "real_2_data_generate.R"),
             envir = functions)
  captured <- capture.output(expanded <- functions$expand_covariates_with_mediation(
    data = df, X_conti = X_conti, X_binar = X_binar, M_conti = M_vars,
    add_conti_binar_interactions = TRUE, remove_zero_var = TRUE))
  df_log <- df
  df_log[, M_vars] <- M_all
  captured <- capture.output(qr_expanded <- functions$expand_covariates_with_mediation(
    data = df_log, X_conti = X_conti, X_binar = X_binar, M_conti = M_vars,
    add_conti_binar_interactions = TRUE, remove_zero_var = TRUE))
  stopifnot(identical(expanded$X_expanded, qr_expanded$X_expanded))
  analysis <- data.frame(source_row = source_rows, A = A_all, X_all,
    M1 = M_all$cd420, M2 = M_all$cd820, Y = Y_all, check.names = FALSE)
  .write_data_csv(analysis, file.path(output_directory, "ACTG175_analysis.csv"))
  prepared <- list(source_rows = source_rows, N = nrow(df), A_all = A_all,
    X_all = X_all, M_all = M_all, Y_all = Y_all, S_all = cbind(X_all, M_all),
    X_expanded = expanded$X_expanded, M_expanded = expanded$M_expanded,
    S_expanded = cbind(expanded$X_expanded, expanded$M_expanded),
    S_sdr_log_expanded = cbind(expanded$X_expanded, qr_expanded$M_expanded))
  saveRDS(prepared, file.path(output_directory, "ACTG175_analysis.rds"))
  message("Prepared ACTG175: 780 participants; 253 control and 527 treated.")
  message("Analysis data: ", normalizePath(output_directory))
  invisible(prepared)
}

prepared_actg175 <- prepare_actg175(.data_root())
