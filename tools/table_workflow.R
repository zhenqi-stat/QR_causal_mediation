# Table generation and comparison for archived results and completed new runs.
.workflow_load <- function(path) {
  env <- new.env(parent = baseenv())
  load(path, envir = env)
  as.list(env)
}

.workflow_support <- function(root) {
  env <- new.env(parent = .GlobalEnv)
  sys.source(file.path(root, "tools", "verify_results.R"), envir = env)
  env
}

.workflow_simulation_summary <- function(record) {
  error <- record$pred_theta - record$theta
  se <- sqrt(record$pred_sig2 / record$N)
  summary <- data.frame(method = colnames(error),
    bias = unname(apply(error, 2L, median)), RMSE = unname(sqrt(colMeans(error^2))),
    AL = unname(apply(2 * qnorm(.975) * se, 2L, median)),
    AC = unname(colMeans(abs(error) <= qnorm(.975) * se)),
    ESD = unname(apply(record$pred_theta, 2L,
      function(x) 1.4826 * median(abs(x - median(x))))),
    ASD = unname(apply(se, 2L, median)), stringsAsFactors = FALSE)
  order <- c("Oracle", "TR", "O-TR", "QR1", "QR2", "O-TR'", "MQR1", "MQR2")
  summary <- summary[match(order, summary$method), , drop = FALSE]
  rownames(summary) <- NULL
  summary
}

.workflow_check_new_result <- function(record, id, case, reference) {
  fail <- function(text) stop(id, ": ", text, call. = FALSE)
  required <- c("source_id", "N", "seed", "tmax", "replication_ids", "method_names",
                "pred_theta", "pred_sig2")
  if (!all(required %in% names(record))) fail("missing result fields")
  if (!identical(record$source_id, id) || record$tmax != 100 ||
      record$N != reference$N || record$seed != reference$seed) fail("configuration metadata disagree")
  ids <- record$replication_ids
  if (!is.numeric(ids) || length(ids) != 100L || anyNA(ids) || anyDuplicated(ids) ||
      !setequal(ids, 1:100)) fail("expected exactly replications 1:100")
  methods <- if (id == "real") c("MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR") else
    c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")
  if (!identical(record$method_names, methods)) fail("incomplete or unexpected methods")
  fields <- if (id == "real") c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
    "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE") else c("pred_theta", "pred_sig2")
  if (!all(fields %in% names(record))) fail("missing estimate/variance matrices")
  for (field in fields) {
    value <- record[[field]]
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(dim(value), c(100L, length(methods))) ||
        !identical(colnames(value), methods) || !all(is.finite(value))) fail(paste("invalid", field))
    if (!is.null(rownames(value)) && !identical(rownames(value), as.character(ids))) {
      fail(paste("row labels disagree with replication IDs in", field))
    }
    if (startsWith(field, "pred_sig2") && any(value < 0)) fail(paste("negative variance in", field))
    record[[field]] <- value[order(ids), , drop = FALSE]
  }
  record$replication_ids <- 1:100
  if (id == "real") {
    for (field in c("pred_theta0", "pred_theta1", "pred_sig2_TE")) {
      value <- record[[field]]
      if (max(abs(sweep(value, 1L, value[, 1L], "-"))) > 1e-12) {
        fail(paste("shared endpoint calculations disagree in", field))
      }
    }
  }
  if (id != "real") {
    if (length(record$theta) != 1L || !is.finite(record$theta) ||
        abs(record$theta - case$theta) > 1e-12) fail("finite-pool target disagrees with DGP")
    oracle_error <- max(abs(record$pred_theta[, "Oracle"] - reference$pred_theta[, "Oracle"]))
    oracle_variance_error <- max(abs(record$pred_sig2[, "Oracle"] - reference$pred_sig2[, "Oracle"]))
    if (oracle_error > 1e-10 || oracle_variance_error > 1e-8) fail("DGP/Oracle reference mismatch")
    record$table_summary <- .workflow_simulation_summary(record)
  }
  record
}

.workflow_read_manifest <- function(path, root, cases, references) {
  path <- normalizePath(path, mustWork = TRUE)
  directory <- dirname(path)
  manifest <- read.csv(path, stringsAsFactors = FALSE, na.strings = "")
  required <- c("source_id", "result_file", "result_md5", "status")
  expected <- c(names(cases), "real")
  if (!all(required %in% names(manifest)) || nrow(manifest) != 25L ||
      anyDuplicated(manifest$source_id) || !setequal(manifest$source_id, expected)) {
    stop("The run manifest must list all 25 configurations exactly once.", call. = FALSE)
  }
  manifest <- manifest[match(expected, manifest$source_id), , drop = FALSE]
  incomplete <- which(is.na(manifest$status) | manifest$status != "complete" |
                        is.na(manifest$result_file) | is.na(manifest$result_md5))
  if (length(incomplete)) stop("Incomplete run: ",
    paste(manifest$source_id[incomplete], collapse = ", "), call. = FALSE)
  records <- setNames(vector("list", 25L), expected)
  for (i in seq_len(nrow(manifest))) {
    file <- normalizePath(file.path(directory, manifest$result_file[i]), mustWork = TRUE)
    if (!startsWith(file, paste0(directory, .Platform$file.sep))) {
      stop("Result path is outside the run directory: ", manifest$source_id[i], call. = FALSE)
    }
    if (unname(tools::md5sum(file)) != manifest$result_md5[i]) {
      stop("Run-result MD5 mismatch: ", manifest$source_id[i], call. = FALSE)
    }
    id <- manifest$source_id[i]
    records[[id]] <- .workflow_check_new_result(.workflow_load(file), id, cases[[id]], references[[id]])
  }
  records
}

# Read both the printed strings and unrounded values from the exported tables.
.workflow_cells <- function(directory, table_map) {
  rows <- list()
  add <- function(table, configuration, method, metric, value, display) {
    rows[[length(rows) + 1L]] <<- data.frame(table = table, configuration = configuration,
      method = method, metric = metric, value = as.numeric(value), display = as.character(display),
      stringsAsFactors = FALSE)
  }
  for (i in seq_len(nrow(table_map))) {
    stem <- paste0(tolower(table_map$document[i]), "_table_", table_map$table[i])
    numeric <- read.csv(file.path(directory, paste0(stem, ".csv")), stringsAsFactors = FALSE)
    display <- read.csv(file.path(directory, paste0(stem, "_display.csv")),
                        colClasses = "character", check.names = FALSE)
    real <- table_map$entrypoint[i] == "real.R"
    for (j in seq_len(nrow(display))) {
      method <- display$method[j]
      matching <- numeric$method == method
      if (!real) matching <- matching & numeric$configuration == display$configuration[j]
      row <- numeric[matching, , drop = FALSE]
      if (nrow(row) != 1L) stop("Ambiguous table row: ", stem, " / ", method)
      metrics <- if (real) setdiff(names(display), "method") else c("Bias", "RMSE", "AL", "AC")
      for (metric in metrics) {
        field <- if (metric == "Bias") "bias" else metric
        value <- row[[field]]
        if (field == "proportion") value <- 100 * value
        add(stem, if (real) "real" else display$configuration[j], method, metric, value, display[[metric]][j])
      }
    }
    if (real) {
      common <- read.csv(file.path(directory, paste0(stem, "_common_ATE.csv")), colClasses = "character")
      for (metric in names(common)) add(stem, "real", "Common ATE", metric, numeric[[metric]][1], common[[metric]][1])
    }
  }
  cells <- do.call(rbind, rows)
  if (nrow(cells) != 816L || any(!is.finite(cells$value))) stop("Expected 816 finite manuscript cells.")
  cells
}

.workflow_compare <- function(reference, rerun, directory) {
  key <- function(x) paste(x$table, x$configuration, x$method, x$metric, sep = "|")
  if (anyDuplicated(key(reference)) || anyDuplicated(key(rerun)) ||
      !setequal(key(reference), key(rerun))) stop("Table comparison keys disagree.")
  rerun <- rerun[match(key(reference), key(rerun)), , drop = FALSE]
  comparison <- reference[, c("table", "configuration", "method", "metric")]
  comparison$reference_value <- reference$value
  comparison$rerun_value <- rerun$value
  comparison$difference <- rerun$value - reference$value
  comparison$absolute_difference <- abs(comparison$difference)
  comparison$reference_display <- reference$display
  comparison$rerun_display <- rerun$display
  comparison$same_display <- reference$display == rerun$display
  summary <- do.call(rbind, lapply(unique(comparison$table), function(table) {
    x <- comparison[comparison$table == table, , drop = FALSE]
    data.frame(table = table, cells = nrow(x), matching_cells = sum(x$same_display),
      differing_cells = sum(!x$same_display), maximum_absolute_difference = max(x$absolute_difference),
      status = if (all(x$same_display)) "MATCH_AT_PAPER_PRECISION" else "NUMERICAL_DIFFERENCES",
      stringsAsFactors = FALSE)
  }))
  write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
  write.csv(summary, file.path(directory, "comparison_summary.csv"), row.names = FALSE)
  text <- c("Result structure: PASS (25 complete files, 100 repetitions/splits each).",
    paste0("Table cells matching at paper precision: ", sum(comparison$same_display), "/", nrow(comparison), "."),
    paste0("Numerical comparison: ", if (all(comparison$same_display)) "MATCH_AT_PAPER_PRECISION" else "NUMERICAL_DIFFERENCES", "."),
    "comparison.csv lists every reference value, new value, difference, and printed value.")
  writeLines(text, file.path(directory, "comparison_summary.txt"))
  message(paste(text[1:3], collapse = "\n"))
  message("Comparison report: ", file.path(directory, "comparison.csv"))
  invisible(comparison)
}

.workflow_render <- function(root, manifest = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  support <- .workflow_support(root)
  support$verify_archive_results(root)
  cases <- readRDS(file.path(root, "configuration_index.rds"))$cases
  ids <- c(names(cases), "real")
  references <- setNames(lapply(ids, function(id) .workflow_load(file.path(root, "results", paste0(id, ".RData")))), ids)
  table_map <- read.csv(file.path(root, "table_map.csv"), stringsAsFactors = FALSE)
  paper_cells <- .workflow_cells(file.path(root, "tables"), table_map)
  if (is.null(manifest)) {
    records <- references
    directory <- file.path(root, "tables")
  } else {
    manifest <- normalizePath(manifest, mustWork = TRUE)
    records <- .workflow_read_manifest(manifest, root, cases, references)
    directory <- file.path(dirname(manifest), "tables")
  }
  .workflow_export(root, records, directory, support$recompute_real_results, archived = is.null(manifest))
  .workflow_compare(paper_cells, .workflow_cells(directory, table_map), directory)
  message("Tables: ", directory)
  invisible(directory)
}

.workflow_write_manifest <- function(manifest, path) {
  temporary <- paste0(path, ".tmp")
  write.csv(manifest, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE)) stop("Cannot save run manifest.")
    unlink(temporary)
  }
}

.workflow_full_fit <- function(root, runner = NULL, output_root = file.path(root, "output")) {
  root <- normalizePath(root, mustWork = TRUE)
  support <- .workflow_support(root)
  support$verify_archive_results(root)
  if (is.null(runner)) {
    execution <- new.env(parent = .GlobalEnv)
    sys.source(file.path(root, "tools", "archive_runtime.R"), envir = execution)
    runner <- execution$run_archive
  }
  parent <- file.path(output_root, "full_run")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent)) stop("Cannot create full-run output directory.")
  stem <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_pid", Sys.getpid())
  created <- FALSE
  for (i in 0:999) {
    directory <- file.path(parent, paste0(stem, if (i == 0L) "" else sprintf("_%02d", i + 1L)))
    if (dir.create(directory, showWarnings = FALSE)) { created <- TRUE; break }
    if (!dir.exists(directory)) stop("Cannot create full-run directory.")
  }
  if (!created) stop("Cannot allocate a full-run directory.")
  directory <- normalizePath(directory)
  dir.create(file.path(directory, "logs"))
  log_file <- file.path(directory, "logs", "progress.log")
  announce <- function(text) {
    cat(text, "\n", sep = ""); flush.console(); flush(stdout())
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", text, "\n",
        sep = "", file = log_file, append = TRUE)
  }
  announce(paste0("Full run: ", directory))
  announce(paste0("Log file: ", log_file))
  manifest_path <- file.path(directory, "run_manifest.csv")
  announce(paste0("Run manifest: ", manifest_path))
  cases <- readRDS(file.path(root, "configuration_index.rds"))$cases
  ids <- c(names(cases), "real")
  manifest <- data.frame(source_id = ids,
    group = c(vapply(cases, function(x) as.character(x$group), character(1)), "real"),
    selector = c(vapply(cases, function(x) as.character(x$config), character(1)), "ACTG175"),
    result_file = "", result_md5 = "", status = "pending", stringsAsFactors = FALSE)
  .workflow_write_manifest(manifest, manifest_path)
  for (i in seq_along(ids)) {
    id <- ids[i]
    group <- if (id == "real") "real" else cases[[id]]$group
    config <- if (id == "real") NULL else cases[[id]]$config
    announce(paste0("Starting ", id))
    tryCatch({
      file <- runner(group, config, replication_ids = 1:100,
                     archive_root = root, output_dir = directory)
      file <- normalizePath(file, mustWork = TRUE)
      if (!startsWith(file, paste0(directory, .Platform$file.sep))) stop("Result outside this full-run directory.")
      reference <- .workflow_load(file.path(root, "results", paste0(id, ".RData")))
      .workflow_check_new_result(.workflow_load(file), id, cases[[id]], reference)
      manifest$result_file[i] <- substring(file, nchar(directory) + 2L)
      manifest$result_md5[i] <- unname(tools::md5sum(file))
      manifest$status[i] <- "complete"
      .workflow_write_manifest(manifest, manifest_path)
      announce(paste0("Completed ", id))
    }, error = function(error) {
      manifest$status[i] <- "failed"
      .workflow_write_manifest(manifest, manifest_path)
      announce(paste0("Failed ", id, ": ", conditionMessage(error)))
      stop(error)
    })
    gc(verbose = FALSE)
  }
  announce("All 25 result files complete. Generating tables and comparison reports.")
  .workflow_render(root, manifest_path)
  announce(paste0("Completed full run. Tables and reports: ", file.path(directory, "tables")))
  invisible(manifest_path)
}

# Export the same manuscript layout for either source of result records.
.workflow_export <- function(archive_root, records, table_directory, recompute_real_results, archived) {
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  load_result <- function(id) records[[id]]
  cases <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
  table_map <- utils::read.csv(file.path(archive_root, "table_map.csv"),
                               stringsAsFactors = FALSE)
  table_stem <- function(entrypoint) {
    mapping <- table_map[table_map$entrypoint == entrypoint, , drop = FALSE]
    stopifnot(nrow(mapping) == 1L)
    paste0(tolower(mapping$document), "_table_", mapping$table)
  }
  table_stems <- vapply(paste0("simu_", 1:6, ".R"), table_stem, character(1))
  for (group in 1:6) {
    selected <- Filter(function(x) identical(as.integer(x$group), group), cases)
    summary_rows <- list()
    display_rows <- list()
    latex_rows <- character()
    for (case in selected) {
      result <- load_result(case$id)
      summary <- result$table_summary
      summary_rows[[case$id]] <- cbind(
        configuration = case$id, N = case$N, selector = case$config, summary)
      display <- data.frame(
        configuration = case$id, method = summary$method,
        Bias = sprintf("%.3f", summary$bias), RMSE = sprintf("%.3f", summary$RMSE),
        AL = sprintf("%.3f", summary$AL), AC = sprintf("%.2f", summary$AC),
        stringsAsFactors = FALSE)
      display_rows[[case$id]] <- display
      latex_rows <- c(latex_rows, paste0("% ", case$id),
        vapply(seq_len(nrow(display)), function(i) {
          paste0(paste(unlist(display[i, c("method", "Bias", "RMSE", "AL", "AC")]),
                       collapse = " & "), " \\\\")
        }, character(1)), "")
    }
    stem <- file.path(table_directory, table_stems[group])
    utils::write.csv(do.call(rbind, summary_rows), paste0(stem, ".csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, display_rows), paste0(stem, "_display.csv"), row.names = FALSE)
    writeLines(latex_rows, paste0(stem, "_rows.tex"))
  }
  real <- load_result("real")
  aggregate <- recompute_real_results(real)
  # Keep all six saved methods intact; export the five rows in Main Table 6.
  real_table_methods <- c("TR", "O-TR", "QR2", "O-TR'", "MQR2")
  stopifnot(all(real_table_methods %in% real$method_names))
  summary <- aggregate$result_summary[
    match(real_table_methods, aggregate$result_summary$method), , drop = FALSE]
  real_stem <- file.path(table_directory, table_stem("real.R"))
  utils::write.csv(summary, paste0(real_stem, ".csv"), row.names = FALSE)
  format_p <- function(x) ifelse(x < .001, sprintf("%.1e", x), sprintf("%.3f", x))
  display <- data.frame(
    method = summary$method,
    theta = sprintf("%.3f", summary$theta), theta_length = sprintf("%.3f", summary$theta_length),
    NDE = sprintf("%.3f", summary$NDE), NDE_length = sprintf("%.3f", summary$NDE_length),
    NDE_p = format_p(summary$NDE_p), NIE = sprintf("%.3f", summary$NIE),
    NIE_length = sprintf("%.3f", summary$NIE_length), NIE_p = format_p(summary$NIE_p),
    proportion = sprintf("%.2f", 100 * summary$proportion), stringsAsFactors = FALSE)
  # Check every printed string against the saved manuscript display, not only rounded numerics.
  if (archived) {
  for (i in seq_len(nrow(real$manuscript_display))) {
    row <- real$manuscript_display[i, ]
    if (!row$method %in% real_table_methods) next
    actual <- display[match(row$method, display$method), row$metric]
    stopifnot(identical(unname(actual), row$table_display))
  }
  }
  utils::write.csv(display, paste0(real_stem, "_display.csv"), row.names = FALSE)
  latex_rows <- vapply(seq_len(nrow(display)), function(i) {
    row <- unlist(display[i, ])
    row[length(row)] <- paste0(row[length(row)], "\\%")
    paste0(paste(row, collapse = " & "), " \\\\")
  }, character(1))
  common <- data.frame(ATE = sprintf("%.3f", summary$ATE[1]),
    ATE_length = sprintf("%.3f", summary$ATE_length[1]), ATE_p = format_p(summary$ATE_p[1]))
  utils::write.csv(common, paste0(real_stem, "_common_ATE.csv"), row.names = FALSE)
  writeLines(c(latex_rows, "", paste0("% Common ATE: ", paste(unlist(common), collapse = ", "))),
             paste0(real_stem, "_rows.tex"))
  message("Regenerated main Tables 4-6 and supplementary Tables 7-10 in: ", table_directory)
}
