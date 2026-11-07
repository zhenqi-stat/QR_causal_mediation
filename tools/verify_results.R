# Recompute the ACTG175 table from the seven saved split-level matrices.
# No patient-level data or nuisance estimation is required.
recompute_real_results <- function(record) {
  metrics <- c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
               "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")
  if (!is.list(record) || !all(c("N", metrics) %in% names(record))) {
    stop("The real-data result must contain N and all seven result matrices.")
  }
  if (length(record$N) != 1L || !is.finite(record$N) || record$N <= 0) {
    stop("N must be a positive finite scalar.")
  }
  reference_dim <- dim(record$pred_theta)
  methods <- colnames(record$pred_theta)
  if (length(reference_dim) != 2L || reference_dim[1] != 100L ||
      reference_dim[2] < 1L || is.null(methods) || anyDuplicated(methods)) {
    stop("The result matrices must contain 100 splits and named method columns.")
  }
  for (metric in metrics) {
    x <- record[[metric]]
    if (!is.matrix(x) || !is.numeric(x) ||
        !identical(dim(x), reference_dim) ||
        !identical(colnames(x), methods) || !all(is.finite(x))) {
      stop("Invalid or incomplete result matrix: ", metric)
    }
    if (grepl("^pred_sig2", metric) && any(x < 0)) {
      stop("Negative variance in ", metric)
    }
  }
  out <- list()
  out$DE_split <- record$pred_theta - record$pred_theta0
  out$IDE_split <- record$pred_theta1 - record$pred_theta
  out$TE_split <- record$pred_theta1 - record$pred_theta0
  for (mapping in list(c("theta", "pred_theta"), c("theta0", "pred_theta0"),
                       c("theta1", "pred_theta1"))) {
    out[[mapping[1]]] <- apply(record[[mapping[2]]], 2, median)
  }
  for (mapping in list(c("DE", "DE_split"), c("IDE", "IDE_split"), c("TE", "TE_split"))) {
    out[[mapping[1]]] <- apply(out[[mapping[2]]], 2, median)
  }
  # Apply the saved analysis's scalar split-aggregation formula unchanged.
  out$sig2 <- apply(record$pred_sig2 + sweep(record$pred_theta, 2, out$theta, "-")^2, 2, median)
  out$sig2_DE <- apply(record$pred_sig2_DE + sweep(out$DE_split, 2, out$DE, "-")^2, 2, median)
  out$sig2_IDE <- apply(record$pred_sig2_IDE + sweep(out$IDE_split, 2, out$IDE, "-")^2, 2, median)
  out$sig2_TE <- apply(record$pred_sig2_TE + sweep(out$TE_split, 2, out$TE, "-")^2, 2, median)
  out$p_theta <- 2 * pnorm(-abs(out$theta / sqrt(out$sig2 / record$N)))
  out$p_de <- 2 * pnorm(-abs(out$DE / sqrt(out$sig2_DE / record$N)))
  out$p_ide <- 2 * pnorm(-abs(out$IDE / sqrt(out$sig2_IDE / record$N)))
  out$p_te <- 2 * pnorm(-abs(out$TE / sqrt(out$sig2_TE / record$N)))
  out$proportion <- out$IDE / out$TE
  out$result_summary <- data.frame(
    method = methods,
    theta = unname(out$theta),
    theta_length = 2 * qnorm(.975) * sqrt(out$sig2 / record$N),
    NDE = unname(out$DE),
    NDE_length = 2 * qnorm(.975) * sqrt(out$sig2_DE / record$N),
    NDE_p = unname(out$p_de),
    NIE = unname(out$IDE),
    NIE_length = 2 * qnorm(.975) * sqrt(out$sig2_IDE / record$N),
    NIE_p = unname(out$p_ide),
    proportion = unname(out$proportion),
    ATE = unname(out$TE),
    ATE_length = 2 * qnorm(.975) * sqrt(out$sig2_TE / record$N),
    ATE_p = unname(out$p_te),
    row.names = methods
  )
  out$confidence_intervals <- lapply(list(
    theta10 = c("theta", "sig2"), NDE = c("DE", "sig2_DE"),
    NIE = c("IDE", "sig2_IDE"), ATE = c("TE", "sig2_TE")
  ), function(mapping) {
    half_length <- qnorm(.975) * sqrt(out[[mapping[2]]] / record$N)
    cbind(lower = out[[mapping[1]]] - half_length,
          upper = out[[mapping[1]]] + half_length)
  })
  out
}

# Check the bundled final tables without fitting any nuisance model.
verify_archive_results <- function(archive_root) {
  archive_root <- normalizePath(archive_root, mustWork = TRUE)
  directory <- file.path(archive_root, "results")
  cases <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
  expected_files <- c(paste0(names(cases), ".RData"), "real.RData")
  files <- list.files(directory, pattern = "[.]RData$")
  stopifnot(length(files) == 25L, setequal(files, expected_files))
  methods <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")
  equal <- function(x, y) {
    stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-12)))
  }
  loaded <- function(name) {
    env <- new.env(parent = baseenv())
    load(file.path(directory, name), envir = env)
    as.list(env)
  }
  checked <- 0L
  for (id in names(cases)) {
    result <- loaded(paste0(id, ".RData"))
    case <- cases[[id]]
    stopifnot(identical(result$source_id, id),
              identical(result$replication_ids, 1:100),
              identical(result$method_names, methods),
              identical(dim(result$pred_theta), c(100L, 8L)),
              identical(dim(result$pred_sig2), c(100L, 8L)),
              identical(colnames(result$pred_theta), methods),
              identical(colnames(result$pred_sig2), methods),
              identical(rownames(result$pred_theta), as.character(1:100)),
              identical(rownames(result$pred_sig2), as.character(1:100)),
              all(is.finite(result$pred_theta)), all(is.finite(result$pred_sig2)),
              all(result$pred_sig2 >= 0), result$N == case$N,
              result$seed == case$seed, result$tmax == 100L)
    equal(result$theta, case$theta)
    error <- result$pred_theta - result$theta
    se <- sqrt(result$pred_sig2 / result$N)
    summary <- data.frame(
      method = methods, bias = unname(apply(error, 2, median)),
      RMSE = unname(sqrt(colMeans(error^2))),
      AL = unname(apply(2 * qnorm(.975) * se, 2, median)),
      AC = unname(colMeans(abs(error) <= qnorm(.975) * se)),
      ESD = unname(apply(result$pred_theta, 2,
        function(x) 1.4826 * median(abs(x - median(x))))),
      ASD = unname(apply(se, 2, median)), row.names = NULL)
    summary <- summary[match(result$paper_method_names, summary$method), , drop = FALSE]
    rownames(summary) <- NULL
    equal(result$table_summary, summary)
    for (field in c("bias", "RMSE", "AL", "AC", "ESD", "ASD")) {
      equal(unname(result[[field]]), summary[[field]][match(methods, summary$method)])
    }
    equal(unname(result$sig2), unname(apply(result$pred_sig2, 2, median)))
    equal(unname(result$results_par), unname(cbind(result$pred_theta, result$pred_sig2)))
    stopifnot(identical(dim(result$results_par), c(100L, 16L)),
              identical(result$manuscript_table$method, summary$method))
    for (field in c("bias", "RMSE", "AL", "AC")) {
      digits <- if (field == "AC") 2L else 3L
      stopifnot(all(round(summary[[field]], digits) == result$manuscript_table[[field]]))
    }
    map <- result$source_map
    stopifnot(nrow(map) == 800L,
              !anyDuplicated(paste(map$replication_id, map$method)),
              setequal(map$method, methods), all(table(map$method) == 100L),
              all(map$replication_id %in% 1:100), all(map$source_row >= 1L),
              all(map$source_file == basename(map$source_file)),
              all(grepl("^[0-9a-f]{32}$", map$source_md5)))
    checked <- checked + 32L
  }
  real <- loaded("real.RData")
  stopifnot(identical(real$source_id, "real"),
            identical(real$replication_ids, 1:100),
            identical(real$method_names, c("MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR")),
            identical(real$paper_method_names, c("TR", "O-TR", "QR1", "QR2", "O-TR'", "MQR2")),
            identical(dim(real$results_par), c(100L, 42L)))
  recomputed <- recompute_real_results(real)
  for (field in names(recomputed)) equal(real[[field]], recomputed[[field]])
  matrices <- c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
                "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")
  equal(unname(real$results_par), unname(do.call(cbind, real[matrices])))
  equal(real$table_summary, {
    x <- recomputed$result_summary[match(real$paper_method_names,
      recomputed$result_summary$method), , drop = FALSE]
    rownames(x) <- NULL
    x
  })
  recipe <- real$source_recipe
  stopifnot(nrow(recipe) == 6L, identical(recipe$method, real$method_names),
            all(recipe$source_row_start == 1L), all(recipe$source_row_end == 100L),
            all(recipe$replication_id_start == 1L), all(recipe$replication_id_end == 100L),
            all(recipe$source_file == basename(recipe$source_file)),
            all(grepl("^[0-9a-f]{32}$", recipe$source_md5)))
  display <- real$manuscript_display
  stopifnot(nrow(display) == 57L, is.character(display$table_display))
  for (i in seq_len(nrow(display))) {
    row <- display[i, ]
    method <- if (row$method == "Common ATE") real$method_names[[1L]] else row$method
    value <- recomputed$result_summary[method, row$metric]
    if (row$metric == "proportion") value <- 100 * value
    expected <- row$table_display
    scientific <- grepl("e", expected, fixed = TRUE)
    mantissa <- if (scientific) strsplit(expected, "e", fixed = TRUE)[[1L]][[1L]] else expected
    fractional <- strsplit(mantissa, ".", fixed = TRUE)[[1L]]
    digits <- if (length(fractional) > 1L) nchar(fractional[[2L]]) else 0L
    actual <- sprintf(paste0("%.", digits, if (scientific) "e" else "f"), value)
    stopifnot(identical(actual, expected))
  }
  checked <- checked + nrow(display)
  # Validate the final source mapping. The raw author-side inputs are kept outside
  # the self-contained reviewer archive and receive additional checks when present.
  # QR1 has historical method column 5, but local source column 1 in its saved file.
  selection <- read.csv(file.path(archive_root, "provenance", "source_selection.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  required_qscale_cases <- paste0("simu_4_N", c(600L, 1000L, 1400L, 1800L))
  allowed_qscale_cases <- c(required_qscale_cases, paste0("simu_6_", c(1L, 3L, 5L, 6L)))
  qscale_algorithm <- "TR/O-TR only; q-scale type-7 truncation; outer-training A=1 reference"
  qscale_cases <- unique(selection$configuration[selection$source_algorithm == qscale_algorithm])
  stopifnot(all(required_qscale_cases %in% qscale_cases), all(qscale_cases %in% allowed_qscale_cases),
            length(qscale_cases) >= 4L, length(qscale_cases) <= 8L,
            nrow(selection) == 78L + 2L * length(qscale_cases),
            setequal(selection$configuration, c(names(cases), "real")),
            !anyDuplicated(paste(selection$configuration, selection$method)),
            all(selection$source_file == basename(selection$source_file)),
            all(selection$source_row_start == 1L), all(selection$source_row_end == 100L))
  for (id in names(cases)) {
    selected <- selection[selection$configuration == id, , drop = FALSE]
    expected_methods <- c("MQR1", "MQR2", "QR1")
    expected_local_columns <- c(1L, 2L, 1L)
    expected_legacy_columns <- c(1L, 2L, 5L)
    if (id %in% qscale_cases) {
      expected_methods <- c(expected_methods, "TR", "O-TR")
      expected_local_columns <- c(expected_local_columns, 1L, 2L)
      expected_legacy_columns <- c(expected_legacy_columns, 7L, 8L)
    }
    stopifnot(nrow(selected) == length(expected_methods),
              setequal(selected$method, expected_methods),
              all(selected$source_local_column[match(expected_methods, selected$method)] == expected_local_columns),
              all(selected$source_legacy_method_column[match(expected_methods, selected$method)] == expected_legacy_columns))
  }
  real_selected <- selection[selection$configuration == "real", , drop = FALSE]
  real_methods <- c("MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR")
  at <- match(real_methods, real_selected$method)
  stopifnot(nrow(real_selected) == 6L, setequal(real_selected$method, real_methods),
            identical(real_selected$source_method[at], real_methods),
            identical(as.integer(real_selected$source_local_column[at]), c(1L,2L,3L,4L,1L,2L)),
            identical(as.integer(real_selected$source_legacy_method_column[at]), c(2L,3L,5L,6L,7L,8L)),
            length(unique(real_selected$source_file[at[1:4]])) == 1L,
            length(unique(real_selected$source_file[at[5:6]])) == 1L,
            length(unique(real_selected$source_file)) == 2L)
  merge_audit <- file.path(dirname(archive_root), "author_only", "build_audit")
  has_author_inputs <- dir.exists(merge_audit)
  if (has_author_inputs) {
    manifest <- read.csv(file.path(merge_audit, "rerun_source_manifest.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(nrow(manifest) == 50L + length(qscale_cases),
              !anyDuplicated(manifest$source_file),
              setequal(manifest$source_file, selection$source_file),
              all(basename(manifest$source_path) == manifest$source_file),
              sum(startsWith(manifest$source_path, "rerun_results/")) == 48L,
              sum(startsWith(manifest$source_path, "../real_final_update/raw_sources/")) == 2L,
              sum(startsWith(manifest$source_path, "../qscale_update/raw_sources/")) == length(qscale_cases))
    source_paths <- file.path(merge_audit, manifest$source_path)
    stopifnot(all(file.exists(source_paths)),
              identical(unname(tools::md5sum(source_paths)), manifest$source_md5))
  }
  for (i in seq_len(nrow(selection))) {
    row <- selection[i, ]
    result <- loaded(paste0(row$configuration, ".RData"))
    if (has_author_inputs) {
      source_index <- match(row$source_file, manifest$source_file)
      stopifnot(!is.na(source_index), manifest$source_md5[source_index] == row$source_md5)
      source_file <- source_paths[source_index]
      stopifnot(identical(unname(tools::md5sum(source_file)), row$source_md5))
      e <- new.env(parent = baseenv())
      load(source_file, envir = e)
      input <- as.list(e)
      stopifnot(identical(input$replication_ids, 1:100),
                input$method_names[row$source_local_column] == row$source_method,
                input$N == result$N, input$seed == result$seed)
      if (row$configuration %in% qscale_cases && row$method %in% c("TR", "O-TR")) {
        stopifnot(identical(input$method_names, c("TR", "O-TR")),
                  identical(as.integer(input$method_columns), c(1L, 2L)),
                  identical(as.integer(input$historical_method_columns), c(7L, 8L)),
                  input$historical_method_columns[row$source_local_column] == row$source_legacy_method_column)
      } else {
        stopifnot(input$method_columns[row$source_local_column] == row$source_legacy_method_column)
      }
      fields <- if (row$configuration == "real") matrices else c("pred_theta", "pred_sig2")
      for (field in fields) stopifnot(identical(
        as.numeric(result[[field]][, row$method]),
        as.numeric(input[[field]][, row$source_local_column])))
    }
    if (row$configuration == "real") {
      source_map <- result$source_recipe[result$source_recipe$method == row$method, , drop = FALSE]
      stopifnot(nrow(source_map) == 1L)
    } else {
      source_map <- result$source_map[result$source_map$method == row$method, , drop = FALSE]
      stopifnot(nrow(source_map) == 100L, identical(source_map$source_row, 1:100))
    }
    stopifnot(all(source_map$source_file == row$source_file),
              all(source_map$source_md5 == row$source_md5),
              all(source_map$source_column == row$source_local_column))
  }
  if (has_author_inputs) {
    preservation <- read.csv(file.path(merge_audit, "merge_verification.csv"), stringsAsFactors = FALSE)
    qscale_preservation <- preservation[preservation$action == "replace_from_qscale", , drop = FALSE]
    expected_qscale_preservation <- expand.grid(
      configuration = qscale_cases, method = c("TR", "O-TR"),
      field = c("pred_theta", "pred_sig2"), stringsAsFactors = FALSE)
    preservation_key <- function(x) paste(x$configuration, x$method, x$field, sep = ":")
    stopifnot(!anyDuplicated(preservation_key(qscale_preservation)),
              setequal(preservation_key(qscale_preservation), preservation_key(expected_qscale_preservation)),
              all(qscale_preservation$n_values == 100L))
    stopifnot(nrow(preservation) == 450L, all(preservation$exact_values),
              all(preservation$max_abs_difference == 0),
              sum(preservation$action == "preserve_original" & preservation$configuration != "real" & preservation$field != "theta") == 240L - 4L * length(qscale_cases),
              sum(preservation$action == "replace_from_qscale" & preservation$configuration %in% qscale_cases) == 4L * length(qscale_cases),
              sum(preservation$action == "replace_from_qscale" & !preservation$configuration %in% qscale_cases) == 0L,
              sum(preservation$action == "replace_from_real_final" & preservation$configuration == "real") == 42L,
              sum(preservation$action != "replace_from_real_final" & preservation$configuration == "real") == 0L)
  }
  merged_manifest <- read.csv(file.path(archive_root, "provenance", "merged_results_md5.csv"), stringsAsFactors = FALSE)
  stopifnot(nrow(merged_manifest) == 25L, setequal(merged_manifest$file, paste0("results/", expected_files)),
            identical(unname(tools::md5sum(file.path(archive_root, merged_manifest$file))), merged_manifest$md5))
  message("PASS: 25 final RData files, 100 repetitions/splits each, and ",
          checked, " displayed table values recomputed.")
  invisible(TRUE)
}
