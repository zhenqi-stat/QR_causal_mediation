# Check code and saved results by default; add --dgp to regenerate all simulation DGP pools.
.verification_root <- function() {
  frames <- sys.frames()
  paths <- vapply(frames, function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L) path else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  if (length(paths)) {
    path <- paths[[length(paths)]]
  } else {
    path <- commandArgs(trailingOnly = FALSE)
    path <- path[startsWith(path, "--file=")]
    if (!length(path)) stop("Run verify.R with Rscript or source().")
    path <- sub("^--file=", "", path[[1L]])
  }
  if (!file.exists(path) && file.exists(basename(path))) path <- basename(path)
  dirname(normalizePath(path, mustWork = TRUE))
}

.has_constant_true_if <- function(expression) {
  if (!is.call(expression)) return(FALSE)
  if (identical(expression[[1L]], as.name("if")) &&
      identical(expression[[2L]], TRUE)) return(TRUE)
  any(vapply(as.list(expression), .has_constant_true_if, logical(1)))
}

archive_root <- .verification_root()
arguments <- commandArgs(trailingOnly = TRUE)
manifest <- NULL
manifest_option <- which(arguments == "--manifest")
if (length(manifest_option)) {
  if (length(manifest_option) != 1L || manifest_option == length(arguments)) {
    stop("Usage: Rscript verify.R [--dgp] [--manifest /path/to/run_manifest.csv]")
  }
  manifest <- arguments[manifest_option + 1L]
  arguments <- arguments[-c(manifest_option, manifest_option + 1L)]
}
if (length(setdiff(arguments, "--dgp"))) {
  stop("Usage: Rscript verify.R [--dgp] [--manifest /path/to/run_manifest.csv]")
}
checksum_file <- file.path(archive_root, "file_checksums.csv")
checksums <- utils::read.csv(checksum_file, stringsAsFactors = FALSE)
stopifnot(nrow(checksums) > 0L, !anyDuplicated(checksums$path))
paths <- file.path(archive_root, checksums$path)
if (any(!file.exists(paths))) stop("An archived file is missing.")
actual <- unname(tools::md5sum(paths))
if (!identical(actual, checksums$md5)) {
  stop("Checksum mismatch: ", paste(checksums$path[actual != checksums$md5], collapse = ", "))
}

source(file.path(archive_root, "tools", "archive_runtime.R"), local = TRUE)
archive_runtime_self_test()
index <- readRDS(file.path(archive_root, "configuration_index.rds"))$cases
stopifnot(length(index) > 0L, !anyDuplicated(names(index)),
          setequal(vapply(index, function(case) as.integer(case$group), integer(1)), 1:7))
seventh <- Filter(function(case) as.integer(case$group) == 7L, index)
stopifnot(setequal(vapply(seventh, function(case) as.integer(case$N), integer(1)),
                   c(1000L, 1800L, 2500L, 5000L)))
simu7 <- new.env(parent = .GlobalEnv)
old_autorun <- getOption("qr.archive.no_autorun")
options(qr.archive.no_autorun = TRUE)
tryCatch(sys.source(file.path(archive_root, "simu_7.R"), envir = simu7),
         finally = options(qr.archive.no_autorun = old_autorun))
stopifnot(is.function(simu7$run_experiment), is.function(simu7$.simu7_prepare_dgp),
  identical(names(formals(simu7$run_experiment)),
    c("rerun_config", "rerun_replication_ids", "rerun_workers", "rerun_output_dir")),
  identical(simu7$.simu7_method_names,
    c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR", "Oracle")),
  identical(as.integer(simu7$N), c(1000L, 1800L, 2500L, 5000L)),
  identical(simu7$replication_ids, 1:100), simu7$workers == 30L)
simu7_expressions <- parse(file.path(archive_root, "simu_7.R"), keep.source = FALSE)
stopifnot(!any(vapply(simu7_expressions, .has_constant_true_if, logical(1))),
  sum(vapply(simu7_expressions, .archive_lhs_is, logical(1), name = "replication_ids")) == 1L)
for (case in seventh) {
  # Each DGP configuration must work without inheriting helpers from simu_7.R.
  standalone <- new.env(parent = .GlobalEnv)
  sys.source(.archive_source_path(archive_root, case$path), envir = standalone)
  stopifnot(exists("shift_truth", envir = standalone, inherits = FALSE))
  configuration <- standalone$.simu7_config(case$N)
  parameters <- standalone$.simu7_parameters(configuration)
  actual_data <- standalone$.simu7_generate(configuration, parameters, 1L)
  reference_config <- simu7$.simu7_config(case$N)
  reference_data <- simu7$.simu7_generate(reference_config,
    simu7$.simu7_parameters(reference_config), 1L)
  stopifnot(identical(actual_data, reference_data))
}
for (case in index) {
  selected <- .archive_select_case(archive_root, case$group, case$config)
  stopifnot(identical(selected$id, case$id), case$tmax == 100)
  path <- .archive_source_path(archive_root, case$path)
  expressions <- parse(path, keep.source = FALSE)
  stopifnot(!any(vapply(expressions, .has_constant_true_if, logical(1))))
  active_text <- paste(deparse(expressions, width.cutoff = 500L), collapse = "\n")
  completions <- c(" | MQR1 | complete", " | MQR2 | complete",
    " | QR1 | complete", " | QR2 | complete", " | O-TR' | complete",
    " | TR & O-TR | complete")
  if (as.integer(case$group) <= 6L) {
    stopifnot(all(vapply(completions, grepl, logical(1), x = active_text, fixed = TRUE)))
  } else {
    stopifnot(identical(as.integer(case$group), 7L), case$config == case$N,
              case$seed == 2026L)
  }
  stopifnot(sum(vapply(expressions, .archive_lhs_is, logical(1), name = "replication_ids")) == 1L)
}
r_files <- list.files(archive_root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
r_files <- r_files[!grepl("/output/", r_files, fixed = TRUE)]
for (path in r_files) parse(path, keep.source = FALSE)
message("PASS: checksums, ", length(r_files), " R source files, ", length(index), " configurations, method interfaces, and runtime self-tests.")

source(file.path(archive_root, "tools", "verify_results.R"), local = TRUE)
verify_archive_results(archive_root)

if ("--dgp" %in% arguments) {
  source(file.path(archive_root, "tools", "verify_dgp.R"), local = TRUE)
  result <- verify_archive_dgp(archive_root)
  print(result, row.names = FALSE)
  message("PASS: ", length(index), " finite-pool targets and ", 100L * length(index), " Oracle estimates/variances; no nuisance-model fitting.")
}

if (!is.null(manifest)) {
  source(file.path(archive_root, "tools", "table_workflow.R"), local = TRUE)
  .workflow_render(archive_root, manifest)
}
