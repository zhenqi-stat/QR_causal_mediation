# Master workflow. Default: regenerate tables. --fit: rerun all 25 analyses.
.master_paths <- vapply(sys.frames(), function(x) {
  if (is.character(x$ofile) && length(x$ofile) == 1L) x$ofile else ""
}, character(1))
.master_paths <- .master_paths[nzchar(.master_paths)]
if (!length(.master_paths)) {
  .master_paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
}
if (!length(.master_paths)) stop("Run with source() or Rscript.")
archive_root <- dirname(normalizePath(tail(.master_paths, 1L), mustWork = TRUE))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(arguments, "--fit"))) stop("Usage: Rscript run_all.R [--fit]")
source(file.path(archive_root, "tools", "table_workflow.R"), local = TRUE)
if (!"--fit" %in% arguments) {
  .workflow_render(archive_root)
} else {
  .workflow_full_fit(archive_root)
}
