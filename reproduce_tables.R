# Regenerate the seven numerical tables from the bundled replication-level results.
# This script uses base R only and does not fit nuisance models.
.table_root <- function() {
  paths <- vapply(sys.frames(), function(x) {
    if (is.character(x$ofile) && length(x$ofile) == 1L) x$ofile else ""
  }, character(1))
  paths <- paths[nzchar(paths)]
  if (!length(paths)) {
    paths <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  }
  if (!length(paths)) stop("Run with source() or Rscript.")
  dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
}
archive_root <- .table_root()
source(file.path(archive_root, "tools", "table_workflow.R"), local = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  .workflow_render(archive_root)
} else if (length(arguments) == 2L && arguments[[1L]] == "--manifest") {
  .workflow_render(archive_root, arguments[[2L]])
} else {
  stop("Usage: Rscript reproduce_tables.R [--manifest /path/to/run_manifest.csv]")
}
