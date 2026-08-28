#!/usr/bin/env Rscript
# Mega-batch Design Module ANOVA without-replication perf-eval.
# Sets ANOVA_NOREP=1 and delegates to perf_eval_anova_megabatch.R.

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", ca, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[[1]]), winslash = "/")))
  }
  getwd()
}

Sys.setenv(ANOVA_NOREP = "1")
code_dir <- Sys.getenv("QST_CODE_DIR", unset = "")
if (!nzchar(code_dir)) code_dir <- get_script_dir()
Sys.setenv(PERF_EVAL_ANOVA_MEGABATCH_SOURCE_ONLY = "1")
source(file.path(code_dir, "perf_eval_anova_megabatch.R"), local = FALSE)
Sys.unsetenv("PERF_EVAL_ANOVA_MEGABATCH_SOURCE_ONLY")

if (!interactive()) {
  run_cli()
}
