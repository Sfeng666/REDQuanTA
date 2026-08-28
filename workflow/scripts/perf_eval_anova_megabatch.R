#!/usr/bin/env Rscript
# Mega-batch Design Module ANOVA perf-eval (with or without replication).
# Writes RData batches compatible with aggregate_perf_eval_multicombo_fast.R.

get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", ca, value = TRUE)
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[[1]]), winslash = "/")))
  }
  getwd()
}

source_core <- function() {
  code_dir <- Sys.getenv("QST_CODE_DIR", unset = "")
  if (!nzchar(code_dir)) code_dir <- get_script_dir()
  core <- file.path(code_dir, "qst_abc_sim.R")
  if (!file.exists(core)) stop("qst_abc_sim.R not found: ", core)
  old <- Sys.getenv("QST_ABC_SOURCE_ONLY", unset = NA)
  Sys.setenv(QST_ABC_SOURCE_ONLY = "1")
  source(core, local = FALSE)
  if (is.na(old)) Sys.unsetenv("QST_ABC_SOURCE_ONLY") else Sys.setenv(QST_ABC_SOURCE_ONLY = old)
}

source_core()

DEFAULT_ADAPTIVE_QST <- seq(0.50, 1.00, by = 0.05)
DEFAULT_VE_RATIOS <- c(0.01, 0.1, 1.0, 10.0, 100.0)
DEFAULT_NUM_TRIALS <- 10000L
DEFAULT_BATCH_SIZE <- 1000L
ANOVA_REQUIRED_STATS <- c("QST", "both_var_negative", "among_pop_sd")

get_mc_cores <- function() {
  for (env_var in c("POC_P2_MC_CORES", "_CONDOR_NPROCS")) {
    val <- Sys.getenv(env_var, unset = "")
    if (nzchar(val)) {
      cores_val <- suppressWarnings(as.integer(val))
      if (!is.na(cores_val) && cores_val > 1L) return(cores_val)
    }
  }
  1L
}

use_norep_mode <- function() {
  nzchar(Sys.getenv("ANOVA_NOREP", "")) && Sys.getenv("ANOVA_NOREP") != "0"
}

combo_label <- function() {
  if (use_norep_mode()) "ANOVA_norep" else "ANOVA"
}

generate_neutral <- function(fst, ve, seed, required_stats) {
  if (use_norep_mode()) {
    generate_neutral_obs_stats_norep(fst, ve, seed, required_stats)
  } else {
    generate_neutral_obs_stats(fst, ve, seed, required_stats)
  }
}

generate_adaptive <- function(qst, ve, seed, required_stats) {
  if (use_norep_mode()) {
    generate_adaptive_obs_stats_norep(qst, ve, seed, required_stats)
  } else {
    generate_adaptive_obs_stats(qst, ve, seed, required_stats)
  }
}

save_neutral_batches <- function(df, ratio_idx, output_root, batch_size, combo) {
  ratio_dir <- file.path(output_root, paste0("neutral_ratio_", ratio_idx))
  dir.create(ratio_dir, recursive = TRUE, showWarnings = FALSE)
  n <- nrow(df)
  n_batches <- ceiling(n / batch_size)
  for (b in seq_len(n_batches)) {
    start <- (b - 1L) * batch_size + 1L
    end <- min(b * batch_size, n)
    batch_df <- df[start:end, , drop = FALSE]
    result <- list(
      mode = "batch_neutral",
      n_fst = nrow(batch_df),
      is_multi_combo = TRUE,
      results = batch_df
    )
    save(result, file = file.path(ratio_dir, paste0("neutral_batch_", b, ".RData")))
  }
}

save_adaptive_batches <- function(df, qst, ratio_idx, output_root, batch_size, combo) {
  qst_str <- gsub("\\.", "_", sprintf("%.2f", qst))
  cond_dir <- file.path(output_root, paste0("adaptive_q", qst_str, "_r", ratio_idx))
  dir.create(cond_dir, recursive = TRUE, showWarnings = FALSE)
  n <- nrow(df)
  n_batches <- ceiling(n / batch_size)
  for (b in seq_len(n_batches)) {
    start <- (b - 1L) * batch_size + 1L
    end <- min(b * batch_size, n)
    batch_df <- df[start:end, , drop = FALSE]
    result <- list(
      mode = "batch_evaluate",
      true_qst = qst,
      ve_ratio = unique(batch_df$ve_ratio)[1],
      n_repeats = nrow(batch_df),
      start_id = batch_df$repeat_id[1],
      is_multi_combo = TRUE,
      results = batch_df
    )
    save(result, file = file.path(cond_dir, paste0("adaptive_batch_", b, ".RData")))
  }
}

run_neutral_chr <- function(fst_file, output_root, num_trials = DEFAULT_NUM_TRIALS,
                            batch_size = DEFAULT_BATCH_SIZE, ve_ratios = DEFAULT_VE_RATIOS) {
  fst_values <- scan(fst_file, what = numeric(), quiet = TRUE)
  if (length(fst_values) < num_trials) {
    stop("FST file has ", length(fst_values), " values; need ", num_trials)
  }
  fst_values <- fst_values[seq_len(num_trials)]
  combo <- combo_label()
  cores <- get_mc_cores()
  cat("ANOVA neutral_chr: n_fst=", num_trials, " cores=", cores, " combo=", combo, "\n", sep = "")

  for (ri in seq_along(ve_ratios)) {
    ve <- ve_ratios[ri]
    cat("  neutral_ratio_", ri, " VE=", ve, "\n", sep = "")
    run_fst <- function(i) {
      fst <- fst_values[i]
      seed <- (round(fst * 1e6) + i * 10000) %% .Machine$integer.max
      obs <- generate_neutral(fst, ve, seed, ANOVA_REQUIRED_STATS)
      qst <- if (is.null(obs)) NA_real_ else estimate_qst_anova(obs)
      data.frame(fst = fst, combo = combo, qst = qst, stringsAsFactors = FALSE)
    }
    parts <- mclapply(seq_along(fst_values), run_fst, mc.cores = cores)
    df <- do.call(rbind, parts)
    save_neutral_batches(df, ri, output_root, batch_size, combo)
  }
}

run_adaptive_conditions <- function(output_root, qst_levels, num_trials = DEFAULT_NUM_TRIALS,
                                    batch_size = DEFAULT_BATCH_SIZE, ve_ratios = DEFAULT_VE_RATIOS) {
  combo <- combo_label()
  cores <- get_mc_cores()
  cat("ANOVA adaptive: qst_levels=", length(qst_levels), " trials=", num_trials,
      " cores=", cores, " combo=", combo, "\n", sep = "")

  for (qst in qst_levels) {
    for (ri in seq_along(ve_ratios)) {
      ve <- ve_ratios[ri]
      cat("  adaptive QST=", qst, " ratio_", ri, " VE=", ve, "\n", sep = "")
      run_trial <- function(i) {
        repeat_id <- i
        seed <- (round(qst * 100) * 10000 + round(ve * 1000) * 100 + repeat_id) %%
          .Machine$integer.max
        obs <- generate_adaptive(qst, ve, seed, ANOVA_REQUIRED_STATS)
        est <- if (is.null(obs)) NA_real_ else estimate_qst_anova(obs)
        data.frame(
          repeat_id = repeat_id,
          true_qst = qst,
          ve_ratio = ve,
          combo = combo,
          estimated_qst = est,
          stringsAsFactors = FALSE
        )
      }
      parts <- mclapply(seq_len(num_trials), run_trial, mc.cores = cores)
      df <- do.call(rbind, parts)
      save_adaptive_batches(df, qst, ri, output_root, batch_size, combo)
    }
  }
}

run_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("Usage:\n")
    cat("  neutral_chr: Rscript perf_eval_anova_megabatch.R neutral_chr FST_FILE OUTPUT_ROOT [num_trials] [batch_size]\n")
    cat("  adaptive_chr: Rscript perf_eval_anova_megabatch.R adaptive_chr OUTPUT_ROOT [num_trials] [batch_size]\n")
    cat("  adaptive_q: Rscript perf_eval_anova_megabatch.R adaptive_q QST OUTPUT_ROOT [num_trials] [batch_size]\n")
    quit(status = 1)
  }

  mode <- args[[1]]
  if (mode == "neutral_chr") {
    fst_file <- args[[2]]
    output_root <- args[[3]]
    num_trials <- if (length(args) >= 4) as.integer(args[[4]]) else DEFAULT_NUM_TRIALS
    batch_size <- if (length(args) >= 5) as.integer(args[[5]]) else DEFAULT_BATCH_SIZE
    run_neutral_chr(fst_file, output_root, num_trials, batch_size)
  } else if (mode == "adaptive_chr") {
    output_root <- args[[2]]
    num_trials <- if (length(args) >= 3) as.integer(args[[3]]) else DEFAULT_NUM_TRIALS
    batch_size <- if (length(args) >= 4) as.integer(args[[4]]) else DEFAULT_BATCH_SIZE
    run_adaptive_conditions(output_root, DEFAULT_ADAPTIVE_QST, num_trials, batch_size)
  } else if (mode == "adaptive_q") {
    qst <- as.numeric(args[[2]])
    output_root <- args[[3]]
    num_trials <- if (length(args) >= 4) as.integer(args[[4]]) else DEFAULT_NUM_TRIALS
    batch_size <- if (length(args) >= 5) as.integer(args[[5]]) else DEFAULT_BATCH_SIZE
    run_adaptive_conditions(output_root, qst, num_trials, batch_size)
  } else {
    stop("Unknown mode: ", mode)
  }
  cat("Done:", mode, "->", output_root, "\n")
}

if (!interactive() && Sys.getenv("PERF_EVAL_ANOVA_MEGABATCH_SOURCE_ONLY", "") != "1") {
  run_cli()
}
