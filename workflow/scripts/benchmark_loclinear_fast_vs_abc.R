#!/usr/bin/env Rscript
# Compare perf_eval_all_combos_fast_cell_v2 (vectorized loclinear) vs abc::abc
# on the *same* reference pool (build once). Two separate full runs with the same
# set.seed() are not bit-identical here because mcmapply + abc() advance RNG
# differently than a second isolated simulation run.
#
# Usage: Rscript benchmark_loclinear_fast_vs_abc.R [num_sim] [n_combos]
# Run from any cwd; resolves sibling scripts in this directory.

args_ca <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", args_ca, value = TRUE)
CODE_DIR <- if (length(f)) {
  dirname(normalizePath(sub("^--file=", "", f[[1]]), winslash = "/"))
} else {
  getwd()
}

trailing <- commandArgs(trailingOnly = TRUE)
num_sim <- if (length(trailing) >= 1L) as.integer(trailing[[1]]) else 12000L
n_combos <- if (length(trailing) >= 2L) as.integer(trailing[[2]]) else 28L
if (is.na(num_sim) || num_sim < 500L) stop("num_sim too small")
if (is.na(n_combos) || n_combos < 1L) stop("n_combos invalid")

Sys.setenv(
  PERF_EVAL_ALL_COMBOS_FAST_SOURCE_ONLY = "1",
  QST_ABC_METHOD = "loclinear"
)
source(file.path(CODE_DIR, "perf_eval_all_combos_fast_cell_v2.R"))

stats_for_pairs <- c(
  "QST", "ratioVbetweenVtotal", "F_within_pop", "F_among_pop",
  "among_pop_sd", "within_pop_sd", "ratioVwithinVtotal", "ratioVextVtotal"
)
pair_mat <- combn(stats_for_pairs, 2, simplify = FALSE)
if (n_combos > length(pair_mat)) {
  stop("n_combos too large; max unique pairs is ", length(pair_mat))
}
combos <- pair_mat[seq_len(n_combos)]

required <- get_required_stats(combos)
obs_stats <- generate_adaptive_obs_stats(0.75, 0.1, 424242L, required)
if (is.null(obs_stats)) stop("generate_adaptive_obs_stats returned NULL")

cat("num_sim =", num_sim, " n_combos =", n_combos, " stat cols =",
    length(required), "\n", sep = "")

set.seed(10001L)
t_sim <- system.time({
  pool <- build_fast_reference_pool(obs_stats, num_sim, combos)
})
n_valid <- nrow(pool$sim_stats_valid)
tol <- max(0.001, get_abc_tol_numerator() / n_valid)
P <- pool$prior_params_valid
S <- pool$sim_stats_valid
os <- obs_stats[pool$required_stats]

t_abc_only <- system.time({
  abc_vec <- vapply(combos, function(cb) {
    res <- run_abc_qst(
      target = os[cb],
      param = P,
      sumstat = S[, cb, drop = FALSE],
      tol = tol,
      transf = rep("none", 3L)
    )
    qst_mean_from_abc(res)
  }, numeric(1))
})
names(abc_vec) <- vapply(combos, function(cb) paste(cb, collapse = ","), character(1))

t_fast <- system.time({
  fast_vec <- fast_abc_estimate_multi(
    os, P, S, combos,
    tol = tol,
    estimator = "loclinear",
    hcorr = TRUE
  )
})

diffs <- abc_vec[names(fast_vec)] - fast_vec
ok <- is.finite(diffs)

t_full <- system.time({
  set.seed(10001L)
  estimate_qst_abc_multi(obs_stats, num_sim, combos)
})

cat("\n--- Timing (seconds, this host) ---\n")
cat(sprintf("Full estimate_qst_abc_multi (sim + abc each combo): %8.3f elapsed\n", t_full["elapsed"]))
cat(sprintf("Sim only (build_fast_reference_pool):              %8.3f elapsed\n", t_sim["elapsed"]))
cat(sprintf("abc::abc loclinear only (same pool, loop combos):   %8.3f elapsed\n", t_abc_only["elapsed"]))
cat(sprintf("fast_abc_estimate_multi only (same pool):          %8.3f elapsed\n", t_fast["elapsed"]))
cat(sprintf("Speedup abc-loop vs fast (post-sim):                 %8.3fx\n",
            t_abc_only["elapsed"] / t_fast["elapsed"]))
cat(sprintf("Speedup full multi vs sim+fast:                     %8.3fx\n",
            t_full["elapsed"] / (t_sim["elapsed"] + t_fast["elapsed"])))

cat("\n--- Agreement abc vs fast (same pool, ", sum(ok), "/", length(ok), " finite) ---\n", sep = "")
cat(sprintf("max |difference|:              %g\n", max(abs(diffs[ok]), na.rm = TRUE)))
cat(sprintf("max |diff| / max(|abc|,1e-12): %g\n",
            max(abs(diffs[ok]) / pmax(abs(abc_vec[names(fast_vec)][ok]), 1e-12), na.rm = TRUE)))
