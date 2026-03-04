#!/usr/bin/env Rscript
# ABC-based QST estimation - adapted from run_abc_sim.R
# MEMORY-OPTIMIZED VERSION: Processes in batches to reduce peak memory usage
#
# This script estimates QST using ABC with 100,000 simulations
# Can be used for both trait QST and neutral QST estimation
#
# Usage:
#   Rscript qst_abc_sim.R trait obs_stats.RData ignored output.RData [num_sim] [summary_stats]
#   Rscript qst_abc_sim.R neutral 0.05 0.1 output.RData [num_sim] [summary_stats]
#
# Arguments:
#   1. mode: "trait" or "neutral"
#   2. input_file: For trait mode: path to observed summary stats
#                  For neutral mode: FST value
#   3. ext_sd: Extrinsic SD (from observed trait data)
#   4. output_file: Path to save results
#   5. num_sim: Number of ABC simulations (default: 100000)
#   6. summary_stats: Comma-separated summary statistics (default: "QST,F_within_pop")

suppressPackageStartupMessages({
  library(abc)
  library(parallel)
})

# Load e1071 only if needed (for skewness and kurtosis)
.e1071_loaded <- FALSE
load_e1071_if_needed <- function() {
  if (!.e1071_loaded) {
    suppressPackageStartupMessages(library(e1071))
    .e1071_loaded <<- TRUE
  }
}

# Use fewer cores to limit memory usage from forking
# Each forked process duplicates memory, so fewer cores = lower peak memory
# 4 cores is a good balance between speed and memory for CHTC jobs with 4GB limit
num_cores <- min(4, max(1, detectCores() - 1))

# All possible summary statistics
ALL_SUMMARY_STATS <- c("among_pop_sd", "within_pop_sd", "ext_sd", "QST", "ratioVext", 
                       "F_among_pop", "F_within_pop", "skewness_data", "kurtosis_data")

# Basic stats always calculated (needed for variance components)
BASIC_STATS <- c("among_pop_sd", "within_pop_sd", "ext_sd")

#' Determine which summary statistics need to be calculated
#' Based on the combinations that will be used for ABC
#' 
#' @param summary_stat_combos List of character vectors specifying combinations
#' @return Character vector of unique stats that need to be calculated
get_required_stats <- function(summary_stat_combos) {
  # Get all unique stats from all combinations
  all_stats <- unique(unlist(summary_stat_combos))
  
  # Always include basic stats (needed for prior generation and variance calculation)
  required <- union(BASIC_STATS, all_stats)
  
  # Order according to ALL_SUMMARY_STATS for consistent column ordering
  required <- intersect(ALL_SUMMARY_STATS, required)
  
  return(required)
}

#' Parse summary stats argument
#' Supports: comma-separated string, "all", or file path
#' 
#' @param summary_stats_arg The argument value (string, "all", or file path)
#' @return List of character vectors, each specifying a combination
parse_summary_stats <- function(summary_stats_arg) {
  # Check if it's a file path
  if (file.exists(summary_stats_arg)) {
    cat("Loading summary stat combinations from file:", summary_stats_arg, "\n")
    lines <- readLines(summary_stats_arg)
    lines <- lines[lines != ""]  # Remove empty lines
    
    # Parse tab-separated combinations
    combos <- lapply(lines, function(line) {
      stats <- strsplit(line, "\t")[[1]]
      # Handle ratioVenv -> ratioVext conversion
      stats <- gsub("ratioVenv", "ratioVext", stats)
      stats <- gsub("env_sd", "ext_sd", stats)
      return(stats)
    })
    
    cat("Loaded", length(combos), "summary stat combinations\n")
    return(combos)
  }
  
  # If "all" is specified, return NULL (will be generated from simulation)
  if (summary_stats_arg == "all") {
    return(NULL)  # Signal to generate from simulation
  }
  
  # Otherwise, parse as comma-separated single combination
  stats <- strsplit(summary_stats_arg, ",")[[1]]
  return(list(stats))
}

#' Calculate correlation p-value matrix
#' Simple implementation without corrplot dependency
cor_pvalue_matrix <- function(mat) {
  n <- ncol(mat)
  pmat <- matrix(1, nrow = n, ncol = n)
  colnames(pmat) <- rownames(pmat) <- colnames(mat)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      test <- cor.test(mat[, i], mat[, j])
      pmat[i, j] <- pmat[j, i] <- test$p.value
    }
  }
  return(pmat)
}

#' Generate uncorrelated summary stat combinations from simulation data
#' Based on legacy script's approach
#' 
#' @param sim_stats_matrix Matrix of simulated summary statistics
#' @param pvalue_threshold P-value threshold for correlation (default: 0.01)
#' @return List of character vectors, each specifying an uncorrelated combination
generate_uncorrelated_combos <- function(sim_stats_matrix, pvalue_threshold = 0.01) {
  # Calculate correlation matrix and p-values
  correlation_matrix <- cor(sim_stats_matrix, method = "pearson")
  pvalue_matrix <- cor_pvalue_matrix(sim_stats_matrix)
  
  # Set diagonal to 1 and 0 to avoid self-exclusion
  diag(pvalue_matrix) <- 1
  diag(correlation_matrix) <- 0
  
  # Function to calculate weighted score
  calculate_weighted_score <- function(corr_mat) {
    apply(corr_mat, 1, function(x) sum(abs(x)))
  }
  
  # Function to exclude correlated stats
  exclude_correlated <- function(corr_mat, pval_mat, threshold) {
    while (any(pval_mat < threshold, na.rm = TRUE)) {
      weighted_scores <- calculate_weighted_score(corr_mat)
      stat_to_exclude <- which.max(weighted_scores)
      corr_mat <- corr_mat[-stat_to_exclude, -stat_to_exclude, drop = FALSE]
      pval_mat <- pval_mat[-stat_to_exclude, -stat_to_exclude, drop = FALSE]
    }
    return(list(correlation_matrix = corr_mat, pvalue_matrix = pval_mat))
  }
  
  # Generate all possible combinations
  stat_names <- colnames(correlation_matrix)
  all_combos <- list()
  for (size in 1:length(stat_names)) {
    combos <- combn(stat_names, size, simplify = FALSE)
    all_combos <- c(all_combos, combos)
  }
  
  # Filter to keep only uncorrelated combinations
  filtered_combos <- lapply(all_combos, function(combo) {
    if (length(combo) == 1) return(combo)
    
    # Get submatrix for this combination
    sub_corr <- correlation_matrix[combo, combo, drop = FALSE]
    sub_pval <- pvalue_matrix[combo, combo, drop = FALSE]
    
    result <- exclude_correlated(sub_corr, sub_pval, pvalue_threshold)
    return(rownames(result$correlation_matrix))
  })
  
  # Remove NULL and duplicate combinations
  filtered_combos <- Filter(Negate(is.null), filtered_combos)
  filtered_combos <- unique(filtered_combos)
  
  cat("Generated", length(filtered_combos), "uncorrelated combinations\n")
  return(filtered_combos)
}

# Fixed parameters matching sample structure - load from environment or use defaults
# Defaults match the original script: 2 populations, 6 individuals, 3 replicates
num_pop <- as.integer(Sys.getenv("SIM_NUM_POP", "2"))
num_ind <- as.integer(Sys.getenv("SIM_NUM_IND", "6"))
num_rep <- as.integer(Sys.getenv("SIM_NUM_REP", "3"))

cat("Sample structure configuration:\n")
cat("  Populations:", num_pop, "\n")
cat("  Individuals per pop:", num_ind, "\n")
cat("  Replicates per ind:", num_rep, "\n")

mean_trait <- 100
var_additive <- 1

# Pre-compute indices for variance calculation (MEMORY EFFICIENT)
# These are computed once at script load, not per simulation
.pop_idx <- rep(1:num_pop, each = num_ind * num_rep)
.line_idx <- rep(1:(num_pop * num_ind), each = num_rep)
.n_total <- num_pop * num_ind * num_rep
.n_lines <- num_pop * num_ind

# Pre-compute degrees of freedom
.DF_among <- num_pop - 1
.DF_within <- num_pop * (num_ind - 1)
.DF_residual <- num_pop * num_ind * (num_rep - 1)

#' Calculate variance components using vectorized operations (NO data.frame)
#' This is much more memory efficient than the tapply-based version
calc_variance_fast <- function(trait_values) {
  # Calculate means using matrix operations instead of tapply
  # Pop means: average of all values in each population
  # Use rowSums on matrix reshape for arbitrary num_pop
  pop_sums <- rowSums(matrix(trait_values, nrow = num_pop, ncol = num_ind * num_rep, byrow = TRUE))
  pop_means <- pop_sums / (num_ind * num_rep)
  overall_mean <- sum(trait_values) / .n_total
  
  # Line means using vectorized reshape
  line_sums <- rowSums(matrix(trait_values, nrow = .n_lines, ncol = num_rep, byrow = TRUE))
  line_means <- line_sums / num_rep
  
  # Expand means for SS calculations
  pop_means_expanded <- rep(pop_means, each = num_ind)
  line_means_expanded <- rep(line_means, each = num_rep)
  
  # Sum of squares
  SS_among <- sum((pop_means - overall_mean)^2) * num_ind * num_rep
  SS_within <- sum((line_means - pop_means_expanded)^2) * num_rep
  SS_residual <- sum((trait_values - line_means_expanded)^2)
  
  # Mean squares and variance components
  MS_among <- SS_among / .DF_among
  MS_within <- SS_within / .DF_within
  MS_residual <- SS_residual / .DF_residual
  
  var_among <- max((MS_among - MS_within) / (num_ind * num_rep), 0)
  var_within <- max((MS_within - MS_residual) / num_rep, 0)
  var_residual <- max(MS_residual, 0)
  
  return(c(var_among, var_within, var_residual))
}

#' Generate simulated trait data and return summary statistics
#' MEMORY EFFICIENT: Uses vectorized operations, no data.frame creation
#' SELECTIVE: Only calculates the summary statistics specified in required_stats
#' 
#' @param required_stats Character vector of stats to calculate (default: all)
generate_sim_data_summarystats <- function(
  num_pop = 2, num_ind = 6, num_rep = 3,
  mean_trait = 100, sd_between_pop = 20, sd_within_pop = 10, sd_ext = 1,
  required_stats = NULL
) {
  # Default to all stats if not specified
  if (is.null(required_stats)) {
    required_stats <- ALL_SUMMARY_STATS
  }
  
  # Generate trait values directly as a vector (no data.frame)
  mu_between_pop <- sqrt(2 * sd_between_pop^2)
  means_pop <- c(mean_trait, mean_trait + mu_between_pop)
  means_ind <- rep(means_pop, each = num_ind) + rnorm(num_pop * num_ind, 0, sd_within_pop)
  trait_values <- rep(means_ind, each = num_rep) + rnorm(.n_total, 0, sd_ext)
  
  # Calculate variance using fast vectorized method (always needed)
  var_components <- calc_variance_fast(trait_values)
  among_pop_var <- var_components[1]
  within_pop_var <- var_components[2]
  ext_var <- var_components[3]
  
  # Initialize result vector with only required stats
  result <- numeric(length(required_stats))
  names(result) <- required_stats
  
  # Basic stats (always calculated as part of variance components)
  if ("among_pop_sd" %in% required_stats) result["among_pop_sd"] <- sqrt(among_pop_var)
  if ("within_pop_sd" %in% required_stats) result["within_pop_sd"] <- sqrt(within_pop_var)
  if ("ext_sd" %in% required_stats) result["ext_sd"] <- sqrt(ext_var)
  
  # Derived stats (only calculate if needed)
  total_genetic_var <- among_pop_var + within_pop_var
  
  if ("QST" %in% required_stats) {
    result["QST"] <- if (total_genetic_var == 0) 0 else among_pop_var / (among_pop_var + 2 * within_pop_var)
  }
  
  if ("ratioVext" %in% required_stats) {
    result["ratioVext"] <- if (total_genetic_var == 0) 0 else ext_var / total_genetic_var
  }
  
  if ("F_among_pop" %in% required_stats) {
    result["F_among_pop"] <- if (total_genetic_var == 0 || within_pop_var == 0) 0 else among_pop_var / within_pop_var
  }
  
  if ("F_within_pop" %in% required_stats) {
    result["F_within_pop"] <- if (total_genetic_var == 0 || ext_var == 0) 0 else within_pop_var / ext_var
  }
  
  # Higher-order moments (expensive - only calculate if needed)
  if ("skewness_data" %in% required_stats) {
    load_e1071_if_needed()
    result["skewness_data"] <- skewness(trait_values)
  }
  
  if ("kurtosis_data" %in% required_stats) {
    load_e1071_if_needed()
    result["kurtosis_data"] <- kurtosis(trait_values)
  }
  
  return(result)
}

#' Run batch simulations - MEMORY OPTIMIZED
#' @param required_stats Character vector of stats to calculate (passed to generate_sim_data_summarystats)
run_batch_simulations <- function(batch_size, sd_between_pop_batch, sd_within_pop_batch, sd_ext_batch, 
                                   required_stats = NULL) {
  results <- mcmapply(generate_sim_data_summarystats,
                      MoreArgs = list(num_pop = num_pop, num_ind = num_ind, num_rep = num_rep,
                                     mean_trait = mean_trait, required_stats = required_stats),
                      sd_between_pop = sd_between_pop_batch,
                      sd_within_pop = sd_within_pop_batch,
                      sd_ext = sd_ext_batch,
                      SIMPLIFY = TRUE,
                      mc.cores = num_cores,
                      mc.preschedule = TRUE,
                      mc.set.seed = TRUE)
  
  if (is.matrix(results)) {
    return(t(results))
  } else {
    return(do.call(rbind, results))
  }
}

#' Estimate QST using ABC with MULTIPLE summary stat combinations
#' Runs simulations once and applies ABC with each combination
#' OPTIMIZED: Only calculates summary statistics that are needed
#' MEMORY EFFICIENT: Aggressive cleanup after each ABC run
#' 
#' @param obs_stats Named vector of observed summary statistics
#' @param num_sim Number of ABC simulations (default: 100000)
#' @param summary_stat_combos List of character vectors, each specifying a combination
#' @return Named list with QST estimates for each combination
estimate_qst_abc_multi <- function(obs_stats, num_sim = num_sim, summary_stat_combos) {
  
  # Ensure cleanup happens even on error
  on.exit(gc(verbose = FALSE, full = TRUE, reset = TRUE), add = TRUE)
  
  # OPTIMIZATION: Determine which stats are actually needed
  required_stats <- get_required_stats(summary_stat_combos)
  n_stats <- length(required_stats)
  cat("Required stats (", n_stats, "/9):", paste(required_stats, collapse = ", "), "\n")
  
  batch_size <- 5000  # Smaller batches reduce peak memory from forking
  num_batches <- ceiling(num_sim / batch_size)
  
  # Generate prior parameters based on observed stats
  sd_between_pop_prior <- runif(num_sim, 0, 2 * obs_stats['among_pop_sd'])
  sd_within_pop_prior <- runif(num_sim, sqrt(0.001) * sd_between_pop_prior, 
                               sqrt(10) * sd_between_pop_prior)
  sd_ext_prior <- runif(num_sim, 0, 2 * obs_stats['ext_sd'])
  
  # Pre-allocate result matrix - ONLY for required stats (memory optimization)
  sim_stats_matrix <- matrix(NA, nrow = num_sim, ncol = n_stats)
  prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
  
  rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  gc(verbose = FALSE, reset = TRUE)
  
  # Process in batches
  cat("Simulating", num_sim, "for", length(summary_stat_combos), "combinations...\n")
  for (b in 1:num_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, num_sim)
    
    batch_results <- run_batch_simulations(
      end_idx - start_idx + 1,
      prior_params[start_idx:end_idx, 1],
      prior_params[start_idx:end_idx, 2],
      prior_params[start_idx:end_idx, 3],
      required_stats = required_stats  # Pass required stats
    )
    
    sim_stats_matrix[start_idx:end_idx, ] <- batch_results
    rm(batch_results)
    gc(verbose = FALSE, reset = TRUE)
  }
  
  colnames(sim_stats_matrix) <- required_stats
  
  # Clean up invalid rows
  valid_rows <- complete.cases(sim_stats_matrix) & 
                (rowSums(is.infinite(sim_stats_matrix)) == 0)
  
  sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
  prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
  n_valid <- nrow(sim_stats_valid)
  
  rm(sim_stats_matrix, prior_params, valid_rows)
  gc(verbose = FALSE, reset = TRUE)
  
  cat("Valid simulations:", n_valid, "/", num_sim, "\n")
  
  # Prepare obs_stats with only required stats for ABC
  obs_stats_subset <- obs_stats[required_stats]
  
  # Run ABC for EACH combination on the SAME simulated data
  # MEMORY OPTIMIZATION: Aggressive cleanup after each combo
  tol <- max(0.001, 50 / n_valid)
  results <- list()
  
  for (i in seq_along(summary_stat_combos)) {
    combo <- summary_stat_combos[[i]]
    combo_name <- paste(combo, collapse = ",")
    
    qst_estimate <- tryCatch({
      # Extract subset for this combo
      combo_sumstat <- sim_stats_valid[, combo, drop = FALSE]
      combo_target <- obs_stats_subset[combo]
      
      res_abc <- abc(
        target = combo_target,
        param = prior_params_valid,
        sumstat = combo_sumstat,
        sizenet = 10,
        tol = tol,
        transf = rep("none", ncol(prior_params_valid)),
        method = "neuralnet"
      )
      
      result <- mean(res_abc$adj.values[, 'sd_between_pop']^2 / 
                     (res_abc$adj.values[, 'sd_between_pop']^2 + 
                      2 * res_abc$adj.values[, 'sd_within_pop']^2))
      
      # AGGRESSIVE CLEANUP after each ABC run
      rm(res_abc, combo_sumstat, combo_target)
      gc(verbose = FALSE, reset = TRUE)
      
      result
      
    }, error = function(e) {
      warning(paste("ABC failed for", combo_name, ":", e$message))
      gc(verbose = FALSE, reset = TRUE)
      NA
    })
    
    results[[combo_name]] <- qst_estimate
    
    # Progress indicator
    if (i %% 10 == 0 || i == length(summary_stat_combos)) {
      cat("  ABC combo", i, "/", length(summary_stat_combos), "\n")
    }
  }
  
  rm(sim_stats_valid, prior_params_valid, obs_stats_subset)
  gc(verbose = FALSE, full = TRUE, reset = TRUE)
  
  return(results)
}

#' Estimate QST using ABC - disk-efficient version
#' Runs 100,000 simulations in memory, returns only QST estimate
#' 
#' @param obs_stats Named vector of observed summary statistics
#' @param num_sim Number of ABC simulations (default: 100000)
#' @param summary_stat_names Vector of summary statistic names to use for ABC
#' @return Estimated QST value
estimate_qst_abc <- function(obs_stats, num_sim = num_sim, summary_stat_names = c("QST", "ratioVext")) {
  
  # Ensure cleanup happens even on error - use reset=TRUE to return memory to OS
  on.exit(gc(verbose = FALSE, full = TRUE, reset = TRUE), add = TRUE)
  
  batch_size <- 5000  # Smaller batches reduce peak memory from forking
  num_batches <- ceiling(num_sim / batch_size)
  
  # Generate prior parameters based on observed stats
  sd_between_pop_prior <- runif(num_sim, 0, 2 * obs_stats['among_pop_sd'])
  sd_within_pop_prior <- runif(num_sim, sqrt(0.001) * sd_between_pop_prior, 
                                sqrt(10) * sd_between_pop_prior)
  sd_ext_prior <- runif(num_sim, 0, 2 * obs_stats['ext_sd'])
  
  # Pre-allocate result matrices
  # Dimensions determined by ALL_SUMMARY_STATS (defined at top of script)
  sim_stats_matrix <- matrix(NA, nrow = num_sim, ncol = length(ALL_SUMMARY_STATS))
  prior_params <- cbind(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  colnames(prior_params) <- c("sd_between_pop", "sd_within_pop", "sd_ext")
  
  # Free the individual prior vectors early - they're now in prior_params
  rm(sd_between_pop_prior, sd_within_pop_prior, sd_ext_prior)
  gc(verbose = FALSE)
  
  # Process in batches
  cat("Simulating", num_sim, "in", num_batches, "batches...\n")
  for (b in 1:num_batches) {
    start_idx <- (b - 1) * batch_size + 1
    end_idx <- min(b * batch_size, num_sim)
    
    batch_results <- run_batch_simulations(
      end_idx - start_idx + 1,
      prior_params[start_idx:end_idx, 1],
      prior_params[start_idx:end_idx, 2],
      prior_params[start_idx:end_idx, 3]
    )
    
    sim_stats_matrix[start_idx:end_idx, ] <- batch_results
    rm(batch_results)
    # Use reset=TRUE to help return memory to OS after forked processes
    gc(verbose = FALSE, reset = TRUE)
    
    if (b %% 5 == 0 || b == num_batches) {
      cat("  Batch", b, "/", num_batches, "complete\n")
    }
  }
  
  colnames(sim_stats_matrix) <- ALL_SUMMARY_STATS
  
  # Clean up invalid rows - use vectorized check (faster and less memory)
  valid_rows <- complete.cases(sim_stats_matrix) & 
                (rowSums(is.infinite(sim_stats_matrix)) == 0)
  
  sim_stats_valid <- sim_stats_matrix[valid_rows, , drop = FALSE]
  prior_params_valid <- prior_params[valid_rows, , drop = FALSE]
  n_valid <- nrow(sim_stats_valid)
  
  # Free original matrices immediately after subsetting
  rm(sim_stats_matrix, prior_params, valid_rows)
  gc(verbose = FALSE)
  
  cat("Valid simulations:", n_valid, "/", num_sim, "\n")
  
  # Run ABC
  tol <- max(0.001, 50 / n_valid)
  
  qst_estimate <- tryCatch({
    res_abc <- abc(
      target = obs_stats[summary_stat_names],
      param = prior_params_valid,
      sumstat = sim_stats_valid[, summary_stat_names],
      sizenet = 10,
      tol = tol,
      transf = rep("none", ncol(prior_params_valid)),
      method = "neuralnet"
    )
    
    # Extract result immediately
    result <- mean(res_abc$adj.values[, 'sd_between_pop']^2 / 
                   (res_abc$adj.values[, 'sd_between_pop']^2 + 
                    2 * res_abc$adj.values[, 'sd_within_pop']^2))
    
    # Clean up ABC result object immediately
    rm(res_abc)
    gc(verbose = FALSE)
    
    result
    
  }, error = function(e) {
    warning(paste("ABC failed:", e$message))
    NA
  })
  
  # Final cleanup of remaining objects
  rm(sim_stats_valid, prior_params_valid)
  # Use reset=TRUE to return memory to OS (important for forked processes in cgroups)
  gc(verbose = FALSE, full = TRUE, reset = TRUE)
  
  return(qst_estimate)
}

#' Generate observed summary statistics for neutral simulation
#' Uses FST to derive sd_between_pop and sd_within_pop
#' OPTIMIZED: Only calculates required stats
#' 
#' @param fst_value FST value to simulate
#' @param ratioVext V_E / V_G ratio observed from real traits
#' @param seed Random seed for reproducibility
#' @param required_stats Character vector of stats to calculate (default: all)
generate_neutral_obs_stats <- function(fst_value, ratioVext, seed = NULL, required_stats = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Derive variance parameters from FST
  sd_within_pop <- sqrt(abs((1 - fst_value) * var_additive))
  sd_between_pop <- sqrt(abs(2 * fst_value * var_additive))
  
  # Derive ext_sd properly mapped via V_E / V_G ratio to match simulated additive variance limits
  ext_sd <- sqrt(ratioVext * var_additive)
  
  # Generate "observed" data for this neutral scenario
  obs_valid <- FALSE
  attempts <- 0
  max_attempts <- 100
  
  while (!obs_valid && attempts < max_attempts) {
    obs <- generate_sim_data_summarystats(
      num_pop = num_pop, num_ind = num_ind, num_rep = num_rep,
      mean_trait = mean_trait,
      sd_between_pop = sd_between_pop,
      sd_within_pop = sd_within_pop,
      sd_ext = ext_sd,
      required_stats = required_stats
    )
    
    if (obs['among_pop_sd'] != 0 && !any(is.nan(obs)) && !any(is.infinite(obs))) {
      obs_valid <- TRUE
    }
    attempts <- attempts + 1
  }
  
  if (!obs_valid) {
    warning("Could not generate valid neutral observation")
    return(NULL)
  }
  
  return(obs)
}

#' Generate observed summary statistics for adaptive trait simulation (Module 2)
#' Uses pre-defined QST and V_E/V_G ratio to derive variance parameters
#' OPTIMIZED: Only calculates required stats
#' 
#' @param qst_value Pre-defined QST value (0-1)
#' @param ve_ratio V_E / V_G ratio (extrinsic variance / genetic variance)
#' @param seed Random seed for reproducibility
#' @param required_stats Character vector of stats to calculate (default: all)
#' @return Named vector of summary statistics, or NULL if generation fails
generate_adaptive_obs_stats <- function(qst_value, ve_ratio, seed = NULL, required_stats = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Derive variance parameters from QST
  # QST = V_GB / (V_GB + 2*V_GW), so if we fix total genetic variance = var_additive
  # V_GB = QST * (V_GB + 2*V_GW) => V_GB = QST * V_total_genetic / (1 - QST + 2*QST) 
  # Simpler: use same formula as neutral (FST = QST under neutrality)
  sd_within_pop <- sqrt(abs((1 - qst_value) * var_additive))
  sd_between_pop <- sqrt(abs(2 * qst_value * var_additive))
  
  # Derive ext_sd from V_E ratio
  # ve_ratio = V_E / V_G, where V_G = V_GB + V_GW = var_additive
  ext_sd <- sqrt(ve_ratio * var_additive)
  
  # Generate "observed" data for this adaptive trait
  obs_valid <- FALSE
  attempts <- 0
  max_attempts <- 100
  
  while (!obs_valid && attempts < max_attempts) {
    obs <- generate_sim_data_summarystats(
      num_pop = num_pop, num_ind = num_ind, num_rep = num_rep,
      mean_trait = mean_trait,
      sd_between_pop = sd_between_pop,
      sd_within_pop = sd_within_pop,
      sd_ext = ext_sd,
      required_stats = required_stats
    )
    
    if (obs['among_pop_sd'] != 0 && !any(is.nan(obs)) && !any(is.infinite(obs))) {
      obs_valid <- TRUE
    }
    attempts <- attempts + 1
  }
  
  if (!obs_valid) {
    warning(paste("Could not generate valid adaptive observation for QST =", qst_value))
    return(NULL)
  }
  
  return(obs)
}

#' Process a batch of repeats for adaptive QST estimation (Module 2)
#' Simulates traits with pre-defined QST and V_E ratio, estimates QST via ABC
#' 
#' @param n_repeats Number of repeats to process
#' @param qst_value Pre-defined QST value
#' @param ve_ratio V_E / V_G ratio
#' @param num_sim Number of ABC simulations per repeat
#' @param summary_stat_names Summary statistics to use for ABC
#' @param start_id Starting ID for reproducible seeds
#' @return Data frame with repeat_id, true_qst, ve_ratio, and estimated_qst columns
process_evaluate_batch <- function(n_repeats, qst_value, ve_ratio, num_sim, 
                                   summary_stat_names = c("QST", "ratioVext"),
                                   start_id = 1) {
  results <- data.frame(
    repeat_id = seq(start_id, start_id + n_repeats - 1),
    true_qst = rep(qst_value, n_repeats),
    ve_ratio = rep(ve_ratio, n_repeats),
    estimated_qst = rep(NA_real_, n_repeats)
  )
  
  for (i in 1:n_repeats) {
    repeat_id <- results$repeat_id[i]
    
    # Set seed for reproducibility (unique per repeat/qst/ratio combination)
    seed <- (round(qst_value * 100) * 10000 + round(ve_ratio * 1000) * 100 + repeat_id) %% .Machine$integer.max
    
    # Generate adaptive observed stats
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed)
    
    if (!is.null(obs_stats)) {
      # Run ABC estimation
      qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stat_names)
      results$estimated_qst[i] <- qst_estimate
    }
    
    # Force garbage collection after each repeat to prevent memory accumulation
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
    
    # Progress update every 10 repeats
    if (i %% 10 == 0 || i == n_repeats) {
      cat("  Processed", i, "/", n_repeats, "repeats\n")
    }
  }
  
  return(results)
}

#' Process a batch of repeats with MULTIPLE summary stat combinations (Module 2)
#' Runs ABC with all combinations on the SAME simulated data for efficiency
#' OPTIMIZED: Only calculates required summary statistics
#' 
#' @param n_repeats Number of repeats to process
#' @param qst_value Pre-defined QST value
#' @param ve_ratio V_E / V_G ratio
#' @param num_sim Number of ABC simulations per repeat
#' @param summary_stat_combos List of character vectors, each specifying a combination
#' @param start_id Starting ID for reproducible seeds
#' @return Data frame with repeat_id, true_qst, ve_ratio, combo, and estimated_qst columns
process_evaluate_batch_multi <- function(n_repeats, qst_value, ve_ratio, num_sim, 
                                         summary_stat_combos,
                                         start_id = 1) {
  n_combos <- length(summary_stat_combos)
  combo_names <- sapply(summary_stat_combos, function(x) paste(x, collapse = ","))
  
  # OPTIMIZATION: Determine required stats upfront
  required_stats <- get_required_stats(summary_stat_combos)
  cat("Required stats for batch:", paste(required_stats, collapse = ", "), "\n")
  
  # Create result data frame with all combinations
  results <- data.frame(
    repeat_id = rep(seq(start_id, start_id + n_repeats - 1), each = n_combos),
    true_qst = rep(qst_value, n_repeats * n_combos),
    ve_ratio = rep(ve_ratio, n_repeats * n_combos),
    combo = rep(combo_names, n_repeats),
    estimated_qst = rep(NA_real_, n_repeats * n_combos),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_repeats) {
    repeat_id <- start_id + i - 1
    
    # Set seed for reproducibility
    seed <- (round(qst_value * 100) * 10000 + round(ve_ratio * 1000) * 100 + repeat_id) %% .Machine$integer.max
    
    # Generate adaptive observed stats (only required ones)
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed, required_stats)
    
    if (!is.null(obs_stats)) {
      # Run ABC with ALL combinations on the SAME simulated data
      qst_estimates <- estimate_qst_abc_multi(obs_stats, num_sim, summary_stat_combos)
      
      # Store results for each combination
      for (j in seq_along(combo_names)) {
        idx <- (i - 1) * n_combos + j
        results$estimated_qst[idx] <- qst_estimates[[combo_names[j]]]
      }
    }
    
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
    
    if (i %% 10 == 0 || i == n_repeats) {
      cat("  Processed", i, "/", n_repeats, "repeats\n")
    }
  }
  
  return(results)
}

#' Process a batch of FST values with MULTIPLE summary stat combinations
#' Runs ABC with all combinations on the SAME simulated data for efficiency
#' OPTIMIZED: Only calculates required summary statistics
#' 
#' @param fst_values Vector of FST values to process
#' @param ratioVext V_E / V_G ratio observed from real traits
#' @param num_sim Number of ABC simulations per FST
#' @param summary_stat_combos List of character vectors, each specifying a combination
#' @return Data frame with fst, combo, and qst columns
process_fst_batch_multi <- function(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx = 1) {
  n_fst <- length(fst_values)
  n_combos <- length(summary_stat_combos)
  combo_names <- sapply(summary_stat_combos, function(x) paste(x, collapse = ","))
  
  # OPTIMIZATION: Determine required stats upfront
  required_stats <- get_required_stats(summary_stat_combos)
  cat("Required stats for batch:", paste(required_stats, collapse = ", "), "\n")
  
  # Create result data frame with all combinations
  results <- data.frame(
    fst = rep(fst_values, each = n_combos),
    combo = rep(combo_names, n_fst),
    qst = rep(NA_real_, n_fst * n_combos),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(fst_values)) {
    fst_value <- fst_values[i]
    
    # Set seed for reproducibility - Use start_idx to ensure unique seeds for identical FST values
    global_idx <- start_idx + i - 1
    seed <- (round(fst_value * 1e6) + global_idx * 10000) %% .Machine$integer.max
    
    # Generate neutral observed stats (only required ones)
    obs_stats <- generate_neutral_obs_stats(fst_value, ratioVext, seed, required_stats)
    
    if (!is.null(obs_stats)) {
      # Run ABC with ALL combinations on the SAME simulated data
      qst_estimates <- estimate_qst_abc_multi(obs_stats, num_sim, summary_stat_combos)
      
      # Store results for each combination
      for (j in seq_along(combo_names)) {
        idx <- (i - 1) * n_combos + j
        results$qst[idx] <- qst_estimates[[combo_names[j]]]
      }
    }
    
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
    
    if (i %% 10 == 0 || i == n_fst) {
      cat("  Processed", i, "/", n_fst, "FST values\n")
    }
  }
  
  return(results)
}

#' Process a batch of FST values for neutral QST estimation
#' Processes FST values sequentially (ABC already uses parallelism internally)
#' 
#' @param fst_values Vector of FST values to process
#' @param ratioVext V_E / V_G ratio observed from real traits
#' @param num_sim Number of ABC simulations per FST
#' @param summary_stat_names Summary statistics to use for ABC
#' @return Data frame with fst and qst columns
process_fst_batch <- function(fst_values, ratioVext, num_sim, summary_stat_names = c("QST", "ratioVext"), start_idx = 1) {
  n_fst <- length(fst_values)
  results <- data.frame(fst = fst_values, qst = rep(NA_real_, n_fst))
  
  for (i in seq_along(fst_values)) {
    fst_value <- fst_values[i]
    
    # Set seed for reproducibility - Use start_idx to ensure unique seeds for identical FST values
    global_idx <- start_idx + i - 1
    seed <- (round(fst_value * 1e6) + global_idx * 10000) %% .Machine$integer.max
    
    # Generate neutral observed stats
    obs_stats <- generate_neutral_obs_stats(fst_value, ratioVext, seed)
    
    if (!is.null(obs_stats)) {
      # Run ABC estimation
      qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stat_names)
      results$qst[i] <- qst_estimate
    }
    
    # Force garbage collection after each FST to prevent memory accumulation
    # Use reset=TRUE to return memory to OS (critical for cgroup memory limits)
    gc(verbose = FALSE, full = TRUE, reset = TRUE)
    
    # Progress update every 10 FST values
    if (i %% 10 == 0 || i == n_fst) {
      cat("  Processed", i, "/", n_fst, "FST values\n")
    }
  }
  
  return(results)
}

# Command-line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 4) {
    cat("Usage: Rscript qst_abc_sim.R <mode> <input> <ratioVext_or_ve_ratio> <output_file> [num_sim] [summary_stats]\n")
    cat("\nModes:\n")
    cat("  Module 1 (Detection):\n")
    cat("    trait        - Estimate QST for observed trait data\n")
    cat("    neutral      - Estimate QST for single neutral FST value\n")
    cat("    batch_neutral - Estimate QST for batch of neutral FST values\n")
    cat("  Module 2 (Performance Evaluation):\n")
    cat("    evaluate     - Estimate QST for single adaptive QST/VE_ratio\n")
    cat("    batch_evaluate - Estimate QST for batch of repeats with same QST/VE_ratio\n")
    cat("\nArguments:\n")
    cat("  input: For trait mode: path to observed stats RData\n")
    cat("         For neutral mode: FST value\n")
    cat("         For batch_neutral mode: path to FST file\n")
    cat("         For evaluate mode: QST value (e.g., 0.75)\n")
    cat("         For batch_evaluate mode: 'start_id_n_repeats_qst_value'\n")
    cat("  ratioVext_or_ve_ratio: For trait mode: ignored (read from input file)\n")
    cat("                      For neutral/batch_neutral mode: ratioVext extracted from true trait observations\n")
    cat("                      For evaluate/batch_evaluate mode: V_E/V_G ratio\n")
    cat("  output_file: Path to save result\n")
    cat("  num_sim: Number of ABC simulations (default: 100000)\n")
    cat("  summary_stats: Can be one of:\n")
    cat("                 - Comma-separated stats: 'QST,F_within_pop' (default)\n")
    cat("                 - Path to file with tab-separated combinations\n")
    cat("                 - 'all' to generate uncorrelated combinations from simulation\n")
    quit(status = 1)
  }
  
  mode <- args[1]
  input <- args[2]
  ratioVext_arg <- args[3]
  output_file <- args[4]
  num_sim <- if (length(args) >= 5) as.numeric(args[5]) else 100000
  summary_stats_arg <- if (length(args) >= 6) args[6] else "QST,F_within_pop"
  
  # Parse summary stats argument - can be comma-string, file path, or "all"
  summary_stat_combos <- parse_summary_stats(summary_stats_arg)
  is_multi_combo <- !is.null(summary_stat_combos) && length(summary_stat_combos) > 1
  
  # For single combination, extract the vector
  if (!is.null(summary_stat_combos) && length(summary_stat_combos) == 1) {
    summary_stats <- summary_stat_combos[[1]]
  } else if (is.null(summary_stat_combos)) {
    # "all" was specified - will generate from simulation
    summary_stats <- c("QST", "ratioVext")  # Default for now
  } else {
    summary_stats <- summary_stat_combos[[1]]  # Use first for display
  }
  
  cat("Mode:", mode, "\n")
  cat("Input:", input, "\n")
  cat("Num simulations:", num_sim, "\n")
  if (is_multi_combo) {
    cat("Summary stats:", length(summary_stat_combos), "combinations\n")
  } else {
    cat("Summary stats:", paste(summary_stats, collapse = ", "), "\n")
  }
  cat("Using", num_cores, "cores\n")
  
  start_time <- Sys.time()
  
  if (mode == "trait") {
    # Load observed summary statistics
    load(input)  # Expects 'obs_stats' variable
    cat("Ext SD from obs_stats:", obs_stats['ext_sd'], "\n")
    
    # Run ABC - only save QST to minimize disk usage
    qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
    
    result <- list(
      mode = "trait",
      qst = qst_estimate
    )
    
  } else if (mode == "neutral") {
    fst_value <- as.numeric(input)
    
    # Get ext_sd - either a numeric value or path to obs_stats file
    if (file.exists(ratioVext_arg)) {
      # Load from obs_stats file
      load(ratioVext_arg)  # Loads obs_stats
      ratioVext <- as.numeric(obs_stats["ratioVext"])
      cat("RatioVext from obs_stats file:", ratioVext, "\n")
    } else {
      ratioVext <- as.numeric(ratioVext_arg)
      cat("RatioVext from argument:", ratioVext, "\n")
    }
    
    # Set seed based on FST value for reproducibility
    seed <- round(fst_value * 1e6) %% .Machine$integer.max
    
    # Generate neutral observed stats
    obs_stats <- generate_neutral_obs_stats(fst_value, ratioVext, seed)
    
    if (is.null(obs_stats)) {
      result <- list(mode = "neutral", fst = fst_value, qst = NA)
    } else {
      # Run ABC - only save QST to minimize disk usage
      qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
      result <- list(
        mode = "neutral",
        fst = fst_value,
        qst = qst_estimate
      )
    }
    
  } else if (mode == "batch_neutral") {
    # Batch mode: process multiple FST values in one job
    # Input format: "start_idx:end_idx:fst_file" OR just a file path
    
    ratioVext <- as.numeric(ratioVext_arg)
    cat("RatioVext:", ratioVext, "\n")
    
    # Parse input - either "start:end:file" or just file path
    if (grepl(":", input)) {
      parts <- strsplit(input, ":")[[1]]
      start_idx <- as.integer(parts[1])
      end_idx <- as.integer(parts[2])
      fst_file <- parts[3]
      
      # Read FST values from file
      all_fst <- scan(fst_file, what = numeric(), quiet = TRUE)
      fst_values <- all_fst[start_idx:end_idx]
      cat("Processing FST values", start_idx, "to", end_idx, "from", fst_file, "\n")
    } else {
      # Input is a file with FST values
      fst_values <- scan(input, what = numeric(), quiet = TRUE)
      cat("Processing", length(fst_values), "FST values from", input, "\n")
      
      # Extract batch number from filename to ensure unique seeds across batches
      batch_num <- suppressWarnings(as.integer(gsub(".*fst_batch_([0-9]+)\\.txt.*", "\\1", input)))
      if (is.na(batch_num)) batch_num <- 1
      start_idx <- (batch_num - 1) * 10000 + 1
    }
    
    cat("Number of FST values:", length(fst_values), "\n")
    
    # Process all FST values - use multi-combo if multiple combinations
    if (is_multi_combo) {
      cat("Running ABC with", length(summary_stat_combos), "summary stat combinations\n")
      batch_results <- process_fst_batch_multi(fst_values, ratioVext, num_sim, summary_stat_combos, start_idx)
    } else {
      batch_results <- process_fst_batch(fst_values, ratioVext, num_sim, summary_stats, start_idx)
    }
    
    result <- list(
      mode = "batch_neutral",
      n_fst = length(fst_values),
      is_multi_combo = is_multi_combo,
      results = batch_results
    )
    
  } else if (mode == "evaluate") {
    # Module 2: Single adaptive QST evaluation
    # Input: QST value, ratioVext_arg: V_E/V_G ratio
    
    qst_value <- as.numeric(input)
    ve_ratio <- as.numeric(ratioVext_arg)
    
    cat("Adaptive QST:", qst_value, "\n")
    cat("V_E/V_G ratio:", ve_ratio, "\n")
    
    # Set seed for reproducibility
    seed <- (round(qst_value * 100) * 1000 + round(ve_ratio * 1000)) %% .Machine$integer.max
    
    # Generate adaptive observed stats
    obs_stats <- generate_adaptive_obs_stats(qst_value, ve_ratio, seed)
    
    if (is.null(obs_stats)) {
      result <- list(mode = "evaluate", true_qst = qst_value, ve_ratio = ve_ratio, 
                     estimated_qst = NA)
    } else {
      # Run ABC estimation
      qst_estimate <- estimate_qst_abc(obs_stats, num_sim, summary_stats)
      result <- list(
        mode = "evaluate",
        true_qst = qst_value,
        ve_ratio = ve_ratio,
        estimated_qst = qst_estimate
      )
    }
    
  } else if (mode == "batch_evaluate") {
    # Module 2: Batch adaptive QST evaluation
    # Input format: "n_repeats_qst_value" OR "start_id_n_repeats_qst_value" (using underscore separator)
    # ratioVext_arg: V_E/V_G ratio
    
    ve_ratio <- as.numeric(ratioVext_arg)
    cat("V_E/V_G ratio:", ve_ratio, "\n")
    
    # Parse input (use underscore as separator to avoid HTCondor colon issues)
    parts <- strsplit(input, "_")[[1]]
    if (length(parts) == 2) {
      n_repeats <- as.integer(parts[1])
      qst_value <- as.numeric(parts[2])
      start_id <- 1
    } else if (length(parts) == 3) {
      start_id <- as.integer(parts[1])
      n_repeats <- as.integer(parts[2])
      qst_value <- as.numeric(parts[3])
    } else {
      stop("Invalid input format for batch_evaluate. Expected 'n_repeats_qst' or 'start_id_n_repeats_qst'")
    }
    
    cat("Adaptive QST:", qst_value, "\n")
    cat("Number of repeats:", n_repeats, "\n")
    cat("Starting ID:", start_id, "\n")
    
    # Process all repeats - use multi-combo if multiple combinations
    if (is_multi_combo) {
      cat("Running ABC with", length(summary_stat_combos), "summary stat combinations\n")
      batch_results <- process_evaluate_batch_multi(n_repeats, qst_value, ve_ratio, num_sim, 
                                                    summary_stat_combos, start_id)
    } else {
      batch_results <- process_evaluate_batch(n_repeats, qst_value, ve_ratio, num_sim, 
                                              summary_stats, start_id)
    }
    
    result <- list(
      mode = "batch_evaluate",
      true_qst = qst_value,
      ve_ratio = ve_ratio,
      n_repeats = n_repeats,
      start_id = start_id,
      is_multi_combo = is_multi_combo,
      results = batch_results
    )
    
  } else {
    stop("Unknown mode: ", mode)
  }
  
  save(result, file = output_file)
  
  end_time <- Sys.time()
  elapsed_mins <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  if (mode == "batch_neutral") {
    n_valid <- sum(!is.na(result$results$qst))
    cat("Batch complete:", n_valid, "/", result$n_fst, "valid QST estimates\n")
    cat("Average time per FST:", round(elapsed_mins * 60 / result$n_fst, 2), "seconds\n")
  } else if (mode == "batch_evaluate") {
    n_valid <- sum(!is.na(result$results$estimated_qst))
    cat("Batch complete:", n_valid, "/", result$n_repeats, "valid QST estimates\n")
    cat("Average time per repeat:", round(elapsed_mins * 60 / result$n_repeats, 2), "seconds\n")
    cat("True QST:", result$true_qst, "| Mean estimated QST:", 
        round(mean(result$results$estimated_qst, na.rm = TRUE), 4), "\n")
  } else if (mode == "evaluate") {
    cat("True QST:", result$true_qst, "| Estimated QST:", result$estimated_qst, "\n")
  } else {
    cat("QST estimate:", result$qst, "\n")
  }
  cat("Completed in", round(elapsed_mins, 2), "minutes\n")
  cat("Saved to:", output_file, "\n")
}
