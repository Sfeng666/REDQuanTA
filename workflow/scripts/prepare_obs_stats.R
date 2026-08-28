#!/usr/bin/env Rscript
# Prepare observed summary statistics for a trait (Batch & Single trait optimized version)
#
# Usage:
#   Rscript prepare_obs_stats.R --batch <sample_structure_csv> <trait_values_csv> <list_file> <results_dir> <summary_tsv>
#   Rscript prepare_obs_stats.R <sample_structure_csv> <trait_values_csv> <trait_id> <output_file> [ext_sd_file] [ratioVext_file]

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage:\n")
  cat("  Rscript prepare_obs_stats.R --batch <sample_structure_csv> <trait_values_csv> <list_file> <results_dir> <summary_tsv>\n")
  cat("  Rscript prepare_obs_stats.R <sample_structure_csv> <trait_values_csv> <trait_id> <output_file> [ext_sd_file] [ratioVext_file]\n")
  quit(status = 1)
}

is_batch <- args[1] == "--batch"

if (is_batch) {
  if (length(args) < 6) {
    cat("Error: Batch mode requires 5 arguments after --batch\n")
    quit(status = 1)
  }
  sample_structure_path <- args[2]
  trait_values_path <- args[3]
  list_file_path <- args[4]
  results_dir <- args[5]
  summary_tsv_path <- args[6]
} else {
  if (length(args) < 4) {
    cat("Error: Single mode requires at least 4 arguments\n")
    quit(status = 1)
  }
  sample_structure_path <- args[1]
  trait_values_path <- args[2]
  trait_id_single <- args[3]
  output_file_single <- args[4]
  ext_sd_file_single <- if (length(args) >= 5) args[5] else NULL
  ratioVext_file_single <- if (length(args) >= 6) args[6] else NULL
}

# 1. Read sample structure
sample_structure <- read.csv(sample_structure_path)

# Determine sample structure parameters
num_pop <- length(unique(sample_structure$population))
num_ind <- length(unique(sample_structure$strain[sample_structure$population == 1]))
num_rep <- length(unique(sample_structure$replicate[sample_structure$population == 1 & 
                                                      sample_structure$strain == 1]))

# Factors for grouping
pop_factor <- factor(sample_structure$population)
line_factor <- factor(paste0(sample_structure$population, "_", sample_structure$strain))
rep_factor <- factor(sample_structure$replicate)

pop_levels <- levels(pop_factor)
line_levels <- levels(line_factor)

# Degrees of freedom
DF_among <- num_pop - 1
DF_within <- num_pop * (num_ind - 1)
DF_residual <- num_pop * num_ind * (num_rep - 1)

# Pre-map population to each line (since each line belongs to exactly one pop)
line_to_pop <- sapply(line_levels, function(l) {
  as.character(pop_factor[line_factor == l][1])
})

# 2. Read trait values
trait_values <- read.csv(trait_values_path, check.names = FALSE)

# Filter trait list to process
if (is_batch) {
  # Read list file
  trait_ids <- readLines(list_file_path)
  trait_ids <- trait_ids[nzchar(trait_ids)]
  
  # Filter rows in trait_values
  trait_rows <- trait_values[trait_values$trait_id %in% trait_ids, , drop = FALSE]
} else {
  trait_rows <- trait_values[trait_values$trait_id == trait_id_single, , drop = FALSE]
  if (nrow(trait_rows) == 0) {
    stop(paste("Trait ID not found:", trait_id_single))
  }
}

if (nrow(trait_rows) == 0) {
  cat("No traits found to process.\n")
  if (is_batch) {
    # Write empty summary tsv
    writeLines(character(0), summary_tsv_path)
  }
  quit(status = 0)
}

# Extract chromosomes and values matrix
chrs <- as.character(trait_rows$chr)
trait_ids_processed <- as.character(trait_rows$trait_id)
Y <- as.matrix(trait_rows[, 3:ncol(trait_rows), drop = FALSE])
storage.mode(Y) <- "double"

# 3. Vectorized ANOVA calculations for all genes in Y
overall_mean <- rowMeans(Y)

pop_means <- matrix(0, nrow = nrow(Y), ncol = length(pop_levels))
colnames(pop_means) <- pop_levels
for (p in pop_levels) {
  pop_means[, p] <- rowMeans(Y[, pop_factor == p, drop = FALSE])
}

line_means <- matrix(0, nrow = nrow(Y), ncol = length(line_levels))
colnames(line_means) <- line_levels
for (l in line_levels) {
  line_means[, l] <- rowMeans(Y[, line_factor == l, drop = FALSE])
}

# Expand to full dimensions for sum of squares
pop_means_exp_samples <- pop_means[, as.character(pop_factor), drop = FALSE]
pop_means_exp_lines <- pop_means[, line_to_pop, drop = FALSE]
line_means_exp_samples <- line_means[, as.character(line_factor), drop = FALSE]

SS_among <- rowSums((pop_means_exp_samples - overall_mean)^2)
SS_within <- rowSums((line_means_exp_samples - pop_means_exp_samples)^2)
SS_residual <- rowSums((Y - line_means_exp_samples)^2)

SS_total <- rowSums((Y - overall_mean)^2)
MS_total <- SS_total / max(num_pop * num_ind * num_rep - 1L, 1L)
# Ridge floor (alpha = 0.1), same formula as qst_abc_sim.R FLOOR_POLICY=ridge_floor.
tf <- ifelse(MS_total > 0, 1e-8 * MS_total, 1e-8)
alpha_floor <- 0.1

MS_among <- SS_among / DF_among
MS_within <- SS_within / DF_within
MS_residual <- SS_residual / DF_residual

var_among_raw <- (MS_among - MS_within) / (num_ind * num_rep)
var_within_raw <- (MS_within - MS_residual) / num_rep
var_residual <- pmax(MS_residual, 0)

noise_among <- MS_within / (num_ind * num_rep)
noise_within <- MS_residual / num_rep
lam_a <- pmax(alpha_floor * noise_among, tf)
lam_w <- pmax(alpha_floor * noise_within, tf)
var_among <- sqrt(pmax(var_among_raw, 0)^2 + lam_a^2)
var_within <- sqrt(pmax(var_within_raw, 0)^2 + lam_w^2)
var_residual <- pmax(MS_residual, 0)

both_var_negative <- rep(FALSE, length(var_among_raw))

# Scale statistics
total_genetic <- (var_among + 2 * var_within) / 2
total_var <- total_genetic + var_residual

qst_denom <- var_among + 2 * var_within
QST <- ifelse(qst_denom == 0, 0, var_among / qst_denom)
ratioVext <- ifelse(total_genetic == 0, 0, var_residual / total_genetic)

F_among_pop <- ifelse(total_genetic == 0 | var_within == 0, 0, var_among / var_within)
F_within_pop <- ifelse(total_genetic == 0 | var_residual == 0, 0, var_within / var_residual)
ratioVbetweenVext <- ifelse(var_residual == 0, 0, var_among / var_residual)
ratioVbetweenVtotal <- ifelse(total_var == 0, 0, var_among / total_var)
ratioVwithinVtotal <- ifelse(total_var == 0, 0, var_within / total_var)
ratioVextVtotal <- ifelse(total_var == 0, 0, var_residual / total_var)

# Pre-allocate output variables/structures
summary_rows <- character(nrow(Y))

# 4. Save individual RData files
for (i in seq_len(nrow(Y))) {
  trait_id <- trait_ids_processed[i]
  chr <- chrs[i]
  
  obs_stats <- c(
    among_pop_sd = as.numeric(sqrt(var_among[i])),
    within_pop_sd = as.numeric(sqrt(var_within[i])),
    ext_sd = as.numeric(sqrt(var_residual[i])),
    QST = as.numeric(QST[i]),
    ratioVext = as.numeric(ratioVext[i]),
    F_among_pop = as.numeric(F_among_pop[i]),
    F_within_pop = as.numeric(F_within_pop[i]),
    ratioVbetweenVext = as.numeric(ratioVbetweenVext[i]),
    ratioVbetweenVtotal = as.numeric(ratioVbetweenVtotal[i]),
    ratioVwithinVtotal = as.numeric(ratioVwithinVtotal[i]),
    ratioVextVtotal = as.numeric(ratioVextVtotal[i]),
    both_var_negative = as.numeric(both_var_negative[i])
  )
  
  trait_meta <- list(
    trait_id = trait_id,
    chr = chr,
    num_pop = num_pop,
    num_ind = num_ind,
    num_rep = num_rep
  )
  
  if (is_batch) {
    trait_dir <- file.path(results_dir, paste0("trait_", trait_id))
    dir.create(trait_dir, showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(trait_dir, paste0(trait_id, "_obs_stats.RData"))
    
    summary_rows[i] <- paste(trait_id, ratioVext[i], sep = "\t")
  } else {
    output_file <- output_file_single
  }
  
  save(obs_stats, trait_meta, file = output_file)
}

# 5. Output summary TSV or auxiliary files
if (is_batch) {
  writeLines(summary_rows, summary_tsv_path)
  cat("Batch processing complete. Summary saved to:", summary_tsv_path, "\n")
} else {
  cat("Trait:", trait_id_single, "\n")
  cat("Chr:", chrs[1], "\n")
  cat("Obs stats saved to:", output_file_single, "\n")
  cat("ratioVext:", ratioVext[1], "\n")
  cat("among_pop_sd:", sqrt(var_among[1]), "\n")
  
  if (!is.null(ext_sd_file_single)) {
    writeLines(as.character(sqrt(var_residual[1])), ext_sd_file_single)
    cat("ext_sd saved to:", ext_sd_file_single, "\n")
  }
  
  if (!is.null(ratioVext_file_single)) {
    writeLines(as.character(ratioVext[1]), ratioVext_file_single)
    cat("ratioVext saved to:", ratioVext_file_single, "\n")
  }
}
