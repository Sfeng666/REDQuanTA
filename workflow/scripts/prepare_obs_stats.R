#!/usr/bin/env Rscript
# Prepare observed summary statistics for a trait
# This generates the obs_stats needed for ABC estimation
#
# Usage:
#   Rscript prepare_obs_stats.R <sample_structure_csv> <trait_values_csv> <trait_id> <output_file> [ext_sd_file]
#
# Arguments:
#   1. sample_structure_csv: Path to sample structure CSV
#   2. trait_values_csv: Path to trait values CSV
#   3. trait_id: Trait ID to process
#   4. output_file: Path to save obs_stats RData
#   5. ext_sd_file: (optional) Path to save ext_sd value as text

# No additional libraries needed - only base R functions used

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript prepare_obs_stats.R <sample_structure_csv> <trait_values_csv> <trait_id> <output_file> [ext_sd_file]\n")
  quit(status = 1)
}

sample_structure_path <- args[1]
trait_values_path <- args[2]
trait_id <- args[3]
output_file <- args[4]
ext_sd_file <- if (length(args) >= 5) args[5] else NULL

# Read sample structure
sample_structure <- read.csv(sample_structure_path)

# Determine sample structure parameters
num_pop <- length(unique(sample_structure$population))
num_ind <- length(unique(sample_structure$strain[sample_structure$population == 1]))
num_rep <- length(unique(sample_structure$replicate[sample_structure$population == 1 & 
                                                      sample_structure$strain == 1]))

# Read trait values
trait_values <- read.csv(trait_values_path)

# Find the trait row
trait_row <- trait_values[trait_values$trait_id == trait_id, ]
if (nrow(trait_row) == 0) {
  stop(paste("Trait ID not found:", trait_id))
}

# Extract chromosome info
chr <- as.character(trait_row$chr)

# Extract trait values (columns 3 onwards are sample values)
values <- as.numeric(trait_row[, 3:ncol(trait_row)])

# Build data frame matching sample structure
data <- data.frame(
  trait = values,
  pop = factor(sample_structure$population),
  line = factor(paste0(sample_structure$population, "_", sample_structure$strain)),
  rep = factor(sample_structure$replicate)
)

# Calculate variance components
calc_variance_components_mom <- function(data, num_pop, num_ind, num_rep) {
  overall_mean <- mean(data$trait)
  pop_means <- tapply(data$trait, data$pop, mean)
  line_means <- tapply(data$trait, data$line, mean)
  
  SS_among <- sum((pop_means - overall_mean)^2) * num_ind * num_rep
  SS_within <- sum((line_means - rep(pop_means, each = num_ind))^2) * num_rep
  SS_residual <- sum((data$trait - rep(line_means, each = num_rep))^2)
  
  DF_among <- num_pop - 1
  DF_within <- num_pop * (num_ind - 1)
  DF_residual <- num_pop * num_ind * (num_rep - 1)
  
  MS_among <- SS_among / DF_among
  MS_within <- SS_within / DF_within
  MS_residual <- SS_residual / DF_residual
  
  var_among <- max((MS_among - MS_within) / (num_ind * num_rep), 0)
  var_within <- max((MS_within - MS_residual) / num_rep, 0)
  var_residual <- max(MS_residual, 0)
  
  return(list(var_among = var_among, var_within = var_within, var_residual = var_residual))
}

# Calculate summary statistics
variance_components <- calc_variance_components_mom(data, num_pop, num_ind, num_rep)
among_pop_var <- variance_components$var_among
within_pop_var <- variance_components$var_within
ext_var <- variance_components$var_residual

total_genetic_var <- among_pop_var + within_pop_var

if (total_genetic_var == 0) {
  QST <- 0
  ratioVext <- 0
  F_among_pop <- 0
  F_within_pop <- 0
} else {
  QST <- among_pop_var / (among_pop_var + 2 * within_pop_var)
  ratioVext <- ext_var / total_genetic_var
  F_among_pop <- if (within_pop_var == 0) 0 else among_pop_var / within_pop_var
  F_within_pop <- if (ext_var == 0) 0 else within_pop_var / ext_var
}

obs_stats <- c(
  among_pop_sd = sqrt(among_pop_var),
  within_pop_sd = sqrt(within_pop_var),
  ext_sd = sqrt(ext_var),
  QST = QST,
  ratioVext = ratioVext,
  F_among_pop = F_among_pop,
  F_within_pop = F_within_pop
)

# Also save metadata
trait_meta <- list(
  trait_id = trait_id,
  chr = chr,
  num_pop = num_pop,
  num_ind = num_ind,
  num_rep = num_rep
)

save(obs_stats, trait_meta, file = output_file)
cat("Trait:", trait_id, "\n")
cat("Chr:", chr, "\n")
cat("Obs stats saved to:", output_file, "\n")
cat("ratioVext:", obs_stats['ratioVext'], "\n")
cat("among_pop_sd:", obs_stats['among_pop_sd'], "\n")

if (!is.null(ext_sd_file)) {
  writeLines(as.character(obs_stats['ext_sd']), ext_sd_file)
  cat("ext_sd saved to:", ext_sd_file, "\n")
}

