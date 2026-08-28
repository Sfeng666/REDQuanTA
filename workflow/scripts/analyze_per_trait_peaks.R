#!/usr/bin/env Rscript
# Per-trait neutral QST peak analysis (memory-safe; no full-pool reload).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop("Usage: Rscript analyze_per_trait_peaks.R <results_dir> [out_dir] [pdf_only]")
}
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
results_dir <- args[1]
out_dir <- if (length(args) >= 2) args[2] else file.path(results_dir, "plots")
pdf_only <- length(args) >= 3 && tolower(args[3]) %in% c("pdf_only", "pdf-only", "pdf")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "neutral_qst_memory_safe.R"))

find_density_peaks <- function(x, from = 0, to = 1, n = 512, adjust = 1.05,
                               min_prominence_frac = 0.15) {
  x <- x[is.finite(x)]
  if (length(x) < 10) {
    return(data.table(peak_x = numeric(), peak_y = numeric()))
  }
  d <- density(x, from = from, to = to, n = n, adjust = adjust)
  idx <- which(diff(sign(diff(d$y))) == -2) + 1
  if (length(idx) == 0) {
    return(data.table(peak_x = numeric(), peak_y = numeric()))
  }
  keep <- (d$y[idx] / max(d$y)) >= min_prominence_frac
  data.table(peak_x = d$x[idx[keep]], peak_y = d$y[idx[keep]])
}

has_peak_near <- function(peaks, center, halfwidth) {
  nrow(peaks) > 0 && any(abs(peaks$peak_x - center) <= halfwidth)
}

load_trait_qst_fst <- function(trait_dir) {
  files <- sort(list.files(trait_dir, "neutral_batch_.*RData", full.names = TRUE))
  chunks <- lapply(files, function(f) {
    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    e$result$results
  })
  rbindlist(chunks)
}

qst_path <- file.path(results_dir, "summary", "qst_results_abc.csv")
if (!file.exists(qst_path)) {
  qst_path <- file.path(results_dir, "qst_results_abc.csv")
}
traits <- fread(qst_path)
n_traits <- nrow(traits)
res_pooled <- numeric()
res_with_peak <- numeric()
res_without_peak <- numeric()

# Determine how many cores to use
condor_cores <- suppressWarnings(as.integer(Sys.getenv("_CONDOR_NPROCS", "4")))
num_cores <- suppressWarnings(as.integer(Sys.getenv("QST_ABC_CORES", as.character(condor_cores))))
if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
num_cores <- min(num_cores, parallel::detectCores(), 8L)

cat("Analyzing per-trait neutral peaks in parallel using", num_cores, "cores...\n")

worker_results <- parallel::mclapply(seq_len(n_traits), function(i) {
  tid <- traits$trait_id[i]
  trait_dir <- file.path(results_dir, paste0("trait_", tid))
  
  t <- tryCatch(load_trait_qst_fst(trait_dir), error = function(e) NULL)
  if (is.null(t) || nrow(t) == 0) {
    return(NULL)
  }
  
  # Resolve ratioVext (rv) from traits table if possible to avoid loading obs_stats file
  rv <- NA_real_
  if ("anova_ratioVext" %in% names(traits)) {
    rv <- as.numeric(traits$anova_ratioVext[i])
  } else if ("ratioVext" %in% names(traits)) {
    rv <- as.numeric(traits$ratioVext[i])
  }
  
  # Fallback to loading RData only if not found in CSV
  if (is.na(rv)) {
    obs_path <- file.path(trait_dir, paste0(tid, "_obs_stats.RData"))
    if (file.exists(obs_path)) {
      e <- new.env(parent = emptyenv())
      load(obs_path, envir = e)
      rv <- as.numeric(e$obs_stats["ratioVext"])
    }
  }

  qst <- t$qst
  peaks <- find_density_peaks(qst)
  near_04 <- has_peak_near(peaks, 0.4, 0.08)

  row_dt <- data.table(
    trait_id = tid,
    ratioVext = rv,
    prior_QST = traits$prior_QST[i],
    n_qst = nrow(t),
    median_qst = median(qst),
    median_fst = median(t$fst),
    mean_qst_low_fst = mean(t[fst >= 0 & fst <= 0.1, qst]),
    n_peaks = nrow(peaks),
    peak_x = if (nrow(peaks)) paste(sprintf("%.3f", peaks$peak_x), collapse = ";") else NA_character_,
    peak_y = if (nrow(peaks)) paste(sprintf("%.3f", peaks$peak_y), collapse = ";") else NA_character_,
    peak_near_0.4 = near_04,
    peak_near_0.0 = has_peak_near(peaks, 0.05, 0.05),
    frac_qst_0.3_0.5 = mean(qst >= 0.3 & qst <= 0.5)
  )

  list(
    row = row_dt,
    qst = qst,
    near_04 = near_04
  )
}, mc.cores = num_cores, mc.preschedule = TRUE)

# Filter out NULL results
worker_results <- Filter(Negate(is.null), worker_results)

if (length(worker_results) == 0) {
  stop("No valid neutral peak analysis results compiled.")
}

# Bind rows to construct per_trait table
rows <- lapply(worker_results, function(x) x$row)
per_trait <- rbindlist(rows)
per_trait[, mass_0.3_0.5 := frac_qst_0.3_0.5 * n_qst]
per_trait[, pct_of_pooled_middle := 100 * mass_0.3_0.5 / sum(mass_0.3_0.5)]
per_trait[, rv_tert := cut(
  ratioVext,
  quantile(ratioVext, c(0, 1 / 3, 2 / 3, 1)),
  include.lowest = TRUE,
  labels = c("low", "mid", "high")
)]

fwrite(per_trait, file.path(out_dir, "per_trait_neutral_qst_peaks.csv"))

cat("\nTraits with ~0.4 peak:", sum(per_trait$peak_near_0.4), "/", nrow(per_trait), "\n")
cat("Traits without ~0.4 peak:", sum(!per_trait$peak_near_0.4), "\n")

p2 <- ggplot(per_trait, aes(x = frac_qst_0.3_0.5, fill = peak_near_0.4)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.15, position = "identity", alpha = 0.6) +
  scale_fill_manual(
    values = c("TRUE" = "#E15759", "FALSE" = "#4E79A7"),
    labels = c("TRUE" = "Has ~0.4 peak", "FALSE" = "No ~0.4 peak"),
    name = NULL
  ) +
  labs(
    title = "Share of each trait neutral QST in [0.3, 0.5]",
    subtitle = sprintf("%d/%d traits (%.1f%%) show a secondary density peak near 0.4",
                       sum(per_trait$peak_near_0.4), nrow(per_trait),
                       100 * mean(per_trait$peak_near_0.4)),
    x = "Fraction of neutral QST values in [0.3, 0.5]",
    y = "Number of traits"
  ) +
  theme_classic(base_size = 11)

ggsave(file.path(out_dir, "frac_qst_0.3_0.5_by_trait.pdf"), p2, width = 6.5, height = 4.2)
if (!pdf_only) {
  ggsave(file.path(out_dir, "frac_qst_0.3_0.5_by_trait.png"), p2, width = 6.5, height = 4.2, dpi = 600)
}

if (n_traits <= 500) {
  pick_traits_by <- function(cond, col, decreasing = FALSE, n = 3) {
    sub <- per_trait[eval(substitute(cond))]
    ord <- if (decreasing) order(-sub[[col]]) else order(sub[[col]])
    sub[ord][seq_len(min(n, nrow(sub))), trait_id]
  }
  example_ids <- c(
    pick_traits_by(peak_near_0.4 == TRUE, "ratioVext", decreasing = TRUE),
    pick_traits_by(peak_near_0.4 == FALSE, "ratioVext", decreasing = FALSE)
  )
  ex_dt <- rbindlist(lapply(example_ids, function(tid) {
    qst <- load_trait_qst_fst(file.path(results_dir, paste0("trait_", tid)))$qst
    rv <- per_trait[trait_id == tid, ratioVext]
    has <- per_trait[trait_id == tid, peak_near_0.4]
    data.table(
      qst = qst,
      label = sprintf("%s (rv=%.1f%s)", tid, rv, ifelse(has, ", ~0.4 peak", ", no ~0.4 peak"))
    )
  }))
  ex_dt[, label := factor(label, levels = unique(label))]

  p3 <- ggplot(ex_dt, aes(x = qst)) +
    geom_density(fill = "#4E79A7", colour = "#2F5597", alpha = 0.4, adjust = 1.05) +
    facet_wrap(~label, ncol = 2, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.5)) +
    labs(
      title = "Per-trait neutral QST: examples with and without a ~0.4 density peak",
      x = expression(italic(Q)[ST]), y = "Density"
    ) +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(size = 7.5))

  ggsave(file.path(out_dir, "per_trait_qst_examples.pdf"), p3, width = 7.5, height = 7)
  if (!pdf_only) {
    ggsave(file.path(out_dir, "per_trait_qst_examples.png"), p3, width = 7.5, height = 7, dpi = 600)
  }
}

for (pat in c("^peak_decomposition_by_trait_group\\.", "^neutral_qst_distribution_combined\\.",
              "^qst_distribution_publication(_pooled|_by_class)?\\.", "^qst_distribution_adaptive_only\\.")) {
  for (f in list.files(out_dir, pattern = pat, full.names = TRUE)) {
    unlink(f)
  }
}

cat("Per-trait analysis written to", out_dir, "\n")
