#!/usr/bin/env Rscript
# Compare per-trait QST adaptive calls to a chromosome-matched FST threshold.
#
# Usage:
#   Rscript compare_adaptive_fst_vs_qst.R \
#     --results-root <dir with trait-type subdirs> \
#     --output-dir <out> \
#     --fst-autosomes data/example/qst_neutral_autosomes.txt \
#     --fst-chrx data/example/qst_neutral_chrX.txt \
#     [--trait-types protein,mRNA_expression] \
#     [--labels-tsv labels.tsv] \
#     [--threshold 0.95]
#
# Each trait-type directory must contain summary/qst_results_abc.csv or
# qst_results_abc.csv. If --trait-types is omitted, every subdirectory with
# that file is used. Optional labels.tsv has columns trait_type, label.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function(argv) {
  out <- list(
    results_root = NA_character_,
    output_dir = NA_character_,
    fst_autosomes = NA_character_,
    fst_chrx = NA_character_,
    trait_types = NA_character_,
    labels_tsv = NA_character_,
    threshold = 0.95
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[i]
    val <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (key == "--results-root") out$results_root <- val
    else if (key == "--output-dir") out$output_dir <- val
    else if (key == "--fst-autosomes") out$fst_autosomes <- val
    else if (key == "--fst-chrx") out$fst_chrx <- val
    else if (key == "--trait-types") out$trait_types <- val
    else if (key == "--labels-tsv") out$labels_tsv <- val
    else if (key == "--threshold") out$threshold <- as.numeric(val)
    else if (key %in% c("-h", "--help")) {
      cat("See header of compare_adaptive_fst_vs_qst.R for usage.\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 2L
  }
  if (is.na(out$results_root) || is.na(out$output_dir) ||
      is.na(out$fst_autosomes) || is.na(out$fst_chrx)) {
    stop("Required: --results-root --output-dir --fst-autosomes --fst-chrx")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
results_root <- args$results_root
out_compare <- args$output_dir
threshold_p <- args$threshold
dir.create(out_compare, recursive = TRUE, showWarnings = FALSE)

discover_types <- function(root) {
  kids <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  keep <- vapply(kids, function(d) {
    file.exists(file.path(d, "summary", "qst_results_abc.csv")) ||
      file.exists(file.path(d, "qst_results_abc.csv"))
  }, logical(1))
  basename(kids[keep])
}

trait_types <- if (!is.na(args$trait_types)) {
  strsplit(args$trait_types, ",", fixed = TRUE)[[1]]
} else {
  discover_types(results_root)
}
trait_types <- trimws(trait_types)
if (!length(trait_types)) stop("No trait-type directories with qst_results_abc.csv under ", results_root)

type_labels <- setNames(gsub("_", " ", trait_types), trait_types)
if (!is.na(args$labels_tsv) && file.exists(args$labels_tsv)) {
  lab <- fread(args$labels_tsv)
  if (!all(c("trait_type", "label") %in% names(lab))) {
    stop("labels.tsv must have columns trait_type and label")
  }
  type_labels[lab$trait_type] <- lab$label
}

normalize_chr <- function(chr) {
  x <- tolower(trimws(as.character(chr)))
  fifelse(x %in% c("chrx", "x", "chr_x"), "chrX",
          fifelse(x %in% c("autosomes", "auto", "chr"), "autosomes", x))
}

load_reference_fst_pool <- function(path) {
  v <- scan(path, quiet = TRUE)
  v <- v[is.finite(v) & v >= 0]
  if (!length(v)) stop("Empty or invalid FST pool: ", path)
  v
}

fst_auto <- load_reference_fst_pool(args$fst_autosomes)
fst_chrx <- load_reference_fst_pool(args$fst_chrx)

global_thr <- data.table(
  chr = c("autosomes", "chrX"),
  fst_pool_file = c(basename(args$fst_autosomes), basename(args$fst_chrx)),
  n_fst_values = c(length(fst_auto), length(fst_chrx)),
  fst_threshold_95 = c(
    as.numeric(quantile(fst_auto, threshold_p, na.rm = TRUE)),
    as.numeric(quantile(fst_chrx, threshold_p, na.rm = TRUE))
  ),
  threshold_percentile = threshold_p
)
fwrite(global_thr, file.path(out_compare, "global_fst_threshold_95_by_chr.csv"))
thr_lookup <- setNames(global_thr$fst_threshold_95, global_thr$chr)

find_qst_results <- function(dir_path) {
  candidates <- c(
    file.path(dir_path, "summary", "qst_results_abc.csv"),
    file.path(dir_path, "qst_results_abc.csv")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) NA_character_ else hit
}

summarize_adaptive <- function(dt, col) {
  ok <- dt[!is.na(get(col)) & get(col) %in% c("yes", "no")]
  n <- nrow(ok)
  n_yes <- sum(ok[[col]] == "yes", na.rm = TRUE)
  data.table(
    n_classifiable = n,
    n_adaptive = n_yes,
    n_non_adaptive = sum(ok[[col]] == "no", na.rm = TRUE),
    pct_adaptive = if (n) round(100 * n_yes / n, 3) else NA_real_
  )
}

per_type_summaries <- list()

for (tt in trait_types) {
  dir_path <- file.path(results_root, tt)
  qst_path <- find_qst_results(dir_path)
  if (is.na(qst_path)) {
    warning("Skipping ", tt, ": no qst_results_abc.csv")
    next
  }

  summary_dir <- file.path(dir_path, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== ", type_labels[tt], " (", tt, ") ===\n", sep = "")
  dt <- fread(qst_path)
  dt[, chr := normalize_chr(chr)]
  dt[, fst_threshold_consistent := thr_lookup[chr]]
  dt[is.na(fst_threshold_consistent), adaptive_fst_based := NA_character_]
  dt[!is.na(fst_threshold_consistent), adaptive_fst_based := fifelse(
    is.na(QST), NA_character_,
    fifelse(QST > fst_threshold_consistent, "yes", "no")
  )]
  dt[, trait_type := tt]
  dt[, fst_threshold_percentile := threshold_p]
  dt[, adaptive_qst_based := adaptive]

  fwrite(dt, file.path(summary_dir, "qst_results_adaptive_fst_vs_qst.csv"))

  s_fst <- summarize_adaptive(dt, "adaptive_fst_based")
  s_qst <- summarize_adaptive(dt, "adaptive_qst_based")
  agree <- dt[
    !is.na(adaptive_fst_based) & !is.na(adaptive_qst_based) &
      adaptive_fst_based %in% c("yes", "no") & adaptive_qst_based %in% c("yes", "no"),
    .(agree = adaptive_fst_based == adaptive_qst_based)
  ]
  pct_agree <- if (nrow(agree)) round(100 * mean(agree$agree), 3) else NA_real_

  by_chr <- dt[
    !is.na(adaptive_fst_based) & adaptive_fst_based %in% c("yes", "no"),
    .(
      n_classifiable_fst = .N,
      n_adaptive_fst = sum(adaptive_fst_based == "yes"),
      n_classifiable_qst = sum(!is.na(adaptive_qst_based) & adaptive_qst_based %in% c("yes", "no")),
      n_adaptive_qst = sum(adaptive_qst_based == "yes", na.rm = TRUE)
    ),
    by = chr
  ]
  by_chr[, `:=`(
    pct_adaptive_fst = round(100 * n_adaptive_fst / n_classifiable_fst, 3),
    pct_adaptive_qst = fifelse(
      n_classifiable_qst > 0,
      round(100 * n_adaptive_qst / n_classifiable_qst, 3),
      NA_real_
    )
  )]
  by_chr[, trait_type := tt]
  fwrite(by_chr, file.path(summary_dir, "adaptive_fst_vs_qst_by_chr.csv"))

  type_sum <- data.table(
    trait_type = tt,
    trait_label = unname(type_labels[tt]),
    threshold_percentile = threshold_p,
    fst_threshold_autosomes = thr_lookup["autosomes"],
    fst_threshold_chrX = thr_lookup["chrX"],
    n_traits_autosomes = sum(dt$chr == "autosomes", na.rm = TRUE),
    n_traits_chrX = sum(dt$chr == "chrX", na.rm = TRUE),
    n_traits_total = nrow(dt),
    n_classifiable_fst = s_fst$n_classifiable,
    n_adaptive_fst = s_fst$n_adaptive,
    pct_adaptive_fst = s_fst$pct_adaptive,
    n_classifiable_qst = s_qst$n_classifiable,
    n_adaptive_qst = s_qst$n_adaptive,
    pct_adaptive_qst = s_qst$pct_adaptive,
    pct_both_methods_agree = pct_agree,
    pct_adaptive_delta = s_fst$pct_adaptive - s_qst$pct_adaptive
  )
  fwrite(type_sum, file.path(summary_dir, "adaptive_fst_vs_qst_summary.csv"))
  per_type_summaries[[tt]] <- type_sum

  cat("  QST-based adaptive:", s_qst$n_adaptive, "/", s_qst$n_classifiable,
      " (", s_qst$pct_adaptive, "%)\n", sep = "")
}

all_sum <- rbindlist(per_type_summaries)
fwrite(all_sum, file.path(out_compare, "adaptive_fst_vs_qst_by_trait_type.csv"))
fwrite(
  all_sum[, .(
    trait_type, trait_label,
    fst_threshold_autosomes, fst_threshold_chrX,
    pct_adaptive_fst, pct_adaptive_qst, pct_adaptive_delta, pct_both_methods_agree
  )],
  file.path(out_compare, "adaptive_fst_vs_qst_comparison_table.csv")
)

theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

thr_label <- sprintf(
  "autosomes Fst %s = %.3f; chrX Fst %s = %.3f",
  paste0(round(threshold_p * 100), "%"), thr_lookup["autosomes"],
  paste0(round(threshold_p * 100), "%"), thr_lookup["chrX"]
)

plot_long <- melt(
  all_sum,
  id.vars = c("trait_type", "trait_label"),
  measure.vars = c("pct_adaptive_fst", "pct_adaptive_qst"),
  variable.name = "method",
  value.name = "pct_adaptive"
)
plot_long[, method := fifelse(
  method == "pct_adaptive_fst",
  "FST-based (QST > reference neutral Fst 95%ile by chr)",
  "QST-based (QST > per-trait neutral Qst 95%ile)"
)]

p_bar <- ggplot(plot_long, aes(x = trait_label, y = pct_adaptive, fill = method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f%%", pct_adaptive)),
    position = position_dodge(width = 0.75), vjust = -0.3, size = 2.6
  ) +
  scale_fill_manual(values = c("#4E79A7", "#E15759"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Adaptive trait proportion by threshold method",
    subtitle = paste0(thr_label, " (classifiable traits only)"),
    x = NULL, y = "Percent adaptive"
  ) +
  theme_pub()

ggsave(file.path(out_compare, "adaptive_proportion_fst_vs_qst_methods.pdf"), p_bar, width = 12, height = 5)

p_delta <- ggplot(all_sum, aes(x = trait_label, y = pct_adaptive_delta)) +
  geom_col(fill = "#59A14F", width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(
    title = "Change in adaptive % (FST-based minus QST-based)",
    subtitle = thr_label, x = NULL, y = "Percentage-point difference"
  ) +
  theme_pub()

ggsave(file.path(out_compare, "adaptive_proportion_fst_minus_qst_delta.pdf"), p_delta, width = 11, height = 4.5)

cat("Wrote FST vs QST comparison to", out_compare, "\n")
