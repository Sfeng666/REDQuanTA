#!/usr/bin/env Rscript
# Compare QST distributions and adaptive proportions across trait-type collections.
#
# Usage:
#   Rscript plot_comparison_across_trait_types.R \
#     --results-root <dir with trait-type subdirs> \
#     --output-dir <out> \
#     [--trait-types protein,mRNA_expression] \
#     [--labels-tsv labels.tsv]
#
# Each trait-type directory must contain summary/qst_results_abc.csv.
# If --trait-types is omitted, every subdirectory with that file is used.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

parse_args <- function(argv) {
  out <- list(
    results_root = NA_character_,
    output_dir = NA_character_,
    trait_types = NA_character_,
    labels_tsv = NA_character_
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[i]
    val <- if (i < length(argv)) argv[i + 1L] else NA_character_
    if (key == "--results-root") out$results_root <- val
    else if (key == "--output-dir") out$output_dir <- val
    else if (key == "--trait-types") out$trait_types <- val
    else if (key == "--labels-tsv") out$labels_tsv <- val
    else if (key %in% c("-h", "--help")) {
      cat("See header of plot_comparison_across_trait_types.R for usage.\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 2L
  }
  if (is.na(out$results_root) || is.na(out$output_dir)) {
    stop("Required: --results-root --output-dir")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
results_root <- args$results_root
out_dir <- args$output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

discover_types <- function(root) {
  kids <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  keep <- vapply(kids, function(d) {
    file.exists(file.path(d, "summary", "qst_results_abc.csv"))
  }, logical(1))
  basename(kids[keep])
}

trait_types <- if (!is.na(args$trait_types)) {
  trimws(strsplit(args$trait_types, ",", fixed = TRUE)[[1]])
} else {
  discover_types(results_root)
}
if (!length(trait_types)) stop("No trait-type summaries found under ", results_root)

type_labels <- setNames(gsub("_", " ", trait_types), trait_types)
if (!is.na(args$labels_tsv) && file.exists(args$labels_tsv)) {
  lab <- fread(args$labels_tsv)
  if (!all(c("trait_type", "label") %in% names(lab))) {
    stop("labels.tsv must have columns trait_type and label")
  }
  type_labels[lab$trait_type] <- lab$label
}

load_summary <- function(name) {
  dir_path <- file.path(results_root, name)
  qst_path <- file.path(dir_path, "summary", "qst_results_abc.csv")
  sum_path <- file.path(dir_path, "summary", "adaptive_trait_summary.csv")
  if (!file.exists(qst_path)) {
    warning("Skipping ", name, ": missing ", qst_path)
    return(NULL)
  }
  qst <- fread(qst_path)
  qst[, trait_type := name]
  summ <- if (file.exists(sum_path)) fread(sum_path) else data.table()
  if (nrow(summ)) summ[, trait_type := name]
  list(qst = qst, summ = summ)
}

chunks <- lapply(trait_types, load_summary)
chunks <- chunks[!vapply(chunks, is.null, logical(1))]
if (!length(chunks)) stop("No trait-type summaries found under ", results_root)

all_qst <- rbindlist(lapply(chunks, `[[`, "qst"), use.names = TRUE, fill = TRUE)
all_qst[, adaptive := factor(adaptive, levels = c("no", "yes"))]
all_qst[, trait_type := factor(trait_type, levels = trait_types)]

all_summ <- rbindlist(lapply(chunks, `[[`, "summ"), fill = TRUE)
fwrite(all_summ, file.path(out_dir, "adaptive_summary_by_trait_type.csv"))
all_qst[, trait_label := unname(type_labels[as.character(trait_type)])]

if (nrow(all_summ) && "n_classifiable" %in% names(all_summ)) {
  n_by_type <- all_summ[, .(
    trait_type,
    n = n_traits,
    n_adaptive,
    n_classifiable,
    pct_adaptive
  )]
} else {
  n_by_type <- all_qst[, .(
    n = .N,
    n_adaptive = sum(adaptive == "yes", na.rm = TRUE),
    n_classifiable = sum(adaptive %in% c("yes", "no"), na.rm = TRUE)
  ), by = trait_type]
  n_by_type[, pct_adaptive := fifelse(
    n_classifiable > 0,
    round(100 * n_adaptive / n_classifiable, 3),
    NA_real_
  )]
}
n_by_type[, trait_label := unname(type_labels[as.character(trait_type)])]
fwrite(n_by_type, file.path(out_dir, "trait_counts_by_type.csv"))

theme_pub <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 10, colour = "grey30", hjust = 0),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
      axis.line = element_line(linewidth = 0.5, colour = "black")
    )
}

p_ad <- ggplot(n_by_type, aes(x = trait_label, y = pct_adaptive, fill = trait_type)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f%%\n(n=%d)", pct_adaptive, n_adaptive)),
            vjust = -0.3, size = 3) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Adaptive proportion by trait type",
    subtitle = "Classifiable traits only (excludes NA adaptive call)",
    x = NULL,
    y = "Percent adaptive"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(out_dir, "adaptive_proportion_by_trait_type.pdf"), p_ad, width = 11, height = 4.5)

p_qst_facet <- ggplot(all_qst, aes(x = QST, fill = adaptive, colour = adaptive)) +
  geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1.05) +
  facet_wrap(~trait_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("no" = "#4E79A7", "yes" = "#E15759"), name = NULL) +
  scale_color_manual(values = c("no" = "#2F5597", "yes" = "#B40426"), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "ABC QST distribution by trait type",
    x = expression(italic(Q)[ST] ~ "(ABC estimate)"),
    y = "Density"
  ) +
  theme_pub()

ggsave(
  file.path(out_dir, "qst_distribution_by_trait_type_faceted.pdf"),
  p_qst_facet, width = 9, height = max(6, 3 * ceiling(length(trait_types) / 2))
)

p_qst_overlay <- ggplot(all_qst, aes(x = QST, colour = trait_label, fill = trait_label)) +
  geom_density(alpha = 0.2, linewidth = 0.95, adjust = 1.05) +
  scale_colour_brewer(palette = "Dark2", name = "Trait type") +
  scale_fill_brewer(palette = "Dark2", name = "Trait type") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "ABC QST distribution overlay across trait types",
    x = expression(italic(Q)[ST] ~ "(ABC estimate)"),
    y = "Density"
  ) +
  theme_pub()

ggsave(
  file.path(out_dir, "qst_distribution_overlay_across_trait_types.pdf"),
  p_qst_overlay, width = 8, height = 4.5
)

ad_only <- all_qst[adaptive == "yes"]
if (nrow(ad_only) >= 5) {
  p_ad_qst <- ggplot(ad_only, aes(x = QST, colour = trait_label, fill = trait_label)) +
    geom_density(alpha = 0.25, linewidth = 0.9, adjust = 1.05) +
    scale_colour_brewer(palette = "Dark2", name = "Trait type") +
    scale_fill_brewer(palette = "Dark2", name = "Trait type") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = "Adaptive-trait QST distributions",
      x = expression(italic(Q)[ST]),
      y = "Density"
    ) +
    theme_pub()
  ggsave(
    file.path(out_dir, "qst_distribution_adaptive_only_overlay.pdf"),
    p_ad_qst, width = 8, height = 4.5
  )
}

cat("Cross-trait-type comparison written to", out_dir, "\n")
