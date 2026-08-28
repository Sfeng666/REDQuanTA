# Memory-safe helpers for pooled neutral QST/FST plots (sourced by other scripts).

suppressPackageStartupMessages(library(data.table))

NEUTRAL_MAX_RESERVOIR <- as.integer(Sys.getenv("NEUTRAL_PLOT_MAX_POINTS", "500000"))
NEUTRAL_N_BINS <- as.integer(Sys.getenv("NEUTRAL_PLOT_N_BINS", "512"))

reservoir_update <- function(reservoir, new_x, k = NEUTRAL_MAX_RESERVOIR) {
  new_x <- new_x[is.finite(new_x)]
  if (!length(new_x)) {
    return(reservoir)
  }
  combined <- c(reservoir, new_x)
  if (length(combined) <= k) {
    return(combined)
  }
  sample(combined, k)
}

NEUTRAL_FST_LO <- as.numeric(Sys.getenv("NEUTRAL_PLOT_FST_LO", "-0.1"))

init_neutral_accumulator <- function(n_bins = NEUTRAL_N_BINS, fst_lo = NEUTRAL_FST_LO) {
  edges <- seq(0, 1, length.out = n_bins + 1)
  fst_edges <- seq(fst_lo, 1, length.out = n_bins + 1)
  list(
    n_bins = n_bins,
    edges = edges,
    fst_edges = fst_edges,
    qst_counts = numeric(n_bins),
    fst_counts = numeric(n_bins),
    n_pairs = 0L,
    n_fst_lt0 = 0L,
    qst_reservoir = numeric(),
    fst_reservoir = numeric()
  )
}

accumulate_neutral_batch <- function(acc, fst, qst) {
  valid <- is.finite(qst)
  if (!any(valid)) {
    return(acc)
  }
  fst <- fst[valid]
  qst <- qst[valid]
  n <- length(qst)
  acc$n_pairs <- acc$n_pairs + n
  acc$n_fst_lt0 <- acc$n_fst_lt0 + sum(fst < 0)
  acc$qst_counts <- acc$qst_counts + as.numeric(
    hist(qst, breaks = acc$edges, plot = FALSE)$counts
  )
  acc$fst_counts <- acc$fst_counts + as.numeric(
    hist(fst, breaks = acc$fst_edges, plot = FALSE)$counts
  )
  acc$qst_reservoir <- reservoir_update(acc$qst_reservoir, qst)
  acc$fst_reservoir <- reservoir_update(acc$fst_reservoir, fst)
  acc
}

weighted_quantile_from_hist <- function(counts, edges, p = 0.5) {
  if (!sum(counts)) {
    return(NA_real_)
  }
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  cw <- cumsum(counts) / sum(counts)
  mids[which(cw >= p)[1]]
}

hist_to_density_curve <- function(counts, edges, n_grid = 512, adjust = 1.05, clip_eval = FALSE) {
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  if (!sum(counts)) {
    return(data.table(x = numeric(), y = numeric()))
  }
  w <- counts / sum(counts)
  d <- if (isTRUE(clip_eval)) {
    density(mids, weights = w, from = 0, to = 1, n = n_grid, adjust = adjust)
  } else {
    density(mids, weights = w, n = n_grid, adjust = adjust)
  }
  data.table(x = d$x, y = d$y)
}

# Normalized histogram density at bin midpoints (peak x matches binned data exactly).
hist_to_density_histogram <- function(counts, edges) {
  n <- sum(counts)
  if (!n) {
    return(data.table(x = numeric(), y = numeric()))
  }
  bw <- diff(edges[1:2])
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  data.table(x = mids, y = counts / (n * bw))
}

hist_mode_x <- function(counts, edges) {
  if (!sum(counts)) {
    return(NA_real_)
  }
  mids <- (edges[-1] + edges[-length(edges)]) / 2
  mids[which.max(counts)]
}

kde_peak_x <- function(curve) {
  if (!nrow(curve)) {
    return(NA_real_)
  }
  curve$x[which.max(curve$y)]
}

curve_peak_x <- kde_peak_x

normalize_density_auc <- function(curve, x_col = "x", y_col = "y") {
  out <- copy(curve)
  area <- trapz_area_xy(out[[x_col]], out[[y_col]])
  if (is.finite(area) && area > 0) {
    out[[y_col]] <- out[[y_col]] / area
  }
  out
}

# Full unclipped KDE, then counter-shift x so the mode aligns with the histogram peak.
hist_to_density_smooth <- function(
    counts,
    edges,
    n_grid = 512,
    clip_eval = FALSE,
    adjust = as.numeric(Sys.getenv("NEUTRAL_PLOT_KDE_ADJUST", "1.05"))) {
  if (!sum(counts)) {
    return(list(
      curve = data.table(x = numeric(), y = numeric()),
      mode_x = NA_real_,
      peak_x = NA_real_,
      kde_peak_x = NA_real_,
      peak_shift = NA_real_,
      kde_adjust = NA_real_
    ))
  }
  mode_x <- hist_mode_x(counts, edges)
  curve <- hist_to_density_curve(counts, edges, n_grid, adjust, clip_eval)
  kde_peak_x <- curve_peak_x(curve)
  peak_shift <- mode_x - kde_peak_x
  curve[, x := x + peak_shift]

  list(
    curve = curve,
    mode_x = mode_x,
    peak_x = curve_peak_x(curve),
    kde_peak_x = kde_peak_x,
    peak_shift = peak_shift,
    kde_adjust = adjust
  )
}

trapz_area_xy <- function(x, y) {
  if (length(x) < 2L) {
    return(NA_real_)
  }
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

hist_counts_in_edges <- function(x, edges) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(integer(length(edges) - 1L))
  }
  as.integer(hist(x, breaks = edges, plot = FALSE)$counts)
}

sample_to_density_curve <- function(x, n_grid = 512, adjust = 1.05) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(data.table(x = numeric(), y = numeric()))
  }
  d <- density(x, from = 0, to = 1, n = n_grid, adjust = adjust)
  data.table(x = d$x, y = d$y)
}

read_neutral_batch_vectors <- function(f) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  r <- env$result$results
  if (is.null(r)) {
    return(NULL)
  }
  valid <- !is.na(r$qst)
  if (!any(valid)) {
    return(NULL)
  }
  list(fst = r$fst[valid], qst = r$qst[valid])
}

stream_neutral_batches <- function(batch_files, progress_every = 500L) {
  condor_cores <- suppressWarnings(as.integer(Sys.getenv("_CONDOR_NPROCS", "4")))
  num_cores <- suppressWarnings(as.integer(Sys.getenv("QST_ABC_CORES", as.character(condor_cores))))
  if (is.na(num_cores) || num_cores < 1L) num_cores <- 1L
  num_cores <- min(num_cores, parallel::detectCores(), 8L)

  if (num_cores <= 1L || length(batch_files) < 20L) {
    acc <- init_neutral_accumulator()
    for (i in seq_along(batch_files)) {
      vecs <- tryCatch(
        read_neutral_batch_vectors(batch_files[i]),
        error = function(e) NULL
      )
      if (!is.null(vecs)) {
        acc <- accumulate_neutral_batch(acc, vecs$fst, vecs$qst)
      }
      if (i %% progress_every == 0) {
        cat("  streamed", i, "/", length(batch_files), "batches (",
            format(acc$n_pairs, big.mark = ","), "pairs in memory-safe accum)\n", sep = "")
      }
    }
    return(acc)
  }

  cat("Streaming ", length(batch_files), " neutral batch files in parallel using ", num_cores, " cores...\n", sep = "")

  # Split batch_files into list of chunks for workers
  chunks <- split(batch_files, cut(seq_along(batch_files), num_cores, labels = FALSE))

  worker_results <- parallel::mclapply(seq_along(chunks), function(chunk_idx) {
    files <- chunks[[chunk_idx]]
    acc <- init_neutral_accumulator()
    
    # Pre-allocate to avoid repeated memory growth and copying
    all_fst <- numeric(length(files) * 1000L)
    all_qst <- numeric(length(files) * 1000L)
    fill_idx <- 0L

    for (f in files) {
      vecs <- tryCatch(
        read_neutral_batch_vectors(f),
        error = function(e) NULL
      )
      if (!is.null(vecs)) {
        n <- length(vecs$qst)
        if (n > 0L) {
          # Accumulate counts
          acc$n_pairs <- acc$n_pairs + n
          acc$n_fst_lt0 <- acc$n_fst_lt0 + sum(vecs$fst < 0)
          
          acc$qst_counts <- acc$qst_counts + as.numeric(
            hist(vecs$qst, breaks = acc$edges, plot = FALSE)$counts
          )
          acc$fst_counts <- acc$fst_counts + as.numeric(
            hist(vecs$fst, breaks = acc$fst_edges, plot = FALSE)$counts
          )
          
          # Store vectors
          start_pos <- fill_idx + 1L
          end_pos <- fill_idx + n
          if (end_pos > length(all_qst)) {
            length(all_fst) <- length(all_fst) * 2L
            length(all_qst) <- length(all_qst) * 2L
          }
          all_fst[start_pos:end_pos] <- vecs$fst
          all_qst[start_pos:end_pos] <- vecs$qst
          fill_idx <- end_pos
        }
      }
    }

    # Trim vectors to actual filled size
    if (fill_idx > 0L) {
      all_fst <- all_fst[1:fill_idx]
      all_qst <- all_qst[1:fill_idx]
      
      # Sample down locally to NEUTRAL_MAX_RESERVOIR to keep transfer light
      if (fill_idx > NEUTRAL_MAX_RESERVOIR) {
        keep_idx <- sample.int(fill_idx, NEUTRAL_MAX_RESERVOIR)
        acc$fst_reservoir <- all_fst[keep_idx]
        acc$qst_reservoir <- all_qst[keep_idx]
      } else {
        acc$fst_reservoir <- all_fst
        acc$qst_reservoir <- all_qst
      }
    }
    
    acc
  }, mc.cores = num_cores)

  # Merge accumulator results
  final_acc <- init_neutral_accumulator()
  
  pooled_fst <- numeric()
  pooled_qst <- numeric()

  for (acc in worker_results) {
    if (is.null(acc) || inherits(acc, "try-error")) next
    final_acc$n_pairs <- final_acc$n_pairs + acc$n_pairs
    final_acc$n_fst_lt0 <- final_acc$n_fst_lt0 + acc$n_fst_lt0
    final_acc$qst_counts <- final_acc$qst_counts + acc$qst_counts
    final_acc$fst_counts <- final_acc$fst_counts + acc$fst_counts
    
    pooled_fst <- c(pooled_fst, acc$fst_reservoir)
    pooled_qst <- c(pooled_qst, acc$qst_reservoir)
  }
  
  # Final sample for reservoir
  if (length(pooled_qst) > NEUTRAL_MAX_RESERVOIR) {
    keep_idx <- sample.int(length(pooled_qst), NEUTRAL_MAX_RESERVOIR)
    final_acc$qst_reservoir <- pooled_qst[keep_idx]
    final_acc$fst_reservoir <- pooled_fst[keep_idx]
  } else {
    final_acc$qst_reservoir <- pooled_qst
    final_acc$fst_reservoir <- pooled_fst
  }

  final_acc
}

save_neutral_accumulator <- function(acc, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- data.table(
    n_pairs = acc$n_pairs,
    n_bins = acc$n_bins,
    reservoir_qst = length(acc$qst_reservoir),
    reservoir_fst = length(acc$fst_reservoir)
  )
  saveRDS(
    list(
      edges = acc$edges,
      fst_edges = acc$fst_edges,
      qst_counts = acc$qst_counts,
      fst_counts = acc$fst_counts,
      n_pairs = acc$n_pairs,
      n_fst_lt0 = acc$n_fst_lt0,
      qst_reservoir = acc$qst_reservoir,
      fst_reservoir = acc$fst_reservoir
    ),
    file.path(out_dir, "combined_neutral_hist.rds")
  )
  fwrite(meta, file.path(out_dir, "combined_neutral_hist_meta.csv"))
  invisible(acc)
}
