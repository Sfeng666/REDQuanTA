#!/usr/bin/env bash
# Optional local post-processing: fast aggregation for multi-combo perf-eval results.
#
# Use when jobs wrote many summary-stat combinations into each neutral_ratio_* /
# adaptive_q*_* directory (e.g. fast all-combo HTCondor runs). Loads each condition
# directory once and writes per-combo tpr_fpr_matrix_*_combo_XXXX.csv files, then
# combined_model_ranking.csv.
#
# Usage:
#   bash run_aggregate_perf_eval_multicombo_fast.sh <perf_eval_root> [combo_list_file]
#
# <perf_eval_root> must contain:
#   combinations*.txt (or pass combo file as 2nd argument)
#   autosomes/neutral_ratio_*/  and  autosomes/adaptive_q*/
#   chrX/...   (same layout)
#
# Environment overrides (optional):
#   ADAPTIVE_QST   default: 0.5,0.55,...,1.0
#   VE_RATIOS      default: 0.01,0.1,1.0,10.0,100.0
#   THRESHOLD      default: 0.95
#   RSCRIPT        path to Rscript (default: Rscript on PATH)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT_BIN="${RSCRIPT:-$(command -v Rscript)}"
FAST_R="$SCRIPT_DIR/aggregate_perf_eval_multicombo_fast.R"
RANK_R="$SCRIPT_DIR/generate_combined_ranking.R"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <perf_eval_root> [combo_list_file]" >&2
  exit 1
fi

OUT_ROOT="$(cd "$1" && pwd)"
COMBO_FILE="${2:-}"

if [[ -z "$COMBO_FILE" ]]; then
  for cand in "$OUT_ROOT"/combinations_all_nonempty_11stats.txt \
              "$OUT_ROOT"/combinations*.txt; do
    if [[ -f "$cand" ]]; then
      COMBO_FILE="$cand"
      break
    fi
  done
fi

if [[ -z "$COMBO_FILE" ]] || [[ ! -f "$COMBO_FILE" ]]; then
  echo "ERROR: combo list file not found under $OUT_ROOT (pass as 2nd argument)" >&2
  exit 1
fi
COMBO_FILE="$(cd "$(dirname "$COMBO_FILE")" && pwd)/$(basename "$COMBO_FILE")"

ADAPTIVE_QST="${ADAPTIVE_QST:-0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0}"
VE_RATIOS="${VE_RATIOS:-0.01,0.1,1.0,10.0,100.0}"
THRESH="${THRESHOLD:-0.95}"

for chr in autosomes chrX; do
  if [[ ! -d "$OUT_ROOT/$chr" ]]; then
    echo "WARNING: skipping missing chromosome dir: $OUT_ROOT/$chr" >&2
    continue
  fi
  echo "[$chr] fast aggregate..."
  "$RSCRIPT_BIN" "$FAST_R" \
    "$OUT_ROOT/$chr" "$chr" "$ADAPTIVE_QST" "$VE_RATIOS" "$THRESH" \
    "$OUT_ROOT/$chr/tpr_fpr_matrix_${chr}.csv" \
    "$COMBO_FILE"
done

echo "Building combined_model_ranking.csv..."
"$RSCRIPT_BIN" "$RANK_R" "$OUT_ROOT"

if [[ "${RUN_PUBLICATION_RANKING:-0}" == "1" ]]; then
  echo "Publication ranking (RUN_PUBLICATION_RANKING=1)..."
  bash "$SCRIPT_DIR/run_perf_eval_publication_ranking.sh" "$OUT_ROOT"
else
  echo "Optional next step: publication table"
  echo "  bash workflow/scripts/run_perf_eval_publication_ranking.sh $OUT_ROOT"
fi

echo "Done. Outputs under $OUT_ROOT (per-combo CSVs in autosomes/ and chrX/)."
