#!/bin/bash
# Aggregate QST results for all traits after DAG completion.
#
# Usage: ./aggregate_all_results.sh <traits_dir> [threshold_percentile] [sanity_check]
#
#   traits_dir           - Directory containing per-trait subdirectories (trait_*/),
#                          e.g. .../results/traits_100
#   threshold_percentile - Neutral QST percentile used as detection threshold (default: 0.95)
#   sanity_check         - Set to TRUE to include extra diagnostic columns (default: FALSE)
#
# Output: qst_results_abc.csv is written directly into <traits_dir>.

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <traits_dir> [threshold_percentile] [sanity_check]" >&2
    exit 1
fi

TRAITS_DIR="$(cd "$1" && pwd)"
THRESHOLD="${2:-0.95}"
SANITY_CHECK="${3:-FALSE}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$(cd "$SCRIPT_DIR/../.." && pwd)/aggregate_qst.R"

OUTPUT_FILE="$TRAITS_DIR/qst_results_abc.csv"

# Write header
if [ "$SANITY_CHECK" = "TRUE" ]; then
    echo "trait_id,chr,among_pop_sd,within_pop_sd,ext_sd,prior_QST,QST,threshold_percentile,threshold_value,n_neutral,adaptive" > "$OUTPUT_FILE"
else
    echo "trait_id,chr,QST,threshold_percentile,threshold_value,adaptive" > "$OUTPUT_FILE"
fi

n_ok=0
n_warn=0

for trait_dir in "$TRAITS_DIR"/trait_*/; do
    [ -d "$trait_dir" ] || continue
    trait_id="$(basename "$trait_dir" | sed 's/trait_//')"
    result_file="$trait_dir/${trait_id}_result.csv"

    # Fast path: per-trait CSV already produced by the DAG POST script
    if [ -f "$result_file" ]; then
        tail -n +2 "$result_file" >> "$OUTPUT_FILE"
        echo "Appended: $trait_id"
        (( n_ok++ )) || true
        continue
    fi

    # Slow path: re-run Rscript aggregation
    trait_qst_file="$trait_dir/${trait_id}_trait_qst.RData"
    if [ ! -f "$trait_qst_file" ]; then
        echo "WARNING: Missing trait QST file for $trait_id" >&2
        (( n_warn++ )) || true
        continue
    fi

    n_neutral=$(ls -1 "$trait_dir"/neutral_*.RData 2>/dev/null | wc -l)
    if [ "$n_neutral" -eq 0 ]; then
        echo "WARNING: No neutral QST files for $trait_id" >&2
        (( n_warn++ )) || true
        continue
    fi

    echo "Processing $trait_id ($n_neutral neutral files)..."
    Rscript "$R_SCRIPT" "$trait_qst_file" "$trait_dir" "$THRESHOLD" "$result_file" "$SANITY_CHECK"

    if [ -f "$result_file" ]; then
        tail -n +2 "$result_file" >> "$OUTPUT_FILE"
        (( n_ok++ )) || true
    else
        echo "ERROR: Rscript failed for $trait_id" >&2
        (( n_warn++ )) || true
    fi
done

echo ""
echo "Done. $n_ok traits aggregated, $n_warn warnings."
echo "Results written to: $OUTPUT_FILE"
wc -l "$OUTPUT_FILE"

