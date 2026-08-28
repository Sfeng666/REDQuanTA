#!/usr/bin/env bash
# Detection Module: generate fused HTCondor DAGs and optionally submit them.
#
# Required:
#   TRAIT_VALUES   trait_values.csv
#   RESULTS_DIR    per-trait output directory
#   OUTPUT_DAG     path for the generated .dag (shards use the same stem)
#
# Optional:
#   SAMPLE_STRUCTURE  default: $REPO/data/example/sample_structure.csv
#   FST_AUTOSOMES     default: $REPO/data/example/qst_neutral_autosomes.txt
#   FST_CHRX          default: $REPO/data/example/qst_neutral_chrX.txt
#   NUM_NEUTRAL       default: 1000
#   NUM_SIM           default: 100000
#   BATCH_SIZE        default: 1000
#   SUMMARY_STATS     default: QST,ratioVbetweenVtotal
#   FLOOR_POLICY      default: ridge_floor
#   FLOOR_ALPHA       default: 0.1
#   JOB_PRIORITY      default: 1000
#   DRY_RUN           1 = generate DAGs only
#   R_ENV             osdf|staging|home
#   SANITY_CHECK      1 = extra diagnostic columns in result CSV
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../.." && pwd)"
PREP="$SCRIPT_DIR/prepare_detection_dag.py"

TRAIT_VALUES="${TRAIT_VALUES:?set TRAIT_VALUES to a trait_values.csv}"
RESULTS_DIR="${RESULTS_DIR:?set RESULTS_DIR}"
OUTPUT_DAG="${OUTPUT_DAG:?set OUTPUT_DAG}"
SAMPLE_STRUCTURE="${SAMPLE_STRUCTURE:-$REPO/data/example/sample_structure.csv}"
FST_AUTOSOMES="${FST_AUTOSOMES:-$REPO/data/example/qst_neutral_autosomes.txt}"
FST_CHRX="${FST_CHRX:-$REPO/data/example/qst_neutral_chrX.txt}"
NUM_NEUTRAL="${NUM_NEUTRAL:-1000}"
NUM_SIM="${NUM_SIM:-100000}"
BATCH_SIZE="${BATCH_SIZE:-1000}"
SUMMARY_STATS="${SUMMARY_STATS:-QST,ratioVbetweenVtotal}"
FLOOR_POLICY="${FLOOR_POLICY:-ridge_floor}"
FLOOR_ALPHA="${FLOOR_ALPHA:-0.1}"
JOB_PRIORITY="${JOB_PRIORITY:-1000}"
DRY_RUN="${DRY_RUN:-0}"
R_ENV="${R_ENV:-osdf}"
SANITY_CHECK="${SANITY_CHECK:-0}"
CODE_DIR="${CODE_DIR:-$REPO/workflow/scripts}"
FUSED_SUB="${FUSED_SUB:-$SCRIPT_DIR/abc_fused_qst.sub}"
MAX_TRAITS="${MAX_TRAITS:-}"
MAX_JOBS_PER_DAG="${MAX_JOBS_PER_DAG:-5000}"
PREP_WORKERS="${PREP_WORKERS:-8}"

mkdir -p "$(dirname "$OUTPUT_DAG")" "$RESULTS_DIR"

prep_args=(
  --trait-values "$TRAIT_VALUES"
  --sample-structure "$SAMPLE_STRUCTURE"
  --fst-autosomes "$FST_AUTOSOMES"
  --fst-chrx "$FST_CHRX"
  --results-dir "$RESULTS_DIR"
  --output-dag "$OUTPUT_DAG"
  --num-neutral "$NUM_NEUTRAL"
  --num-sim "$NUM_SIM"
  --batch-size "$BATCH_SIZE"
  --summary-stats "$SUMMARY_STATS"
  --job-mode fused
  --max-jobs-per-dag "$MAX_JOBS_PER_DAG"
  --base-priority "$JOB_PRIORITY"
  --r-env "$R_ENV"
  --fused-sub "$FUSED_SUB"
  --code-dir "$CODE_DIR"
  --environment "FLOOR_POLICY=${FLOOR_POLICY} FLOOR_ALPHA=${FLOOR_ALPHA}"
  --prep-workers "$PREP_WORKERS"
)
if [[ -n "$MAX_TRAITS" ]]; then
  prep_args+=(--max-traits "$MAX_TRAITS")
fi
if [[ "$SANITY_CHECK" == "1" ]]; then
  prep_args+=(--sanity-check)
fi

echo "Generating Detection Module DAGs..."
python3 "$PREP" "${prep_args[@]}"

stem="${OUTPUT_DAG%.dag}"
shards=()
if ls "${stem}_shard"*.dag >/dev/null 2>&1; then
  mapfile -t shards < <(ls "${stem}_shard"*.dag | sort)
elif [[ -f "$OUTPUT_DAG" ]]; then
  shards=("$OUTPUT_DAG")
else
  echo "ERROR: no DAG written at $OUTPUT_DAG" >&2
  exit 1
fi

if [[ "$DRY_RUN" == "1" ]]; then
  echo "DRY_RUN=1; DAGs:"
  printf '  %s\n' "${shards[@]}"
  exit 0
fi

clusters=()
for dag in "${shards[@]}"; do
  echo "Submitting $(basename "$dag") at priority $JOB_PRIORITY"
  out=$(condor_submit_dag -Force -Priority "$JOB_PRIORITY" "$dag")
  echo "$out"
  cluster=$(echo "$out" | sed -n 's/.*cluster \([0-9][0-9]*\).*/\1/p' | head -1)
  clusters+=("${cluster:-}")
done
echo "Submitted clusters: ${clusters[*]}"
