#!/usr/bin/env bash
# REDQuanTA HTCondor validation bundle (three batches, descending Condor priority).
#
# Prerequisites (clone REDQuanTA, conda env optional for local DAG generation):
#   - Run from an HTCondor submit host (e.g. CHTC login node)
#   - htcondor/env/r_env.tar.gz present (conda-pack; see README_details.md)
#
# Priority convention (condor_submit_dag -Priority): higher runs sooner relative to lower.
#   Batch 1 (sample-structure perf-eval replica): Priority 3000
#   Batch 2 (trait QST detection, sanity-check aggregation): Priority 2000
#   Batch 3 (2047-combo fast perf-eval): Priority 1000
#
# Trait input: populate REDQuanTA/data/example/trait_values.csv (see harmonizr_to_trait_values.py
# + workflow/scripts/add_chromosome_info.py from repo root paths in docs below).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REDQUANTA="$(cd "$SCRIPT_DIR/../.." && pwd)"
PY="${PYTHON:-python3}"

cd "$REDQUANTA"

P1="${PRIORITY_BATCH1:-3000}"
P2="${PRIORITY_BATCH2:-2000}"
P3="${PRIORITY_BATCH3:-1000}"

SUMMARY_STATS="QST,ratioVbetweenVtotal"
NUM_NEUTRAL="${NUM_NEUTRAL:-10000}"
NUM_SIM="${NUM_SIM:-100000}"
NUM_REPEATS_SAMPLE_STRUCT="${NUM_REPEATS_SAMPLE_STRUCT:-10000}"
BATCH_SIZE="${BATCH_SIZE:-1000}"

TRAIT_MAX="${TRAIT_MAX:-1040}"
COLLECTION_NAME="${COLLECTION_NAME:-traits_harmonizr_validation}"
TRAIT_VALUES="${TRAIT_VALUES_CSV:-$REDQUANTA/data/example/trait_values.csv}"

STRUCTURES=(
  "4 3"
  "6 3"
  "8 3"
  "10 3"
  "12 3"
  "14 3"
  "6 2"
  "9 2"
  "12 2"
  "15 2"
  "18 2"
  "21 2"
)

SAMPLE_STRUCT_ROOT="$REDQUANTA/htcondor/results/validation_sample_struct_QST_ratioVbetweenVtotal"
PERF_FAST_ROOT="$REDQUANTA/htcondor/results/validation_perf_fast_2047"

usage() {
  echo "Usage: $0 [--dag-only] [--skip-batch1] [--skip-batch2] [--skip-batch3]"
  echo "  --dag-only    Generate DAGs only; do not condor_submit_dag"
  echo "  Env: TRAIT_MAX (default 1040), TRAIT_VALUES_CSV, PRIORITY_BATCH1/2/3"
}

DAG_ONLY=false
SKIP1=false
SKIP2=false
SKIP3=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dag-only) DAG_ONLY=true ;;
    --skip-batch1) SKIP1=true ;;
    --skip-batch2) SKIP2=true ;;
    --skip-batch3) SKIP3=true ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
  shift
done

submit_dag() {
  local dag="$1"
  local pri="$2"
  if [[ ! -f "$dag" ]]; then
    echo "ERROR missing DAG: $dag"
    return 1
  fi
  if [[ "$DAG_ONLY" == true ]]; then
    echo "[dry-run] condor_submit_dag -Force -Priority $pri $dag"
    return 0
  fi
  echo "Submitting $dag at Priority $pri"
  condor_submit_dag -Force "-Priority" "$pri" "$dag"
}

echo "REDQuanTA root: $REDQUANTA"

# ----- Batch 1: sample-structure comparison (matches code/run_sample_structure_comparison_QST_ratioVbetweenVtotal.sh) -----
if [[ "$SKIP1" != true ]]; then
  mkdir -p "$SAMPLE_STRUCT_ROOT"
  for struct in "${STRUCTURES[@]}"; do
    read -r n_ind n_rep <<< "$struct"
    struct_id="n2_i${n_ind}_r${n_rep}"
    out="$SAMPLE_STRUCT_ROOT/$struct_id"
    mkdir -p "$out"
    gen_mode="both"
    need_auto=true
    need_chrx=true
    res_auto="$out/autosomes/tpr_fpr_matrix_autosomes.csv"
    res_x="$out/chrX/tpr_fpr_matrix_chrX.csv"
    if [[ -f "$res_auto" ]]; then need_auto=false; fi
    if [[ -f "$res_x" ]]; then need_chrx=false; fi
    if [[ "$need_auto" == false ]] && [[ "$need_chrx" == false ]]; then
      echo "[$struct_id] already complete"
      continue
    fi
    if [[ "$need_auto" == true ]] && [[ "$need_chrx" == true ]]; then gen_mode="both"
    elif [[ "$need_auto" == true ]]; then gen_mode="autosomes"
    else gen_mode="chrX"; fi

    "$PY" "$SCRIPT_DIR/prepare_perf_eval_dag.py" \
      --chr "$gen_mode" \
      --output-dir "$out" \
      --num-pop 2 \
      --num-ind "$n_ind" \
      --num-rep "$n_rep" \
      --adaptive-qst "0.8" \
      --ve-ratios "0.1,1.0,10.0" \
      --num-neutral "$NUM_NEUTRAL" \
      --num-sim "$NUM_SIM" \
      --num-repeats "$NUM_REPEATS_SAMPLE_STRUCT" \
      --batch-size "$BATCH_SIZE" \
      --summary-stats "$SUMMARY_STATS"

    if [[ "$need_auto" == true ]] && [[ -f "$out/perf_eval_autosomes.dag" ]]; then
      submit_dag "$out/perf_eval_autosomes.dag" "$P1"
    fi
    if [[ "$need_chrx" == true ]] && [[ -f "$out/perf_eval_chrX.dag" ]]; then
      submit_dag "$out/perf_eval_chrX.dag" "$P1"
    fi
  done
fi

# ----- Batch 2: trait DAGs + master SUBDAG (sanity-check mode) -----
if [[ "$SKIP2" != true ]]; then
  if [[ ! -f "$TRAIT_VALUES" ]]; then
    echo "ERROR: trait_values not found: $TRAIT_VALUES"
    exit 1
  fi
  "$PY" "$SCRIPT_DIR/generate_all_dags.py" \
    --max-traits "$TRAIT_MAX" \
    --sanity-check \
    --collection-name "$COLLECTION_NAME" \
    --trait-values "$TRAIT_VALUES" \
    --num-neutral "$NUM_NEUTRAL" \
    --batch-size "$BATCH_SIZE" \
    --num-sim "$NUM_SIM"

  MASTER="$REDQUANTA/htcondor/dags/$COLLECTION_NAME/all_traits.dag"
  submit_dag "$MASTER" "$P2"
fi

# ----- Batch 3: fast all-combo 2047 (REDQuanTA worker scripts) -----
if [[ "$SKIP3" != true ]]; then
  COMBO="$PERF_FAST_ROOT/combinations_all_nonempty_11stats.txt"
  if [[ ! -f "$COMBO" ]]; then
    echo "ERROR: copy combos file first, expected: $COMBO"
    echo "  cp test_abc_method/results/perf_eval_fast_allcombo_2047_v2_loclinear_4gb_softclip_fixed/combinations_all_nonempty_11stats.txt \\"
    echo "    $COMBO"
    exit 1
  fi
  "$PY" "$SCRIPT_DIR/prepare_perf_eval_all_combos_fast_dag.py" \
    --summary-stats "$COMBO" \
    --output-dir "$PERF_FAST_ROOT" \
    --chr both \
    --num-repeats 10000 \
    --num-neutral "$NUM_NEUTRAL" \
    --num-sim "$NUM_SIM" \
    --batch-size "$BATCH_SIZE" \
    --request-memory 4GB

  DAG_AUTO="$PERF_FAST_ROOT/perf_eval_autosomes.dag"
  DAG_X="$PERF_FAST_ROOT/perf_eval_chrX.dag"
  if [[ ! -f "$DAG_AUTO" ]]; then
    echo "ERROR: missing autosomes DAG: $DAG_AUTO"
    exit 1
  fi
  if [[ ! -f "$DAG_X" ]]; then
    echo "ERROR: missing chrX DAG: $DAG_X"
    exit 1
  fi
  submit_dag "$DAG_AUTO" "$P3"
  submit_dag "$DAG_X" "$P3"
fi

if [[ "$DAG_ONLY" != true ]]; then
  sleep 3
  echo "--- condor_q (batch user) ---"
  condor_q -nobatch || true
fi

echo "Done."
