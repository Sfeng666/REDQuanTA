#!/usr/bin/env bash
# HTCondor/local wrapper for perf_eval_all_combos_fast_cell_v2.R (REDQuanTEA workflow/scripts).
set -euo pipefail

MODE="${1:?need mode adaptive|neutral}"
shift
wd="$(pwd)"

if [[ -f "r_env.tar.gz" ]]; then
  mkdir -p r_env
  tar -xzf r_env.tar.gz -C r_env
  R_ENV_PATH="$wd/r_env"
  if [[ -f "$R_ENV_PATH/bin/conda-unpack" ]]; then
    "$R_ENV_PATH/bin/conda-unpack" 2>&1 || true
  fi
elif [[ -d "lib/R" ]]; then
  R_ENV_PATH="$wd"
  if [[ -f "$R_ENV_PATH/bin/conda-unpack" ]]; then
    "$R_ENV_PATH/bin/python" "$R_ENV_PATH/bin/conda-unpack" 2>&1 || true
  fi
else
  R_ENV_PATH=""
fi

if [[ -n "${R_ENV_PATH:-}" ]]; then
  export PATH="$R_ENV_PATH/bin:$PATH"
  export R_HOME="$R_ENV_PATH/lib/R"
  export R_LIBS="$R_ENV_PATH/lib/R/library"
  export R_LIBS_USER="$R_ENV_PATH/lib/R/library"
  export LD_LIBRARY_PATH="$R_ENV_PATH/lib:$R_ENV_PATH/lib/R/lib:${LD_LIBRARY_PATH:-}"
  RSCRIPT_BIN="$R_ENV_PATH/bin/Rscript"
else
  RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
fi

export QST_CODE_DIR="${QST_CODE_DIR:-$wd}"
echo "REDQuanTEA fast all-combo perf-eval (v2)"
echo "Mode: $MODE"
echo "QST_CODE_DIR: $QST_CODE_DIR"
echo "Estimator: ${FAST_ABC_ESTIMATOR:-loclinear}"
echo "Rscript: $RSCRIPT_BIN"

exec "$RSCRIPT_BIN" "$QST_CODE_DIR/perf_eval_all_combos_fast_cell_v2.R" "$MODE" "$@"
