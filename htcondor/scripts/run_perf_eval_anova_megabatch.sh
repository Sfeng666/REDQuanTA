#!/usr/bin/env bash
# HTCondor/local wrapper for perf_eval_anova_megabatch.R (with-rep or norep via ANOVA_NOREP).
set -euo pipefail

MODE="${1:?need mode neutral_chr|adaptive_chr|adaptive_q}"
shift
wd="$(pwd)"

ENV_ARCHIVE=""
for dir in "." "${_CONDOR_SCRATCH_DIR:-}"; do
  if [[ -n "$dir" && -d "$dir" ]]; then
    for candidate in "$dir"/r_env*.pkg "$dir"/r_env*.tar.gz "$dir"/*.pkg; do
      if [[ -f "$candidate" ]]; then
        ENV_ARCHIVE="$candidate"
        break 2
      fi
    done
  fi
done

if [[ -n "$ENV_ARCHIVE" ]]; then
  mkdir -p r_env
  tar -xzf "$ENV_ARCHIVE" -C r_env
  R_ENV_PATH="$wd/r_env"
  if [[ -f "$R_ENV_PATH/bin/conda-unpack" ]]; then
    "$R_ENV_PATH/bin/python" "$R_ENV_PATH/bin/conda-unpack" 2>&1 || true
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

CODE_DIR=""
for d in "." ".." "${_CONDOR_SCRATCH_DIR:-}"; do
  if [[ -n "$d" && -f "$d/qst_abc_sim.R" ]]; then
    CODE_DIR="$(cd "$d" && pwd)"
    break
  fi
done

if [[ -n "$CODE_DIR" ]]; then
  cd "$CODE_DIR"
  export QST_CODE_DIR="$CODE_DIR"
else
  export QST_CODE_DIR="${QST_CODE_DIR:-$wd}"
fi

WORKER="${ANOVA_WORKER_SCRIPT:-perf_eval_anova_megabatch.R}"
if [[ "${ANOVA_NOREP:-0}" == "1" ]]; then
  WORKER="perf_eval_anova_norep_megabatch.R"
fi

echo "REDQuanTEA ANOVA mega-batch perf-eval"
echo "Mode: $MODE"
echo "Worker: $WORKER"
echo "QST_CODE_DIR: $QST_CODE_DIR"
echo "ANOVA_NOREP: ${ANOVA_NOREP:-0}"
echo "Rscript: $RSCRIPT_BIN"

exec "$RSCRIPT_BIN" "$QST_CODE_DIR/$WORKER" "$MODE" "$@"
