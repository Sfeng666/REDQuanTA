#!/bin/bash
# Run ABC QST estimation (trait, neutral, batch_neutral, full_trait)
# Arguments: mode input ext_sd_or_file output_file num_sim summary_stats
#
# R environment modes:
# 1. osdf_unpack: lib/R present (legacy ?pack=auto; unused for conda env)
# 2. osdf_tarball / tarball: r_env.tar.gz transferred, extracted locally

mode=$1
input=$2
ext_sd_or_file=$3
output_file=$4
num_sim=${5:-100000}
summary_stats=${6:-"QST,ratioVbetweenVtotal"}

wd=$(pwd)
START_TS=$(date +%s.%N)
STAGING_PKG="${R_ENV_STAGING_PKG:-/staging/${USER:-sfeng77}/qst_fst/r_env_conda_v2.pkg}"
EXPECTED_R_ENV_BYTES="${R_ENV_EXPECTED_BYTES:-519625507}"

echo "Running ABC QST estimation"
echo "Mode: $mode"
echo "Input: $input"
echo "Ext SD or file: $ext_sd_or_file"
echo "Num sim: $num_sim"
echo "Working directory: $wd"

R_ENV_MODE=""
R_ENV_PATH=""

extract_r_env() {
    local archive="$1"
    local label="$2"
    local size
    size=$(stat -c '%s' "$archive" 2>/dev/null || wc -c < "$archive")
    if [ "$size" -lt 10000000 ]; then
        echo "ERROR: $label size $size too small (expected > 10MB)"
        return 1
    fi
    mkdir -p r_env
    local attempt
    for attempt in 1 2 3; do
        rm -rf r_env/*
        if tar -xzf "$archive" -C r_env 2>/dev/null && [ -f "r_env/bin/Rscript" ]; then
            return 0
        fi
        echo "WARN: tar extract attempt $attempt failed from $label"
        sleep 2
    done
    return 1
}

# Find local/transferred R environment archive
ENV_ARCHIVE=""
for dir in "." "${_CONDOR_SCRATCH_DIR:-}"; do
    if [ -n "$dir" ] && [ -d "$dir" ]; then
        for candidate in "$dir"/r_env*.pkg "$dir"/r_env*.tar.gz "$dir"/*.pkg; do
            if [ -f "$candidate" ]; then
                ENV_ARCHIVE="$candidate"
                break 2
            fi
        done
    fi
done

# Find staging package dynamically if hardcoded STAGING_PKG is not readable
if [ ! -r "$STAGING_PKG" ]; then
    for candidate in /staging/${USER:-sfeng77}/qst_fst/r_env*.pkg /staging/${USER:-sfeng77}/qst_fst/r_env*.tar.gz; do
        if [ -r "$candidate" ]; then
            STAGING_PKG="$candidate"
            break
        fi
    done
fi

if [ -d "lib/R" ]; then
    R_ENV_MODE="osdf_unpack"
    R_ENV_PATH="$wd"
    echo "R_ENV_MODE=$R_ENV_MODE"
    if [ -f "$R_ENV_PATH/bin/conda-unpack" ]; then
        "$R_ENV_PATH/bin/python" "$R_ENV_PATH/bin/conda-unpack" 2>&1 || true
    fi
elif [ -r "$STAGING_PKG" ]; then
    R_ENV_MODE="staging_direct"
    echo "R_ENV_MODE=$R_ENV_MODE source=$STAGING_PKG"
    if ! extract_r_env "$STAGING_PKG" "staging $(basename "$STAGING_PKG")"; then
        echo "ERROR: staging direct extract failed from $STAGING_PKG"
        exit 1
    fi
    R_ENV_PATH="$wd/r_env"
elif [ -n "$ENV_ARCHIVE" ]; then
    R_ENV_MODE="tarball"
    echo "R_ENV_MODE=$R_ENV_MODE source=$ENV_ARCHIVE"
    if ! extract_r_env "$ENV_ARCHIVE" "$(basename "$ENV_ARCHIVE")"; then
        echo "ERROR: tar extract failed from transferred $ENV_ARCHIVE"
        echo "HINT: truncated HTCondor transfer; use staging_direct (remove r_env from transfer_input_files)"
        exit 1
    fi
    R_ENV_PATH="$wd/r_env"
else
    echo "ERROR: No R environment (expected $STAGING_PKG, lib/R, or local/transferred archive)"
    ls -la
    exit 1
fi

if [ -f "$R_ENV_PATH/bin/conda-unpack" ]; then
    source "$R_ENV_PATH/bin/activate" 2>/dev/null || true
    "$R_ENV_PATH/bin/python" "$R_ENV_PATH/bin/conda-unpack" 2>&1 || true
fi

ENV_READY_TS=$(date +%s.%N)
echo "R env ready in $(echo "$ENV_READY_TS - $START_TS" | bc)s mode=$R_ENV_MODE"

export PATH=$R_ENV_PATH/bin:$PATH
export R_HOME=$R_ENV_PATH/lib/R
export R_LIBS=$R_ENV_PATH/lib/R/library
export R_LIBS_USER=$R_ENV_PATH/lib/R/library
export LD_LIBRARY_PATH=$R_ENV_PATH/lib:$R_ENV_PATH/lib/R/lib:$LD_LIBRARY_PATH

if [ ! -f "$R_ENV_PATH/bin/Rscript" ]; then
    echo "ERROR: Rscript not found"
    exit 1
fi
CODE_DIR=""
for d in "." ".." "${_CONDOR_SCRATCH_DIR:-}"; do
    if [ -n "$d" ] && [ -f "$d/qst_abc_sim.R" ]; then
        CODE_DIR="$(cd "$d" && pwd)"
        break
    fi
done

if [ -z "$CODE_DIR" ]; then
    echo "ERROR: R script qst_abc_sim.R not found in ., .., or ${_CONDOR_SCRATCH_DIR:-}"
    echo "Directory listing of pwd ($(pwd)):"
    ls -la
    echo "Directory listing of ..:"
    ls -la ..
    exit 1
fi

cd "$CODE_DIR"
wd="$CODE_DIR"

echo "Running Rscript..."
set +e
"$R_ENV_PATH/bin/Rscript" --vanilla "$wd/qst_abc_sim.R" \
    "$mode" "$input" "$ext_sd_or_file" "$output_file" "$num_sim" "$summary_stats" 2>&1
exit_code=$?
set -e

END_TS=$(date +%s.%N)
echo "R script finished exit_code=$exit_code wall_sec=$(echo "$END_TS - $START_TS" | bc)"

if [ "$exit_code" -ne 0 ]; then
    exit "$exit_code"
fi

if [ "$mode" = "full_trait" ]; then
    if [ -f "$output_file" ]; then
        echo "Result CSV created: $output_file"
    else
        echo "ERROR: Result CSV not created: $output_file"
        exit 1
    fi
elif [ -f "$output_file" ]; then
    echo "Output file created: $output_file"
else
    echo "ERROR: Output file not created: $output_file"
    exit 1
fi
