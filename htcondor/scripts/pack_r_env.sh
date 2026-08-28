#!/usr/bin/env bash
# Build htcondor/env/r_env.tar.gz from htcondor/env/r_qst.yml using conda-pack.
# Run from the REDQuanTEA repository root:
#   bash htcondor/scripts/pack_r_env.sh
#
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENV_FILE="${REPO_ROOT}/htcondor/env/r_qst.yml"
OUT="${REPO_ROOT}/htcondor/env/r_env.tar.gz"
ENV_NAME="r_qst"

CONDA="${CONDA_EXE:-$(command -v conda || true)}"
if [[ -z "${CONDA}" ]]; then
  echo "conda not found. Set CONDA_EXE or install Miniforge/Anaconda." >&2
  exit 1
fi

echo "Using conda: ${CONDA}"
"${CONDA}" env remove -n "${ENV_NAME}" -y 2>/dev/null || true
"${CONDA}" env create -f "${ENV_FILE}"

PREFIX="$("${CONDA}" info --envs | awk -v n="${ENV_NAME}" '$1==n {print $NF}')"
if [[ -z "${PREFIX}" || ! -d "${PREFIX}" ]]; then
  echo "Could not resolve prefix for env ${ENV_NAME}" >&2
  exit 1
fi

"${CONDA}" install -y -n base -c conda-forge conda-pack
CONDA_PACK="$(dirname "${CONDA}")/conda-pack"
if [[ ! -x "${CONDA_PACK}" ]]; then
  echo "conda-pack not found next to conda; try: conda install -n base -c conda-forge conda-pack" >&2
  exit 1
fi

TMP_OUT="$(mktemp -u /tmp/r_env_pack_XXXXXX.tar.gz)"
"${CONDA_PACK}" -p "${PREFIX}" -o "${TMP_OUT}" --compress-level 6 --ignore-editable-packages
mv -f "${TMP_OUT}" "${OUT}"
echo "Wrote ${OUT}"
ls -lh "${OUT}"
