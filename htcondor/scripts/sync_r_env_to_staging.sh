#!/usr/bin/env bash
# Copy r_env.tar.gz to CHTC staging and write the OSDF URL for DAG generators.
#
# Environment:
#   R_ENV_STAGING_DIR   staging directory (default: /staging/$USER/qst_fst)
#   R_ENV_OSDF_URL      OSDF URL written to htcondor/env/r_env_osdf.url
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HTCONDOR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC="$HTCONDOR/env/r_env.tar.gz"
USER_NAME="${USER:-$(id -un)}"
DEST_DIR="${R_ENV_STAGING_DIR:-/staging/${USER_NAME}/qst_fst}"
DEST="$DEST_DIR/r_env.tar.gz"
PKG="$DEST_DIR/r_env_conda_v2.pkg"
URL_FILE="$HTCONDOR/env/r_env_osdf.url"
OSDF_URL="${R_ENV_OSDF_URL:-osdf:///chtc/staging/${USER_NAME}/qst_fst/r_env_conda_v2.pkg}"

log() { echo "$(date '+%F %T') $*"; }

[ -f "$SRC" ] || { log "ERROR: missing $SRC (run pack_r_env.sh first)"; exit 1; }
mkdir -p "$DEST_DIR"

src_hash=$(sha256sum "$SRC" | awk '{print $1}')
if [ -f "$DEST" ]; then
  dest_hash=$(sha256sum "$DEST" | awk '{print $1}')
else
  dest_hash=""
fi

if [ "$src_hash" = "$dest_hash" ]; then
  log "Staging unchanged (sha256 $src_hash)"
else
  log "Copying $SRC -> $DEST"
  cp -f "$SRC" "$DEST"
  log "Synced $(ls -lh "$DEST" | awk '{print $5}')"
fi
cp -f "$SRC" "$PKG"
log "Wrote $PKG"

printf '%s\n' "$OSDF_URL" > "$URL_FILE"
log "Wrote $URL_FILE"
cat "$URL_FILE"
