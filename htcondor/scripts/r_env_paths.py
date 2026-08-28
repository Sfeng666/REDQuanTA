"""Resolve r_env tarball URL for HTCondor transfer_input_files.

Override defaults with environment variables:
  R_ENV_OSDF_URL       OSDF URL of the packed R environment
  R_ENV_STAGING_PKG    Path on CHTC staging (file:// used when R_ENV_TRANSFER=staging)
  R_ENV_TRANSFER       osdf (default), staging, or home
"""

from __future__ import annotations

import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
HTCONDOR_DIR = SCRIPT_DIR.parent
R_ENV_HOME = HTCONDOR_DIR / "env" / "r_env.tar.gz"
R_ENV_OSDF_URL_FILE = HTCONDOR_DIR / "env" / "r_env_osdf.url"

_USER = os.environ.get("USER") or os.environ.get("LOGNAME") or "sfeng77"
DEFAULT_STAGING_PKG = os.environ.get(
    "R_ENV_STAGING_PKG", f"/staging/{_USER}/qst_fst/r_env_conda_v2.pkg"
)
DEFAULT_OSDF_URL = os.environ.get(
    "R_ENV_OSDF_URL", f"osdf:///chtc/staging/{_USER}/qst_fst/r_env_conda_v2.pkg"
)
DEFAULT_STAGING_FILE_URL = f"file://{DEFAULT_STAGING_PKG}"


def read_osdf_url() -> str:
    if R_ENV_OSDF_URL_FILE.is_file():
        line = R_ENV_OSDF_URL_FILE.read_text().strip().splitlines()[0].strip()
        if line:
            return line
    env_url = os.environ.get("R_ENV_OSDF_URL", "").strip()
    if env_url:
        return env_url
    return DEFAULT_OSDF_URL


def resolve_r_env_tarball(transfer: str | None = None) -> str:
    """Return path/URL for r_env in DAG VARS r_env_tarball.

    transfer: 'osdf' (default), 'staging' (file:///staging, CHTC-only),
              'home', or None (uses R_ENV_TRANSFER env).
    """
    mode = (transfer or os.environ.get("R_ENV_TRANSFER", "osdf")).strip().lower()
    if mode == "home":
        return str(R_ENV_HOME.resolve())
    if mode in ("staging", "file", "osdf_chtc"):
        return DEFAULT_STAGING_FILE_URL
    return read_osdf_url()
