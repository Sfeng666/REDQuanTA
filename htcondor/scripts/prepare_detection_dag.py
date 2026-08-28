#!/usr/bin/env python3
"""Generate HTCondor DAGs for the Detection Module (all traits in a collection).

Usage:
    python htcondor/scripts/prepare_detection_dag.py \\
        --trait-values data/example/trait_values.csv \\
        --sample-structure data/example/sample_structure.csv \\
        --results-dir results/detect_chtc \\
        --output-dag results/dags/detect.dag \\
        --job-mode fused --environment 'FLOOR_POLICY=F3 FLOOR_ALPHA=0.1'
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from prepare_trait_dag import (
    DEFAULT_BATCH_SIZE,
    HTCONDOR_DIR,
    INPUT_DIR,
    RSCRIPT,
    SCRIPT_DIR,
    WORKFLOW_SCRIPT_DIR,
    create_fst_batch_files,
    read_fst_values,
    read_trait_ids,
    resolve_r_env_tarball_for_dag,
    run_prep_locally,
)
from r_env_paths import resolve_r_env_tarball

DEFAULT_COLLECTION = "traits_harmonizr_validation"
DAG_DIR = HTCONDOR_DIR / "dags" / DEFAULT_COLLECTION


def safe_job_prefix(trait_id: str) -> str:
    return trait_id.replace("-", "_").replace(".", "_").replace(":", "_").replace("+", "_")


def copy_obs_stats_from_baseline(
    trait_id: str,
    trait_outdir: Path,
    baseline_dir: Path,
) -> Path | None:
    """Copy precomputed obs_stats from another repo tree (read-only)."""
    obs_name = f"{trait_id}_obs_stats.RData"
    src = baseline_dir / f"trait_{trait_id}" / obs_name
    if not src.is_file():
        return None
    dest = trait_outdir / obs_name
    shutil.copy2(src, dest)
    return dest


def existing_obs_stats(trait_id: str, trait_outdir: Path) -> Path | None:
    """Return path when obs_stats already exist in trait output dir."""
    dest = trait_outdir / f"{trait_id}_obs_stats.RData"
    return dest if dest.is_file() else None


def ratioVext_from_obs_stats_file(obs_path: Path) -> float | None:
    """Read ratioVext from an obs_stats RData file."""
    rscript = os.environ.get("RSCRIPT", "Rscript")
    cmd = [
        rscript,
        "-e",
        f'load("{obs_path}"); v <- as.numeric(obs_stats["ratioVext"]); '
        f'if (length(v) != 1L || is.na(v)) quit(status=1); cat(v)',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    try:
        return float(result.stdout.strip())
    except ValueError:
        return None


def batch_ratioVext_from_obs(results_dir: Path, trait_ids: list[str]) -> dict[str, float]:
    """Read ratioVext for many traits in a single R invocation."""
    rscript = os.environ.get("RSCRIPT", "Rscript")
    out: dict[str, float] = {}

    # 1. Read consolidated TSV if it exists
    summary_tsv = results_dir / "batch_ratioVext.tsv"
    if summary_tsv.is_file():
        try:
            with open(summary_tsv, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        out[parts[0]] = float(parts[1])
        except Exception as e:
            print(f"Warning: failed to read {summary_tsv}: {e}", file=sys.stderr)

    # 2. Find missing ones that have existing obs_stats
    missing_rows = []
    for trait_id in trait_ids:
        if trait_id in out:
            continue
        obs = results_dir / f"trait_{trait_id}" / f"{trait_id}_obs_stats.RData"
        if obs.is_file():
            missing_rows.append(f"{trait_id}\t{obs.resolve()}")

    if missing_rows:
        print(f"Loading ratioVext for {len(missing_rows)} pre-existing/copied traits...", file=sys.stderr)
        out.update(_ratioVext_chunk(rscript, missing_rows))

    return out


def _ratioVext_chunk(rscript: str, rows: list[str]) -> dict[str, float]:
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("\n".join(rows))
        f.write("\n")
        list_path = f.name

    r_code = (
        f'lines <- readLines("{list_path}"); '
        'for (line in lines) { '
        'if (!nzchar(line)) next; '
        'parts <- strsplit(line, "\\t")[[1]]; '
        'if (length(parts) != 2) next; '
        'tryCatch({ '
        '  load(parts[2]); '
        '  if (exists("obs_stats") && "ratioVext" %in% names(obs_stats)) { '
        '    cat(parts[1], "\\t", as.numeric(obs_stats["ratioVext"]), "\\n", sep=""); '
        '  } '
        '}, error = function(e) {}) '
        '}'
    )
    result = subprocess.run([rscript, "-e", r_code], capture_output=True, text=True)
    os.unlink(list_path)
    if result.returncode != 0:
        print(result.stderr[:500], file=sys.stderr)
        return {}

    chunk_out: dict[str, float] = {}
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) == 2:
            try:
                chunk_out[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return chunk_out


def ratioVext_from_baseline_dag(trait_id: str, baseline_dag_dir: Path) -> float | None:
    """Read ratioVext from existing per-trait DAG header or VARS (read-only)."""
    dag = baseline_dag_dir / f"trait_{trait_id}.dag"
    if not dag.is_file():
        return None
    text = dag.read_text()
    for line in text.splitlines():
        if "ratioVext:" in line and line.strip().startswith("#"):
            part = line.split("ratioVext:")[1].split(",")[0].strip()
            return float(part)
        if 'ratioVext="' in line:
            start = line.index('ratioVext="') + len('ratioVext="')
            end = line.index('"', start)
            return float(line[start:end])
    return None


def write_trait_section(
    f,
    trait_info: dict,
    params: dict,
    ratioVext: float,
    obs_stats_file: Path,
    batch_size: int | None,
    job_prefix: str,
    sanity_check: bool,
    r_env: str | None = None,
    batch_neutral_sub: str | None = None,
    trait_sub: str | None = None,
) -> tuple[list[str], list[str], int]:
    """Write trait/batch/aggregate jobs; return (parent_jobs, all_jobs, n_execute)."""
    trait_id = trait_info["trait_id"]
    chr_type = trait_info["chr"]
    fst_file = params["fst_chrx"] if chr_type == "chrX" else params["fst_autosomes"]
    fst_values = read_fst_values(fst_file, params["num_neutral"])
    trait_outdir = Path(params["results_dir"]) / f"trait_{trait_id}"
    os.makedirs(trait_outdir, exist_ok=True)

    trait_qst_file = f"{trait_id}_trait_qst.RData"
    result_file = f"{trait_id}_result.csv"
    use_batch = batch_size is not None and batch_size > 1
    parent_jobs: list[str] = []
    neutral_jobs: list[str] = []
    all_jobs: list[str] = []
    neutral_rdata_names: list[str] = []

    r_env = r_env or resolve_r_env_tarball_for_dag()
    batch_neutral_sub = batch_neutral_sub or f"{SCRIPT_DIR}/abc_batch_neutral.sub"
    trait_sub = trait_sub or f"{SCRIPT_DIR}/abc_qst.sub"

    job_trait = f"trait_{job_prefix}"
    f.write(f"JOB {job_trait} {trait_sub}\n")
    f.write(f'VARS {job_trait} mode="trait" ')
    f.write(f'input="{obs_stats_file.name}" ')
    f.write(f'ext_sd_or_file="ignored" ')
    f.write(f'output_file="{trait_qst_file}" ')
    f.write(f'num_sim="{params["num_sim"]}" ')
    f.write(f'summary_stats="{params["summary_stats"]}" ')
    f.write(f'extra_input_comma=", " ')
    f.write(f'extra_input_files="{obs_stats_file}" ')
    f.write(f'outdir="{trait_outdir}" ')
    f.write(f'dir_scripts="{SCRIPT_DIR}" ')
    f.write(f'dir_code="{WORKFLOW_SCRIPT_DIR}" ')
    f.write(f'r_env_tarball="{r_env}" ')
    f.write(f'job="{job_trait}"\n')
    parent_jobs.append(job_trait)
    all_jobs.append(job_trait)

    # Neutral jobs use ABC VE/VG from trait_qst.RData (child of trait ABC job).
    neutral_ratioVext_arg = trait_qst_file

    if use_batch:
        batch_files = create_fst_batch_files(fst_values, trait_outdir, batch_size)
        for i, (batch_file, _n) in enumerate(batch_files, 1):
            job_neutral = f"batch_{job_prefix}_{i}"
            neutral_file = f"neutral_batch_{i}.RData"
            neutral_rdata_names.append(neutral_file)
            f.write(f"JOB {job_neutral} {batch_neutral_sub}\n")
            f.write(f'VARS {job_neutral} fst_input="{batch_file.name}" ')
            f.write(f'ratioVext="{neutral_ratioVext_arg}" ')
            f.write(f'output_file="{neutral_file}" ')
            f.write(f'num_sim="{params["num_sim"]}" ')
            f.write(f'summary_stats="{params["summary_stats"]}" ')
            f.write(f'fst_file="{batch_file.resolve()}" ')
            f.write(f'outdir="{trait_outdir}" ')
            f.write(f'dir_scripts="{SCRIPT_DIR}" ')
            f.write(f'dir_code="{WORKFLOW_SCRIPT_DIR}" ')
            f.write(f'r_env_tarball="{r_env}" ')
            extra_files = f", {trait_outdir}/{trait_qst_file}"
            f.write(f'extra_input_files="{extra_files}" ')
            f.write(f'job="{job_neutral}"\n')
            f.write(f"PARENT {job_trait} CHILD {job_neutral}\n")
            neutral_jobs.append(job_neutral)
            all_jobs.append(job_neutral)
        num_neutral_jobs = len(batch_files)
    else:
        for i, fst_val in enumerate(fst_values, 1):
            job_neutral = f"neutral_{job_prefix}_{i}"
            neutral_file = f"neutral_{i}.RData"
            neutral_rdata_names.append(neutral_file)
            f.write(f"JOB {job_neutral} {SCRIPT_DIR}/abc_neutral.sub\n")
            f.write(f'VARS {job_neutral} fst_value="{fst_val}" ')
            f.write(f'ratioVext="{neutral_ratioVext_arg}" ')
            f.write(f'output_file="{neutral_file}" ')
            f.write(f'num_sim="{params["num_sim"]}" ')
            f.write(f'summary_stats="{params["summary_stats"]}" ')
            f.write(f'outdir="{trait_outdir}" ')
            f.write(f'dir_scripts="{SCRIPT_DIR}" ')
            f.write(f'dir_code="{WORKFLOW_SCRIPT_DIR}" ')
            f.write(f'r_env_tarball="{r_env}" ')
            extra_files = f", {trait_outdir}/{trait_qst_file}"
            f.write(f'extra_input_files="{extra_files}" ')
            f.write(f'job="{job_neutral}"\n')
            f.write(f"PARENT {job_trait} CHILD {job_neutral}\n")
            neutral_jobs.append(job_neutral)
            all_jobs.append(job_neutral)
        num_neutral_jobs = len(fst_values)

    # NOOP + POST aggregation on submit node (obs_stats local; no r_env transfer).
    noop_job = f"noop_{job_prefix}"
    agg_script = trait_outdir / f"aggregate_{trait_id}.sh"
    result_path = trait_outdir / result_file
    rscript = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")
    tmpdir = HTCONDOR_DIR / "tmp"

    with open(agg_script, "w") as agg_f:
        agg_f.write("#!/bin/bash\n")
        agg_f.write(f'export TMPDIR="{tmpdir}"\n')
        agg_f.write('mkdir -p "$TMPDIR"\n')
        agg_f.write(f'"{rscript}" {WORKFLOW_SCRIPT_DIR}/aggregate_qst.R \\\n')
        agg_f.write(f'    "{trait_outdir}/{trait_qst_file}" \\\n')
        agg_f.write(f'    "{trait_outdir}" \\\n')
        agg_f.write(f'    "{params["threshold_percentile"]}" \\\n')
        agg_f.write(f'    "{result_path}" \\\n')
        agg_f.write(f'    "{"TRUE" if sanity_check else "FALSE"}"\n')
    os.chmod(agg_script, 0o755)

    f.write(f"JOB {noop_job} {SCRIPT_DIR}/noop.sub NOOP\n")
    noop_parents = [job_trait] + neutral_jobs
    f.write(f"PARENT {' '.join(noop_parents)} CHILD {noop_job}\n")
    f.write(f"SCRIPT POST {noop_job} {agg_script}\n\n")

    execute_jobs = 1 + num_neutral_jobs
    return parent_jobs, all_jobs + [noop_job], execute_jobs


def write_fused_trait_section(
    f,
    trait_info: dict,
    params: dict,
    obs_stats_file: Path,
    batch_size: int | None,
    job_prefix: str,
    sanity_check: bool,
    r_env: str | None = None,
    fused_sub: str | None = None,
    environment: str = "",
) -> tuple[list[str], list[str], int]:
    """Single fused job: trait ABC + neutral batch + aggregate."""
    trait_id = trait_info["trait_id"]
    chr_type = trait_info["chr"]
    fst_file = params["fst_chrx"] if chr_type == "chrX" else params["fst_autosomes"]
    fst_values = read_fst_values(fst_file, params["num_neutral"])
    trait_outdir = Path(params["results_dir"]) / f"trait_{trait_id}"
    os.makedirs(trait_outdir, exist_ok=True)

    result_file = f"{trait_id}_result.csv"
    neutral_qst_file = result_file.replace(".csv", "_neutral_qst.txt")
    use_batch = batch_size is not None and batch_size > 1
    if not use_batch:
        raise ValueError("fused mode requires batch_size > 1")

    batch_files = create_fst_batch_files(fst_values, trait_outdir, batch_size)
    batch_file, _n = batch_files[0]
    job_fused = f"fused_{job_prefix}"
    r_env = r_env or resolve_r_env_tarball_for_dag()
    fused_sub = fused_sub or f"{SCRIPT_DIR}/abc_fused_qst.sub"

    f.write(f"JOB {job_fused} {fused_sub}\n")
    f.write(f'VARS {job_fused} mode="full_trait" ')
    f.write(f'input="{obs_stats_file.name}" ')
    f.write(f'fst_input="{batch_file.name}" ')
    f.write(f'output_file="{result_file}" ')
    f.write(f'extra_output_files=", neutral_batch_1.RData, {neutral_qst_file}" ')
    f.write(f'num_sim="{params["num_sim"]}" ')
    f.write(f'summary_stats="{params["summary_stats"]}" ')
    f.write(f'sanity_check="{"TRUE" if sanity_check else "FALSE"}" ')
    f.write(f'extra_input_comma=", " ')
    f.write(f'extra_input_files="{obs_stats_file}, {batch_file.resolve()}" ')
    f.write(f'outdir="{trait_outdir}" ')
    f.write(f'dir_scripts="{SCRIPT_DIR}" ')
    f.write(f'dir_code="{WORKFLOW_SCRIPT_DIR}" ')
    f.write(f'r_env_tarball="{r_env}" ')
    if environment:
        f.write(f'environment="{environment}" ')
    f.write(f'job="{job_fused}"\n\n')

    return [job_fused], [job_fused], 1


def trait_is_complete(trait_outdir: Path, trait_id: str) -> bool:
    result = trait_outdir / f"{trait_id}_result.csv"
    neutral = trait_outdir / "neutral_batch_1.RData"
    if not result.is_file() or not neutral.is_file():
        return False
    return ",unknown,NA," not in result.read_text()


def _prep_trait_worker(args: tuple[str, str, str, str]) -> tuple[str, bool]:
    trait_id, trait_values_s, sample_structure_s, trait_outdir_s = args
    trait_outdir = Path(trait_outdir_s)
    if existing_obs_stats(trait_id, trait_outdir) is not None:
        return trait_id, True
    obs_stats, ratioVext = run_prep_locally(
        trait_id, Path(trait_values_s), Path(sample_structure_s), trait_outdir
    )
    return trait_id, obs_stats is not None and ratioVext is not None


def parallel_prep_traits(
    traits: list[dict],
    trait_values: Path,
    sample_structure: Path,
    results_dir: Path,
    baseline: Path | None,
    workers: int,
) -> list[str]:
    """Prep obs_stats in batch mode using prepare_obs_stats.R."""
    todo_traits = []
    for trait_info in traits:
        trait_id = trait_info["trait_id"]
        trait_outdir = results_dir / f"trait_{trait_id}"
        trait_outdir.mkdir(parents=True, exist_ok=True)
        if existing_obs_stats(trait_id, trait_outdir) is not None:
            continue
        if baseline is not None and copy_obs_stats_from_baseline(trait_id, trait_outdir, baseline):
            continue
        todo_traits.append(trait_id)
        
    if not todo_traits:
        return []
        
    print(f"Batch obs prep: {len(todo_traits)} traits to process", file=sys.stderr)
    
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("\n".join(todo_traits))
        f.write("\n")
        list_file_path = f.name
        
    summary_tsv_path = results_dir / "batch_ratioVext.tsv"
    
    rscript_path = RSCRIPT
    cmd = [
        rscript_path,
        str(WORKFLOW_SCRIPT_DIR / "prepare_obs_stats.R"),
        "--batch",
        str(sample_structure),
        str(trait_values),
        list_file_path,
        str(results_dir),
        str(summary_tsv_path)
    ]
    
    print(f"Running command: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(list_file_path)
    
    if result.returncode != 0:
        print(f"ERROR running batch prep: {result.stderr}", file=sys.stderr)
        return todo_traits
        
    return []


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate HTCondor DAGs for the Detection Module.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--max-traits", type=int, default=None,
                        help="Cap on number of traits (default: all)")
    parser.add_argument("--num-neutral", type=int, default=1000)
    parser.add_argument("--num-sim", type=int, default=100000)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--sanity-check", action="store_true")
    parser.add_argument("--base-priority", type=int, default=1000)
    parser.add_argument("--collection-name", default=DEFAULT_COLLECTION,
                        help="Used as results subdirectory when --results-dir is omitted")
    parser.add_argument("--trait-values", type=Path, default=None)
    parser.add_argument("--sample-structure", type=Path, default=None)
    parser.add_argument("--fst-autosomes", type=Path, default=None)
    parser.add_argument("--fst-chrx", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None,
                        help="Per-trait output directory (default: htcondor/results/<collection-name>)")
    parser.add_argument("--summary-stats", default="QST,ratioVbetweenVtotal")
    parser.add_argument("--threshold-percentile", type=float, default=0.95)
    parser.add_argument("--baseline-obs-dir", type=Path, default=None,
                        help="Optional dir with trait_*/{id}_obs_stats.RData to skip local prep")
    parser.add_argument("--output-dag", type=Path, default=None)
    parser.add_argument(
        "--skip-complete",
        action="store_true",
        help="Skip traits with neutral_batch + valid result.csv (not unknown/NA)",
    )
    parser.add_argument(
        "--prep-workers",
        type=int,
        default=1,
        help="Parallel workers for local obs_stats prep (default: 1)",
    )
    parser.add_argument(
        "--job-mode",
        choices=("split", "fused"),
        default="fused",
        help="fused=one job per trait (default); split=trait+neutral+POST",
    )
    parser.add_argument(
        "--max-jobs-per-dag",
        type=int,
        default=5000,
        help="Max execute jobs per DAG shard (0=unlimited)",
    )
    parser.add_argument(
        "--r-env",
        choices=("osdf", "staging", "home"),
        default=None,
        help="r_env transfer mode (default: R_ENV_TRANSFER or osdf)",
    )
    parser.add_argument(
        "--fused-sub",
        type=Path,
        default=None,
        help="HTCondor submit file for fused jobs (default: abc_fused_qst.sub)",
    )
    parser.add_argument(
        "--code-dir",
        type=Path,
        default=None,
        help="Directory with qst_abc_sim.R (default: workflow/scripts)",
    )
    parser.add_argument(
        "--environment",
        default="FLOOR_POLICY=F3 FLOOR_ALPHA=0.1",
        help="HTCondor job environment (default: F3 ridge floor)",
    )
    args = parser.parse_args()

    global WORKFLOW_SCRIPT_DIR
    if args.code_dir:
        WORKFLOW_SCRIPT_DIR = args.code_dir.resolve()

    fused_r_env = resolve_r_env_tarball(args.r_env) if args.r_env else resolve_r_env_tarball_for_dag()
    fused_sub = str(args.fused_sub) if args.fused_sub else None

    results_dir = (args.results_dir or (HTCONDOR_DIR / "results" / args.collection_name)).resolve()
    dag_path = args.output_dag or (DAG_DIR / "all_traits_flat.dag")
    dag_path.parent.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    trait_values = args.trait_values or (INPUT_DIR / "trait_values.csv")
    sample_structure = args.sample_structure or (INPUT_DIR / "sample_structure.csv")
    fst_auto = args.fst_autosomes or (INPUT_DIR / "qst_neutral_autosomes.txt")
    fst_chrx = args.fst_chrx or (INPUT_DIR / "qst_neutral_chrX.txt")
    traits = read_trait_ids(trait_values, max_traits=args.max_traits)

    baseline = args.baseline_obs_dir
    baseline_dag_dir = None

    params = {
        "fst_autosomes": str(fst_auto),
        "fst_chrx": str(fst_chrx),
        "num_neutral": args.num_neutral,
        "num_sim": args.num_sim,
        "threshold_percentile": args.threshold_percentile,
        "summary_stats": args.summary_stats,
        "results_dir": str(results_dir),
    }
    effective_batch = min(args.batch_size, args.num_neutral) if args.batch_size > 1 else None

    total_execute = 0
    n_copied_obs = 0
    n_skipped = 0

    trait_ids_all = [t["trait_id"] for t in traits]
    prep_failed = parallel_prep_traits(
        traits, trait_values, sample_structure, results_dir, baseline, args.prep_workers
    )
    failed: list[str] = list(prep_failed)
    ratioVext_cache = batch_ratioVext_from_obs(results_dir, trait_ids_all)
    print(f"Preloaded ratioVext for {len(ratioVext_cache)} traits from existing obs_stats", file=sys.stderr)

    max_per_dag = args.max_jobs_per_dag if args.max_jobs_per_dag > 0 else None
    shard_idx = 0
    shard_execute = 0
    dag_f = None
    dag_paths: list[Path] = []

    def start_new_shard() -> None:
        nonlocal shard_idx, shard_execute, dag_f
        if dag_f is not None:
            dag_f.close()
        shard_idx += 1
        shard_execute = 0
        if max_per_dag:
            shard_path = dag_path.parent / f"{dag_path.stem}_shard{shard_idx:03d}{dag_path.suffix}"
        else:
            shard_path = dag_path
        dag_paths.append(shard_path)
        dag_f = open(shard_path, "w")
        dag_f.write(
            f"# Flat harmonizr DAG shard {shard_idx}: {len(traits)} traits, "
            f"num_neutral={args.num_neutral}, job_mode={args.job_mode}\n"
            f"# Aggregation: {'inline fused' if args.job_mode == 'fused' else 'DAG POST on submit node'}\n\n"
        )
        dag_f.write(f"CONFIG {SCRIPT_DIR / 'unlimited.config'}\n\n")

    def ensure_shard_capacity(n_exec: int) -> None:
        if dag_f is None:
            start_new_shard()
        elif max_per_dag and shard_execute > 0 and shard_execute + n_exec > max_per_dag:
            start_new_shard()

    for i, trait_info in enumerate(traits):
        trait_id = trait_info["trait_id"]
        trait_outdir = results_dir / f"trait_{trait_id}"
        trait_outdir.mkdir(parents=True, exist_ok=True)

        if args.skip_complete and trait_is_complete(trait_outdir, trait_id):
            n_skipped += 1
            continue

        obs_stats = existing_obs_stats(trait_id, trait_outdir)
        if obs_stats is None and baseline is not None:
            obs_stats = copy_obs_stats_from_baseline(trait_id, trait_outdir, baseline)
            if obs_stats is not None:
                n_copied_obs += 1

        ratioVext = None
        if args.job_mode == "fused":
            if obs_stats is None:
                obs_stats, ratioVext = run_prep_locally(
                    trait_id, trait_values, sample_structure, trait_outdir
                )
                if ratioVext is None:
                    failed.append(trait_id)
                    continue
        elif obs_stats is None:
            obs_stats, ratioVext = run_prep_locally(
                trait_id, trait_values, sample_structure, trait_outdir
            )
            if ratioVext is None:
                failed.append(trait_id)
                continue
        else:
            ratioVext = ratioVext_cache.get(trait_id)
            if ratioVext is None:
                ratioVext = ratioVext_from_obs_stats_file(obs_stats)
            if ratioVext is None and baseline_dag_dir is not None:
                ratioVext = ratioVext_from_baseline_dag(trait_id, baseline_dag_dir)
            if ratioVext is None:
                failed.append(trait_id)
                continue

        prefix = safe_job_prefix(trait_id)
        priority = args.base_priority
        if args.job_mode == "fused":
            n_exec = 1
            ensure_shard_capacity(n_exec)
            rv_note = f" ratioVext={ratioVext:.5g}" if ratioVext is not None else ""
            dag_f.write(f"# --- trait {trait_id}{rv_note} fused P={priority} ---\n")
            _parents, all_jobs, n_exec = write_fused_trait_section(
                dag_f, trait_info, params, obs_stats,
                effective_batch, prefix, args.sanity_check,
                r_env=fused_r_env, fused_sub=fused_sub,
                environment=args.environment,
            )
        else:
            n_exec = 1 + (1 if effective_batch else params["num_neutral"])
            ensure_shard_capacity(n_exec)
            dag_f.write(f"# --- trait {trait_id} anova_ratioVext={ratioVext:.5g} abcVEVG P={priority} ---\n")
            _parents, all_jobs, n_exec = write_trait_section(
                dag_f, trait_info, params, ratioVext, obs_stats,
                effective_batch, prefix, args.sanity_check,
            )
        for job in all_jobs:
            dag_f.write(f"PRIORITY {job} {priority}\n")
        dag_f.write("\n")
        total_execute += n_exec
        shard_execute += n_exec

        if (i + 1) % 100 == 0:
            print(f"  [{i + 1}/{len(traits)}] {trait_id}", file=sys.stderr)

    if dag_f is not None:
        dag_f.close()

    print(f"DAG: {dag_paths[0] if dag_paths else dag_path}")
    if len(dag_paths) > 1:
        print(f"Shards: {len(dag_paths)} files")
        for p in dag_paths:
            print(f"  {p}")
    print(f"Results: {results_dir}")
    print(f"Traits: {len(traits) - len(failed)} ok, {len(failed)} failed prep")
    print(f"Traits skipped (complete): {n_skipped}")
    print(f"Obs stats copied from baseline: {n_copied_obs}")
    print(f"Execute jobs (approx): {total_execute}")
    if failed:
        print(f"Failed: {', '.join(failed[:10])}{'...' if len(failed) > 10 else ''}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
