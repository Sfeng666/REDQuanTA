#!/usr/bin/env python3
"""
Generate HTCondor DAG for QST detection with proper ABC estimation.

Supports two modes:
1. Single FST per job (original mode) - for backward compatibility
2. Batch FST per job (new mode) - for efficiency and meeting 5-min policy

Job structure per trait:
1. (Local) prepare_obs: Calculate observed statistics and ratioVext
2. trait_qst: ABC estimation of trait QST (1 job)
3. neutral_qst: ABC estimation for neutral FST (N jobs, batch or single)
4. (Local) aggregation: Compare trait QST to neutral distribution

Usage:
    python prepare_trait_dag.py --trait-id L0MQ04
    python prepare_trait_dag.py --trait-id L0MQ04 --num-neutral 10000 --batch-size 1000 --num-sim 100000
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Derive paths relative to this script's location
SCRIPT_DIR = Path(__file__).resolve().parent
HTCONDOR_DIR = SCRIPT_DIR.parent
BASE_DIR = HTCONDOR_DIR.parent
WORKFLOW_DIR = BASE_DIR / "workflow"
WORKFLOW_SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
DATA_DIR = BASE_DIR / "data"
INPUT_DIR = DATA_DIR / "example"  # Default to example data

RSCRIPT = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")

# Default batch size for neutral jobs (FST values per job)
# Set to 1000 for optimal efficiency (based on benchmark results)
# Actual batch size will be min(batch_size, num_neutral) to handle small FST counts
DEFAULT_BATCH_SIZE = 1000


def read_trait_ids(trait_values_path, max_traits=None):
    """Read trait IDs and chromosome info from trait values CSV."""
    traits = []
    with open(trait_values_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            traits.append({'trait_id': row['trait_id'], 'chr': row['chr']})
            if max_traits and len(traits) >= max_traits:
                break
    return traits


def read_fst_values(fst_path, num_values):
    """Read first N FST values from file."""
    values = []
    with open(fst_path, 'r') as f:
        for i, line in enumerate(f):
            if i >= num_values:
                break
            val = line.strip()
            if val:
                values.append(float(val))
    return values


def create_fst_batch_files(fst_values, output_dir, batch_size):
    """Create FST batch files for batch mode processing."""
    num_batches = (len(fst_values) + batch_size - 1) // batch_size
    batch_files = []
    
    for i in range(num_batches):
        start = i * batch_size
        end = min((i + 1) * batch_size, len(fst_values))
        batch_fst = fst_values[start:end]
        
        batch_file = output_dir / f"fst_batch_{i+1}.txt"
        with open(batch_file, 'w') as f:
            for fst in batch_fst:
                f.write(f"{fst}\n")
        batch_files.append((batch_file, len(batch_fst)))
    
    return batch_files


def run_prep_locally(trait_id, trait_values, sample_structure, output_dir):
    """Run prep job locally to get obs_stats and ratioVext."""
    import subprocess
    
    obs_stats_file = output_dir / f"{trait_id}_obs_stats.RData"
    
    # Use Rscript from conda environment
    rscript_path = RSCRIPT
    
    # Run R script to prepare obs stats
    cmd = [
        rscript_path, str(WORKFLOW_SCRIPTS_DIR / "prepare_obs_stats.R"),
        str(trait_values),
        str(sample_structure),
        trait_id,
        str(obs_stats_file)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running prep: {result.stderr}")
        return None, None
    
    # Parse ratioVext from output
    ratioVext = None
    for line in result.stdout.split('\n'):
        if 'ratioVext:' in line:
            ratioVext = float(line.split(':')[1].strip())
            break
    
    return obs_stats_file, ratioVext


def generate_trait_dag(dag_path, trait_info, params, ratioVext, obs_stats_file, 
                       sanity_check=False, batch_size=None, priority=None):
    """Generate DAG for a single trait.
    
    Args:
        dag_path: Path to output DAG file
        trait_info: Dict with trait_id and chr
        params: Dict with workflow parameters
        ratioVext: Extrinsic standard deviation
        obs_stats_file: Path to observed stats RData file
        sanity_check: Enable sanity check mode
        batch_size: FST values per job (None = single FST per job)
        priority: Job priority for this trait (None = no priority set)
    """
    trait_id = trait_info['trait_id']
    chr_type = trait_info['chr']
    
    # Select FST file based on chromosome
    if chr_type == "chrX":
        fst_file = params['fst_chrx']
    else:
        fst_file = params['fst_autosomes']
    
    # Read FST values
    fst_values = read_fst_values(fst_file, params['num_neutral'])
    
    # Output directory
    trait_outdir = Path(params['results_dir']) / f"trait_{trait_id}"
    os.makedirs(trait_outdir, exist_ok=True)
    
    trait_qst_file = f"{trait_id}_trait_qst.RData"
    use_batch = batch_size is not None and batch_size > 1
    
    with open(dag_path, 'w') as f:
        f.write("# HTCondor DAG for QST detection with ABC estimation\n")
        f.write(f"# Trait: {trait_id}, ratioVext: {ratioVext}")
        if use_batch:
            f.write(f", batch_size: {batch_size}")
        if priority is not None:
            f.write(f", priority: {priority}")
        f.write("\n\n")
        f.write(f"CONFIG {SCRIPTS_DIR / 'unlimited.config'}\n\n")
        
        all_jobs = []
        
        # Job 1: Trait QST estimation
        job_trait = f"trait_{trait_id}"
        f.write(f"JOB {job_trait} {SCRIPTS_DIR}/abc_qst.sub\n")
        f.write(f'VARS {job_trait} mode="trait" ')
        f.write(f'input="{obs_stats_file.name}" ')
        f.write(f'ratioVext_or_file="ignored" ')
        f.write(f'output_file="{trait_qst_file}" ')
        f.write(f'num_sim="{params["num_sim"]}" ')
        f.write(f'summary_stats="{params["summary_stats"]}" ')
        f.write(f'extra_input_comma=", " ')
        f.write(f'extra_input_files="{obs_stats_file}" ')
        f.write(f'outdir="{trait_outdir}" ')
        f.write(f'dir_scripts="{SCRIPTS_DIR}" ')
        f.write(f'dir_code="{WORKFLOW_SCRIPTS_DIR}" ')
        f.write(f'r_env_tarball="{HTCONDOR_DIR}/env/r_env.tar.gz" ')
        f.write(f'job="{job_trait}"\n')
        all_jobs.append(job_trait)
        
        if use_batch:
            # Batch mode: multiple FST values per job
            batch_files = create_fst_batch_files(fst_values, trait_outdir, batch_size)
            
            f.write(f"\n# Batch neutral jobs: {len(batch_files)} jobs, {batch_size} FSTs each\n")
            for i, (batch_file, n_fst) in enumerate(batch_files, 1):
                job_neutral = f"batch_{trait_id}_{i}"
                neutral_file = f"neutral_batch_{i}.RData"
                
                f.write(f"JOB {job_neutral} {SCRIPTS_DIR}/abc_batch_neutral.sub\n")
                f.write(f'VARS {job_neutral} fst_input="{batch_file.name}" ')
                f.write(f'ratioVext="{ratioVext}" ')
                f.write(f'output_file="{neutral_file}" ')
                f.write(f'num_sim="{params["num_sim"]}" ')
                f.write(f'summary_stats="{params["summary_stats"]}" ')
                f.write(f'fst_file="{batch_file}" ')
                f.write(f'outdir="{trait_outdir}" ')
                f.write(f'dir_scripts="{SCRIPTS_DIR}" ')
                f.write(f'dir_code="{WORKFLOW_SCRIPTS_DIR}" ')
                f.write(f'r_env_tarball="{HTCONDOR_DIR}/env/r_env.tar.gz" ')
                f.write(f'job="{job_neutral}"\n')
                all_jobs.append(job_neutral)
            
            num_neutral_jobs = len(batch_files)
        else:
            # Single FST per job (original mode)
            f.write("\n# Neutral jobs: 1 FST per job\n")
            for i, fst_val in enumerate(fst_values, 1):
                job_neutral = f"neutral_{trait_id}_{i}"
                neutral_file = f"neutral_{i}.RData"
                
                f.write(f"JOB {job_neutral} {SCRIPTS_DIR}/abc_neutral.sub\n")
                f.write(f'VARS {job_neutral} fst_value="{fst_val}" ')
                f.write(f'ratioVext="{ratioVext}" ')
                f.write(f'output_file="{neutral_file}" ')
                f.write(f'num_sim="{params["num_sim"]}" ')
                f.write(f'summary_stats="{params["summary_stats"]}" ')
                f.write(f'outdir="{trait_outdir}" ')
                f.write(f'dir_scripts="{SCRIPTS_DIR}" ')
                f.write(f'dir_code="{WORKFLOW_SCRIPTS_DIR}" ')
                f.write(f'r_env_tarball="{HTCONDOR_DIR}/env/r_env.tar.gz" ')
                f.write(f'job="{job_neutral}"\n')
                all_jobs.append(job_neutral)
            
            num_neutral_jobs = len(fst_values)
        
        f.write("\n")
        
        # Add PRIORITY for all jobs if specified
        if priority is not None:
            f.write(f"# Priority for sequential trait completion\n")
            for job in all_jobs:
                f.write(f"PRIORITY {job} {priority}\n")
            f.write("\n")
        
        # No dependencies between trait and neutral - all run in parallel
        f.write("# All ABC jobs run independently (prep was done locally)\n\n")
        
        # Create aggregation script that runs on submit node
        agg_script = trait_outdir / f"aggregate_{trait_id}.sh"
        result_file = trait_outdir / f"{trait_id}_result.csv"
        rscript_path = RSCRIPT
        
        with open(agg_script, 'w') as agg_f:
            agg_f.write("#!/bin/bash\n")
            agg_f.write(f'# Aggregation script for trait {trait_id}\n')
            agg_f.write(f'# Runs on submit node as POST script after all jobs complete\n')
            agg_f.write(f'{rscript_path} {WORKFLOW_SCRIPTS_DIR}/aggregate_qst.R \\\n')
            agg_f.write(f'    "{trait_outdir}/{trait_qst_file}" \\\n')
            agg_f.write(f'    "{trait_outdir}" \\\n')
            agg_f.write(f'    "{params["threshold_percentile"]}" \\\n')
            agg_f.write(f'    "{result_file}" \\\n')
            agg_f.write(f'    "{"TRUE" if sanity_check else "FALSE"}"\n')
            agg_f.write(f'echo "Result saved to: {result_file}"\n')
        os.chmod(agg_script, 0o755)
        
        # Add NOOP job that depends on all other jobs, with POST script for aggregation
        noop_job = f"noop_{trait_id}"
        f.write(f"# NOOP job that triggers aggregation after all jobs complete\n")
        f.write(f"JOB {noop_job} {SCRIPTS_DIR}/noop.sub NOOP\n")
        
        # NOOP depends on all jobs
        f.write(f"PARENT {' '.join(all_jobs)} CHILD {noop_job}\n")
        
        # POST script runs on submit node after NOOP completes
        f.write(f"SCRIPT POST {noop_job} {agg_script}\n")
    
    total_jobs = 1 + num_neutral_jobs  # trait + neutral jobs (NOOP doesn't count)
    return total_jobs, agg_script, all_jobs


def main():
    parser = argparse.ArgumentParser(description="Generate HTCondor DAG for QST detection")
    parser.add_argument("--trait-id", required=True, help="Trait ID to process")
    parser.add_argument("--num-neutral", type=int, default=10000,
                        help="Number of neutral FST values (default: 10000)")
    parser.add_argument("--num-sim", type=int, default=100000,
                        help="Number of ABC simulations (default: 100000)")
    parser.add_argument("--threshold-percentile", type=float, default=0.95,
                        help="Threshold percentile (default: 0.95)")
    parser.add_argument("--summary-stats", default="QST,F_within_pop",
                        help="Summary statistics for ABC (default: QST,F_within_pop)")
    parser.add_argument("--sanity-check", action="store_true",
                        help="Enable sanity check mode")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                        help=f"FST values per neutral job (default: {DEFAULT_BATCH_SIZE}, use 1 for single FST per job)")
    parser.add_argument("--priority", type=int, default=None,
                        help="Job priority for this trait (higher = scheduled first)")
    parser.add_argument("--results-dir",
                        default=str(HTCONDOR_DIR / "results" / "traits_100"),
                        help="Directory for per-trait result outputs (default: results/traits_100)")
    parser.add_argument("--output-dir", default=str(HTCONDOR_DIR / "dags"),
                        help="Output directory for DAG files")
    
    args = parser.parse_args()
    
    # Input file paths
    trait_values = INPUT_DIR / "trait_values.csv"
    sample_structure = INPUT_DIR / "sample_structure.csv"
    fst_autosomes = INPUT_DIR / "qst_neutral_autosomes.txt"
    fst_chrx = INPUT_DIR / "qst_neutral_chrX.txt"
    
    # Read trait info
    traits = read_trait_ids(trait_values)
    trait_info = None
    for t in traits:
        if t['trait_id'] == args.trait_id:
            trait_info = t
            break
    
    if trait_info is None:
        print(f"ERROR: Trait ID not found: {args.trait_id}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory for this trait's results
    trait_outdir = Path(args.results_dir) / f"trait_{args.trait_id}"
    os.makedirs(trait_outdir, exist_ok=True)
    
    # Step 1: Run prep LOCALLY to get ratioVext
    print(f"Step 1: Running prep locally for trait {args.trait_id}...")
    obs_stats_file, ratioVext = run_prep_locally(
        args.trait_id, trait_values, sample_structure, trait_outdir
    )
    
    if ratioVext is None:
        print("ERROR: Failed to get ratioVext from prep job")
        sys.exit(1)
    
    print(f"  ratioVext = {ratioVext}")
    print(f"  obs_stats = {obs_stats_file}")
    
    # Parameters
    params = {
        "trait_values": str(trait_values),
        "sample_structure": str(sample_structure),
        "fst_autosomes": str(fst_autosomes),
        "fst_chrx": str(fst_chrx),
        "num_neutral": args.num_neutral,
        "num_sim": args.num_sim,
        "threshold_percentile": args.threshold_percentile,
        "summary_stats": args.summary_stats,
        "results_dir": args.results_dir,
    }
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Step 2: Generate DAG
    print(f"\nStep 2: Generating DAG...")
    # Actual batch size is min of user-defined batch size and total neutral FST count
    # This handles cases where num_neutral < batch_size
    effective_batch_size = min(args.batch_size, args.num_neutral) if args.batch_size > 1 else None
    dag_path = os.path.join(args.output_dir, f"trait_{args.trait_id}.dag")
    total_jobs, agg_script, all_jobs = generate_trait_dag(
        dag_path, trait_info, params, ratioVext, obs_stats_file, 
        args.sanity_check, effective_batch_size, args.priority
    )
    
    print(f"\nDAG written to: {dag_path}")
    print(f"Total jobs: {total_jobs}")
    print(f"  - 1 trait QST job")
    if effective_batch_size:
        num_batches = (args.num_neutral + effective_batch_size - 1) // effective_batch_size
        print(f"  - {num_batches} batch neutral jobs ({effective_batch_size} FSTs each)")
    else:
        print(f"  - {args.num_neutral} neutral QST jobs")
    if args.priority is not None:
        print(f"  - Priority: {args.priority}")
    print(f"\nSubmit with: condor_submit_dag {dag_path}")
    print(f"\nAfter DAG completes, run aggregation:")
    print(f"  {agg_script}")


if __name__ == "__main__":
    main()

