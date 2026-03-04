#!/usr/bin/env python3
"""
Generate HTCondor DAG for QST Performance Evaluation (Module 2).

This script generates jobs to evaluate the performance of ABC-based QST estimation
by simulating traits with known adaptive QST values and calculating TPR/FPR.

Job structure:
1. Neutral QST jobs: Estimate QST for neutral FST values (for FPR calculation)
2. Adaptive QST jobs: Estimate QST for simulated adaptive traits (for TPR calculation)
3. (Local) Aggregation: Calculate TPR/FPR matrices and generate heatmaps

Usage:
    python prepare_perf_eval_dag.py --chr autosomes --output-dir results/perf_eval
    python prepare_perf_eval_dag.py --chr both --num-repeats 100 --num-sim 10000 --batch-size 100
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

# Derive paths relative to this script's location
# Script is at: htcondor/scripts/prepare_perf_eval_dag.py
SCRIPT_DIR = Path(__file__).resolve().parent
HTCONDOR_DIR = SCRIPT_DIR.parent
BASE_DIR = HTCONDOR_DIR.parent
WORKFLOW_DIR = BASE_DIR / "workflow"
WORKFLOW_SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
DATA_DIR = BASE_DIR / "data"
INPUT_DIR = DATA_DIR / "example"

RSCRIPT = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")

# Default parameters
DEFAULT_BATCH_SIZE = 1000  # Same as Module 1
DEFAULT_NUM_NEUTRAL = 10000  # Same as Module 1
DEFAULT_NUM_REPEATS = 10000  # Repeats per adaptive QST level
DEFAULT_NUM_SIM = 100000  # ABC simulations
DEFAULT_THRESHOLD = 0.95

# Default adaptive QST levels (0.50 to 1.00 in 0.05 steps)
DEFAULT_ADAPTIVE_QST = [round(0.50 + i * 0.05, 2) for i in range(11)]

# Default V_E/V_G ratios (5 levels from 0.01 to 100, geometric scale)
# Based on reference codebase: 0.01, 0.1, 1, 10, 100
DEFAULT_VE_RATIOS = [0.01, 0.1, 1.0, 10.0, 100.0]

# Default sample structure (matching original script defaults)
DEFAULT_NUM_POP = 2
DEFAULT_NUM_IND = 6
DEFAULT_NUM_REP = 3

# Default combinations file path
DEFAULT_COMBOS_FILE = BASE_DIR / "data" / "input" / "combinations_sumstats_corrfree.txt"


def parse_summary_stats_arg(summary_stats_arg):
    """
    Parse --summary-stats argument.
    
    Returns:
        tuple: (summary_stats_value, combos_file_path)
        - For file path: (file_path, absolute_file_path)
        - For comma-separated: (comma_string, None)
        - For "all": use default combos file
    """
    if summary_stats_arg == "all":
        # Use default combinations file
        if DEFAULT_COMBOS_FILE.exists():
            return str(DEFAULT_COMBOS_FILE.resolve()), str(DEFAULT_COMBOS_FILE.resolve())
        else:
            raise ValueError(f"'all' specified but default combos file not found: {DEFAULT_COMBOS_FILE}")
    
    # Check if it's a file path
    path = Path(summary_stats_arg)
    if path.exists():
        # Return absolute path for transfer_input_files
        return str(path.resolve()), str(path.resolve())
    
    # Otherwise, treat as comma-separated string
    return summary_stats_arg, None


def read_combos_from_file(combos_file):
    """Read summary stat combinations from file (tab-separated lines)."""
    combos = []
    with open(combos_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                # Tab-separated stats -> comma-separated for R script
                stats = line.split('\t')
                # Convert env_sd to ext_sd, ratioVenv to ratioVext
                stats = [s.replace('env_sd', 'ext_sd').replace('ratioVenv', 'ratioVext') for s in stats]
                combos.append(','.join(stats))
    return combos


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


def generate_perf_eval_dag(dag_path, params, chr_type, output_dir):
    """Generate DAG for performance evaluation.
    
    Args:
        dag_path: Path to output DAG file
        params: Dict with workflow parameters
        chr_type: 'autosomes' or 'chrX'
        output_dir: Output directory for results
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Read neutral FST values
    fst_file = params['fst_chrx'] if chr_type == 'chrX' else params['fst_autosomes']
    fst_values = read_fst_values(fst_file, params['num_neutral'])
    
    all_jobs = []
    neutral_jobs = []
    adaptive_jobs = []
    
    # Environment variables for sample structure
    env_vars = f'SIM_NUM_POP={params["num_pop"]} SIM_NUM_IND={params["num_ind"]} SIM_NUM_REP={params["num_rep"]}'
    
    with open(dag_path, 'w') as f:
        f.write("# HTCondor DAG for QST Performance Evaluation (Module 2)\n")
        f.write(f"# Chromosome: {chr_type}\n")
        f.write(f"# Sample Structure: {params['num_pop']} pops, {params['num_ind']} inds, {params['num_rep']} reps\n")
        f.write(f"# Adaptive QST levels: {params['adaptive_qst']}\n")
        f.write(f"# V_E/V_G ratios: {params['ve_ratios']}\n")
        f.write(f"# Neutral FSTs: {len(fst_values)}, Repeats per condition: {params['num_repeats']}\n\n")
        f.write(f"CONFIG {SCRIPTS_DIR / 'unlimited.config'}\n\n")
        
        # --- Neutral QST Jobs ---
        # For each V_E ratio, estimate QST for all neutral FST values
        f.write("# ===== Neutral QST Jobs (for FPR calculation) =====\n\n")
        
        batch_size = params['batch_size']
        
        for ratio_idx, ve_ratio in enumerate(params['ve_ratios'], 1):
            ratio_dir = Path(output_dir) / f"neutral_ratio_{ratio_idx}"
            ratio_dir_abs = ratio_dir.resolve()
            os.makedirs(ratio_dir_abs, exist_ok=True)
            
            # Compute ext_sd from V_E ratio
            ratioVext = ve_ratio # ratioVext IS ve_ratio * 1.0) ** 0.5  # sqrt(ve_ratio * var_additive)
            
            # Create FST batch files
            batch_files = create_fst_batch_files(fst_values, ratio_dir_abs, batch_size)
            
            f.write(f"# V_E/V_G ratio {ratio_idx}: {ve_ratio} (ratioVext = {ratioVext:.6f})\n")
            for batch_idx, (batch_file, n_fst) in enumerate(batch_files, 1):
                job_name = f"neutral_r{ratio_idx}_b{batch_idx}"
                output_file = f"neutral_batch_{batch_idx}.RData"
                batch_file_abs = batch_file.resolve()
                
                # Handle extra input files (combos file)
                if params.get('combos_file'):
                    extra_input = f', {params["combos_file"]}'
                    # Use basename for summary_stats argument (file will be in working dir)
                    combos_basename = Path(params['combos_file']).name
                    summary_stats_arg = combos_basename
                else:
                    extra_input = ''
                    summary_stats_arg = params["summary_stats"]
                
                f.write(f"JOB {job_name} {SCRIPTS_DIR}/abc_batch_neutral.sub\n")
                f.write(f'VARS {job_name} fst_input="{batch_file.name}" ')
                f.write(f'ratioVext="{ratioVext}" ')
                f.write(f'output_file="{output_file}" ')
                f.write(f'num_sim="{params["num_sim"]}" ')
                f.write(f'summary_stats="{summary_stats_arg}" ')
                f.write(f'fst_file="{batch_file_abs}" ')
                f.write(f'outdir="{ratio_dir_abs}" ')
                f.write(f'dir_scripts="{SCRIPTS_DIR}" ')
                f.write(f'dir_code="{WORKFLOW_SCRIPTS_DIR}" ')
                f.write(f'r_env_tarball="{HTCONDOR_DIR}/env/r_env.tar.gz" ')
                f.write(f'extra_input_files="{extra_input}" ')
                f.write(f'job="{job_name}"\n')
                # Pass sample structure via environment variables
                f.write(f'VARS {job_name} environment="{env_vars}"\n')
                
                all_jobs.append(job_name)
                neutral_jobs.append(job_name)
            
            f.write("\n")
        
        # --- Adaptive QST Jobs ---
        f.write("# ===== Adaptive QST Jobs (for TPR calculation) =====\n\n")
        
        for qst_value in params['adaptive_qst']:
            qst_str = f"{qst_value:.2f}".replace(".", "_")
            
            for ratio_idx, ve_ratio in enumerate(params['ve_ratios'], 1):
                condition_dir = Path(output_dir) / f"adaptive_q{qst_str}_r{ratio_idx}"
                condition_dir_abs = condition_dir.resolve()
                os.makedirs(condition_dir_abs, exist_ok=True)
                
                # Calculate number of batches for repeats
                num_repeats = params['num_repeats']
                num_batches = (num_repeats + batch_size - 1) // batch_size
                
                f.write(f"# QST={qst_value}, V_E/V_G ratio {ratio_idx}: {ve_ratio}\n")
                
                # Handle extra input files (combos file) - same as neutral jobs
                if params.get('combos_file'):
                    extra_input = f', {params["combos_file"]}'
                    combos_basename = Path(params['combos_file']).name
                    summary_stats_arg = combos_basename
                else:
                    extra_input = ''
                    summary_stats_arg = params["summary_stats"]
                
                for batch_idx in range(num_batches):
                    start_id = batch_idx * batch_size + 1
                    n_in_batch = min(batch_size, num_repeats - batch_idx * batch_size)
                    
                    job_name = f"adaptive_q{qst_str}_r{ratio_idx}_b{batch_idx + 1}"
                    output_file = f"adaptive_batch_{batch_idx + 1}.RData"
                    
                    # Input format: "start_id_n_repeats_qst_value" (underscore separator to avoid HTCondor colon issues)
                    input_str = f"{start_id}_{n_in_batch}_{qst_value}"
                    
                    f.write(f"JOB {job_name} {SCRIPTS_DIR}/abc_batch_evaluate.sub\n")
                    f.write(f'VARS {job_name} eval_params="{input_str}" ')
                    f.write(f've_ratio="{ve_ratio}" ')
                    f.write(f'output_file="{output_file}" ')
                    f.write(f'num_sim="{params["num_sim"]}" ')
                    f.write(f'summary_stats="{summary_stats_arg}" ')
                    f.write(f'outdir="{condition_dir_abs}" ')
                    f.write(f'dir_scripts="{SCRIPTS_DIR}" ')
                    f.write(f'dir_code="{WORKFLOW_SCRIPTS_DIR}" ')
                    f.write(f'r_env_tarball="{HTCONDOR_DIR}/env/r_env.tar.gz" ')
                    f.write(f'extra_input_files="{extra_input}" ')
                    f.write(f'job="{job_name}"\n')
                    # Pass sample structure via environment variables
                    f.write(f'VARS {job_name} environment="{env_vars}"\n')
                    
                    all_jobs.append(job_name)
                    adaptive_jobs.append(job_name)
                
                f.write("\n")
        
        # Add NOOP job and aggregation POST script
        noop_job = f"noop_perf_eval_{chr_type}"
        f.write("# ===== Aggregation =====\n\n")
        f.write(f"JOB {noop_job} {SCRIPTS_DIR}/noop.sub NOOP\n")
        f.write(f"PARENT {' '.join(all_jobs)} CHILD {noop_job}\n")
        
        # Create aggregation script - use absolute paths
        output_dir_abs = Path(output_dir).resolve()
        agg_script = output_dir_abs / f"aggregate_perf_eval_{chr_type}.sh"
        rscript_path = RSCRIPT
        
        result_file = output_dir_abs / f"tpr_fpr_matrix_{chr_type}.csv"
        
        with open(agg_script, 'w') as agg_f:
            agg_f.write("#!/bin/bash\n")
            agg_f.write(f'# Aggregation script for performance evaluation ({chr_type})\n')
            agg_f.write(f'# Runs on submit node as POST script after all jobs complete\n\n')
            agg_f.write(f'{rscript_path} {WORKFLOW_SCRIPTS_DIR}/aggregate_perf_eval.R \\\n')
            agg_f.write(f'    "{output_dir_abs}" \\\n')
            agg_f.write(f'    "{chr_type}" \\\n')
            agg_f.write(f'    "{",".join(map(str, params["adaptive_qst"]))}" \\\n')
            agg_f.write(f'    "{",".join(map(str, params["ve_ratios"]))}" \\\n')
            agg_f.write(f'    "{params["threshold_percentile"]}" \\\n')
            agg_f.write(f'    "{result_file}"\n\n')
            agg_f.write(f'echo "Results saved to: {result_file}"\n')
        os.chmod(agg_script, 0o755)
        
        f.write(f"SCRIPT POST {noop_job} {agg_script}\n")
    
    return len(all_jobs), len(neutral_jobs), len(adaptive_jobs), agg_script


def main():
    parser = argparse.ArgumentParser(
        description="Generate HTCondor DAG for QST Performance Evaluation (Module 2)"
    )
    
    # Module 2-specific parameters
    parser.add_argument("--adaptive-qst", type=str, default=None,
                        help=f"Comma-separated adaptive QST levels (default: {DEFAULT_ADAPTIVE_QST})")
    parser.add_argument("--ve-ratios", type=str, default=None,
                        help=f"Comma-separated V_E/V_G ratios (default: {DEFAULT_VE_RATIOS})")
    parser.add_argument("--num-repeats", type=int, default=DEFAULT_NUM_REPEATS,
                        help=f"Repeats per adaptive QST level (default: {DEFAULT_NUM_REPEATS})")
    parser.add_argument("--chr", choices=['autosomes', 'chrX', 'both'], default='autosomes',
                        help="Chromosome type (default: autosomes)")
    
    # Shared parameters from Module 1
    parser.add_argument("--num-neutral", type=int, default=DEFAULT_NUM_NEUTRAL,
                        help=f"Number of neutral FST values (default: {DEFAULT_NUM_NEUTRAL})")
    parser.add_argument("--num-sim", type=int, default=DEFAULT_NUM_SIM,
                        help=f"ABC simulations per estimation (default: {DEFAULT_NUM_SIM})")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                        help=f"FST values per job (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--threshold-percentile", type=float, default=DEFAULT_THRESHOLD,
                        help=f"Threshold for adaptive detection (default: {DEFAULT_THRESHOLD})")
    parser.add_argument("--summary-stats", default="QST,F_within_pop",
                        help="""Summary statistics for ABC:
                        - Comma-separated stats: 'QST,F_within_pop' (default)
                        - Path to file with tab-separated combinations
                        - 'all' to use default combinations file for model ranking""")
    parser.add_argument("--priority", type=int, default=None,
                        help="Job priority for HTCondor scheduling")
    
    # Sample Structure parameters
    parser.add_argument("--num-pop", type=int, default=DEFAULT_NUM_POP,
                        help=f"Number of populations (default: {DEFAULT_NUM_POP})")
    parser.add_argument("--num-ind", type=int, default=DEFAULT_NUM_IND,
                        help=f"Number of individuals per population (default: {DEFAULT_NUM_IND})")
    parser.add_argument("--num-rep", type=int, default=DEFAULT_NUM_REP,
                        help=f"Number of replicates per individual (default: {DEFAULT_NUM_REP})")
    
    # Output
    parser.add_argument("--output-dir", default=str(HTCONDOR_DIR / "results" / "perf_eval"),
                        help="Output directory for results")
    
    args = parser.parse_args()
    
    # Parse adaptive QST levels
    if args.adaptive_qst:
        adaptive_qst = [float(x) for x in args.adaptive_qst.split(',')]
    else:
        adaptive_qst = DEFAULT_ADAPTIVE_QST
    
    # Parse V_E/V_G ratios
    if args.ve_ratios:
        ve_ratios = [float(x) for x in args.ve_ratios.split(',')]
    else:
        ve_ratios = DEFAULT_VE_RATIOS
    
    # Effective batch size
    effective_batch_size = min(args.batch_size, args.num_neutral, args.num_repeats)
    
    # Parameters
    params = {
        'adaptive_qst': adaptive_qst,
        've_ratios': ve_ratios,
        'num_repeats': args.num_repeats,
        'num_neutral': args.num_neutral,
        'num_sim': args.num_sim,
        'batch_size': effective_batch_size,
        'threshold_percentile': args.threshold_percentile,
        'summary_stats': args.summary_stats,
        'fst_autosomes': str(INPUT_DIR / "qst_neutral_autosomes.txt"),
        'fst_chrx': str(INPUT_DIR / "qst_neutral_chrX.txt"),
        'num_pop': args.num_pop,
        'num_ind': args.num_ind,
        'num_rep': args.num_rep,
    }
    
    print("=" * 60)
    print("QST Performance Evaluation DAG Generator (Module 2)")
    print("=" * 60)
    print(f"Adaptive QST levels: {adaptive_qst}")
    print(f"V_E/V_G ratios: {ve_ratios}")
    print(f"Neutral FSTs: {args.num_neutral}")
    print(f"Repeats per condition: {args.num_repeats}")
    print(f"Batch size: {effective_batch_size}")
    print(f"Batch size: {effective_batch_size}")
    print(f"ABC simulations: {args.num_sim}")
    print(f"Sample structure: {args.num_pop} pops, {args.num_ind} inds, {args.num_rep} reps")
    print()
    
    # Determine chromosomes to process
    chr_types = ['autosomes', 'chrX'] if args.chr == 'both' else [args.chr]
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse summary stats argument
    summary_stats_value, combos_file = parse_summary_stats_arg(args.summary_stats)
    
    if combos_file:
        # File path mode - run all combinations on same simulated data
        combos = read_combos_from_file(combos_file)
        print(f"Multi-combo mode: {len(combos)} combinations from {combos_file}")
        print(f"  Combinations will be run on SAME simulated data (efficient)")
        for i, c in enumerate(combos[:5], 1):
            print(f"    {i}. {c}")
        if len(combos) > 5:
            print(f"    ... and {len(combos) - 5} more")
        print()
        params['summary_stats'] = str(combos_file)  # Pass file path to R script
        params['combos_file'] = str(combos_file)
        params['n_combos'] = len(combos)
    else:
        # Single combination mode
        print(f"Single-combo mode: {summary_stats_value}")
        params['summary_stats'] = summary_stats_value
        params['combos_file'] = None
        params['n_combos'] = 1
    
    all_dag_paths = []
    
    # Generate single DAG (not per-combination)
    for chr_type in chr_types:
        print(f"Generating DAG for {chr_type}...")
        
        output_dir = Path(args.output_dir) / chr_type
        dag_path = Path(args.output_dir) / f"perf_eval_{chr_type}.dag"
        
        total_jobs, n_neutral, n_adaptive, agg_script = generate_perf_eval_dag(
            dag_path, params, chr_type, output_dir
        )
        
        print(f"  DAG: {dag_path}")
        print(f"  Jobs: {total_jobs} ({n_neutral} neutral, {n_adaptive} adaptive)")
        print(f"  Aggregation script: {agg_script}")
        all_dag_paths.append(dag_path)
        print()
    
    # Generate model ranking script if using file with multiple combos
    if params.get('combos_file'):
        ranking_script = Path(args.output_dir) / "rank_models.sh"
        combos = read_combos_from_file(params['combos_file'])
        with open(ranking_script, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# Model ranking script - run after DAG completes\n")
            f.write(f"# Evaluates {len(combos)} summary stat combinations\n\n")
            f.write(f'RESULTS_DIR="{args.output_dir}"\n')
            f.write(f'RSCRIPT="{RSCRIPT}"\n')
            f.write(f'PLOT_SCRIPT="{WORKFLOW_SCRIPTS_DIR}/plot_tpr_performance.R"\n\n')
            f.write("# Aggregate results and generate model ranking\n")
            stats_arg = ";".join(combos)
            f.write(f'$RSCRIPT $PLOT_SCRIPT "$RESULTS_DIR" "$RESULTS_DIR/plots" "{stats_arg}"\n')
        os.chmod(ranking_script, 0o755)
        print(f"Model ranking script: {ranking_script}")
    
    print()
    print("=" * 60)
    print("To submit, run:")
    for dag_path in all_dag_paths:
        print(f"  condor_submit_dag {dag_path}")
    if params.get('combos_file'):
        print()
        print("After DAG completes, run model ranking:")
        print(f"  bash {ranking_script}")
    print("=" * 60)


if __name__ == "__main__":
    main()
