#!/usr/bin/env python3
"""
Generate DAGs for all traits (or specified number).

This script generates individual trait DAGs with:
  - Batch mode support (multiple FST values per job)
  - Priority-based scheduling (traits complete in order)

Each trait DAG includes:
  - 1 trait QST estimation job
  - N batch neutral QST jobs (default: 200 jobs of 50 FSTs each)
  - A NOOP job with POST script for aggregation

Usage:
    python generate_all_dags.py --max-traits 100 --batch-size 50 --sanity-check
    python generate_all_dags.py --max-traits 100 --num-neutral 10000 --batch-size 100
"""

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path
import concurrent.futures

# Derive paths relative to this script's location
SCRIPT_DIR = Path(__file__).resolve().parent
HTCONDOR_DIR = SCRIPT_DIR.parent
BASE_DIR = HTCONDOR_DIR.parent
WORKFLOW_DIR = BASE_DIR / "workflow"
WORKFLOW_SCRIPTS_DIR = WORKFLOW_DIR / "scripts"
DATA_DIR = BASE_DIR / "data" / "example"
SCRIPTS_DIR = SCRIPT_DIR

# Default batch size (FST values per job) - set to meet 5-min policy
DEFAULT_BATCH_SIZE = 1000


def main():
    parser = argparse.ArgumentParser(description="Generate DAGs for multiple traits")
    parser.add_argument("--max-traits", type=int, default=100,
                        help="Maximum number of traits (default: 100)")
    parser.add_argument("--collection-name", default=None,
                        help="Name for this trait collection, sets dags/NAME and results/NAME dirs "
                             "(default: traits_<max-traits>)")
    parser.add_argument("--num-neutral", type=int, default=10000,
                        help="Number of neutral FST values (default: 10000)")
    parser.add_argument("--num-sim", type=int, default=100000,
                        help="Number of ABC simulations (default: 100000)")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                        help=f"FST values per neutral job (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--sanity-check", action="store_true",
                        help="Enable sanity check mode")
    parser.add_argument("--base-priority", type=int, default=1000,
                        help="Base priority for first trait (default: 1000)")

    args = parser.parse_args()

    collection_name = args.collection_name or f"traits_{args.max_traits}"
    DAGS_DIR = HTCONDOR_DIR / "dags" / collection_name
    RESULTS_DIR = HTCONDOR_DIR / "results" / collection_name
    DAGS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Read trait IDs
    trait_file = DATA_DIR / "trait_values.csv"
    traits = []
    with open(trait_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            traits.append(row['trait_id'])
            if len(traits) >= args.max_traits:
                break
    
    num_batches = (args.num_neutral + args.batch_size - 1) // args.batch_size
    jobs_per_trait = 1 + num_batches  # 1 trait + N batch jobs
    
    print(f"Generating DAGs for {len(traits)} traits...")
    print(f"  Neutral FST values: {args.num_neutral}")
    print(f"  Batch size: {args.batch_size} FSTs/job")
    print(f"  Jobs per trait: {jobs_per_trait}")
    print(f"  ABC simulations: {args.num_sim}")
    print(f"  Sanity check mode: {args.sanity_check}")
    print(f"  Priority: {args.base_priority} (first) to {args.base_priority - len(traits) + 1} (last)")
    print()
    
    total_jobs = 0
    dag_files = []
    failed_traits = []
    
    def worker(i, trait_id):
        priority = args.base_priority - i
        cmd = [
            sys.executable, str(SCRIPTS_DIR / "prepare_trait_dag.py"),
            "--trait-id", trait_id,
            "--num-neutral", str(args.num_neutral),
            "--num-sim", str(args.num_sim),
            "--batch-size", str(args.batch_size),
            "--priority", str(priority),
            "--output-dir", str(DAGS_DIR),
            "--results-dir", str(RESULTS_DIR),
        ]
        if args.sanity_check:
            cmd.append("--sanity-check")

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return False, trait_id, result.stderr
        
        jobs = 0
        for line in result.stdout.split('\n'):
            if 'Total jobs:' in line:
                jobs = int(line.split(':')[1].strip())
        
        dag_file = DAGS_DIR / f"trait_{trait_id}.dag"
        return True, trait_id, (dag_file, priority, jobs, i)

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(worker, i, trait_id): trait_id for i, trait_id in enumerate(traits)}
        
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            success, trait_id, result = future.result()
            completed += 1
            if not success:
                print(f"ERROR generating DAG for {trait_id}: {result}")
                failed_traits.append(trait_id)
            else:
                dag_file, priority, jobs, orig_idx = result
                dag_files.append((dag_file, priority, orig_idx))
                total_jobs += jobs
                print(f"  [{completed}/{len(traits)}] Generated: trait_{trait_id}.dag ({jobs} jobs, P={priority})")

    # Sort dag_files back into original order
    dag_files.sort(key=lambda x: x[2])
    dag_files = [(d[0], d[1]) for d in dag_files]
    
    print()
    print(f"Total DAGs: {len(dag_files)}")
    print(f"Total jobs: {total_jobs}")
    if failed_traits:
        print(f"Failed traits: {', '.join(failed_traits)}")
    
    # Create master DAG that submits all trait DAGs with PRIORITY
    master_dag = DAGS_DIR / "all_traits.dag"
    with open(master_dag, 'w') as f:
        f.write("# Master DAG for all traits\n")
        f.write(f"# Generated for {len(dag_files)} traits with priority-based scheduling\n")
        f.write(f"# Higher priority traits complete first (Trait 1 = P{args.base_priority})\n\n")
        f.write(f"CONFIG {SCRIPTS_DIR / 'unlimited.config'}\n\n")
        
        for dag_file, priority in dag_files:
            trait_id = dag_file.stem.replace('trait_', '')
            f.write(f"SUBDAG EXTERNAL trait_{trait_id} {dag_file}\n")
            f.write(f"PRIORITY trait_{trait_id} {priority}\n\n")
    
    # Create aggregation script for all results - delegates to aggregate_all_results.sh
    agg_all_script = DAGS_DIR / "aggregate_all.sh"
    with open(agg_all_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Aggregate all trait results after master DAG completion.\n")
        f.write("# Usage: bash aggregate_all.sh [threshold] [sanity_check]\n")
        f.write(f"# Run this script after the master DAG completes\n\n")
        f.write(f'bash "{SCRIPTS_DIR}/aggregate_all_results.sh" "{RESULTS_DIR}" "${{1:-0.95}}" "${{2:-FALSE}}"\n')
    os.chmod(agg_all_script, 0o755)
    
    print(f"\nMaster DAG written to: {master_dag}")
    print(f"Aggregation script: {agg_all_script}")
    print(f"\nTo submit all traits:")
    print(f"  condor_submit_dag {master_dag}")
    print(f"\nAfter all DAGs complete, run:")
    print(f"  {agg_all_script}")


if __name__ == "__main__":
    main()

