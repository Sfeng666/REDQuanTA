#!/usr/bin/env python3
"""Summarize QST detection results.

Usage:
    python summarize_results.py
    python summarize_results.py --results-dir code/results --output data/output/combined_results.csv
"""

import argparse
import csv
import os
import re
from datetime import datetime, timedelta
from pathlib import Path
from collections import defaultdict

# Derive paths relative to this script's location
SCRIPT_DIR = Path(__file__).resolve().parent
CHTC_DIR = SCRIPT_DIR.parent
CODE_DIR = CHTC_DIR.parent
BASE_DIR = CODE_DIR.parent
RESULTS_DIR = CHTC_DIR / "results"
OUTPUT_DIR = BASE_DIR / "data" / "output"


def combine_results(results_dir, output_file):
    """Combine all per-trait result CSVs into a single file."""
    results_dir = Path(results_dir)
    all_rows = []
    header = None
    
    for trait_dir in sorted(results_dir.glob("trait_*")):
        if not trait_dir.is_dir():
            continue
        result_files = list(trait_dir.glob("*_result.csv"))
        if not result_files:
            continue
        
        with open(result_files[0], 'r') as f:
            reader = csv.reader(f)
            file_header = next(reader)
            if header is None:
                header = file_header
            for row in reader:
                all_rows.append(row)
    
    if not all_rows:
        print("No results found!")
        return None
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(all_rows)
    
    print(f"Combined {len(all_rows)} trait results to: {output_file}")
    return output_file


def summarize_outcomes(combined_file):
    """Summarize outcomes from combined results."""
    with open(combined_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        print("No data to summarize!")
        return
    
    total = len(rows)
    adaptive = sum(1 for r in rows if r['adaptive'].lower() == 'yes')
    
    qst_values = [float(r['QST']) for r in rows if r['QST'] and r['QST'] != 'NA']
    threshold_values = [float(r['threshold_value']) for r in rows 
                        if r['threshold_value'] and r['threshold_value'] != 'NA']
    
    chr_counts = defaultdict(lambda: {'total': 0, 'adaptive': 0})
    for r in rows:
        chr_counts[r['chr']]['total'] += 1
        if r['adaptive'].lower() == 'yes':
            chr_counts[r['chr']]['adaptive'] += 1
    
    print("\n" + "=" * 60)
    print("QST DETECTION RESULTS SUMMARY")
    print("=" * 60)
    print(f"\nTotal traits analyzed: {total}")
    print(f"  Adaptive (QST > threshold): {adaptive} ({100*adaptive/total:.1f}%)")
    print(f"  Non-adaptive: {total - adaptive} ({100*(total-adaptive)/total:.1f}%)")
    
    if qst_values:
        print(f"\nQST Distribution:")
        print(f"  Min: {min(qst_values):.4f}, Max: {max(qst_values):.4f}")
        print(f"  Mean: {sum(qst_values)/len(qst_values):.4f}")
    
    if threshold_values:
        print(f"\nThreshold Distribution:")
        print(f"  Min: {min(threshold_values):.4f}, Max: {max(threshold_values):.4f}")
        print(f"  Mean: {sum(threshold_values)/len(threshold_values):.4f}")
    
    print(f"\nBreakdown by Chromosome:")
    for chr_type in sorted(chr_counts.keys()):
        c = chr_counts[chr_type]
        pct = 100 * c['adaptive'] / c['total'] if c['total'] > 0 else 0
        print(f"  {chr_type}: {c['adaptive']}/{c['total']} adaptive ({pct:.1f}%)")
    print("=" * 60)


def estimate_time_costs(results_dir):
    """Estimate time costs from job logs."""
    results_dir = Path(results_dir)
    trait_times = []
    
    for trait_dir in sorted(results_dir.glob("trait_*")):
        if not trait_dir.is_dir():
            continue
        trait_id = trait_dir.name.replace("trait_", "")
        job_logs = list(trait_dir.glob("*.log"))
        job_times = []
        
        for log_file in job_logs:
            try:
                with open(log_file, 'r') as f:
                    content = f.read()
                submit = re.search(r'Job submitted.*?(\d{2}/\d{2} \d{2}:\d{2}:\d{2})', content)
                term = re.search(r'Job terminated.*?(\d{2}/\d{2} \d{2}:\d{2}:\d{2})', content)
                if submit and term:
                    year = datetime.now().year
                    s = datetime.strptime(f"{year}/{submit.group(1)}", '%Y/%m/%d %H:%M:%S')
                    t = datetime.strptime(f"{year}/{term.group(1)}", '%Y/%m/%d %H:%M:%S')
                    job_times.append((t - s).total_seconds())
            except:
                continue
        
        if job_times:
            trait_times.append({
                'trait_id': trait_id,
                'n_jobs': len(job_times),
                'total_seconds': sum(job_times),
                'avg_seconds': sum(job_times) / len(job_times)
            })
    
    if not trait_times:
        print("\nNo timing information available.")
        return
    
    print("\n" + "=" * 60)
    print("TIME COST ESTIMATES")
    print("=" * 60)
    
    total_traits = len(trait_times)
    total_jobs = sum(t['n_jobs'] for t in trait_times)
    total_seconds = sum(t['total_seconds'] for t in trait_times)
    avg_job = total_seconds / total_jobs if total_jobs > 0 else 0
    avg_per_trait = total_seconds / total_traits if total_traits > 0 else 0
    
    print(f"\nTraits processed: {total_traits}")
    print(f"Total jobs: {total_jobs}")
    print(f"Total CPU time: {timedelta(seconds=int(total_seconds))}")
    print(f"Average job time: {avg_job:.1f} seconds")
    print(f"Average per trait: {timedelta(seconds=int(avg_per_trait))}")
    
    print(f"\nProjections:")
    for n in [100, 1000, 10000]:
        cpu = avg_per_trait * n
        print(f"  {n} traits: ~{timedelta(seconds=int(cpu))} total CPU time")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Summarize QST results")
    parser.add_argument("--results-dir", default=str(RESULTS_DIR))
    parser.add_argument("--output", default=str(OUTPUT_DIR / "combined_results.csv"))
    args = parser.parse_args()
    
    output_file = combine_results(args.results_dir, args.output)
    if output_file:
        summarize_outcomes(output_file)
        estimate_time_costs(args.results_dir)


if __name__ == "__main__":
    main()
