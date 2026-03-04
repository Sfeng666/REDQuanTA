#!/usr/bin/env python3
"""
Generate neutral FST values from MS coalescent simulation.

This script wraps the existing MS simulation pipeline to:
1. Run MS with user-provided demographic parameters
2. Calculate Reynolds' FST per SNP from simulation output
3. Output FST values for use in QST detection

Usage:
    python fst_from_simulation.py --params demo_params.txt --output-autosomes fst_auto.txt --output-chrx fst_chrx.txt
    python fst_from_simulation.py --params demo_params.txt --output-autosomes fst_auto.txt --output-chrx fst_chrx.txt --num-reps 50000

Reference: code/reference/estimate_neutral_fst_from_simulation/
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate neutral FST from MS simulation"
    )
    
    parser.add_argument("--params", required=True,
                        help="Path to MS demographic parameters file")
    parser.add_argument("--output-autosomes", required=True,
                        help="Output FST file for autosomes")
    parser.add_argument("--output-chrx", required=True,
                        help="Output FST file for X chromosome")
    parser.add_argument("--sample-size", type=int, default=18,
                        help="Sample size per population (default: 18)")
    parser.add_argument("--num-reps", type=int, default=100000,
                        help="Number of independent replicates (default: 100000)")
    parser.add_argument("--window-size", type=int, default=5000,
                        help="Window size in bp (default: 5000)")
    
    return parser.parse_args()


def read_ms_params(params_file):
    """
    Read MS demographic parameters from file.
    
    Expected format:
    - Lines starting with # are comments
    - Each non-comment line is a parameter set for one chromosome type
    - Format: chr_type<TAB>theta<TAB>demo_params
    
    Example:
    chrX    34.66    -en 0 1 3.17 -en 0 2 1.74 ...
    autosomes    39.13    -en 0 1 1.70 -en 0 2 0.33 ...
    """
    params = {}
    
    with open(params_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                chr_type = parts[0].strip()
                theta = float(parts[1].strip())
                demo_params = parts[2].strip()
                params[chr_type] = {'theta': theta, 'demo_params': demo_params}
            elif len(parts) == 2:
                # Simpler format: chr_type<TAB>demo_params (use default theta)
                chr_type = parts[0].strip()
                demo_params = parts[1].strip()
                params[chr_type] = {'theta': None, 'demo_params': demo_params}
    
    return params


def run_ms_simulation(sample_size, num_reps, theta, demo_params, output_file):
    """
    Run MS simulation.
    
    Args:
        sample_size: Sample size per population
        num_reps: Number of independent replicates
        theta: Mutation rate parameter (4*N*mu*L)
        demo_params: Demographic parameters string
        output_file: Path to output file
    
    Returns:
        Path to MS output file
    """
    number_samples = 2 * sample_size  # Total haploid samples
    
    # Build MS command
    cmd = ["ms", str(number_samples), str(num_reps)]
    
    if theta:
        cmd.extend(["-t", str(theta)])
    
    # Add population structure (2 populations)
    cmd.extend(["-I", "2", str(sample_size), str(sample_size), "0"])
    
    # Add demographic parameters
    cmd.extend(demo_params.split())
    
    print(f"Running MS: {' '.join(cmd[:10])}...")
    
    with open(output_file, 'w') as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print(f"MS error: {result.stderr}", file=sys.stderr)
        return None
    
    return output_file


def parse_ms_output(ms_output_file, sample_size):
    """
    Parse MS output and calculate FST for each segregating site.
    
    Args:
        ms_output_file: Path to MS output file
        sample_size: Sample size per population
    
    Yields:
        FST values
    """
    with open(ms_output_file, 'r') as f:
        in_data = False
        
        for line in f:
            line = line.strip()
            
            if line.startswith('//'):
                in_data = True
                continue
            
            if not in_data:
                continue
            
            if line.startswith('segsites:'):
                segsites = int(line.split(':')[1].strip())
                if segsites == 0:
                    in_data = False
                continue
            
            if line.startswith('positions:'):
                continue
            
            if line and len(line) > 0 and line[0] in '01':
                # This is a haplotype line - collect all haplotypes
                haplotypes = [line]
                
                # Read remaining haplotypes for this replicate
                for _ in range(2 * sample_size - 1):
                    hap = f.readline().strip()
                    if hap:
                        haplotypes.append(hap)
                
                # Calculate FST for each site
                if len(haplotypes) == 2 * sample_size:
                    num_sites = len(haplotypes[0])
                    
                    for site_idx in range(num_sites):
                        # Count alleles in each population
                        pop1_count = [0, 0]  # [ancestral, derived]
                        pop2_count = [0, 0]
                        
                        for i, hap in enumerate(haplotypes):
                            allele = int(hap[site_idx])
                            if i < sample_size:
                                pop1_count[allele] += 1
                            else:
                                pop2_count[allele] += 1
                        
                        # Calculate FST
                        fst = calculate_fst_reynolds(pop1_count, pop2_count)
                        if fst is not None:
                            yield fst
                
                in_data = False


def calculate_fst_reynolds(p1_counts, p2_counts):
    """
    Calculate Reynolds' FST for biallelic site.
    
    Args:
        p1_counts: [ancestral_count, derived_count] for population 1
        p2_counts: [ancestral_count, derived_count] for population 2
    
    Returns:
        FST value or None if invalid
    """
    n1 = sum(p1_counts) / 2  # Diploid sample size
    n2 = sum(p2_counts) / 2
    
    if n1 == 0 or n2 == 0:
        return None
    
    total1 = sum(p1_counts)
    total2 = sum(p2_counts)
    
    if total1 == 0 or total2 == 0:
        return None
    
    # Allele frequencies
    p1 = p1_counts[1] / total1  # Derived allele frequency in pop1
    p2 = p2_counts[1] / total2  # Derived allele frequency in pop2
    
    # Reynolds' FST calculation
    first_term = ((p1 - p2)**2) / 2
    
    h1 = 2 * p1 * (1 - p1)  # Expected heterozygosity in pop1
    h2 = 2 * p2 * (1 - p2)  # Expected heterozygosity in pop2
    second_term_num2 = n1 * h1 + n2 * h2
    
    al_frac_denom = 4 * n1 * n2 * (n1 + n2 - 1)
    
    if al_frac_denom == 0:
        return 0
    
    al_frac_num = (n1 + n2) * second_term_num2
    frac = al_frac_num / al_frac_denom
    al = first_term - frac
    
    albl_frac_num = (4 * n1 * n2 - (n1 + n2)) * second_term_num2
    albl_frac = albl_frac_num / al_frac_denom
    albl = first_term + albl_frac
    
    if albl == 0:
        return 0
    
    fst = al / albl
    
    # Correct floating-point precision
    if abs(fst) < 1e-10:
        fst = 0
    
    return fst


def generate_default_demo_params():
    """Generate default demographic parameters based on Drosophila models."""
    return {
        'autosomes': {
            'theta': 39.13,
            'demo_params': '-en 0 1 1.70 -en 0 2 0.33'
        },
        'chrX': {
            'theta': 34.66,
            'demo_params': '-en 0 1 3.17 -en 0 2 1.74'
        }
    }


def main():
    args = parse_args()
    
    # Read parameters
    if os.path.exists(args.params):
        params = read_ms_params(args.params)
    else:
        print(f"Warning: params file not found, using defaults", file=sys.stderr)
        params = generate_default_demo_params()
    
    # Ensure we have parameters for both chromosome types
    if 'autosomes' not in params:
        params['autosomes'] = generate_default_demo_params()['autosomes']
    if 'chrX' not in params:
        params['chrX'] = generate_default_demo_params()['chrX']
    
    # Create output directories
    os.makedirs(os.path.dirname(os.path.abspath(args.output_autosomes)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_chrx)), exist_ok=True)
    
    # Run simulations
    with tempfile.TemporaryDirectory() as tmpdir:
        for chr_type, output_file in [('autosomes', args.output_autosomes),
                                       ('chrX', args.output_chrx)]:
            print(f"\nSimulating {chr_type}...")
            
            ms_output = os.path.join(tmpdir, f"ms_{chr_type}.txt")
            
            chr_params = params.get(chr_type, params.get('autosomes'))
            theta = chr_params.get('theta')
            if theta is None:
                theta = args.window_size * 0.007  # Default per-site rate * window size
            
            # Run MS
            run_ms_simulation(
                sample_size=args.sample_size,
                num_reps=args.num_reps,
                theta=theta,
                demo_params=chr_params['demo_params'],
                output_file=ms_output
            )
            
            # Parse and calculate FST
            fst_values = list(parse_ms_output(ms_output, args.sample_size))
            
            print(f"  Generated {len(fst_values)} FST values")
            
            # Write output
            with open(output_file, 'w') as f:
                for fst in fst_values:
                    f.write(f"{fst:.6f}\n")
    
    print(f"\nOutput written to:")
    print(f"  Autosomes: {args.output_autosomes}")
    print(f"  X chromosome: {args.output_chrx}")


if __name__ == "__main__":
    main()

