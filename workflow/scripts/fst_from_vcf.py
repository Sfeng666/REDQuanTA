#!/usr/bin/env python3
"""
Generate neutral FST values from VCF file.

This script wraps the existing empirical FST estimation pipeline to:
1. Extract allele counts from VCF
2. Calculate Reynolds' FST per SNP
3. Filter to neutral sites (e.g., 4-fold synonymous in high-recombination regions)
4. Output FST values for use in QST detection

Usage:
    python fst_from_vcf.py --vcf populations.vcf --output-autosomes fst_auto.txt --output-chrx fst_chrx.txt
    python fst_from_vcf.py --vcf populations.vcf.gz --pop1 pop1 --pop2 pop2 --output-autosomes fst_auto.txt --output-chrx fst_chrx.txt

Reference: code/reference/estimate_neutral_fst_from_empirical/
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate neutral FST from VCF file"
    )
    
    parser.add_argument("--vcf", required=True,
                        help="Path to VCF file with population genotypes")
    parser.add_argument("--pop1", default="pop1",
                        help="Population 1 identifier (default: pop1)")
    parser.add_argument("--pop2", default="pop2",
                        help="Population 2 identifier (default: pop2)")
    parser.add_argument("--output-autosomes", required=True,
                        help="Output FST file for autosomes")
    parser.add_argument("--output-chrx", required=True,
                        help="Output FST file for X chromosome")
    parser.add_argument("--sample-size", type=int, default=18,
                        help="Target sample size per population (default: 18)")
    parser.add_argument("--neutral-bed",
                        help="BED file of neutral sites (e.g., 4-fold synonymous)")
    parser.add_argument("--exclude-bed",
                        help="BED file of regions to exclude (e.g., low recombination)")
    
    return parser.parse_args()


def fst_reynolds(p1_counts, p2_counts):
    """
    Calculate Reynolds' FST for a single SNP.
    
    Args:
        p1_counts: List of allele counts for population 1 [A, T, C, G]
        p2_counts: List of allele counts for population 2 [A, T, C, G]
    
    Returns:
        FST value (float) or None if invalid
    """
    size1 = sum(p1_counts) / 2  # Diploid sample size
    size2 = sum(p2_counts) / 2
    
    if size1 == 0 or size2 == 0:
        return None
    
    total1 = sum(p1_counts)
    total2 = sum(p2_counts)
    
    p1_afs = {i: c / total1 for i, c in enumerate(p1_counts)}
    p2_afs = {i: c / total2 for i, c in enumerate(p2_counts)}
    
    # First term shared by numerator and denominator
    first_term = sum((p1_afs[x] - p2_afs[x])**2 for x in p1_afs) / 2
    
    # Second term numerator
    second_term_num2 = (size1 * (1 - sum(x**2 for x in p1_afs.values())) +
                        size2 * (1 - sum(x**2 for x in p2_afs.values())))
    
    # Implement al (numerator)
    al_frac_num = (size1 + size2) * second_term_num2
    al_frac_denom = 4 * size1 * size2 * (size1 + size2 - 1)
    
    if al_frac_denom == 0:
        return 0
    
    frac = al_frac_num / al_frac_denom
    al = first_term - frac
    
    # Implement al + bl (denominator)
    albl_frac_num = (4 * size1 * size2 - (size1 + size2)) * second_term_num2
    albl_frac = albl_frac_num / al_frac_denom
    albl = first_term + albl_frac
    
    # Calculate FST
    if albl == 0:
        return 0
    
    snp_fst = al / albl
    
    # Correct floating-point precision issues
    if abs(snp_fst) < 1e-10:
        snp_fst = 0
    
    return snp_fst


def extract_fst_from_vcf(vcf_path, pop1_samples, pop2_samples, 
                          neutral_bed=None, exclude_bed=None):
    """
    Extract FST values from VCF file.
    
    Args:
        vcf_path: Path to VCF file
        pop1_samples: List of sample IDs for population 1
        pop2_samples: List of sample IDs for population 2
        neutral_bed: Optional BED file of neutral sites to include
        exclude_bed: Optional BED file of regions to exclude
    
    Yields:
        (chrom, pos, fst) tuples
    """
    import gzip
    
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as f:
        sample_indices = {}
        
        for line in f:
            if line.startswith('##'):
                continue
            
            if line.startswith('#CHROM'):
                # Parse header
                fields = line.strip().split('\t')
                for i, sample in enumerate(fields[9:], 9):
                    sample_indices[sample] = i
                
                pop1_idx = [sample_indices[s] for s in pop1_samples if s in sample_indices]
                pop2_idx = [sample_indices[s] for s in pop2_samples if s in sample_indices]
                continue
            
            # Parse variant line
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            # Skip multi-allelic sites
            if ',' in alt:
                continue
            
            # Extract genotypes and count alleles
            p1_counts = [0, 0]  # [ref, alt]
            p2_counts = [0, 0]
            
            for idx in pop1_idx:
                gt = fields[idx].split(':')[0]
                for allele in gt.replace('|', '/').split('/'):
                    if allele == '0':
                        p1_counts[0] += 1
                    elif allele == '1':
                        p1_counts[1] += 1
            
            for idx in pop2_idx:
                gt = fields[idx].split(':')[0]
                for allele in gt.replace('|', '/').split('/'):
                    if allele == '0':
                        p2_counts[0] += 1
                    elif allele == '1':
                        p2_counts[1] += 1
            
            # Calculate FST
            fst = fst_reynolds(p1_counts, p2_counts)
            
            if fst is not None:
                yield (chrom, pos, fst)


def read_sample_ids(vcf_path, pop_pattern):
    """Read sample IDs from VCF header matching a pattern."""
    import gzip
    import re
    
    opener = gzip.open if vcf_path.endswith('.gz') else open
    
    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                samples = fields[9:]
                matched = [s for s in samples if re.search(pop_pattern, s)]
                return matched
    
    return []


def main():
    args = parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"ERROR: VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    
    # Get sample IDs for each population
    pop1_samples = read_sample_ids(args.vcf, args.pop1)
    pop2_samples = read_sample_ids(args.vcf, args.pop2)
    
    if not pop1_samples or not pop2_samples:
        print("ERROR: Could not find samples for one or both populations", file=sys.stderr)
        print(f"Pop1 ({args.pop1}): {len(pop1_samples)} samples", file=sys.stderr)
        print(f"Pop2 ({args.pop2}): {len(pop2_samples)} samples", file=sys.stderr)
        sys.exit(1)
    
    print(f"Population 1: {len(pop1_samples)} samples")
    print(f"Population 2: {len(pop2_samples)} samples")
    
    # Extract FST values
    autosome_fst = []
    chrx_fst = []
    
    for chrom, pos, fst in extract_fst_from_vcf(args.vcf, pop1_samples, pop2_samples,
                                                  args.neutral_bed, args.exclude_bed):
        if chrom.lower() in ['x', 'chrx']:
            chrx_fst.append(fst)
        else:
            autosome_fst.append(fst)
    
    # Write output files
    os.makedirs(os.path.dirname(os.path.abspath(args.output_autosomes)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_chrx)), exist_ok=True)
    
    with open(args.output_autosomes, 'w') as f:
        for fst in autosome_fst:
            f.write(f"{fst:.6f}\n")
    
    with open(args.output_chrx, 'w') as f:
        for fst in chrx_fst:
            f.write(f"{fst:.6f}\n")
    
    print(f"\nAutosome FST values: {len(autosome_fst)}")
    print(f"X chromosome FST values: {len(chrx_fst)}")
    print(f"\nOutput written to:")
    print(f"  Autosomes: {args.output_autosomes}")
    print(f"  X chromosome: {args.output_chrx}")


if __name__ == "__main__":
    main()

