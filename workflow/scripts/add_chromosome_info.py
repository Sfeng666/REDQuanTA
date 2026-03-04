#!/usr/bin/env python3
"""
Add chromosome information ('chr' column) to trait_values.csv by querying UniProt.
Values: 'chrX' for X chromosome, 'autosomes' for others.

Usage:
    python add_chromosome_info.py
    python add_chromosome_info.py --input data/input/trait_values.csv --output data/input/trait_values.csv
"""

import argparse
import csv
import urllib.request
import urllib.parse
import json
import time
import sys

BATCH_SIZE = 100  # UniProt API batch size limit


def get_chromosome_batch(accessions):
    """Query UniProt for chromosome info for a batch of accessions."""
    query = " OR ".join(f"accession:{acc}" for acc in accessions)
    params = {
        "query": query,
        "fields": "accession,xref_proteomes",
        "format": "json",
        "size": len(accessions)
    }
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)
    
    req = urllib.request.Request(url)
    req.add_header("Accept", "application/json")
    
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            data = json.loads(response.read().decode())
            return data.get("results", [])
    except Exception as e:
        print(f"Error querying UniProt: {e}", file=sys.stderr)
        return []


def get_single_entry(accession):
    """Query a single UniProt entry (handles inactive/merged entries)."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            return json.loads(response.read().decode())
    except Exception as e:
        print(f"Error querying {accession}: {e}", file=sys.stderr)
        return None


def get_chromosome_for_merged(accession):
    """Get chromosome info, following merged entry if needed."""
    entry = get_single_entry(accession)
    if not entry:
        return ""
    
    # Check if entry is inactive/merged
    inactive = entry.get("inactiveReason", {})
    if inactive.get("inactiveReasonType") == "MERGED":
        merged_to = inactive.get("mergeDemergeTo", [])
        if merged_to:
            # Follow the merge and get chromosome from new entry
            merged_entry = get_single_entry(merged_to[0])
            if merged_entry:
                return extract_chromosome(merged_entry)
    
    return extract_chromosome(entry)


def extract_chromosome(entry):
    """Extract chromosome from UniProt entry."""
    xrefs = entry.get("uniProtKBCrossReferences", [])
    for xref in xrefs:
        if xref.get("database") == "Proteomes":
            for prop in xref.get("properties", []):
                if prop.get("key") == "Component":
                    return prop.get("value", "")
    return ""


def classify_chromosome(chrom_str):
    """Classify as 'chrX' or 'autosomes' based on chromosome string."""
    if not chrom_str:
        return "unknown"
    chrom_lower = chrom_str.lower()
    if "x" in chrom_lower and "chromosome" in chrom_lower:
        return "chrX"
    return "autosomes"


def main():
    parser = argparse.ArgumentParser(description="Add chromosome info to trait_values.csv via UniProt")
    parser.add_argument("--input", default="data/input/trait_values.csv",
                        help="Input trait values CSV (default: data/input/trait_values.csv)")
    parser.add_argument("--output", default="data/input/trait_values.csv",
                        help="Output CSV (default: same as input, overwrites)")
    args = parser.parse_args()

    # Read input CSV
    with open(args.input, "r") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)
    
    trait_ids = [row["trait_id"] for row in rows]
    print(f"Found {len(trait_ids)} proteins to query")
    
    # Query UniProt in batches
    chr_map = {}
    for i in range(0, len(trait_ids), BATCH_SIZE):
        batch = trait_ids[i:i + BATCH_SIZE]
        print(f"Querying batch {i // BATCH_SIZE + 1}/{(len(trait_ids) + BATCH_SIZE - 1) // BATCH_SIZE}...")
        
        results = get_chromosome_batch(batch)
        for entry in results:
            acc = entry.get("primaryAccession", "")
            chrom = extract_chromosome(entry)
            chr_map[acc] = classify_chromosome(chrom)
        
        time.sleep(0.5)  # Be polite to the API
    
    # Handle unknown entries (likely merged/inactive) by querying individually
    unknown = [tid for tid in trait_ids if chr_map.get(tid) == "unknown"]
    if unknown:
        print(f"\nResolving {len(unknown)} unknown entries (possibly merged)...")
        for acc in unknown:
            chrom = get_chromosome_for_merged(acc)
            chr_map[acc] = classify_chromosome(chrom)
            time.sleep(0.3)
    
    # Add chr column and write output
    new_fieldnames = ["trait_id", "chr"] + [f for f in fieldnames if f not in ("trait_id", "chr")]
    
    for row in rows:
        tid = row["trait_id"]
        row["chr"] = chr_map.get(tid, "unknown")
    
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=new_fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    # Summary
    chrx_count = sum(1 for r in rows if r["chr"] == "chrX")
    auto_count = sum(1 for r in rows if r["chr"] == "autosomes")
    unknown_count = sum(1 for r in rows if r["chr"] == "unknown")
    print(f"\nDone! Results: chrX={chrx_count}, autosomes={auto_count}, unknown={unknown_count}")


if __name__ == "__main__":
    main()

