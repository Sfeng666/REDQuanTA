#!/usr/bin/env python3
"""Convert Harmonizr abundance TSV (FeatureID × samples) to REDQuanTA trait_values.csv.

Output columns: trait_id, then sample columns in ascending sample_id order (default 1..36),
matching data/input/sample_structure.csv sample_id values.

Example:
  python harmonizr_to_trait_values.py \\
    ../../../data/input/harmonizr_combat_grouped_TMT9_0426_filterPSMs.tsv \\
    ../../../data/input/trait_values_from_harmonizr.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("harmonizr_tsv", type=Path, help="Tab-separated matrix with FeatureID column")
    parser.add_argument("out_csv", type=Path, help="Output CSV (trait_id + ordered sample columns)")
    parser.add_argument(
        "--sample-ids",
        default=",".join(str(i) for i in range(1, 37)),
        help="Comma-separated sample_id order (default: 1..36)",
    )
    args = parser.parse_args()

    order = [x.strip() for x in args.sample_ids.split(",") if x.strip()]

    with open(args.harmonizr_tsv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit("Empty or invalid TSV")
        id_col = reader.fieldnames[0]
        rows = list(reader)

    missing = [sid for sid in order if sid not in reader.fieldnames]
    if missing:
        raise SystemExit(f"Sample columns missing from TSV header: {missing[:10]}")

    out_headers = ["trait_id"] + order
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(out_headers)
        for row in rows:
            tid = row[id_col].strip()
            if not tid:
                continue
            w.writerow([tid] + [row[sid] for sid in order])

    print(f"Wrote {len(rows)} traits to {args.out_csv.resolve()}")


if __name__ == "__main__":
    main()
