from __future__ import annotations

# Purpose: Find example pairs of overlapping human intervals in hg38_mm39.synNet.human.bed12.
# Inputs:
#   - --bed: data/external/hg38_mm39.synNet.human.bed12
#   - --max: maximum number of overlapping interval pairs to report
# Outputs:
#   - Printed overlap examples to stdout
#
# Notes:
#   - Records are sorted by chromosome and start position before overlap scanning.
#   - Overlaps are detected only among intervals on the same chromosome.

import argparse
from pathlib import Path
from typing import List, Tuple

def read_bed(bed_path: Path) -> List[Tuple[str, int, int, str]]:
    records = []
    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) > 3 else ""
            records.append((chrom, start, end, name))
    records.sort(key=lambda x: (x[0], x[1], x[2]))
    return records

def find_overlaps(records: List[Tuple[str, int, int, str]], max_examples: int = 10) -> List[Tuple[Tuple[str,int,int,str], Tuple[str,int,int,str]]]:
    overlaps = []
    i = 0
    n = len(records)
    while i < n and len(overlaps) < max_examples:
        chrom = records[i][0]
        chrom_records = []
        while i < n and records[i][0] == chrom:
            chrom_records.append(records[i])
            i += 1
        active: List[Tuple[str,int,int,str]] = []
        for rec in chrom_records:
            _, start, end, _ = rec
            active = [r for r in active if r[2] > start]
            for a in active:
                if a[2] > start:
                    overlaps.append((a, rec))
                    if len(overlaps) >= max_examples:
                        return overlaps
            active.append(rec)
    return overlaps

def main() -> None:
    parser = argparse.ArgumentParser(description="Find example overlapping human intervals in BED12")
    parser.add_argument("--bed", type=Path, default=Path("data/external/hg38_mm39.synNet.human.bed12"))
    parser.add_argument("--max", type=int, default=10, help="Number of overlap examples to report")
    args = parser.parse_args()

    recs = read_bed(args.bed)
    overlaps = find_overlaps(recs, max_examples=args.max)
    for idx, (a, b) in enumerate(overlaps, 1):
        print(f"{idx}: {a[0]}\t{a[1]}-{a[2]} ({a[3]}) overlaps {b[0]}\t{b[1]}-{b[2]} ({b[3]})")
    if not overlaps:
        print("No overlaps found")

if __name__ == "__main__":
    main()
