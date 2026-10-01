#!/usr/bin/env python3
"""
Exploratory utilities for coverage analysis on hg38/mm39 BED outputs.

Task 1: Sum the coverage (last column) of a BED-like file produced from the
MAF (default: results/hg38_chr19_mm39_coverage.bed).

Usage:
  python3 exploratory/coverage_explore.py [PATH_TO_BED]

If PATH_TO_BED is omitted, uses the default path.
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BED = ROOT / "results" / "hg38_chr19_mm39_coverage.bed"


def human_readable_bp(n: int) -> str:
    if n >= 1_000_000_000:
        return f"{n/1e9:.2f} Gb"
    if n >= 1_000_000:
        return f"{n/1e6:.2f} Mb"
    if n >= 1_000:
        return f"{n/1e3:.2f} kb"
    return f"{n} bp"


def sum_coverage(bed_path: Path) -> int:
    total = 0
    with bed_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts or len(parts) < 9:
                continue
            try:
                size = int(parts[8])
            except ValueError:
                continue
            total += size
    return total


def main(argv: list[str]) -> int:
    bed_path = Path(argv[1]) if len(argv) > 1 else DEFAULT_BED
    if not bed_path.exists():
        print(f"BED not found: {bed_path}", file=sys.stderr)
        return 1
    total = sum_coverage(bed_path)
    print(f"File: {bed_path}")
    print(f"Total coverage (sum of last column): {total:,} bp ({human_readable_bp(total)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

