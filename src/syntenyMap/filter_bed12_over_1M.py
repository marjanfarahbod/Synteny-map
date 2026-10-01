#!/usr/bin/env python3
"""
Filter merged BED12 entries by length > THRESHOLD bp.

Input (default): data/external/hg38_mm39.synNet.human.merge.bed12
Output (default): data/external/hg38_mm39.synNet.human.merge.<label>filter.bed12
  where <label> is derived from the threshold (e.g., 1Mfilter, 100kfilter).

Length is computed as (end - start + 1) using columns 2 and 3.

Usage:
  python3 exploratory/filter_bed12_over_1M.py [THRESHOLD_BP] [OUT_PATH]

Examples:
  python3 exploratory/filter_bed12_over_1M.py           # defaults to 1,000,000 bp
  python3 exploratory/filter_bed12_over_1M.py 100000    # writes ...merge.100kfilter.bed12
  python3 exploratory/filter_bed12_over_1M.py 500000 my.out.bed12
"""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
IN_PATH = ROOT / "data" / "external" / "hg38_mm39.synNet.human.merge.bed12"
DEFAULT_OUT = ROOT / "data" / "external" / "hg38_mm39.synNet.human.merge.1Mfilter.bed12"


def filter_over_1m(in_path: Path, out_path: Path, threshold: int = 1_000_000) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with in_path.open("r", encoding="utf-8", errors="ignore") as fh, out_path.open("w") as out:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            length = end - start + 1
            if length > threshold:
                out.write(line)
                kept += 1
    return kept


def make_label(threshold: int) -> str:
    if threshold % 1_000_000 == 0:
        return f"{threshold // 1_000_000}Mfilter"
    if threshold % 1_000 == 0:
        return f"{threshold // 1_000}kfilter"
    return f"{threshold}filter"


def main() -> int:
    in_path = IN_PATH
    if not in_path.exists():
        print(f"Input not found: {in_path}", file=sys.stderr)
        return 1

    # Parse threshold and optional output path
    threshold = 1_000_000
    out_path: Path
    if len(sys.argv) > 1:
        try:
            threshold = int(sys.argv[1])
        except ValueError:
            print("First argument must be an integer threshold in bp.", file=sys.stderr)
            return 2
    if len(sys.argv) > 2:
        out_path = Path(sys.argv[2])
    else:
        label = make_label(threshold)
        out_path = in_path.with_name(in_path.name.replace(".merge.bed12", f".merge.{label}.bed12"))

    kept = filter_over_1m(in_path, out_path, threshold=threshold)
    print(f"Kept {kept:,} rows > {threshold:,} bp in {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
