from __future__ import annotations

# Purpose: Count unique human bp covered by hg38 rows in a synNet MAF file.
# Inputs:
#   - --maf: data/external/hg38.mm39.synNet.maf
# Outputs:
#   - Printed summary statistics to stdout
#
# Notes:
#   - Only lines starting with "s hg38." are counted.
#   - Overlapping and partially overlapping human intervals are merged before
#     counting bp, so each human base is counted at most once.

import argparse
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count unique human bp aligned to mouse from a synNet MAF file."
    )
    parser.add_argument(
        "--maf",
        type=Path,
        default=Path("data/external/hg38.mm39.synNet.maf"),
        help="Input synNet MAF file",
    )
    return parser.parse_args()


def parse_human_interval(parts: list[str]) -> tuple[str, int, int]:
    src = parts[1]
    chrom = src.split(".", 1)[1]
    start = int(parts[2])
    size = int(parts[3])
    strand = parts[4]
    src_size = int(parts[5])

    if strand == "+":
        interval_start = start
        interval_end = start + size
    else:
        interval_start = src_size - (start + size)
        interval_end = src_size - start

    return chrom, interval_start, interval_end


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[int, int]:
    if not intervals:
        return 0, 0

    intervals.sort()
    merged_count = 0
    total_bp = 0
    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= current_end:
            if end > current_end:
                current_end = end
            continue
        total_bp += current_end - current_start
        merged_count += 1
        current_start, current_end = start, end

    total_bp += current_end - current_start
    merged_count += 1
    return total_bp, merged_count


def main() -> None:
    args = parse_args()
    intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    hg38_line_count = 0
    raw_total_bp = 0
    max_size = 0

    with args.maf.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("s hg38."):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            chrom, start, end = parse_human_interval(parts)
            size = end - start
            intervals_by_chrom[chrom].append((start, end))
            hg38_line_count += 1
            raw_total_bp += size
            if size > max_size:
                max_size = size

    unique_total_bp = 0
    merged_interval_count = 0

    print(f"maf\t{args.maf}")
    print(f"hg38_line_count\t{hg38_line_count}")
    print(f"raw_total_bp\t{raw_total_bp}")
    print(f"largest_interval_bp\t{max_size}")

    for chrom in sorted(intervals_by_chrom):
        chrom_total_bp, chrom_merged_count = merge_intervals(intervals_by_chrom[chrom])
        unique_total_bp += chrom_total_bp
        merged_interval_count += chrom_merged_count
        print(
            f"chrom_summary\t{chrom}\tunique_bp={chrom_total_bp}\tmerged_intervals={chrom_merged_count}"
        )

    print(f"merged_interval_count\t{merged_interval_count}")
    print(f"unique_total_bp\t{unique_total_bp}")


if __name__ == "__main__":
    main()
