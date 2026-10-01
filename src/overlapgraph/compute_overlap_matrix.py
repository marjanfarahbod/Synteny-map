from __future__ import annotations

# Purpose: Build a 100x100 overlap matrix between the first 100 Mb of hg38 and mm39 (1 Mb bins)
# using alignments from hg38_mm39.synNet.human.bed12.
# Inputs:
#   - --bed: data/external/hg38_mm39.synNet.human.bed12
#   - --hg-chrom: data/external/hg38.chrom.sizes
#   - --mm-chrom: data/external/mm39.chromsizes
# Outputs:
#   - --output: data/processed/hg38_mm39_overlap_1Mb_matrix.tsv (100x100 table of overlap lengths)
#
# Notes:
#   - Bins are built sequentially across chromosomes in chrom.sizes order until 100 bins (1 Mb each).
#   - Overlap contribution for a bin pair = min(overlap in human bin, overlap in mouse bin) for each fill.

import argparse
from pathlib import Path
from typing import List, Tuple

import numpy as np

Bin = Tuple[str, int, int]


def load_bins(chrom_sizes: Path, bin_size: int = 1_000_000, max_bins: int = 100) -> List[Bin]:
    bins: List[Bin] = []
    with chrom_sizes.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, size_str = line.strip().split()[:2]
            size = int(size_str)
            start = 0
            while start < size and len(bins) < max_bins:
                end = min(start + bin_size, size)
                bins.append((chrom, start, end))
                start += bin_size
            if len(bins) >= max_bins:
                break
    return bins


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def bins_for_interval(chrom: str, start: int, end: int, bins: List[Bin]) -> List[int]:
    return [
        idx
        for idx, (b_chrom, b_start, b_end) in enumerate(bins)
        if b_chrom == chrom and interval_overlap(start, end, b_start, b_end) > 0
    ]


def parse_mouse_interval(name_field: str) -> Tuple[str, int, int] | None:
    try:
        q_chrom, pos = name_field.split(":")
        q_start_str, q_end_str = pos.split("-")
        return q_chrom, int(q_start_str), int(q_end_str)
    except ValueError:
        return None


def build_matrix(bed_path: Path, hg_bins: List[Bin], mm_bins: List[Bin]) -> np.ndarray:
    matrix = np.zeros((len(hg_bins), len(mm_bins)), dtype=np.int64)
    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("	")
            if len(parts) < 4:
                continue
            h_chrom, h_start, h_end = parts[0], int(parts[1]), int(parts[2])
            mouse_iv = parse_mouse_interval(parts[3])
            if mouse_iv is None:
                continue
            m_chrom, m_start, m_end = mouse_iv

            h_indices = bins_for_interval(h_chrom, h_start, h_end, hg_bins)
            m_indices = bins_for_interval(m_chrom, m_start, m_end, mm_bins)
            if not h_indices or not m_indices:
                continue

            for hi in h_indices:
                h_overlap = interval_overlap(h_start, h_end, hg_bins[hi][1], hg_bins[hi][2])
                if h_overlap <= 0:
                    continue
                for mi in m_indices:
                    m_overlap = interval_overlap(m_start, m_end, mm_bins[mi][1], mm_bins[mi][2])
                    if m_overlap <= 0:
                        continue
                    matrix[hi, mi] += min(h_overlap, m_overlap)
    return matrix


def write_matrix(matrix: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        f.write("\t".join([""] + [f"mouse_bin_{i+1}" for i in range(matrix.shape[1])]) + "\n")
        for i in range(matrix.shape[0]):
            row = [f"human_bin_{i+1}"] + [str(val) for val in matrix[i]]
            f.write("\t".join(row) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute 1Mb bin overlap matrix between hg38 and mm39.")
    parser.add_argument("--bed", type=Path, default=Path("data/external/hg38_mm39.synNet.human.bed12"))
    parser.add_argument("--hg-chrom", type=Path, default=Path("data/external/hg38.chrom.sizes"))
    parser.add_argument("--mm-chrom", type=Path, default=Path("data/external/mm39.chrom.sizes"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/hg38_mm39_overlap_1Mb_matrix.tsv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hg_bins = load_bins(args.hg_chrom)
    mm_bins = load_bins(args.mm_chrom)
    matrix = build_matrix(args.bed, hg_bins, mm_bins)
    write_matrix(matrix, args.output)
    print(f"Built matrix of shape {matrix.shape} -> {args.output}")


if __name__ == "__main__":
    main()
