from __future__ import annotations

# Purpose: Build a human x mouse overlap matrix in 1 Mb bins, save it as a pickle, and plot a capped log10 heatmap.
# Inputs:
#   - --bed: data/external/hg38_mm39.synNet.human.bed12
#   - --hg-chrom: data/external/hg38.chrom.sizes
#   - --mm-chrom: data/external/mm39.chrom.sizes
#   - --hg-start: start bin index for human
#   - --mm-start: start bin index for mouse
#   - --regions: number of 1 Mb bins to include per species
# Outputs:
#   - --output: pickle file containing the overlap matrix - "data/processed/hg38_mm39_overlap_1Mb_matrix.pkl
#   - --plot: PDF heatmap of the overlap matrix
#
# Notes:
#   - Bins are built sequentially across chromosomes in chrom.sizes order.
#   - Plot values are capped at 1 Mb, transformed with log10, and zero cells are masked.

import argparse
import pickle
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

Bin = Tuple[str, int, int]


def generate_bins(chrom_sizes: Path, start_bins: int, num_bins: int, bin_size: int = 1_000_000) -> List[Bin]:
    bins: List[Bin] = []
    total_bins_seen = 0
    with chrom_sizes.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, size_str = line.strip().split()[:2]
            size = int(size_str)
            start = 0
            while start < size and len(bins) < num_bins + max(0, start_bins - total_bins_seen):
                end = min(start + bin_size, size)
                if total_bins_seen >= start_bins and len(bins) < num_bins:
                    bins.append((chrom, start, end))
                total_bins_seen += 1
                start += bin_size
            if len(bins) >= num_bins:
                break
    if len(bins) < num_bins:
        raise ValueError(f"Not enough genome length to build {num_bins} bins starting at bin {start_bins}")
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


def write_pickle(matrix: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as f:
        pickle.dump(matrix, f)


def plot_heatmap(matrix: np.ndarray, out_path: Path) -> None:
    clipped = np.clip(matrix, 0, 1_000_000)
    masked = np.ma.masked_less_equal(clipped, 0)
    log_matrix = np.ma.log10(masked)

    plt.figure(figsize=(8, 6))
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color="white")
    im = plt.imshow(log_matrix, aspect="auto", cmap=cmap, vmin=0, vmax=6)
    plt.colorbar(im, label="log10 overlap length (capped at 1 Mb)")
    plt.xlabel("Mouse bins (1 Mb)")
    plt.ylabel("Human bins (1 Mb)")
    plt.title("hg38 vs mm39 overlap (log10, capped at 1 Mb)")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute 1Mb overlap matrix and save to pickle.")
    parser.add_argument("--bed", type=Path, default=Path("data/external/hg38_mm39.synNet.human.bed12"))
    parser.add_argument("--hg-chrom", type=Path, default=Path("data/external/hg38.chrom.sizes"))
    parser.add_argument("--mm-chrom", type=Path, default=Path("data/external/mm39.chrom.sizes"))
    parser.add_argument("--hg-start", type=int, default=0, help="Start bin (1Mb bins) for human")
    parser.add_argument("--mm-start", type=int, default=0, help="Start bin (1Mb bins) for mouse")
    parser.add_argument("--regions", type=int, default=100, help="Number of 1Mb bins to include per species")
    parser.add_argument("--output", type=Path, default=Path("data/processed/hg38_mm39_overlap_1Mb_matrix.pkl"))
    parser.add_argument("--plot", type=Path, default=Path("results/figures/hg38_mm39_overlap_1Mb_matrix.pdf"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hg_bins = generate_bins(args.hg_chrom, start_bins=args.hg_start, num_bins=args.regions)
    mm_bins = generate_bins(args.mm_chrom, start_bins=args.mm_start, num_bins=args.regions)
    matrix = build_matrix(args.bed, hg_bins, mm_bins)
    write_pickle(matrix, args.output)
    plot_heatmap(matrix, args.plot)
    print(
        f"Built matrix of shape {matrix.shape} (hg start bin {args.hg_start}, mm start bin {args.mm_start}) -> {args.output}; plot: {args.plot}"
    )


if __name__ == "__main__":
    main()
