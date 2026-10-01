from __future__ import annotations

# Purpose: Find the earliest contiguous high-overlap block in a pickled overlap matrix and list BED records for each bin pair.
# Inputs:
#   - --matrix: pickle file containing the human x mouse overlap matrix
#   - --block-size: size k of the contiguous k x k block to find
#   - --threshold: minimum matrix value required in every cell of the block
#   - --bed: data/external/hg38_mm39.synNet.human.bed12
#   - --hg-chrom: data/external/hg38.chrom.sizes
#   - --mm-chrom: data/external/mm39.chrom.sizes
#   - --hg-start: start bin index used when the matrix was generated for human
#   - --mm-start: start bin index used when the matrix was generated for mouse
#   - --regions: number of bins per species represented in the matrix
# Outputs:
#   - Printed selected local/absolute indices, genomic bin coordinates, and overlapping BED records to stdout
#
# Notes:
#   - The selected block is the earliest contiguous k x k block whose values are all >= threshold.
#   - Local indices refer to matrix rows/columns; absolute indices add hg-start/mm-start offsets.

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.overlapgraph.compute_overlap_matrix_pickle import generate_bins
from src.overlapgraph.list_bin_overlaps import list_overlaps


def parse_interval(name: str) -> tuple[str, int, int]:
    chrom, coords = name.split(":")
    start_str, end_str = coords.split("-")
    return chrom, int(start_str), int(end_str)


def find_block(matrix: np.ndarray, block_size: int, threshold: int) -> tuple[int, int]:
    n_rows, n_cols = matrix.shape
    for row_start in range(n_rows - block_size + 1):
        for col_start in range(n_cols - block_size + 1):
            block = matrix[row_start : row_start + block_size, col_start : col_start + block_size]
            if np.all(block >= threshold):
                return row_start, col_start
    raise ValueError(f"No contiguous {block_size}x{block_size} block found with all values >= {threshold}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Investigate the earliest contiguous high-overlap block in a pickled overlap matrix."
    )
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("data/processed/hg38_mm39_overlap_1Mb_matrix_hg150_mm100_100bins.pkl"),
        help="Pickled overlap matrix",
    )
    parser.add_argument("--block-size", type=int, default=3, help="Size of contiguous block to find")
    parser.add_argument("--threshold", type=int, default=1_000_000, help="Minimum value for every block cell")
    parser.add_argument("--bed", type=Path, default=Path("data/external/hg38_mm39.synNet.human.bed12"))
    parser.add_argument("--hg-chrom", type=Path, default=Path("data/external/hg38.chrom.sizes"))
    parser.add_argument("--mm-chrom", type=Path, default=Path("data/external/mm39.chrom.sizes"))
    parser.add_argument("--hg-start", type=int, default=150, help="Human start bin used to generate the matrix")
    parser.add_argument("--mm-start", type=int, default=100, help="Mouse start bin used to generate the matrix")
    parser.add_argument("--regions", type=int, default=100, help="Number of bins per species in the matrix")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.matrix.open("rb") as handle:
        matrix = pickle.load(handle)

    expected_shape = (args.regions, args.regions)
    if matrix.shape != expected_shape:
        raise ValueError(f"Matrix shape {matrix.shape} does not match expected shape {expected_shape}")

    row_start, col_start = find_block(matrix, args.block_size, args.threshold)
    human_local = list(range(row_start, row_start + args.block_size))
    mouse_local = list(range(col_start, col_start + args.block_size))
    human_absolute = [args.hg_start + idx for idx in human_local]
    mouse_absolute = [args.mm_start + idx for idx in mouse_local]

    human_bins = generate_bins(args.hg_chrom, start_bins=args.hg_start, num_bins=args.regions)
    mouse_bins = generate_bins(args.mm_chrom, start_bins=args.mm_start, num_bins=args.regions)
    print(f"matrix\t{args.matrix}")
    print(f"block_size\t{args.block_size}")
    print(f"threshold\t{args.threshold}")
    print(f"human_local_indices\t{','.join(map(str, human_local))}")
    print(f"mouse_local_indices\t{','.join(map(str, mouse_local))}")
    print(f"human_absolute_indices\t{','.join(map(str, human_absolute))}")
    print(f"mouse_absolute_indices\t{','.join(map(str, mouse_absolute))}")

    for human_idx in human_local:
        for mouse_idx in mouse_local:
            h_bin = human_bins[human_idx]
            m_bin = mouse_bins[mouse_idx]
            print(f"\nHuman {args.hg_start + human_idx} and Mouse {args.mm_start + mouse_idx}:")
            rows = list_overlaps(args.bed, h_bin, m_bin)
            if not rows:
                print("  No overlapping records found")
                continue
            for _, (hc, hs, he, name, h_ov, m_ov, combined) in enumerate(rows, 1):
                _, ms, me = parse_interval(name)
                human_syntenic_length = he - hs
                mouse_syntenic_length = me - ms
                print(
                    f"  {hc}:{hs}-{he}\t{name}\t"
                    f"human_syntenic_length={human_syntenic_length}\t"
                    f"mouse_syntenic_length={mouse_syntenic_length}"
                )


if __name__ == "__main__":
    main()
