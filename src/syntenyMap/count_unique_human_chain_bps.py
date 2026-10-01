from __future__ import annotations

# Purpose: Count unique human bp covered by alignment blocks and top-level chains.
# Inputs:
#   - --chain: data/external/hg38.mm39.all.chain
# Outputs:
#   - Printed summary statistics to stdout
#
# Notes:
#   - The target side of the chain is treated as the human genome.
#   - Aligned block sizes are counted separately from top-level chain spans.
#   - Chain gaps are excluded from block-based coverage but included in the
#     top-level chain-span coverage.
#   - Overlapping and partially overlapping human target intervals are merged
#     before counting bp, so each human base is counted at most once.

import argparse
from array import array
from collections import defaultdict
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count unique human bp mapped to mouse from a chain file."
    )
    parser.add_argument(
        "--chain",
        type=Path,
        default=Path("data/external/hg38.mm39.all.chain"),
        help="Input chain file",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory for histogram outputs",
    )
    return parser.parse_args()


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


def plot_histogram(
    values: list[float] | array,
    out_path: Path,
    title: str,
    xlabel: str,
    weights: list[int] | None = None,
) -> None:
    plt.figure(figsize=(8, 5))
    plt.hist(values, bins=100, weights=weights, color="#4c72b0", edgecolor="black")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def log10_plus_one(values: list[float] | array) -> list[float]:
    return [math.log10(value + 1.0) for value in values]


def main() -> None:
    args = parse_args()
    block_intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    top_level_intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    chain_count = 0
    block_count = 0
    raw_total_bp = 0
    max_block_bp = 0
    current_t_name: str | None = None
    current_t_pos: int | None = None
    in_chain = False
    mean_gap_sizes_per_chain = array("d")
    chain_sizes = array("I")
    gap_size_counts: dict[int, int] = defaultdict(int)
    current_chain_gap_sum = 0
    current_chain_gap_count = 0

    def finalize_chain_gap_stats() -> None:
        if not in_chain:
            return
        if current_chain_gap_count == 0:
            mean_gap_sizes_per_chain.append(0.0)
        else:
            mean_gap_sizes_per_chain.append(current_chain_gap_sum / current_chain_gap_count)

    with args.chain.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            if line.startswith("chain "):
                finalize_chain_gap_stats()
                parts = line.split()
                if len(parts) < 13:
                    current_t_name = None
                    current_t_pos = None
                    continue
                current_t_name = parts[2]
                current_t_pos = int(parts[5])
                in_chain = True
                current_chain_gap_sum = 0
                current_chain_gap_count = 0
                top_level_intervals_by_chrom[current_t_name].append(
                    (current_t_pos, int(parts[6]))
                )
                chain_sizes.append(int(parts[6]) - current_t_pos)
                chain_count += 1
                continue

            if current_t_name is None or current_t_pos is None:
                continue

            parts = line.split()
            if len(parts) not in (1, 3):
                continue

            size = int(parts[0])
            block_start = current_t_pos
            block_end = current_t_pos + size
            block_intervals_by_chrom[current_t_name].append((block_start, block_end))
            block_count += 1
            raw_total_bp += size
            if size > max_block_bp:
                max_block_bp = size

            if len(parts) == 3:
                dt = int(parts[1])
                gap_size_counts[dt] += 1
                current_chain_gap_sum += dt
                current_chain_gap_count += 1
                current_t_pos = block_end + dt
            else:
                finalize_chain_gap_stats()
                in_chain = False
                current_t_name = None
                current_t_pos = None

    finalize_chain_gap_stats()

    unique_total_bp = 0
    merged_interval_count = 0
    unique_top_level_bp = 0
    merged_top_level_interval_count = 0

    print(f"chain\t{args.chain}")
    print(f"chain_count\t{chain_count}")
    print(f"block_count\t{block_count}")
    print(f"raw_total_bp\t{raw_total_bp}")
    print(f"largest_block_bp\t{max_block_bp}")

    for chrom in sorted(block_intervals_by_chrom):
        chrom_total_bp, chrom_merged_count = merge_intervals(block_intervals_by_chrom[chrom])
        unique_total_bp += chrom_total_bp
        merged_interval_count += chrom_merged_count
        print(
            "chrom_block_summary\t"
            f"{chrom}\tunique_bp={chrom_total_bp}\tmerged_intervals={chrom_merged_count}"
        )

    for chrom in sorted(top_level_intervals_by_chrom):
        chrom_total_bp, chrom_merged_count = merge_intervals(
            top_level_intervals_by_chrom[chrom]
        )
        unique_top_level_bp += chrom_total_bp
        merged_top_level_interval_count += chrom_merged_count
        print(
            "chrom_top_level_summary\t"
            f"{chrom}\tunique_bp={chrom_total_bp}\tmerged_intervals={chrom_merged_count}"
        )

    print(f"merged_interval_count\t{merged_interval_count}")
    print(f"unique_total_bp\t{unique_total_bp}")
    print(f"merged_top_level_interval_count\t{merged_top_level_interval_count}")
    print(f"unique_top_level_bp\t{unique_top_level_bp}")

    chain_stem = args.chain.stem
    mean_gap_plot = args.figures_dir / f"{chain_stem}_mean_gap_per_chain_hist.pdf"
    gap_size_plot = args.figures_dir / f"{chain_stem}_gap_size_hist.pdf"
    chain_size_plot = args.figures_dir / f"{chain_stem}_top_level_chain_size_hist.pdf"

    plot_histogram(
        log10_plus_one(mean_gap_sizes_per_chain),
        mean_gap_plot,
        "log10 Mean Target Gap Size Per Chain",
        "log10(mean target gap size + 1)",
    )
    sorted_gap_sizes = sorted(gap_size_counts)
    plot_histogram(
        log10_plus_one(sorted_gap_sizes),
        gap_size_plot,
        "log10 Target Gap Sizes Across All Chains",
        "log10(target gap size + 1)",
        weights=[gap_size_counts[size] for size in sorted_gap_sizes],
    )
    plot_histogram(
        log10_plus_one(chain_sizes),
        chain_size_plot,
        "log10 Top-Level Chain Sizes",
        "log10(top-level chain span + 1)",
    )

    print(f"mean_gap_histogram\t{mean_gap_plot}")
    print(f"gap_size_histogram\t{gap_size_plot}")
    print(f"chain_size_histogram\t{chain_size_plot}")


if __name__ == "__main__":
    main()
