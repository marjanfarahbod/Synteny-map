import argparse
from collections import defaultdict
from pathlib import Path
import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot chain spans, aligned base pairs, and chromosome lengths."
    )
    parser.add_argument(
        "--chain",
        type=Path,
        default=Path("data/external/hg38.mm39.all.chain"),
        help="Path to the input chain file"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("chain_coverage_by_chromosome.pdf"),
        help="Output path for the plot"
    )
    return parser.parse_args()

def merge_intervals(intervals: list[tuple[int, int]]) -> int:
    """Merges overlapping intervals and returns the total unique length."""
    if not intervals:
        return 0

    intervals.sort()
    total_len = 0
    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= current_end:
            if end > current_end:
                current_end = end
            continue
        total_len += current_end - current_start
        current_start, current_end = start, end

    total_len += current_end - current_start
    return total_len

def main() -> None:
    args = parse_args()
    
    # Store intervals for both top-level spans and actual aligned blocks
    top_level_intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    block_intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    chrom_sizes: dict[str, int] = {}
    
    valid_chrom_pattern = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|[XYM])$")
    
    current_t_name: str | None = None
    current_t_pos: int | None = None

    try:
        with args.chain.open("r", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                
                # 1. Parse Header Lines
                if line.startswith("chain "):
                    parts = line.split()
                    if len(parts) < 13:
                        current_t_name = None
                        current_t_pos = None
                        continue
                    
                    t_name = parts[2]
                    
                    # Filter out non-standard contigs
                    if not valid_chrom_pattern.match(t_name):
                        current_t_name = None
                        current_t_pos = None
                        continue
                        
                    current_t_name = t_name
                    t_size = int(parts[3])
                    current_t_pos = int(parts[5])
                    t_end = int(parts[6])
                    
                    # Log top-level span and true chromosome size
                    top_level_intervals_by_chrom[current_t_name].append((current_t_pos, t_end))
                    if current_t_name not in chrom_sizes:
                        chrom_sizes[current_t_name] = t_size
                    continue
                
                # 2. Parse Data Lines (Blocks and Gaps)
                if current_t_name is None or current_t_pos is None:
                    continue
                    
                parts = line.split()
                if len(parts) not in (1, 3):
                    continue
                    
                size = int(parts[0])
                block_start = current_t_pos
                block_end = current_t_pos + size
                
                # Log the actual aligned block
                block_intervals_by_chrom[current_t_name].append((block_start, block_end))
                
                if len(parts) == 3:
                    dt = int(parts[1])
                    current_t_pos = block_end + dt
                else:
                    # End of chain
                    current_t_name = None
                    current_t_pos = None
                        
    except FileNotFoundError:
        print(f"Error: Could not find the file at {args.chain}")
        return

    if not top_level_intervals_by_chrom:
        print("No valid standard chromosome chains found.")
        return

    # Sort logic (1-22, X, Y, M)
    def chrom_sort_key(c: str) -> int:
        val = c.replace('chr', '')
        if val.isdigit():
            return int(val)
        if val == 'X': return 23
        if val == 'Y': return 24
        if val == 'M': return 25
        return 99

    sorted_chroms = sorted(top_level_intervals_by_chrom.keys(), key=chrom_sort_key)

    chrom_labels = []
    unique_spans_mb = []
    redundant_spans_mb = []
    unique_aligned_bps_mb = []
    actual_sizes_mb = []

    print(f"{'Chrom':<7} | {'Total Span':<12} | {'Unique Span':<12} | {'Unique Aligned BP':<17} | {'Total Size'}")
    print("-" * 75)

    for chrom in sorted_chroms:
        # Calculate span metrics (header)
        span_intervals = top_level_intervals_by_chrom[chrom]
        chrom_raw_span = sum(end - start for start, end in span_intervals)
        chrom_unique_span = merge_intervals(span_intervals)
        chrom_redundant_span = chrom_raw_span - chrom_unique_span
        
        # Calculate actual mapped bases (blocks)
        chrom_unique_bps = merge_intervals(block_intervals_by_chrom[chrom])
        
        chrom_size = chrom_sizes[chrom]
        
        chrom_labels.append(chrom)
        unique_spans_mb.append(chrom_unique_span / 1e6)
        redundant_spans_mb.append(chrom_redundant_span / 1e6)
        unique_aligned_bps_mb.append(chrom_unique_bps / 1e6)
        actual_sizes_mb.append(chrom_size / 1e6)
        
        print(f"{chrom:<7} | {chrom_raw_span:<12,} | {chrom_unique_span:<12,} | {chrom_unique_bps:<17,} | {chrom_size:,}")

    # Generate Plot
    plt.figure(figsize=(15, 8))
    
    # Plot Stacked Bars (Top-Level Spans)
    plt.bar(chrom_labels, unique_spans_mb, color='#4C72B0', edgecolor='black', label='Unique Chain Span (Includes Gaps)')
    plt.bar(chrom_labels, redundant_spans_mb, bottom=unique_spans_mb, color='#DD8452', edgecolor='black', label='Redundant / Overlapping Chain Span')
    
    # Plot Horizontal Markers (Actual Length & Actual Aligned BP)
    for i in range(len(chrom_labels)):
        # Red line for physical chromosome length
        plt.hlines(y=actual_sizes_mb[i], xmin=i - 0.4, xmax=i + 0.4, color='#D62728', linewidth=3, zorder=5)
        # Green line for unique aligned base pairs (will sit inside the blue bar)
        plt.hlines(y=unique_aligned_bps_mb[i], xmin=i - 0.4, xmax=i + 0.4, color='#2CA02C', linewidth=3, zorder=6)

    # Custom Legend Handles
    size_line = mlines.Line2D([], [], color='#D62728', linewidth=3, label='Actual Chromosome Length')
    aligned_bp_line = mlines.Line2D([], [], color='#2CA02C', linewidth=3, label='Unique Aligned BP Coverage (Gap-Free)')
    
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([size_line, aligned_bp_line])

    plt.title("Chain Span & Block Alignment Coverage vs. Actual Chromosome Length", pad=20, fontsize=14, fontweight='bold')
    plt.xlabel("Human Chromosomes (hg38)", fontsize=12)
    plt.ylabel("Length (Megabase Pairs)", fontsize=12)
    
    # Place legend outside the main plot area so it doesn't cover chromosome bars
    plt.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.)
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(args.out, format="pdf")
    print(f"\nPlot saved successfully to {args.out}")

if __name__ == "__main__":
    main()
