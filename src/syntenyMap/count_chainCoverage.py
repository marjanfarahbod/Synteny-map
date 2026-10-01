import argparse
from collections import defaultdict
from pathlib import Path
import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot unique and total chain span coverage per standard chromosome."
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
    """Merges overlapping intervals and returns the total unique span."""
    if not intervals:
        return 0

    intervals.sort()
    total_span = 0
    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= current_end:
            if end > current_end:
                current_end = end
            continue
        total_span += current_end - current_start
        current_start, current_end = start, end

    total_span += current_end - current_start
    return total_span

def main() -> None:
    args = parse_args()
    top_level_intervals_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    chrom_sizes: dict[str, int] = {}
    
    # Regex to match only standard chromosomes (chr1-22, chrX, chrY, chrM)
    valid_chrom_pattern = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|[XYM])$")

    try:
        with args.chain.open("r", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                
                if line.startswith("chain "):
                    parts = line.split()
                    if len(parts) < 13:
                        continue
                    
                    t_name = parts[2]
                    
                    # Skip contigs with extended names
                    if not valid_chrom_pattern.match(t_name):
                        continue
                        
                    t_size = int(parts[3])
                    t_start = int(parts[5])
                    t_end = int(parts[6])
                    
                    top_level_intervals_by_chrom[t_name].append((t_start, t_end))
                    
                    # Store actual chromosome length for the plot line
                    if t_name not in chrom_sizes:
                        chrom_sizes[t_name] = t_size
                        
    except FileNotFoundError:
        print(f"Error: Could not find the file at {args.chain}")
        return

    if not top_level_intervals_by_chrom:
        print("No valid standard chromosome chains found.")
        return

    # Sort chromosomes logically (1-22, X, Y, M) rather than alphabetically
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
    actual_sizes_mb = []

    print(f"{'Chromosome':<15} {'Total Chain Span':<20} {'Unique Chain Span':<20} {'Chromosome Size'}")
    print("-" * 75)

    for chrom in sorted_chroms:
        intervals = top_level_intervals_by_chrom[chrom]
        chrom_size = chrom_sizes[chrom]
        
        chrom_raw_coverage = sum(end - start for start, end in intervals)
        chrom_unique_coverage = merge_intervals(intervals)
        chrom_redundant_coverage = chrom_raw_coverage - chrom_unique_coverage
        
        chrom_labels.append(chrom)
        # Convert base pairs to Megabases (Mb) for cleaner plot axes
        unique_spans_mb.append(chrom_unique_coverage / 1e6)
        redundant_spans_mb.append(chrom_redundant_coverage / 1e6)
        actual_sizes_mb.append(chrom_size / 1e6)
        
        print(f"{chrom:<15} {chrom_raw_coverage:<20,} {chrom_unique_coverage:<20,} {chrom_size:,}")

    # Generate Plot
    plt.figure(figsize=(14, 7))
    
    # 1. Plot the unique coverage (bottom of stack)
    plt.bar(chrom_labels, unique_spans_mb, color='#4C72B0', edgecolor='black', label='Unique Chain Span')
    
    # 2. Plot the redundant coverage on top of the unique coverage
    plt.bar(chrom_labels, redundant_spans_mb, bottom=unique_spans_mb, color='#DD8452', edgecolor='black', label='Overlapping (Redundant) Span')
    
    # 3. Draw horizontal lines representing the actual size of each chromosome
    for i, size in enumerate(actual_sizes_mb):
        plt.hlines(y=size, xmin=i - 0.4, xmax=i + 0.4, color='red', linewidth=3, zorder=5)

    # Create a custom legend handle for the horizontal line
    size_line = mlines.Line2D([], [], color='red', linewidth=3, label='Actual Chromosome Length')
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.append(size_line)

    plt.title("Chain Span Coverage vs. Actual Chromosome Length", pad=20, fontsize=14)
    plt.xlabel("Human Chromosomes (hg38)", fontsize=12)
    plt.ylabel("Length (Megabase Pairs)", fontsize=12)
    plt.legend(handles=handles, loc='upper right')
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(args.out, format="pdf")
    print(f"\nPlot saved successfully to {args.out}")

if __name__ == "__main__":
    main()
