import argparse
import math
from pathlib import Path
import matplotlib.pyplot as plt

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot a histogram of raw top-level chain spans from a chain file."
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
        default=Path("top_level_chain_span_hist.pdf"),
        help="Output path for the histogram"
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=50,
        help="Number of bins for the histogram"
    )
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    chain_spans = []

    try:
        with args.chain.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                # We only need to look at the overarching header for each chain
                if line.startswith("chain "):
                    parts = line.split()
                    if len(parts) < 13:
                        continue
                    
                    t_start = int(parts[5])
                    t_end = int(parts[6])
                    
                    # Calculate raw top-level span (tEnd - tStart)
                    span_length = t_end - t_start
                    chain_spans.append(span_length)
                    
    except FileNotFoundError:
        print(f"Error: Could not find the file at {args.chain}")
        return

    if not chain_spans:
        print("No valid chains found to plot.")
        return

    print(f"Extracted {len(chain_spans)} raw chain spans. Generating histogram...")

    # Convert lengths to log10 space to handle the extreme variance in chain sizes
    log_spans = [math.log10(span) for span in chain_spans if span > 0]

    plt.figure(figsize=(10, 6))
    
    plt.hist(log_spans, bins=args.bins, color='#4C72B0', edgecolor='black')
    
    plt.title("Distribution of Raw Top-Level Chain Spans")
    plt.xlabel("Log10(Chain Span in Base Pairs)")
    plt.ylabel("Frequency")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(args.out, format="pdf")
    print(f"Histogram saved successfully to {args.out}")

if __name__ == "__main__":
    main()
