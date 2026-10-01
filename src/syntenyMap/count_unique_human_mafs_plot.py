import argparse
from pathlib import Path
import matplotlib.pyplot as plt

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract hg38 alignment lengths from a synNet MAF file and plot a histogram."
    )
    parser.add_argument(
        "--maf",
        type=Path,
        default=Path("data/external/hg38.mm39.synNet.maf"),
        help="Input synNet MAF file"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("alignment_lengths_histogram.pdf"),
        help="Output path for the histogram image file"
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
    alignment_lengths = []

    try:
        with args.maf.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("s hg38."):
                    continue
                
                parts = line.split()
                if len(parts) < 6:
                    continue
                
                size = int(parts[3])
                alignment_lengths.append(size)
    except FileNotFoundError:
        print(f"Error: Could not find the file at {args.maf}")
        print("Ensure you are running the script from the correct working directory.")
        return

    if not alignment_lengths:
        print("No hg38 alignments found. Histogram not generated.")
        return

    print(f"Extracted {len(alignment_lengths)} alignment lengths. Generating histogram...")

    plt.figure(figsize=(10, 6))
    
    plt.hist(alignment_lengths, bins=args.bins, color='#4C72B0', edgecolor='black', log=True)
    
    plt.title("Distribution of hg38 Alignment Lengths")
    plt.xlabel("Alignment Length (bp)")
    plt.ylabel("Frequency (Log Scale)")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    # Explicitly saving as PDF
    plt.savefig(args.out, format="pdf")
    print(f"Histogram saved successfully to {args.out}")

if __name__ == "__main__":
    main()
