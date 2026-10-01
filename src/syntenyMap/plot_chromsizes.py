#!/usr/bin/env python3
"""
Generate bar plots of chromosome lengths for human (hg38) and mouse (mm39).

Inputs:
- data/external/hg38.chrom.sizes
- data/external/mm39.chrom.sizes

Outputs (PDF):
- results/figures/hg38_chrom_lengths.pdf
- results/figures/mm39_chrom_lengths.pdf
"""

from pathlib import Path
import csv
import sys
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for saving files
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external"
FIG_DIR = ROOT / "results" / "figures"


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            chrom, size = row[0], row[1]
            try:
                sizes[chrom] = int(size)
            except ValueError:
                # skip malformed lines
                continue
    return sizes


def canonical_order(species: str) -> List[str]:
    if species == "hg38":
        base = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
    elif species == "mm39":
        base = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM", "chrMT"]
    else:
        base = []
    # Deduplicate while preserving order
    seen = set()
    ordered = []
    for c in base:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    return ordered


def filter_canonical(sizes: Dict[str, int], order: List[str]) -> List[Tuple[str, int]]:
    allowed = set(order)
    items = [(c, sizes[c]) for c in order if c in sizes and c in allowed]
    return items


def plot_chrom_lengths(items: List[Tuple[str, int]], title: str, out_path: Path) -> None:
    if not items:
        print(f"No chromosomes to plot for {title}", file=sys.stderr)
        return
    labels = [c for c, _ in items]
    lengths = [l for _, l in items]

    # Convert to megabases for readability
    lengths_mb = [l / 1e6 for l in lengths]

    fig_width = max(8, min(20, 0.35 * len(labels)))
    fig, ax = plt.subplots(figsize=(fig_width, 5))
    ax.bar(labels, lengths_mb, color="#4C78A8")
    ax.set_title(title)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Length (Mb)")
    ax.set_ylim(0, max(lengths_mb) * 1.1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.7)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


def main() -> int:
    hg38_path = DATA_DIR / "hg38.chrom.sizes"
    mm39_path = DATA_DIR / "mm39.chrom.sizes"

    if not hg38_path.exists():
        print(f"Missing file: {hg38_path}", file=sys.stderr)
        return 1
    if not mm39_path.exists():
        print(f"Missing file: {mm39_path}", file=sys.stderr)
        return 1

    hg38_sizes = read_chrom_sizes(hg38_path)
    mm39_sizes = read_chrom_sizes(mm39_path)

    hg38_order = canonical_order("hg38")
    mm39_order = canonical_order("mm39")

    hg38_items = filter_canonical(hg38_sizes, hg38_order)
    mm39_items = filter_canonical(mm39_sizes, mm39_order)

    plot_chrom_lengths(hg38_items, "Human (hg38) Chromosome Lengths", FIG_DIR / "hg38_chrom_lengths.pdf")
    plot_chrom_lengths(mm39_items, "Mouse (mm39) Chromosome Lengths", FIG_DIR / "mm39_chrom_lengths.pdf")

    print("Saved:")
    print(f" - {FIG_DIR / 'hg38_chrom_lengths.pdf'}")
    print(f" - {FIG_DIR / 'mm39_chrom_lengths.pdf'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

