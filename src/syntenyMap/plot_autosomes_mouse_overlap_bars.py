#!/usr/bin/env python3
"""
Plot colored overlap bars for all human autosomes (chr1–chr22) from the
10k-filtered merged BED12, with bar lengths proportional to chromosome sizes.

Design:
- Overall width corresponds to chr1 at 20 cm. All bars share the same x-scale
  (0..len(chr1)), and each chromosome's bar spans its own length, leaving any
  remaining width blank.
- For each human block, color by the mouse chromosome (chr1..chr19) using the
  rainbow colormap; ignore other mouse chromosomes.
- Excludes chrX and chrY.
- Adds a legend mapping mouse chr1..chr19 to colors.

Inputs:
- data/external/hg38_mm39.synNet.human.merge.10kfilter.bed12
- data/external/hg38.chrom.sizes

Output:
- results/figures/hg38_autosomes_mouse_overlap_bars_10kfilter.pdf
"""

from pathlib import Path
from typing import Dict, List, Tuple, DefaultDict
from collections import defaultdict
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import colormaps
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external"
FIG_DIR = ROOT / "results" / "figures"

BED12_PATH = DATA_DIR / "hg38_mm39.synNet.human.merge.10kfilter.bed12"
CHROMSIZES = DATA_DIR / "hg38.chrom.sizes"


def read_chrom_sizes(path: Path) -> Dict[str, int]:
    sizes: Dict[str, int] = {}
    with path.open("r") as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.strip().split("\t")
            try:
                sizes[chrom] = int(size)
            except ValueError:
                continue
    return sizes


def parse_blocks_by_chrom(bed12_path: Path) -> DefaultDict[str, List[Tuple[int, int, str]]]:
    # Returns mapping human chrom -> list of (start, end, mouse_chr)
    by_chrom: DefaultDict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    with bed12_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            hchr = parts[0]
            if not hchr.startswith("chr"):
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            name = parts[3]
            mchr = name.split(":", 1)[0]
            by_chrom[hchr].append((start, end, mchr))
    return by_chrom


def mouse_chr_color_map() -> Dict[str, Tuple[float, float, float, float]]:
    # Map chr1..chr19 along the "terrain" colormap
    cmap = colormaps["terrain"]
    colors: Dict[str, Tuple[float, float, float, float]] = {}
    # Skip the top 5% of the colormap to avoid near-white colors
    vals = np.linspace(0.0, 0.95, 19)
    for i in range(1, 20):
        colors[f"chr{i}"] = cmap(vals[i - 1])
    return colors


def plot_bars(by_chrom: Dict[str, List[Tuple[int, int, str]]], sizes: Dict[str, int], out_path: Path) -> None:
    # Order chromosomes 1..22
    order = [f"chr{i}" for i in range(1, 23)]
    # chr1 sets the width scale
    if "chr1" not in sizes:
        raise RuntimeError("chr1 not found in chrom sizes")
    chr1_len = sizes["chr1"]

    # Figure size: width 20 cm (chr1), height ~ 1 cm per bar + small gaps
    bar_cm = 1.0
    gap_cm = 0.1
    fig_w_in = 20.0 / 2.54
    fig_h_in = ((bar_cm + gap_cm) * len(order) - gap_cm) / 2.54
    fig, axes = plt.subplots(nrows=len(order), ncols=1, figsize=(fig_w_in, fig_h_in))
    if len(order) == 1:
        axes = [axes]

    mcolors = mouse_chr_color_map()

    for ax, chrom in zip(axes, order):
        ax.set_xlim(0, chr1_len)
        ax.set_ylim(0, 1)
        ax.axis("off")

        # Base white bar for this chromosome length
        clen = sizes.get(chrom, 0)
        if clen <= 0:
            continue
        ax.add_patch(Rectangle((0, 0), clen, 1, facecolor="white", edgecolor="black", linewidth=0.5))

        # Title/label at left
        ax.text(-chr1_len * 0.002, 0.5, chrom, va="center", ha="right", fontsize=8)

        # Draw colored blocks for this chromosome
        for start, end, mchr in by_chrom.get(chrom, []):
            if start >= end:
                continue
            color = mcolors.get(mchr)
            if color is None:
                continue
            # Clip to chromosome length just in case
            s = max(0, min(start, clen))
            e = max(0, min(end, clen))
            if e <= s:
                continue
            ax.add_patch(Rectangle((s, 0), e - s + 1, 1, facecolor=color, edgecolor="none"))

    # Build legend handles/labels
    handles = []
    labels = []
    for i in range(1, 20):
        lab = f"chr{i}"
        patch = Rectangle((0, 0), 1, 1, facecolor=mcolors[lab], edgecolor="black", linewidth=0.2)
        handles.append(patch)
        labels.append(lab)

    # Reserve bottom space and add legend in 3 rows: 7, 7, 5 items
    # Try tight_layout to reduce overlaps while reserving space for legend.
    fig.tight_layout(rect=(0.05, 0.25, 0.995, 0.995))

    # Row 1: chr1..chr7
    fig.legend(
        handles[:7],
        labels[:7],
        loc="lower center",
        ncol=7,
        frameon=False,
        bbox_to_anchor=(0.5, 0.10),
        fontsize=7,
    )
    # Row 2: chr8..chr14
    fig.legend(
        handles[7:14],
        labels[7:14],
        loc="lower center",
        ncol=7,
        frameon=False,
        bbox_to_anchor=(0.5, 0.06),
        fontsize=7,
    )
    # Row 3: chr15..chr19
    fig.legend(
        handles[14:],
        labels[14:],
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.02),
        fontsize=7,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf")
    plt.close(fig)


def main() -> int:
    if not BED12_PATH.exists():
        print(f"Missing BED12: {BED12_PATH}", file=sys.stderr)
        return 1
    if not CHROMSIZES.exists():
        print(f"Missing chrom sizes: {CHROMSIZES}", file=sys.stderr)
        return 1

    sizes = read_chrom_sizes(CHROMSIZES)
    by_chrom = parse_blocks_by_chrom(BED12_PATH)
    out = FIG_DIR / "hg38_autosomes_mouse_overlap_bars_10kfilter.pdf"
    plot_bars(by_chrom, sizes, out)
    print(f"Saved: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
