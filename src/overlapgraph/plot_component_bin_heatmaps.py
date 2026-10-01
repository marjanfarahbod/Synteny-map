from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np


def parse_index_list(value: str) -> list[int]:
    if not value.strip():
        return []
    return [int(part) for part in value.split(",") if part]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot one 2xN occupancy heatmap per connected component from a component TSV."
    )
    parser.add_argument(
        "--components",
        type=Path,
        default=Path("results/tables/hg38_mm39_overlap_1Mb_matrix_490bins_components_gt10k.tsv"),
        help="Component TSV file",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=490,
        help="Number of bins per species",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/figures/hg38_mm39_overlap_1Mb_matrix_490bins_component_bin_heatmaps.pdf"),
        help="Output figure path",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    components: list[dict[str, object]] = []

    with args.components.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            human_indices = parse_index_list(row["indices_of_human_bins"])
            mouse_indices = parse_index_list(row["indices_of_mouse_bins"])
            components.append(
                {
                    "component": row["component"],
                    "human_indices": human_indices,
                    "mouse_indices": mouse_indices,
                }
            )

    if not components:
        raise ValueError(f"No components found in {args.components}")

    n_components = len(components)
    ncols = min(3, n_components)
    nrows = math.ceil(n_components / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 1.8 * nrows), squeeze=False)
    cmap = ListedColormap(["white", "#2166ac"])

    for ax, component in zip(axes.flat, components):
        data = np.zeros((2, args.bins), dtype=np.int8)
        data[0, component["human_indices"]] = 1
        data[1, component["mouse_indices"]] = 1

        ax.imshow(data, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
        ax.set_title(f"Component {component['component']}")
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["human", "mouse"])
        ax.set_xlabel("Bin index")
        ax.set_xticks(np.arange(-0.5, args.bins, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, 2, 1), minor=True)
        ax.grid(which="minor", color="#bdbdbd", linestyle="-", linewidth=0.2)
        ax.tick_params(which="minor", bottom=False, left=False)

    for ax in axes.flat[n_components:]:
        ax.axis("off")

    fig.suptitle("Component Bin Occupancy Heatmaps", y=0.995)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print(f"Saved {n_components} component heatmaps to {args.output}")


if __name__ == "__main__":
    main()
