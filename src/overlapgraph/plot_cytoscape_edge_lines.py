from __future__ import annotations

# Purpose: Plot bipartite graph edges as lines between human and mouse absolute bin indices.
# Inputs:
#   - --edges: Cytoscape edge TSV with human_node, mouse_node, and component columns
# Outputs:
#   - --output: PDF line plot
#
# Notes:
#   - Human nodes are drawn at y=0 and mouse nodes at y=100.
#   - Edges from the same component share the same color from the Paired colormap.

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def node_index(node_name: str) -> int:
    return int(node_name.split("_", 1)[1])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Cytoscape edge table as human-to-mouse component lines.")
    parser.add_argument(
        "--edges",
        type=Path,
        default=Path("results/tables/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_cytoscape_edges.tsv"),
        help="Input Cytoscape edge TSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/figures/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_edge_lines.pdf"),
        help="Output PDF path",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    edges: list[tuple[int, int, int]] = []

    with args.edges.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            edges.append(
                (
                    node_index(row["human_node"]),
                    node_index(row["mouse_node"]),
                    int(row["component"]),
                )
            )

    if not edges:
        raise ValueError(f"No edges found in {args.edges}")

    components = sorted({component for _, _, component in edges})
    cmap = plt.get_cmap("Paired", len(components))
    component_colors = {component: cmap(i) for i, component in enumerate(components)}

    fig, ax = plt.subplots(figsize=(12, 6))
    for human_index, mouse_index, component in edges:
        ax.plot(
            [human_index, mouse_index],
            [0, 100],
            color=component_colors[component],
            linewidth=0.8,
            alpha=0.85,
        )

    ax.set_xlim(0, 490)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Absolute bin index")
    ax.set_ylabel("Species axis")
    ax.set_yticks([0, 100])
    ax.set_yticklabels(["human", "mouse"])
    ax.set_title("Human-mouse bin edges by connected component")
    ax.grid(axis="x", color="#dddddd", linewidth=0.4)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output)
    plt.close(fig)
    print(f"edges\t{len(edges)}")
    print(f"components\t{len(components)}")
    print(f"plot\t{args.output}")


if __name__ == "__main__":
    main()
