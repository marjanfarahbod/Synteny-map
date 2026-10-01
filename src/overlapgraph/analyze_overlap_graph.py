from __future__ import annotations

# Purpose: Analyze an overlap-matrix pickle by generating a PDF heatmap and reporting connected components in the induced bipartite graph.
# Inputs:
#   - --matrix: pickle file containing a human x mouse overlap matrix
#   - --hg-start: start bin index used for the human rows
#   - --mm-start: start bin index used for the mouse columns
#   - --plot: output PDF path for the heatmap
# Outputs:
#   - --plot: PDF heatmap of the matrix
#   - Printed connected components with more than one vertex to stdout
#
# Notes:
#   - Each matrix cell greater than --edge-threshold is treated as an edge in the bipartite graph.
#   - Heatmap values are capped at 1 Mb, transformed with log10, and zero cells are masked.

import argparse
import pickle
from collections import deque
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot_heatmap(matrix: np.ndarray, out_path: Path) -> None:
    clipped = np.clip(matrix, 0, 1_000_000)
    masked = np.ma.masked_less_equal(clipped, 0)
    log_matrix = np.ma.log10(masked)

    plt.figure(figsize=(8, 6))
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color="white")
    im = plt.imshow(log_matrix, aspect="auto", cmap=cmap, vmin=0, vmax=6)
    plt.colorbar(im, label="log10 overlap length (capped at 1 Mb)")
    plt.xlabel("Mouse bins (1 Mb)")
    plt.ylabel("Human bins (1 Mb)")
    plt.title("Overlap graph matrix (log10, capped at 1 Mb)")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def connected_components(matrix: np.ndarray, edge_threshold: int) -> list[tuple[list[int], list[int]]]:
    human_to_mouse = [set(np.where(matrix[i] > edge_threshold)[0].tolist()) for i in range(matrix.shape[0])]
    mouse_to_human = [set(np.where(matrix[:, j] > edge_threshold)[0].tolist()) for j in range(matrix.shape[1])]

    visited_human = [False] * matrix.shape[0]
    visited_mouse = [False] * matrix.shape[1]
    components: list[tuple[list[int], list[int]]] = []

    for start_human in range(matrix.shape[0]):
        if visited_human[start_human] or not human_to_mouse[start_human]:
            continue

        queue = deque([("human", start_human)])
        visited_human[start_human] = True
        human_nodes: list[int] = []
        mouse_nodes: list[int] = []

        while queue:
            kind, index = queue.popleft()
            if kind == "human":
                human_nodes.append(index)
                for mouse_idx in human_to_mouse[index]:
                    if not visited_mouse[mouse_idx]:
                        visited_mouse[mouse_idx] = True
                        queue.append(("mouse", mouse_idx))
            else:
                mouse_nodes.append(index)
                for human_idx in mouse_to_human[index]:
                    if not visited_human[human_idx]:
                        visited_human[human_idx] = True
                        queue.append(("human", human_idx))

        if len(human_nodes) + len(mouse_nodes) > 1:
            components.append((sorted(human_nodes), sorted(mouse_nodes)))

    return components


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate an overlap heatmap and report connected components from a matrix pickle."
    )
    parser.add_argument("--matrix", type=Path, required=True, help="Pickled human x mouse overlap matrix")
    parser.add_argument("--hg-start", type=int, default=0, help="Human start bin used to generate the matrix")
    parser.add_argument("--mm-start", type=int, default=0, help="Mouse start bin used to generate the matrix")
    parser.add_argument("--plot", type=Path, required=True, help="Output PDF path for the heatmap")
    parser.add_argument("--edge-threshold", type=int, default=10_000, help="Minimum overlap required to treat a cell as a graph edge")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.matrix.open("rb") as handle:
        matrix = pickle.load(handle)

    plot_heatmap(matrix, args.plot)
    components = connected_components(matrix, args.edge_threshold)

    print(f"matrix\t{args.matrix}")
    print(f"plot\t{args.plot}")
    print(f"edge_threshold\t{args.edge_threshold}")
    print(f"component_count\t{len(components)}")
    for component_id, (human_nodes, mouse_nodes) in enumerate(components, start=1):
        human_absolute = ",".join(str(args.hg_start + idx) for idx in human_nodes)
        mouse_absolute = ",".join(str(args.mm_start + idx) for idx in mouse_nodes)
        print(f"{component_id}:  human {human_absolute}    mouse {mouse_absolute}")


if __name__ == "__main__":
    main()
