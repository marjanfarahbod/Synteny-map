from __future__ import annotations

# Purpose: Export connected-component edge and node tables for Cytoscape from an overlap-matrix pickle.
# Inputs:
#   - --matrix: pickle file containing a human x mouse overlap matrix
#   - --edge-threshold: minimum overlap required to treat a cell as a graph edge
# Outputs:
#   - edge table TSV with hg/mm node names and component ids
#   - node table TSV with node name, absolute value, component id, and species code
#
# Notes:
#   - Human and mouse absolute indices are inferred from the pickle filename.
#   - Only nodes and edges that belong to connected components with more than one vertex are exported.

import argparse
import csv
import pickle
import re
from collections import deque
from pathlib import Path

import numpy as np


def infer_starts(matrix_path: Path) -> tuple[int, int]:
    name = matrix_path.stem

    match = re.search(r"hg(\d+)_(\d+)_mm(\d+)_(\d+)", name)
    if match:
        return int(match.group(1)), int(match.group(3))

    match = re.search(r"hg(\d+)_mm(\d+)_", name)
    if match:
        return int(match.group(1)), int(match.group(2))

    raise ValueError(f"Could not infer hg/mm start bins from filename: {matrix_path.name}")


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
        description="Export Cytoscape node and edge tables from an overlap-matrix pickle."
    )
    parser.add_argument("--matrix", type=Path, required=True, help="Pickled human x mouse overlap matrix")
    parser.add_argument("--edge-threshold", type=int, default=10_000, help="Minimum overlap required to treat a cell as a graph edge")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.matrix.open("rb") as handle:
        matrix = pickle.load(handle)

    hg_start, mm_start = infer_starts(args.matrix)
    components = connected_components(matrix, args.edge_threshold)

    base = args.matrix.stem
    edges_path = Path("results/tables") / f"{base}_cytoscape_edges.tsv"
    nodes_path = Path("results/tables") / f"{base}_cytoscape_nodes.tsv"
    edges_path.parent.mkdir(parents=True, exist_ok=True)

    with edges_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["human_node", "mouse_node", "component", "overlap_value"])
        for component_id, (human_nodes, mouse_nodes) in enumerate(components, start=1):
            for human_idx in human_nodes:
                neighbors = np.where(matrix[human_idx] > args.edge_threshold)[0]
                for mouse_idx in neighbors:
                    if mouse_idx in mouse_nodes:
                        writer.writerow(
                            [
                                f"hg_{hg_start + human_idx}",
                                f"mm_{mm_start + mouse_idx}",
                                component_id,
                                int(matrix[human_idx, mouse_idx]),
                            ]
                        )

    with nodes_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["node", "absolute_value", "component", "species_code"])
        for component_id, (human_nodes, mouse_nodes) in enumerate(components, start=1):
            for human_idx in human_nodes:
                writer.writerow(
                    [
                        f"hg_{hg_start + human_idx}",
                        hg_start + human_idx,
                        component_id,
                        1,
                    ]
                )
            for mouse_idx in mouse_nodes:
                writer.writerow(
                    [
                        f"mm_{mm_start + mouse_idx}",
                        mm_start + mouse_idx,
                        component_id,
                        2,
                    ]
                )

    print(f"matrix\t{args.matrix}")
    print(f"edge_threshold\t{args.edge_threshold}")
    print(f"hg_start\t{hg_start}")
    print(f"mm_start\t{mm_start}")
    print(f"components_exported\t{len(components)}")
    print(f"edges_out\t{edges_path}")
    print(f"nodes_out\t{nodes_path}")


if __name__ == "__main__":
    main()
