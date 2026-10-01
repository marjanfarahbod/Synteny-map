from __future__ import annotations

# Purpose: Threshold a human x mouse overlap matrix into a bipartite graph and report its connected components.
# Inputs:
#   - --matrix: pickled overlap matrix - data/processed/hg38_mm39_overlap_1Mb_matrix_490bins.pkl
#   - --hg-chrom: data/external/hg38.chrom.sizes
#   - --mm-chrom: data/external/mm39.chrom.sizes
#   - --hg-start: start bin index for human
#   - --mm-start: start bin index for mouse
#   - --regions: number of 1 Mb bins per species (must match the matrix shape)
#   - --threshold: minimum overlap in bp for an edge (strictly greater than; default 10,000)
# Outputs:
#   - --table-out: TSV of components with at least one edge - results/tables/hg38_mm39_overlap_1Mb_matrix_490bins_components_gt10k.tsv
#   - --edges-out: TSV edge list (human_index, mouse_index) - results/tables/hg38_mm39_overlap_1Mb_matrix_490bins_edges_gt10k.tsv
#   - stdout: component counts and a per-component summary (including isolated bins)
#
# Notes:
#   - Bins are rebuilt with generate_bins from compute_overlap_matrix_pickle.py.
#   - Components are found by BFS; edges are counted once from the human side.
#   - Component IDs follow discovery order; rows are sorted by edges > 0, then total human + mouse bp (descending).

import argparse
import pickle
import sys
from collections import deque
from pathlib import Path
import csv

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.overlapgraph.compute_overlap_matrix_pickle import Bin, generate_bins


def build_adjacency(matrix: np.ndarray, threshold: int) -> tuple[list[set[int]], list[set[int]]]:
    human_to_mouse = [set() for _ in range(matrix.shape[0])]
    mouse_to_human = [set() for _ in range(matrix.shape[1])]

    for human_idx in range(matrix.shape[0]):
        mouse_indices = np.flatnonzero(matrix[human_idx] > threshold)
        for mouse_idx in mouse_indices:
            human_to_mouse[human_idx].add(int(mouse_idx))
            mouse_to_human[int(mouse_idx)].add(human_idx)

    return human_to_mouse, mouse_to_human


def find_components(
    human_to_mouse: list[set[int]],
    mouse_to_human: list[set[int]],
) -> list[tuple[list[int], list[int], int]]:
    visited_human = [False] * len(human_to_mouse)
    visited_mouse = [False] * len(mouse_to_human)
    components: list[tuple[list[int], list[int], int]] = []

    for start_human in range(len(human_to_mouse)):
        if visited_human[start_human]:
            continue

        queue = deque([("human", start_human)])
        visited_human[start_human] = True
        human_nodes: list[int] = []
        mouse_nodes: list[int] = []
        edge_count = 0

        while queue:
            kind, index = queue.popleft()
            if kind == "human":
                human_nodes.append(index)
                neighbors = human_to_mouse[index]
                edge_count += len(neighbors)
                for mouse_idx in neighbors:
                    if not visited_mouse[mouse_idx]:
                        visited_mouse[mouse_idx] = True
                        queue.append(("mouse", mouse_idx))
            else:
                mouse_nodes.append(index)
                for human_idx in mouse_to_human[index]:
                    if not visited_human[human_idx]:
                        visited_human[human_idx] = True
                        queue.append(("human", human_idx))

        components.append((human_nodes, mouse_nodes, edge_count))

    for start_mouse in range(len(mouse_to_human)):
        if visited_mouse[start_mouse]:
            continue
        visited_mouse[start_mouse] = True
        components.append(([], [start_mouse], 0))

    return components


def component_bp(bin_indices: list[int], bins: list[Bin]) -> int:
    return sum(bins[idx][2] - bins[idx][1] for idx in bin_indices)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze connected components in a thresholded bipartite overlap matrix."
    )
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("data/processed/hg38_mm39_overlap_1Mb_matrix_490bins.pkl"),
        help="Pickled human x mouse overlap matrix",
    )
    parser.add_argument("--hg-chrom", type=Path, default=Path("data/external/hg38.chrom.sizes"))
    parser.add_argument("--mm-chrom", type=Path, default=Path("data/external/mm39.chrom.sizes"))
    parser.add_argument("--hg-start", type=int, default=0, help="Start bin for human")
    parser.add_argument("--mm-start", type=int, default=0, help="Start bin for mouse")
    parser.add_argument("--regions", type=int, default=490, help="Number of bins per species")
    parser.add_argument("--threshold", type=int, default=10_000, help="Edge threshold in bp")
    parser.add_argument(
        "--table-out",
        type=Path,
        default=Path("results/tables/hg38_mm39_overlap_1Mb_matrix_490bins_components_gt10k.tsv"),
        help="Write components with at least one edge to this TSV file",
    )
    parser.add_argument(
        "--edges-out",
        type=Path,
        default=Path("results/tables/hg38_mm39_overlap_1Mb_matrix_490bins_edges_gt10k.tsv"),
        help="Write the sparse bipartite adjacency list to this TSV file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.matrix.open("rb") as handle:
        matrix = pickle.load(handle)

    human_bins = generate_bins(args.hg_chrom, start_bins=args.hg_start, num_bins=args.regions)
    mouse_bins = generate_bins(args.mm_chrom, start_bins=args.mm_start, num_bins=args.regions)

    expected_shape = (len(human_bins), len(mouse_bins))
    if matrix.shape != expected_shape:
        raise ValueError(f"Matrix shape {matrix.shape} does not match expected shape {expected_shape}")

    human_to_mouse, mouse_to_human = build_adjacency(matrix, args.threshold)
    components = find_components(human_to_mouse, mouse_to_human)

    rows = []
    for component_id, (human_nodes, mouse_nodes, edge_count) in enumerate(components, start=1):
        rows.append(
            {
                "component": component_id,
                "human_bins": len(human_nodes),
                "mouse_bins": len(mouse_nodes),
                "human_bp": component_bp(human_nodes, human_bins),
                "mouse_bp": component_bp(mouse_nodes, mouse_bins),
                "edges": edge_count,
                "human_indices": sorted(human_nodes),
                "mouse_indices": sorted(mouse_nodes),
            }
        )

    rows.sort(
        key=lambda row: (row["edges"] > 0, row["human_bp"] + row["mouse_bp"], row["component"]),
        reverse=True,
    )

    edge_components = sum(1 for row in rows if row["edges"] > 0)
    isolated_components = len(rows) - edge_components

    args.table_out.parent.mkdir(parents=True, exist_ok=True)
    with args.table_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "component",
                "number_of_mouse_bins",
                "number_of_human_bins",
                "indices_of_mouse_bins",
                "indices_of_human_bins",
                "total_human_bp",
                "total_mouse_bp",
                "edges",
            ]
        )
        for row in rows:
            if row["edges"] <= 0:
                continue
            mouse_indices = ",".join(map(str, row["mouse_indices"]))
            human_indices = ",".join(map(str, row["human_indices"]))
            writer.writerow(
                [
                    row["component"],
                    row["mouse_bins"],
                    row["human_bins"],
                    mouse_indices,
                    human_indices,
                    row["human_bp"],
                    row["mouse_bp"],
                    row["edges"],
                ]
            )

    args.edges_out.parent.mkdir(parents=True, exist_ok=True)
    with args.edges_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["human_index", "mouse_index"])
        for human_idx, mouse_neighbors in enumerate(human_to_mouse):
            for mouse_idx in sorted(mouse_neighbors):
                writer.writerow([human_idx, mouse_idx])

    print(f"total_components\t{len(rows)}")
    print(f"components_with_edges\t{edge_components}")
    print(f"isolated_components\t{isolated_components}")
    print(f"table_out\t{args.table_out}")
    print(f"edges_out\t{args.edges_out}")
    print("component\thuman_bins\tmouse_bins\thuman_bp\tmouse_bp\tedges\thuman_indices\tmouse_indices")
    for row in rows:
        print(
            f"{row['component']}\t{row['human_bins']}\t{row['mouse_bins']}\t"
            f"{row['human_bp']}\t{row['mouse_bp']}\t{row['edges']}\t"
            f"{','.join(map(str, row['human_indices']))}\t"
            f"{','.join(map(str, row['mouse_indices']))}"
        )


if __name__ == "__main__":
    main()
