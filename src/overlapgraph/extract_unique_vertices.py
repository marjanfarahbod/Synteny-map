from __future__ import annotations

import argparse
import csv
from pathlib import Path


def vertex_type(vertex_name: str) -> int:
    if vertex_name.startswith("hh_"):
        return 1
    if vertex_name.startswith("mm_"):
        return 2
    raise ValueError(f"Unsupported vertex prefix for '{vertex_name}'")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract unique vertices from an hg/mm edge TSV and assign species codes."
    )
    parser.add_argument("input_tsv", type=Path, help="Input TSV with hg_index and mm_index columns")
    parser.add_argument("output_tsv", type=Path, help="Output TSV with unique vertices and type code")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    vertices: set[str] = set()

    with args.input_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "hg_index" not in fieldnames or "mm_index" not in fieldnames:
            raise ValueError("Input TSV must contain 'hg_index' and 'mm_index' columns")

        for row in reader:
            vertices.add(row["hg_index"])
            vertices.add(row["mm_index"])

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["vertex", "value"])
        for vertex in sorted(vertices, key=lambda name: (vertex_type(name), name)):
            writer.writerow([vertex, vertex_type(vertex)])


if __name__ == "__main__":
    main()
