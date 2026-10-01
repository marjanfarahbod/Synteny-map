from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rename bipartite edge TSV columns from human_index/mouse_index to hg_index/mm_index."
    )
    parser.add_argument("input_tsv", type=Path, help="Input edge TSV file")
    parser.add_argument("output_tsv", type=Path, help="Output TSV file with renamed columns")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    with args.input_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "human_index" not in fieldnames or "mouse_index" not in fieldnames:
            raise ValueError("Input TSV must contain 'human_index' and 'mouse_index' columns")

        args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_tsv.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(
                out_handle,
                fieldnames=["hg_index", "mm_index"],
                delimiter="\t",
            )
            writer.writeheader()
            for row in reader:
                writer.writerow(
                    {
                        "hg_index": f"hh_{row['human_index']}",
                        "mm_index": f"mm_{row['mouse_index']}",
                    }
                )


if __name__ == "__main__":
    main()
