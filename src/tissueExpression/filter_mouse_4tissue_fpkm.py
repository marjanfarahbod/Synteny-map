from __future__ import annotations

# Purpose: Filter the mouse FPKM matrix to ortholog genes and keep only the four tissue columns (mLi, mLu, mKi, mMu).
# Inputs:
#   - --expression: data/external/ST6_Lietal.tsv (default)
#   - --ortholog-map: data/processed/mouse_human_ortholog_one2one.tsv (default)
# Outputs:
#   - --output: data/processed/mouse_orthologs_4Tissue_FPKM.tsv

import argparse
import csv
from pathlib import Path


KEYWORDS = ("mLi", "mLu", "mKi", "mMu")


def load_mouse_gene_names(map_path: Path) -> set[str]:
    """Collect mouse gene names from the ortholog map."""
    with map_path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None or "Mouse gene name" not in reader.fieldnames:
            raise KeyError("Column 'Mouse gene name' not found in ortholog map file")
        return {row["Mouse gene name"] for row in reader if row.get("Mouse gene name")}


def filter_expression(expr_path: Path, gene_names: set[str], output_path: Path) -> tuple[int, int]:
    """Filter expression rows to mouse ortholog genes and select four-tissue columns."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with expr_path.open("r", newline="", encoding="utf-8") as infile:
        # First line is descriptive text; second is header
        meta_line = infile.readline()
        header_line = infile.readline()
        if not header_line:
            raise ValueError("Expression file is missing a header line")
        header = header_line.rstrip("\n").split("\t")
        if len(header) < 2:
            raise ValueError("Expression file header must have at least two columns")

        selected_indices = [0, 1] + [
            idx for idx, col in enumerate(header) if any(key in col for key in KEYWORDS)
        ]
        selected_header = [header[idx] for idx in selected_indices]

        reader = csv.reader(infile, delimiter="\t")
        kept_rows = []
        for row in reader:
            if len(row) < 2:
                continue
            gene_symbol = row[1]
            if gene_symbol not in gene_names:
                continue
            kept_rows.append([row[idx] for idx in selected_indices])

    with output_path.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(selected_header)
        writer.writerows(kept_rows)

    return len(kept_rows), len(selected_header)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter mouse FPKM table to genes with mouse-human orthologs and four tissue columns."
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/external/ST6_Lietal.tsv"),
        help="Input expression TSV file (default: data/external/ST6_Lietal.tsv)",
    )
    parser.add_argument(
        "--ortholog-map",
        type=Path,
        default=Path("data/processed/mouse_human_ortholog_one2one.tsv"),
        help="Ortholog map TSV file (default: data/processed/mouse_human_ortholog_one2one.tsv)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/mouse_orthologs_4Tissue_FPKM.tsv"),
        help="Output TSV path (default: data/processed/mouse_orthologs_4Tissue_FPKM.tsv)",
    )
    args = parser.parse_args()

    gene_names = load_mouse_gene_names(args.ortholog_map)
    rows, cols = filter_expression(args.expression, gene_names, args.output)
    print(f"Wrote {rows} rows and {cols} columns to {args.output}")


if __name__ == "__main__":
    main()
