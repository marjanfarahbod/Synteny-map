from __future__ import annotations

# Purpose: Filter GTEx expression rows to genes present in the human/mouse ortholog map.
# Inputs:
#   - --expression: data/processed/GTEx_4Tissues.txt (GCT-style)
#   - --ortholog-map: data/processed/mouse_human_ortholog_one2one.tsv (default)
# Outputs:
#   - --output: data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv (filtered expression table)

import argparse
import csv
from pathlib import Path


def load_gene_names(map_path: Path) -> set[str]:
    """Load human gene names from the last column of the ortholog map."""
    with map_path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError("Ortholog map file is empty")
        target_idx = len(header) - 1
        names = {
            row[target_idx]
            for row in reader
            if len(row) > target_idx and row[target_idx].strip()
        }
    return names


def filter_expression(expr_path: Path, map_names: set[str], output_path: Path) -> tuple[int, int]:
    """Filter expression rows to those whose gene symbol (col 2) is in map_names."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with expr_path.open("r", newline="", encoding="utf-8") as infile:
        version_line = infile.readline()
        if not version_line:
            raise ValueError("Expression file missing version line")
        dims_line = infile.readline()
        if not dims_line:
            raise ValueError("Expression file missing dims line")
        header_line = infile.readline()
        if not header_line:
            raise ValueError("Expression file missing header line")

        header_fields = header_line.rstrip("\n").split("\t")
        if len(header_fields) < 2:
            raise ValueError("Expression file header must contain at least two columns")

        kept_rows: list[list[str]] = []
        total = 0
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            total += 1
            if row[1] not in map_names:
                continue
            kept_rows.append(row)

    dims_parts = dims_line.rstrip("\n").split("\t")
    if dims_parts:
        try:
            dims_parts[0] = str(len(kept_rows))
        except ValueError:
            pass
    new_dims_line = "\t".join(dims_parts) + "\n"

    with output_path.open("w", newline="", encoding="utf-8") as outfile:
        outfile.write(version_line)
        outfile.write(new_dims_line)
        outfile.write("\t".join(header_fields) + "\n")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in kept_rows:
            writer.writerow(row)

    return len(kept_rows), total


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter GTEx expression rows to those with human gene names present in the ortholog map."
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/processed/GTEx_4Tissues.txt"),
        help="Input expression GCT-like file (default: data/processed/GTEx_4Tissues.txt)",
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
        default=Path("data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv"),
        help="Output filtered expression file (default: data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv)",
    )
    args = parser.parse_args()

    gene_names = load_gene_names(args.ortholog_map)
    kept, total = filter_expression(args.expression, gene_names, args.output)
    print(f"Filtered expression rows: kept {kept} of {total} (output: {args.output})")


if __name__ == "__main__":
    main()
