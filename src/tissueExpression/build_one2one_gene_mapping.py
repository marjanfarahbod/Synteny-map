from __future__ import annotations

# Purpose: Build a one-to-one human–mouse gene name mapping using the processed ortholog file.
# Inputs:
#   - --orthologs: first matching data/processed/mouse_human_*one2one.tsv (default)
#   - --human GTEx file: data/processed/GTEx_*mapped.tsv)
# Outputs:
#   - data/interim/one2oneGeneNameMapping.tsv (two columns: human gene name, mouse gene name)
# Usage:
#   python3 src/tissueExpression/build_one2one_gene_mapping.py #       --orthologs data/processed/mouse_human_ortholog_one2one.tsv #       --human data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped.tsv #       --output data/interim/one2oneGeneNameMapping.tsv

import argparse
import csv
import glob
from pathlib import Path


def find_file(pattern: str) -> Path:
    matches = [Path(p) for p in glob.glob(pattern) if Path(p).is_file()]
    if not matches:
        raise FileNotFoundError(f"No file found for pattern: {pattern}")
    if len(matches) > 1:
        # pick first sorted for determinism
        matches.sort()
    return matches[0]


def load_ortholog_map(path: Path) -> dict[str, str]:
    """Map human gene name -> mouse gene name."""
    with path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Ortholog map missing header row")
        if "Gene name" not in reader.fieldnames or "Mouse gene name" not in reader.fieldnames:
            raise KeyError("Required columns 'Gene name' and 'Mouse gene name' not found")
        return {
            row["Gene name"]: row["Mouse gene name"]
            for row in reader
            if row.get("Gene name") and row.get("Mouse gene name")
        }


def load_human_genes(path: Path) -> list[str]:
    """Collect human gene names from column 2 (skip GCT header if present)."""
    genes: list[str] = []
    with path.open("r", newline="", encoding="utf-8") as infile:
        first = infile.readline()
        if first.startswith("#1."):
            infile.readline()  # dims
            header = infile.readline()
        else:
            header = first
        if not header:
            return genes
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                genes.append(row[1])
    return genes


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build one-to-one human->mouse gene name mapping from expression files and ortholog map."
    )
    parser.add_argument(
        "--human",
        type=Path,
        default=None,
        help="Human GTEx mapped file (default: first matching data/processed/GTEx_*mapped.tsv)",
    )
    parser.add_argument(
        "--ortholog-map",
        type=Path,
        default=None,
        help="Ortholog map file (default: first matching data/processed/mouse_human_*one2one.tsv)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/interim/one2oneGeneNameMapping.tsv"),
        help="Output mapping TSV (default: data/interim/one2oneGeneNameMapping.tsv)",
    )
    args = parser.parse_args()

    human_path = args.human or find_file("data/processed/GTEx_*mapped.tsv")
    ortholog_path = args.ortholog_map or find_file("data/processed/mouse_human_*one2one.tsv")

    ortholog_map = load_ortholog_map(ortholog_path)
    human_genes = load_human_genes(human_path)

    mapped_pairs = [
        (h, ortholog_map[h])
        for h in human_genes
        if h in ortholog_map
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerow(["Gene name", "Mouse gene name"])
        writer.writerows(mapped_pairs)

    print(f"Wrote {len(mapped_pairs)} mappings to {args.output}")
    print(f"Human file: {human_path}")
    print(f"Ortholog map: {ortholog_path}")


if __name__ == "__main__":
    main()
