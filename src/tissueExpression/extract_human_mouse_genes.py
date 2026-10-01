from __future__ import annotations

# Purpose: Extract Ensembl ID (col1) and gene symbol (col2) from the mapped/ordered GTEx file and
# store them as a pickle list of tuples for downstream use.
# Inputs:
#   - --input: data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped_ordered.tsv (default, GCT-style)
# Outputs:
#   - --output: data/processed/human15370genes_mouseOrtholog.pkl containing (Ensembl_ID, Gene_symbol) tuples

import argparse
import csv
import pickle
from pathlib import Path


def extract_genes(expr_path: Path) -> list[tuple[str, str]]:
    genes: list[tuple[str, str]] = []
    with expr_path.open("r", newline="", encoding="utf-8") as infile:
        # GCT-style: skip version and dims
        infile.readline()
        infile.readline()
        header = infile.readline()
        if not header:
            raise ValueError("Expression file missing header")
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if len(row) >= 2:
                genes.append((row[0], row[1]))
    return genes


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract Ensembl ID (col1) and gene symbol (col2) from mapped GTEx file."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped_ordered.tsv"),
        help="Input GTEx mapped ordered TSV (GCT-style).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/human15370genes_mouseOrtholog.pkl"),
        help="Output pickle path.",
    )
    args = parser.parse_args()

    genes = extract_genes(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as handle:
        pickle.dump(genes, handle)
    print(f"Wrote {len(genes)} gene entries to {args.output}")


if __name__ == "__main__":
    main()
