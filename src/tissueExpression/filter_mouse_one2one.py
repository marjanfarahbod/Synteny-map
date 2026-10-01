from __future__ import annotations

# Purpose: Filter the mouse/human ortholog table to retain unique one-to-one ortholog pairs.
# Inputs:
#   - data/external/mouse_human_ortholog_ensemble.txt
# Outputs:
#   - data/processed/mouse_human_ortholog_one2one.tsv (unique Mouse gene name / Gene name pairs with homology type ortholog_one2one)

import csv
from pathlib import Path

def main() -> None:
    input_path = Path("data/external/mouse_human_ortholog_ensemble.txt")
    output_path = Path("data/processed/mouse_human_ortholog_one2one.tsv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", newline="", encoding="utf-8") as infile, output_path.open(
        "w", newline="", encoding="utf-8"
    ) as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input file is missing a header row")
        if "Mouse homology type" not in reader.fieldnames:
            raise KeyError("Column 'Mouse homology type' not found in input file")
        if "Mouse gene name" not in reader.fieldnames:
            raise KeyError("Column 'Mouse gene name' not found in input file")
        if "Gene name" not in reader.fieldnames:
            raise KeyError("Column 'Gene name' not found in input file")

        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        seen_pairs = set()
        kept = 0
        for row in reader:
            if row.get("Mouse homology type") != "ortholog_one2one":
                continue
            pair = (row.get("Mouse gene name", ""), row.get("Gene name", ""))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            writer.writerow(row)
            kept += 1

    print(f"Wrote {kept} rows to {output_path}")

if __name__ == "__main__":
    main()
