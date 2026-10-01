from __future__ import annotations

# Purpose: Normalize numeric columns so each column sums to 1e6, preserving any GCT-style prefix lines.
# Inputs:
#   - --input: path to TSV/GCT table (default: data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped_ordered.tsv)
# Outputs:
#   - --output: normalized table with `_columnSumNorm` suffix by default
# Notes:
#   - Only columns from index 2 onward are scaled; first two columns are treated as identifiers.

import argparse
import csv
from pathlib import Path
from typing import List, Tuple

SCALE = 1_000_000.0


def read_table(path: Path) -> Tuple[List[str], List[str], List[List[str]]]:
    """Return (prefix_lines, header, rows). prefix_lines holds GCT meta if present."""
    prefix: List[str] = []
    rows: List[List[str]] = []
    with path.open("r", newline="", encoding="utf-8") as infile:
        first = infile.readline()
        if not first:
            raise ValueError(f"{path} is empty")
        if first.startswith("#1."):
            # GCT-style: version, dims, header
            dims = infile.readline()
            header_line = infile.readline()
            if not header_line:
                raise ValueError(f"{path} missing header line")
            prefix = [first, dims]
            header = header_line.rstrip("\n").split("\t")
        else:
            header = first.rstrip("\n").split("\t")
        reader = csv.reader(infile, delimiter="\t")
        rows = [row for row in reader if row]
    return prefix, header, rows


def normalize_columns(rows: List[List[str]]) -> List[List[str]]:
    """Normalize numeric columns (index >=2) so each column sums to SCALE."""
    if not rows:
        return rows
    max_cols = max(len(r) for r in rows)
    col_sums = [0.0] * max_cols
    for row in rows:
        for idx in range(2, min(len(row), max_cols)):
            try:
                col_sums[idx] += float(row[idx])
            except (ValueError, TypeError):
                continue
    normalized: List[List[str]] = []
    for row in rows:
        out_row = list(row)
        for idx in range(2, min(len(row), max_cols)):
            try:
                val = float(row[idx])
            except (ValueError, TypeError):
                norm = 0.0
            else:
                denom = col_sums[idx]
                norm = (val / denom) * SCALE if denom else 0.0
            out_row[idx] = f"{norm:g}"
        normalized.append(out_row)
    return normalized


def write_table(
    path: Path,
    prefix: List[str],
    header: List[str],
    rows: List[List[str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as outfile:
        if prefix:
            outfile.write(prefix[0])
            dims_parts = prefix[1].rstrip("\n").split("\t")
            if dims_parts and dims_parts[0].isdigit():
                dims_parts[0] = str(len(rows))
            outfile.write("\t".join(dims_parts) + "\n")
        outfile.write("\t".join(header) + "\n")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Column-sum normalize expression table (starting at column 3) to 1e6 per column."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input TSV file (GCT-style or plain TSV).",
    )
    args = parser.parse_args()

    prefix, header, rows = read_table(args.input)
    normalized_rows = normalize_columns(rows)
    output_path = args.input.with_name(args.input.stem + "_columnSumNorm.tsv")
    write_table(output_path, prefix, header, normalized_rows)

    print(f"Wrote normalized file: {output_path}")
    print(f"Rows: {len(normalized_rows)}, Columns: {len(header)}")


if __name__ == "__main__":
    main()
