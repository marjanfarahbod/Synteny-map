from __future__ import annotations

# Purpose: Map shared ortholog expression with pre-filtering of duplicates in column 2 before mapping,
# producing aligned human/mouse outputs matching column 2 value sets.
# Inputs:
#   - --gtex: data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv (default; GCT-style)
#   - --mouse: data/processed/mouse_orthologs_4Tissue_FPKM.tsv (default)
#   - --ortholog-map: data/processed/mouse_human_ortholog_one2one.tsv (default)
# Outputs:
#   - --gtex-out: data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped.tsv
#   - --mouse-out: data/processed/mouse_orthologs_4Tissue_FPKM_mapped.tsv
# Notes:
#   - Drops any gene rows with duplicate column-2 values before mapping.

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Iterable


def load_ortholog_map(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Return mappings human->mouse and mouse->human from the ortholog file."""
    with path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Ortholog map file missing header row")
        if "Gene name" not in reader.fieldnames:
            raise KeyError("Column 'Gene name' not found in ortholog map file")
        if "Mouse gene name" not in reader.fieldnames:
            raise KeyError("Column 'Mouse gene name' not found in ortholog map file")

        h_to_m: dict[str, str] = {}
        m_to_h: dict[str, str] = {}
        for row in reader:
            human = row.get("Gene name", "").strip()
            mouse = row.get("Mouse gene name", "").strip()
            if not human or not mouse:
                continue
            h_to_m[human] = mouse
            m_to_h[mouse] = human
    return h_to_m, m_to_h


def load_table(path: Path, has_gct_header: bool) -> tuple[list[str], list[list[str]], list[str]]:
    """Load table rows; returns header, rows, and prefix lines if GCT-like."""
    rows: list[list[str]] = []
    prefix: list[str] = []

    with path.open("r", newline="", encoding="utf-8") as infile:
        if has_gct_header:
            version = infile.readline()
            dims = infile.readline()
            header_line = infile.readline()
            if not header_line:
                raise ValueError(f"{path} missing header line")
            prefix = [version, dims]
        else:
            header_line = infile.readline()
            if not header_line:
                raise ValueError(f"{path} missing header line")
        header = header_line.rstrip("\n").split("\t")

        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if row:
                rows.append(row)

    return header, rows, prefix


def remove_duplicates_by_col2(rows: Iterable[list[str]]) -> list[list[str]]:
    """Drop all rows whose second-column value appears more than once."""
    counts: Counter[str] = Counter()
    for row in rows:
        if len(row) >= 2:
            counts[row[1]] += 1
    return [row for row in rows if len(row) >= 2 and counts[row[1]] == 1]


def filter_rows(
    header: list[str],
    rows: Iterable[list[str]],
    allowed_names: set[str],
    has_gct_header: bool,
    prefix_lines: list[str],
    output_path: Path,
) -> int:
    """Write filtered rows whose second column is in allowed_names."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    kept_rows = [row for row in rows if len(row) >= 2 and row[1] in allowed_names]

    with output_path.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        if has_gct_header:
            outfile.write(prefix_lines[0])
            dims_parts = prefix_lines[1].rstrip("\n").split("\t")
            if dims_parts and dims_parts[0].isdigit():
                dims_parts[0] = str(len(kept_rows))
            outfile.write("\t".join(dims_parts) + "\n")
        writer.writerow(header)
        writer.writerows(kept_rows)

    return len(kept_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter GTEx and mouse FPKM tables to shared ortholog genes after removing duplicates."
    )
    parser.add_argument(
        "--ortholog-map",
        type=Path,
        default=Path("data/processed/mouse_human_ortholog_one2one.tsv"),
        help="Ortholog map TSV (default: data/processed/mouse_human_ortholog_one2one.tsv)",
    )
    parser.add_argument(
        "--gtex",
        type=Path,
        default=Path("data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv"),
        help="GTEx expression file with GCT header (default: data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv)",
    )
    parser.add_argument(
        "--mouse",
        type=Path,
        default=Path("data/processed/mouse_orthologs_4Tissue_FPKM.tsv"),
        help="Mouse expression file (default: data/processed/mouse_orthologs_4Tissue_FPKM.tsv)",
    )
    args = parser.parse_args()

    h_to_m, m_to_h = load_ortholog_map(args.ortholog_map)

    gtex_header, gtex_rows, gtex_prefix = load_table(args.gtex, has_gct_header=True)
    mouse_header, mouse_rows, _ = load_table(args.mouse, has_gct_header=False)

    # Remove any duplicated names within each file
    gtex_unique_rows = remove_duplicates_by_col2(gtex_rows)
    mouse_unique_rows = remove_duplicates_by_col2(mouse_rows)

    gtex_names = {row[1] for row in gtex_unique_rows if len(row) >= 2}
    mouse_names = {row[1] for row in mouse_unique_rows if len(row) >= 2}

    # Candidate pairs based on presence
    candidate_pairs = [
        (h, h_to_m[h])
        for h in gtex_names
        if h in h_to_m and h_to_m[h] in mouse_names
    ]

    # Enforce one-to-one within the present sets
    human_counts = Counter(h for h, _ in candidate_pairs)
    mouse_counts = Counter(m for _, m in candidate_pairs)
    filtered_pairs = [
        (h, m) for h, m in candidate_pairs if human_counts[h] == 1 and mouse_counts[m] == 1
    ]

    shared_pairs = set(filtered_pairs)
    shared_humans_final = {h for h, _ in shared_pairs}
    shared_mice_final = {m for _, m in shared_pairs}

    gtex_out = args.gtex.with_name(args.gtex.stem + "_mapped.tsv")
    mouse_out = args.mouse.with_name(args.mouse.stem + "_mapped.tsv")

    kept_gtex = filter_rows(gtex_header, gtex_unique_rows, shared_humans_final, True, gtex_prefix, gtex_out)
    kept_mouse = filter_rows(mouse_header, mouse_unique_rows, shared_mice_final, False, [], mouse_out)

    print(f"Wrote GTEx mapped: {kept_gtex} rows -> {gtex_out}")
    print(f"Wrote mouse mapped: {kept_mouse} rows -> {mouse_out}")
    print(f"Shared ortholog gene pairs (after duplicate removal): {len(shared_pairs)}")


if __name__ == "__main__":
    main()
