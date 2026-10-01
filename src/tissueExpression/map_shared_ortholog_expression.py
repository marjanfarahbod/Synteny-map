from __future__ import annotations

# Purpose: Filter two expression tables to shared orthologs using an ortholog map, writing aligned outputs.
# Inputs:
#   - --gtex: data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv (default; GCT-style)
#   - --mouse: data/processed/mouse_orthologs_4Tissue_FPKM.tsv (default)
#   - --ortholog-map: data/processed/mouse_human_ortholog_one2one.tsv (default)
# Outputs:
#   - --gtex-out: data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped.tsv
#   - --mouse-out: data/processed/mouse_orthologs_4Tissue_FPKM_mapped.tsv

import argparse
import csv
from pathlib import Path
from typing import Iterable


def load_ortholog_map(path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Return mappings human->mouse and mouse->human from the ortholog file."""
    with path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Ortholog map file missing header row")
        if "Gene name" not in reader.fieldnames:
            raise KeyError("Column 'Gene name' not found in ortholog map file")
        if "Mouse gene name" not in reader.fieldnames:
            raise KeyError("Column 'Mouse gene name' not found in ortholog map file")

        h_to_m: dict[str, set[str]] = {}
        m_to_h: dict[str, set[str]] = {}
        for row in reader:
            human = row.get("Gene name", "").strip()
            mouse = row.get("Mouse gene name", "").strip()
            if not human or not mouse:
                continue
            h_to_m.setdefault(human, set()).add(mouse)
            m_to_h.setdefault(mouse, set()).add(human)
    return h_to_m, m_to_h


def load_names_from_table(
    path: Path, has_gct_header: bool
) -> tuple[list[str], list[list[str]], list[str]]:
    """Load header and rows; returns header, rows, and prefix lines if GCT-like."""
    rows: list[list[str]] = []
    header: list[str] = []
    with path.open("r", newline="", encoding="utf-8") as infile:
        if has_gct_header:
            # skip version + dims lines
            version = infile.readline()
            dims = infile.readline()
            header_line = infile.readline()
            if not header_line:
                raise ValueError(f"{path} missing header line")
            header = header_line.rstrip("\n").split("\t")
            prefix = [version, dims]
        else:
            prefix = []
            header_line = infile.readline()
            if not header_line:
                raise ValueError(f"{path} missing header line")
            header = header_line.rstrip("\n").split("\t")

        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if row:
                rows.append(row)

    return header, rows, prefix if has_gct_header else []


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
            # update dims line first entry to kept_rows count
            dims_parts = prefix_lines[1].rstrip("\n").split("\t")
            if dims_parts and dims_parts[0].isdigit():
                dims_parts[0] = str(len(kept_rows))
            outfile.write("\t".join(dims_parts) + "\n")
        writer.writerow(header)
        writer.writerows(kept_rows)

    return len(kept_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter GTEx and mouse FPKM tables to shared ortholog genes using ortholog map."
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

    gtex_header, gtex_rows, gtex_prefix = load_names_from_table(args.gtex, has_gct_header=True)
    mouse_header, mouse_rows, _ = load_names_from_table(args.mouse, has_gct_header=False)

    gtex_names = {row[1] for row in gtex_rows if len(row) >= 2}
    mouse_names = {row[1] for row in mouse_rows if len(row) >= 2}

    # Genes that have orthologs present in the other table
    gtex_allowed = {
        gene for gene in gtex_names if any(mouse in mouse_names for mouse in h_to_m.get(gene, []))
    }
    mouse_allowed = {
        mouse for mouse in mouse_names if any(human in gtex_names for human in m_to_h.get(mouse, []))
    }

    shared_pairs = {
        (human, mouse)
        for human in gtex_allowed
        for mouse in h_to_m.get(human, set())
        if mouse in mouse_allowed
    }

    gtex_out = args.gtex.with_name(args.gtex.stem + "_mapped.tsv")
    mouse_out = args.mouse.with_name(args.mouse.stem + "_mapped.tsv")

    kept_gtex = filter_rows(gtex_header, gtex_rows, gtex_allowed, True, gtex_prefix, gtex_out)
    kept_mouse = filter_rows(mouse_header, mouse_rows, mouse_allowed, False, [], mouse_out)

    print(f"Wrote GTEx mapped: {kept_gtex} rows -> {gtex_out}")
    print(f"Wrote mouse mapped: {kept_mouse} rows -> {mouse_out}")
    print(f"Shared ortholog gene pairs: {len(shared_pairs)}")


if __name__ == "__main__":
    main()
