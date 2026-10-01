from __future__ import annotations

# Purpose: Reorder mapped human/mouse expression tables according to a one-to-one gene name mapping,
# producing ordered TSV/GCT outputs aligned to the mapping order.
# Inputs:
#   - --human: data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped.tsv (default)
#   - --mouse: data/processed/mouse_orthologs_4Tissue_FPKM_mapped.tsv (default)
#   - --mapping: data/interim/one2oneGeneNameMapping.tsv (default; human gene, mouse gene, in desired order)
# Outputs:
#   - --human-out: data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped_ordered.tsv
#   - --mouse-out: data/processed/mouse_orthologs_4Tissue_FPKM_mapped_ordered.tsv

import argparse
import csv
from pathlib import Path


def load_mapping(path: Path) -> list[tuple[str, str]]:
    """Load ordered human->mouse gene name pairs."""
    pairs: list[tuple[str, str]] = []
    with path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError("Mapping file missing header or columns")
        for row in reader:
            if len(row) >= 2 and row[0] and row[1]:
                pairs.append((row[0], row[1]))
    return pairs


def read_gct_table(path: Path) -> tuple[list[str], list[list[str]], list[str]]:
    """Read GCT-like table: returns header, rows, and prefix lines [version, dims]."""
    rows: list[list[str]] = []
    with path.open("r", newline="", encoding="utf-8") as infile:
        version = infile.readline()
        dims = infile.readline()
        header_line = infile.readline()
        if not header_line:
            raise ValueError(f"{path} missing header line")
        header = header_line.rstrip("\n").split("\t")
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if row:
                rows.append(row)
    return header, rows, [version, dims]


def read_tsv_table(path: Path) -> tuple[list[str], list[list[str]]]:
    """Read standard TSV table with header and rows."""
    rows: list[list[str]] = []
    with path.open("r", newline="", encoding="utf-8") as infile:
        header_line = infile.readline()
        if not header_line:
            raise ValueError(f"{path} missing header line")
        header = header_line.rstrip("\n").split("\t")
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if row:
                rows.append(row)
    return header, rows


def write_gct_table(path: Path, prefix: list[str], header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as outfile:
        dims_parts = prefix[1].rstrip("\n").split("\t")
        if dims_parts and dims_parts[0].isdigit():
            dims_parts[0] = str(len(rows))
        outfile.write(prefix[0])
        outfile.write("\t".join(dims_parts) + "\n")
        outfile.write("\t".join(header) + "\n")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def write_tsv_table(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as outfile:
        outfile.write("\t".join(header) + "\n")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def reorder_rows(rows: list[list[str]], key_order: dict[str, int]) -> list[list[str]]:
    """Filter to keys present in key_order and reorder by that order."""
    filtered = [row for row in rows if len(row) >= 2 and row[1] in key_order]
    return sorted(filtered, key=lambda r: key_order[r[1]])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reorder human/mouse mapped expression files to align with one-to-one mapping order."
    )
    parser.add_argument(
        "--human",
        type=Path,
        default=Path("data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped.tsv"),
        help="Human mapped expression file (GCT-like).",
    )
    parser.add_argument(
        "--mouse",
        type=Path,
        default=Path("data/processed/mouse_orthologs_4Tissue_FPKM_mapped.tsv"),
        help="Mouse mapped expression file (TSV).",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("data/interim/one2oneGeneNameMapping.tsv"),
        help="Ordered mapping file with columns Gene name and Mouse gene name.",
    )
    parser.add_argument(
        "--human-out",
        type=Path,
        default=Path("data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped_ordered.tsv"),
        help="Output path for reordered human file.",
    )
    parser.add_argument(
        "--mouse-out",
        type=Path,
        default=Path("data/processed/mouse_orthologs_4Tissue_FPKM_mapped_ordered.tsv"),
        help="Output path for reordered mouse file.",
    )
    args = parser.parse_args()

    mapping_pairs = load_mapping(args.mapping)
    human_order = {h: idx for idx, (h, _) in enumerate(mapping_pairs)}
    mouse_order = {m: idx for idx, (_, m) in enumerate(mapping_pairs)}

    human_header, human_rows, human_prefix = read_gct_table(args.human)
    mouse_header, mouse_rows = read_tsv_table(args.mouse)

    human_reordered = reorder_rows(human_rows, human_order)
    mouse_reordered = reorder_rows(mouse_rows, mouse_order)

    write_gct_table(args.human_out, human_prefix, human_header, human_reordered)
    write_tsv_table(args.mouse_out, mouse_header, mouse_reordered)

    print(f"Wrote human ordered: {len(human_reordered)} rows -> {args.human_out}")
    print(f"Wrote mouse ordered: {len(mouse_reordered)} rows -> {args.mouse_out}")


if __name__ == "__main__":
    main()
