from __future__ import annotations

# Purpose: List BED records that overlap a specified human bin and mouse bin pair.
# Inputs:
#   - --bed: data/external/hg38_mm39.synNet.human.bed12
#   - --human: human bin as chrom:start-end
#   - --mouse: mouse bin as chrom:start-end
#   - --max: optional maximum number of overlapping records to print
# Outputs:
#   - Printed overlapping BED records and their human/mouse overlap lengths to stdout
#
# Notes:
#   - A record is reported only if it overlaps both the requested human bin and the requested mouse bin.
#   - The reported combined overlap value is min(human_overlap, mouse_overlap) for that record.

import argparse
from pathlib import Path
from typing import List, Tuple


def interval_overlap(a1: int, a2: int, b1: int, b2: int) -> int:
    return max(0, min(a2, b2) - max(a1, b1))


def parse_mouse(name_field: str) -> Tuple[str, int, int] | None:
    try:
        chrom, pos = name_field.split(":")
        s, e = pos.split("-")
        return chrom, int(s), int(e)
    except Exception:
        return None


def list_overlaps(bed_path: Path, h_bin: Tuple[str, int, int], m_bin: Tuple[str, int, int], max_rows: int = None) -> List[Tuple[str, int, int, str, int, int, int]]:
    rows = []
    h_chrom, hb1, hb2 = h_bin
    m_chrom, mb1, mb2 = m_bin
    with bed_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("	")
            if len(parts) < 4:
                continue
            hc, hs, he = parts[0], int(parts[1]), int(parts[2])
            mi = parse_mouse(parts[3])
            if mi is None:
                continue
            mc, ms, me = mi
            if hc != h_chrom or mc != m_chrom:
                continue
            h_ov = interval_overlap(hs, he, hb1, hb2)
            m_ov = interval_overlap(ms, me, mb1, mb2)
            if h_ov > 0 and m_ov > 0:
                rows.append((hc, hs, he, parts[3], h_ov, m_ov, min(h_ov, m_ov)))
                if max_rows and len(rows) >= max_rows:
                    break
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="List overlaps for a given human/mouse bin pair.")
    parser.add_argument("--bed", type=Path, default=Path("data/external/hg38_mm39.synNet.human.bed12"))
    parser.add_argument("--human", type=str, default="chr1:159000000-160000000", help="Human bin as chrom:start-end")
    parser.add_argument("--mouse", type=str, default="chr1:172000000-173000000", help="Mouse bin as chrom:start-end")
    parser.add_argument("--max", type=int, default=None, help="Optional limit of overlaps to print")
    args = parser.parse_args()

    def parse_bin(spec: str) -> Tuple[str, int, int]:
        chrom, pos = spec.split(":")
        s, e = pos.split("-")
        return chrom, int(s), int(e)

    h_bin = parse_bin(args.human)
    m_bin = parse_bin(args.mouse)
    rows = list_overlaps(args.bed, h_bin, m_bin, max_rows=args.max)
    for i, (hc, hs, he, name, h_ov, m_ov, combined) in enumerate(rows, 1):
        print(f"{i}: {hc} {hs}-{he} {name} | human_ov={h_ov} mouse_ov={m_ov} min={combined}")
    print(f"total overlaps found: {len(rows)}")


if __name__ == "__main__":
    main()
