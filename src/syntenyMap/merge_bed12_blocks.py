#!/usr/bin/env python3
"""
Merge adjacent BED12 records when:
- They are on the same human chromosome (col 1),
- The gap between them is < max_gap bp (next.start - prev.end), and
- They map to the same mouse chromosome (chr in col 4, which is like 'chrN:start-end').

After merging two blocks, continue merging with subsequent blocks while the
criteria hold. For the merged record:
- Human: chrom = first.chrom, start = first.start, end = last.end
- Strand: keep from first record
- Name (col 4): keep 'mouseChr:firstMouseStart-lastMouseEnd' where firstMouseStart
  comes from the first record's name field, and lastMouseEnd from the last
  record's name field. The mouse chromosome must remain the same across merged
  records.
- Set score/itemRgb from the first record; set thickStart/thickEnd to merged
  human start/end.
- BED12 blocks: write a single block (blockCount=1) spanning the merged length,
  with blockSizes = merged_len and blockStarts = 0.

Usage:
  python3 exploratory/merge_bed12_blocks.py [MAX_GAP_BP]

If MAX_GAP_BP is omitted, defaults to 10,000.

Input (default): data/external/hg38_mm39.synNet.human.bed12
Output: data/external/hg38_mm39.synNet.human.merge.bed12
"""

from pathlib import Path
from typing import List, Optional, Tuple
import sys


ROOT = Path(__file__).resolve().parents[1]
IN_PATH = ROOT / "data" / "external" / "hg38_mm39.synNet.human.bed12"
OUT_PATH = ROOT / "data" / "external" / "hg38_mm39.synNet.human.merge.bed12"


def parse_mouse_name(name_field: str) -> Tuple[str, int, int]:
    # Expect: 'chrN:start-end'
    try:
        chrom_part, range_part = name_field.split(":", 1)
        start_s, end_s = range_part.split("-", 1)
        return chrom_part, int(start_s), int(end_s)
    except Exception:
        # Fallback: treat whole as chrom, unknown coords
        return name_field, 0, 0


def merge_bed12(in_path: Path, out_path: Path, max_gap: int = 10_000) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0

    def flush_cluster(cluster: List[List[str]]) -> None:
        nonlocal written
        if not cluster:
            return
        first = cluster[0]
        last = cluster[-1]

        chrom = first[0]
        start = int(first[1])
        end = int(last[2])
        score = first[4] if len(first) > 4 else "0"
        strand = first[5] if len(first) > 5 else "+"
        thickStart = str(start)
        thickEnd = str(end)
        itemRgb = first[8] if len(first) > 8 else "0"

        # Mouse name field
        mchr1, mstart1, _ = parse_mouse_name(first[3])
        _, _, mendN = parse_mouse_name(last[3])
        name = f"{mchr1}:{mstart1}-{mendN}"

        merged_len = max(0, end - start + 1)
        blockCount = "1"
        blockSizes = str(merged_len)
        blockStarts = "0"

        out_fields = [
            chrom,
            str(start),
            str(end),
            name,
            score,
            strand,
            thickStart,
            thickEnd,
            itemRgb,
            blockCount,
            blockSizes,
            blockStarts,
        ]
        out.write("\t".join(out_fields) + "\n")
        written += 1

    with in_path.open("r", encoding="utf-8", errors="ignore") as fh, out_path.open("w") as out:
        cluster: List[List[str]] = []
        cluster_mouse_chr: Optional[str] = None
        cluster_chrom: Optional[str] = None
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                # Not a valid BED12 row
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            mouse_chr, mouse_start, mouse_end = parse_mouse_name(parts[3])

            if not cluster:
                cluster = [parts]
                cluster_mouse_chr = mouse_chr
                cluster_chrom = chrom
                prev_end = end
                continue

            # Ensure same human chrom and mouse chrom
            if chrom != cluster_chrom or mouse_chr != cluster_mouse_chr:
                flush_cluster(cluster)
                cluster = [parts]
                cluster_mouse_chr = mouse_chr
                cluster_chrom = chrom
                prev_end = end
                continue

            # Check gap to previous record in cluster
            gap = start - prev_end
            # Only merge if the gap is non-negative and strictly less than max_gap
            if 0 <= gap < max_gap:
                cluster.append(parts)
                prev_end = end
            else:
                flush_cluster(cluster)
                cluster = [parts]
                cluster_mouse_chr = mouse_chr
                cluster_chrom = chrom
                prev_end = end

        # Flush last cluster
        flush_cluster(cluster)

    return written


def main() -> int:
    in_path = IN_PATH
    out_path = OUT_PATH
    if not in_path.exists():
        print(f"Input not found: {in_path}", file=sys.stderr)
        return 1

    max_gap = 10_000
    if len(sys.argv) > 1:
        try:
            max_gap = int(sys.argv[1])
            if max_gap < 0:
                raise ValueError
        except ValueError:
            print("MAX_GAP_BP must be a non-negative integer.", file=sys.stderr)
            return 2

    n = merge_bed12(in_path, out_path, max_gap=max_gap)
    print(f"Wrote {n:,} merged rows to {out_path} (max_gap={max_gap} bp)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
