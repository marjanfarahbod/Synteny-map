from __future__ import annotations

# Purpose: Convert UCSC net alignment files into BED12 rows by emitting one BED12 line
# per "fill" entry. Useful for regenerating hg38_mm39.synNet.human.bed12 from
# hg38.mm39.syn.net.
#
# Inputs:
#   - --net: path to a .net file (default: data/external/hg38.mm39.syn.net)
# Outputs:
#   - --output: BED12 file with one record per fill (default: data/external/hg38_mm39.synNet.human_generated.bed12)
#
# Notes:
#   - The script tracks the current target chromosome from "net" header lines.
#   - Only "fill" lines are converted; gaps/other records are skipped.

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert UCSC net file fill lines to BED12.")
    parser.add_argument(
        "--net",
        type=Path,
        default=Path("data/external/hg38.mm39.syn.net"),
        help="Input net file (e.g., hg38.mm39.syn.net).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/external/hg38_mm39.synNet.human_generated.bed12"),
        help="Output BED12 path.",
    )
    return parser.parse_args()


def convert_net_to_bed12(net_path: Path, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    current_chrom: str | None = None

    with net_path.open("r", encoding="utf-8") as net_file, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as out_file:
        for line in net_file:
            if not line.strip():
                continue
            if line.startswith("net"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    current_chrom = parts[1]
                continue

            if line.lstrip().startswith("fill"):
                if current_chrom is None:
                    raise ValueError("Encountered fill line before any net header.")
                parts = line.strip().split()
                if len(parts) < 7:
                    continue
                t_start = int(parts[1])
                t_size = int(parts[2])
                q_name = parts[3]
                strand = parts[4]
                q_start = int(parts[5])
                q_size = int(parts[6])

                chrom_start = t_start
                chrom_end = t_start + t_size
                name = f"{q_name}:{q_start}-{q_start + q_size}"
                score = "0"
                thick_start = chrom_start
                thick_end = chrom_end
                item_rgb = "0"
                block_count = "1"
                block_sizes = str(t_size)
                block_starts = "0"

                fields = [
                    current_chrom,
                    str(chrom_start),
                    str(chrom_end),
                    name,
                    score,
                    strand,
                    str(thick_start),
                    str(thick_end),
                    item_rgb,
                    block_count,
                    block_sizes,
                    block_starts,
                ]
                out_file.write("\t".join(fields) + "\n")


def main() -> None:
    args = parse_args()
    convert_net_to_bed12(args.net, args.output)


if __name__ == "__main__":
    main()
