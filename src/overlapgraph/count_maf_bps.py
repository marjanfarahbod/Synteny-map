from __future__ import annotations

# Purpose: Count hg38 sequence rows in a synNet MAF file and summarize their aligned lengths.
# Inputs:
#   - --maf: data/external/hg38.mm39.synNet.maf
# Outputs:
#   - Printed summary statistics to stdout
#
# Notes:
#   - Only lines starting with "s hg38." are counted.
#   - The third numeric field on those lines is added to totalbps.

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count total aligned bp for s hg38.* rows in a synNet MAF file."
    )
    parser.add_argument(
        "--maf",
        type=Path,
        default=Path("data/external/hg38.mm39.synNet.maf"),
        help="Input synNet MAF file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    total_bps = 0
    hg38_line_count = 0
    max_size = 0

    with args.maf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("s hg38."):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            size = int(parts[3])
            total_bps += size
            hg38_line_count += 1
            if size > max_size:
                max_size = size

    print(f"maf\t{args.maf}")
    print(f"hg38_line_count\t{hg38_line_count}")
    print(f"totalbps\t{total_bps}")
    print(f"largest_integer\t{max_size}")


if __name__ == "__main__":
    main()
