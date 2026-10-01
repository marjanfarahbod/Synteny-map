from __future__ import annotations

# Purpose: Build a human x mouse 1 Mb overlap matrix from a UCSC syn.net file using only paired fill segments.
# Inputs:
#   - --net: data/external/hg38.mm39.syn.net
#   - --hg-chrom: data/external/hg38.chrom.sizes
#   - --mm-chrom: data/external/mm39.chrom.sizes
#   - --hg-start: start bin index for human
#   - --mm-start: start bin index for mouse
#   - --hg-regions: number of 1 Mb bins to include for human
#   - --mm-regions: number of 1 Mb bins to include for mouse
# Outputs:
#   - --output-tsv: tab-delimited overlap matrix
#   - --output-pkl: pickle file containing the overlap matrix
#
# Notes:
#   - Only paired sequence is counted: each fill is split into aligned subsegments between immediate gap lines.
#   - Nested fills are also counted, but parent fills contribute only their own paired subsegments, so gaps are not double-counted.
#   - Matrix entries are computed by intersecting the paired-coordinate axis of each human/mouse subsegment with the requested bins.

import argparse
import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import numpy as np

Bin = tuple[str, int, int]


@dataclass
class GapRecord:
    t_start: int
    t_size: int
    q_start: int
    q_size: int


@dataclass
class FillState:
    t_chrom: str
    indent: int
    t_start: int
    t_size: int
    q_name: str
    strand: str
    q_start: int
    q_size: int
    q_chrom_size: int
    gaps: list[GapRecord] = field(default_factory=list)


@dataclass
class PairedSegment:
    t_chrom: str
    t_start: int
    t_end: int
    q_chrom: str
    q_start: int
    q_end: int
    strand: str


def load_chrom_sizes(chrom_sizes: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with chrom_sizes.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, size_str = line.strip().split()[:2]
            sizes[chrom] = int(size_str)
    return sizes


def generate_bins(chrom_sizes: Path, start_bins: int, num_bins: int, bin_size: int = 1_000_000) -> list[Bin]:
    bins: list[Bin] = []
    total_bins_seen = 0
    with chrom_sizes.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, size_str = line.strip().split()[:2]
            size = int(size_str)
            start = 0
            while start < size and len(bins) < num_bins + max(0, start_bins - total_bins_seen):
                end = min(start + bin_size, size)
                if total_bins_seen >= start_bins and len(bins) < num_bins:
                    bins.append((chrom, start, end))
                total_bins_seen += 1
                start += bin_size
            if len(bins) >= num_bins:
                break
    if len(bins) < num_bins:
        raise ValueError(f"Not enough genome length to build {num_bins} bins starting at bin {start_bins}")
    return bins


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def bins_for_interval(chrom: str, start: int, end: int, bins: list[Bin]) -> list[int]:
    return [
        idx
        for idx, (b_chrom, b_start, b_end) in enumerate(bins)
        if b_chrom == chrom and interval_overlap(start, end, b_start, b_end) > 0
    ]


def finalize_fill(fill: FillState) -> Iterator[PairedSegment]:
    t_cursor = fill.t_start
    q_cursor = fill.q_start if fill.strand == "+" else fill.q_start + fill.q_size
    t_end = fill.t_start + fill.t_size
    q_end = fill.q_start + fill.q_size if fill.strand == "+" else fill.q_start

    for gap in sorted(fill.gaps, key=lambda record: record.t_start):
        t_block = gap.t_start - t_cursor
        q_block = gap.q_start - q_cursor if fill.strand == "+" else q_cursor - gap.q_start
        pair_len = min(t_block, q_block)
        if pair_len > 0:
            yield build_segment(fill, t_cursor, q_cursor, pair_len)
        t_cursor = gap.t_start + gap.t_size
        q_cursor = gap.q_start + gap.q_size if fill.strand == "+" else gap.q_start - gap.q_size

    t_block = t_end - t_cursor
    q_block = q_end - q_cursor if fill.strand == "+" else q_cursor - q_end
    pair_len = min(t_block, q_block)
    if pair_len > 0:
        yield build_segment(fill, t_cursor, q_cursor, pair_len)


def build_segment(fill: FillState, t_cursor: int, q_cursor: int, pair_len: int) -> PairedSegment:
    if fill.strand == "+":
        q_start = q_cursor
        q_end = q_cursor + pair_len
    else:
        q_end = q_cursor
        q_start = q_cursor - pair_len
    return PairedSegment(
        t_chrom=fill.t_chrom,
        t_start=t_cursor,
        t_end=t_cursor + pair_len,
        q_chrom=fill.q_name,
        q_start=q_start,
        q_end=q_end,
        strand=fill.strand,
    )


def iter_paired_segments(net_path: Path, q_sizes: dict[str, int]) -> Iterator[PairedSegment]:
    current_t_chrom: str | None = None
    stack: list[FillState] = []

    with net_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue

            stripped = line.strip()
            indent = len(line) - len(line.lstrip(" "))

            if stripped.startswith("net "):
                while stack:
                    yield from finalize_fill(stack.pop())
                parts = stripped.split()
                current_t_chrom = parts[1]
                continue

            while stack and indent <= stack[-1].indent:
                yield from finalize_fill(stack.pop())

            if stripped.startswith("fill "):
                if current_t_chrom is None:
                    raise ValueError("Encountered fill line before any net header.")
                parts = stripped.split()
                if len(parts) < 7:
                    continue
                q_name = parts[3]
                if q_name not in q_sizes:
                    raise ValueError(f"Missing chromosome size for query chromosome {q_name}")
                stack.append(
                    FillState(
                        t_chrom=current_t_chrom,
                        indent=indent,
                        t_start=int(parts[1]),
                        t_size=int(parts[2]),
                        q_name=q_name,
                        strand=parts[4],
                        q_start=int(parts[5]),
                        q_size=int(parts[6]),
                        q_chrom_size=q_sizes[q_name],
                    )
                )
                continue

            if stripped.startswith("gap ") and stack and indent == stack[-1].indent + 1:
                parts = stripped.split()
                if len(parts) >= 7:
                    stack[-1].gaps.append(
                        GapRecord(
                            t_start=int(parts[1]),
                            t_size=int(parts[2]),
                            q_start=int(parts[5]),
                            q_size=int(parts[6]),
                        )
                    )

    while stack:
        yield from finalize_fill(stack.pop())


def paired_axis_overlap(segment: PairedSegment, h_start: int, h_end: int, m_start: int, m_end: int) -> int:
    human_axis_start = h_start - segment.t_start
    human_axis_end = h_end - segment.t_start
    if segment.strand == "+":
        mouse_axis_start = m_start - segment.q_start
        mouse_axis_end = m_end - segment.q_start
    else:
        mouse_axis_start = segment.q_end - m_end
        mouse_axis_end = segment.q_end - m_start
    return interval_overlap(human_axis_start, human_axis_end, mouse_axis_start, mouse_axis_end)


def build_matrix(net_path: Path, hg_bins: list[Bin], mm_bins: list[Bin], mm_sizes: dict[str, int]) -> np.ndarray:
    matrix = np.zeros((len(hg_bins), len(mm_bins)), dtype=np.int64)
    for segment in iter_paired_segments(net_path, mm_sizes):
        h_indices = bins_for_interval(segment.t_chrom, segment.t_start, segment.t_end, hg_bins)
        m_indices = bins_for_interval(segment.q_chrom, segment.q_start, segment.q_end, mm_bins)
        if not h_indices or not m_indices:
            continue

        for hi in h_indices:
            _, hb_start, hb_end = hg_bins[hi]
            h_start = max(segment.t_start, hb_start)
            h_end = min(segment.t_end, hb_end)
            if h_end <= h_start:
                continue
            for mi in m_indices:
                _, mb_start, mb_end = mm_bins[mi]
                m_start = max(segment.q_start, mb_start)
                m_end = min(segment.q_end, mb_end)
                if m_end <= m_start:
                    continue
                matrix[hi, mi] += paired_axis_overlap(segment, h_start, h_end, m_start, m_end)
    return matrix


def write_matrix(matrix: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join([""] + [f"mouse_bin_{i + 1}" for i in range(matrix.shape[1])]) + "\n")
        for i in range(matrix.shape[0]):
            row = [f"human_bin_{i + 1}"] + [str(value) for value in matrix[i]]
            handle.write("\t".join(row) + "\n")


def write_pickle(matrix: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as handle:
        pickle.dump(matrix, handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute a 1 Mb paired overlap matrix between hg38 and mm39 from a syn.net file."
    )
    parser.add_argument("--net", type=Path, default=Path("data/external/hg38.mm39.syn.net"))
    parser.add_argument("--hg-chrom", type=Path, default=Path("data/external/hg38.chrom.sizes"))
    parser.add_argument("--mm-chrom", type=Path, default=Path("data/external/mm39.chrom.sizes"))
    parser.add_argument("--hg-start", type=int, default=0, help="Start bin (1 Mb bins) for human")
    parser.add_argument("--mm-start", type=int, default=0, help="Start bin (1 Mb bins) for mouse")
    parser.add_argument("--regions", type=int, default=None, help="Number of 1 Mb bins to include per species")
    parser.add_argument("--hg-regions", type=int, default=None, help="Number of 1 Mb bins to include for human")
    parser.add_argument("--mm-regions", type=int, default=None, help="Number of 1 Mb bins to include for mouse")
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=Path("data/processed/hg38_mm39_paired_overlap_1Mb_matrix.tsv"),
    )
    parser.add_argument(
        "--output-pkl",
        type=Path,
        default=Path("data/processed/hg38_mm39_paired_overlap_1Mb_matrix.pkl"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hg_regions = args.hg_regions if args.hg_regions is not None else args.regions
    mm_regions = args.mm_regions if args.mm_regions is not None else args.regions
    if hg_regions is None or mm_regions is None:
        raise ValueError("Provide --regions or both --hg-regions and --mm-regions")
    mm_sizes = load_chrom_sizes(args.mm_chrom)
    hg_bins = generate_bins(args.hg_chrom, start_bins=args.hg_start, num_bins=hg_regions)
    mm_bins = generate_bins(args.mm_chrom, start_bins=args.mm_start, num_bins=mm_regions)
    matrix = build_matrix(args.net, hg_bins, mm_bins, mm_sizes)
    write_matrix(matrix, args.output_tsv)
    write_pickle(matrix, args.output_pkl)
    print(
        f"Built paired matrix of shape {matrix.shape} "
        f"(hg start bin {args.hg_start}, mm start bin {args.mm_start}, "
        f"hg regions {hg_regions}, mm regions {mm_regions}) "
        f"-> {args.output_tsv}; pickle: {args.output_pkl}"
    )


if __name__ == "__main__":
    main()
