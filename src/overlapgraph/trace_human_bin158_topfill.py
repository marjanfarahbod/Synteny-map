from __future__ import annotations

# Purpose: Trace how line 1024924 from hg38.mm39.syn.net is split by immediate gap lines for human bin 158.
# Inputs:
#   - data/external/hg38.mm39.syn.net
# Outputs:
#   - results/tables/human_bin158_topfill_gap_trace.tsv
#
# Notes:
#   - This traces the current compute_paired_overlap_matrix_from_net.py gap-walking behavior.
#   - Human bin 158 is chr1:158000000-159000000.

from collections import defaultdict
from pathlib import Path

FILL_LINE = 1_024_924
HUMAN_BIN_START = 158_000_000
HUMAN_BIN_END = 159_000_000
MOUSE_BIN_SIZE = 1_000_000


def main() -> None:
    net_path = Path("data/external/hg38.mm39.syn.net")
    out_path = Path("results/tables/human_bin158_topfill_gap_trace.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fill = None
    gaps = []
    with net_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line_number < FILL_LINE:
                continue
            stripped = line.strip()
            indent = len(line) - len(line.lstrip(" "))
            if line_number == FILL_LINE:
                parts = stripped.split()
                fill = {
                    "t_start": int(parts[1]),
                    "t_size": int(parts[2]),
                    "q_chrom": parts[3],
                    "strand": parts[4],
                    "q_start": int(parts[5]),
                    "q_size": int(parts[6]),
                    "raw": stripped,
                }
                continue
            if indent <= 1:
                break
            if indent == 2 and stripped.startswith("gap "):
                parts = stripped.split()
                gaps.append(
                    {
                        "line_number": line_number,
                        "t_start": int(parts[1]),
                        "t_size": int(parts[2]),
                        "q_start": int(parts[5]),
                        "q_size": int(parts[6]),
                    }
                )

    if fill is None:
        raise ValueError(f"Fill line {FILL_LINE} not found")

    rows = []
    totals = defaultdict(int)

    t_cursor = fill["t_start"]
    q_cursor = fill["q_start"] + fill["q_size"]
    segment_index = 0

    for gap in gaps:
        t_block = gap["t_start"] - t_cursor
        q_block = q_cursor - gap["q_start"]
        pair_len = min(t_block, q_block)
        if pair_len > 0:
            segment_index += 1
            q_segment_start = q_cursor - pair_len
            q_segment_end = q_cursor
            add_segment_rows(
                rows,
                totals,
                segment_index,
                gap["line_number"],
                t_cursor,
                t_cursor + pair_len,
                q_segment_start,
                q_segment_end,
                pair_len,
                t_block,
                q_block,
            )

        t_cursor = gap["t_start"] + gap["t_size"]
        q_cursor = gap["q_start"] - gap["q_size"]

    header = [
        "segment_index",
        "next_gap_line",
        "segment_human_start",
        "segment_human_end",
        "segment_mouse_start",
        "segment_mouse_end",
        "pair_len",
        "t_block_before_gap",
        "q_block_before_gap",
        "human_bin_overlap_start",
        "human_bin_overlap_end",
        "mapped_mouse_start",
        "mapped_mouse_end",
        "mouse_bin",
        "counted_overlap",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[name]) for name in header) + "\n")

    print(f"fill_line\t{FILL_LINE}")
    print(f"fill_raw\t{fill['raw']}")
    print(f"immediate_gap_count\t{len(gaps)}")
    print(f"segments_counted_in_human_bin158\t{len(rows)}")
    print("mouse_bin_totals")
    for mouse_bin in sorted(totals):
        print(f"{mouse_bin}\t{totals[mouse_bin]}")
    print(f"trace_table\t{out_path}")


def add_segment_rows(
    rows: list[dict[str, int]],
    totals: defaultdict[int, int],
    segment_index: int,
    next_gap_line: int,
    t_start: int,
    t_end: int,
    q_start: int,
    q_end: int,
    pair_len: int,
    t_block: int,
    q_block: int,
) -> None:
    human_overlap_start = max(t_start, HUMAN_BIN_START)
    human_overlap_end = min(t_end, HUMAN_BIN_END)
    if human_overlap_end <= human_overlap_start:
        return

    mapped_mouse_start = q_end - (human_overlap_end - t_start)
    mapped_mouse_end = q_end - (human_overlap_start - t_start)
    first_mouse_bin = mapped_mouse_start // MOUSE_BIN_SIZE
    last_mouse_bin = (mapped_mouse_end - 1) // MOUSE_BIN_SIZE

    for mouse_bin in range(first_mouse_bin, last_mouse_bin + 1):
        bin_start = mouse_bin * MOUSE_BIN_SIZE
        bin_end = bin_start + MOUSE_BIN_SIZE
        counted_overlap = max(0, min(mapped_mouse_end, bin_end) - max(mapped_mouse_start, bin_start))
        if counted_overlap <= 0:
            continue
        totals[mouse_bin] += counted_overlap
        rows.append(
            {
                "segment_index": segment_index,
                "next_gap_line": next_gap_line,
                "segment_human_start": t_start,
                "segment_human_end": t_end,
                "segment_mouse_start": q_start,
                "segment_mouse_end": q_end,
                "pair_len": pair_len,
                "t_block_before_gap": t_block,
                "q_block_before_gap": q_block,
                "human_bin_overlap_start": human_overlap_start,
                "human_bin_overlap_end": human_overlap_end,
                "mapped_mouse_start": mapped_mouse_start,
                "mapped_mouse_end": mapped_mouse_end,
                "mouse_bin": mouse_bin,
                "counted_overlap": counted_overlap,
            }
        )


if __name__ == "__main__":
    main()
