#!/usr/bin/env python3
"""
Scan the MAF (hg38.mm39.synNet.maf) and print an example of two
alignment blocks that overlap on hg38 chr1.

Outputs the full text of both blocks when the first overlap is found.
"""

from pathlib import Path
from typing import Optional, Tuple, List
import sys


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external"


def parse_block_chr1_interval(block_lines: List[str], human_src: str = "hg38.chr1") -> Tuple[Optional[Tuple[int, int]], bool]:
    """Return ((start, end), has_other_species) for hg38 chr1 within the block.
    End is exclusive. If no hg38 chr1 component, returns (None, has_other_species).
    """
    human_start: Optional[int] = None
    human_size: Optional[int] = None
    has_other = False
    for line in block_lines:
        if not line or line[0] == '#':
            continue
        if line.startswith('s '):
            parts = line.split()
            if len(parts) < 7:
                continue
            src = parts[1]
            try:
                start = int(parts[2])
                size = int(parts[3])
            except ValueError:
                continue
            if src == human_src:
                human_start = start
                human_size = size
            else:
                has_other = True
    if human_start is not None and human_size is not None:
        return (human_start, human_start + human_size), has_other
    return None, has_other


def find_overlap_example(maf_path: Path) -> Optional[Tuple[List[str], List[str]]]:
    # First pass: collect all chr1 intervals with file offsets
    intervals: List[Tuple[int, int, int, int]] = []  # (start, end, off_start, off_end)

    with maf_path.open('r', encoding='utf-8', errors='ignore') as fh:
        block_lines: List[str] = []
        in_block = False
        block_start_off: Optional[int] = None

        while True:
            pos = fh.tell()
            raw = fh.readline()
            if raw == '':
                # EOF
                if in_block and block_lines and block_start_off is not None:
                    interval, has_other = parse_block_chr1_interval(block_lines)
                    if interval is not None and has_other:
                        s, e = interval
                        intervals.append((s, e, block_start_off, pos))
                break
            line = raw.rstrip('\n')
            if line.startswith('a '):
                # flush previous block
                if in_block and block_lines and block_start_off is not None:
                    interval, has_other = parse_block_chr1_interval(block_lines)
                    if interval is not None and has_other:
                        s, e = interval
                        intervals.append((s, e, block_start_off, pos))
                # start new block
                in_block = True
                block_lines = [line]
                block_start_off = pos
            elif not in_block:
                continue
            elif line.strip() == '':
                # end of block
                if block_lines and block_start_off is not None:
                    interval, has_other = parse_block_chr1_interval(block_lines)
                    if interval is not None and has_other:
                        s, e = interval
                        intervals.append((s, e, block_start_off, fh.tell()))
                in_block = False
                block_lines = []
                block_start_off = None
            else:
                block_lines.append(line)

    if not intervals:
        return None

    # Second pass: sort by start and find first overlap pair
    intervals.sort(key=lambda x: (x[0], x[1]))
    pair: Optional[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]]] = None
    prev = intervals[0]
    for cur in intervals[1:]:
        if cur[0] < prev[1]:  # overlap
            pair = (prev, cur) if prev[0] <= cur[0] else (cur, prev)
            break
        prev = cur

    if pair is None:
        return None

    # Read and return the two blocks by offsets
    with maf_path.open('r', encoding='utf-8', errors='ignore') as fh:
        (s1, e1, o1s, o1e), (s2, e2, o2s, o2e) = pair
        fh.seek(o1s)
        b1 = fh.read(o1e - o1s).splitlines()
        fh.seek(o2s)
        b2 = fh.read(o2e - o2s).splitlines()
    return b1, b2


def main() -> int:
    maf_path = DATA_DIR / 'hg38.mm39.synNet.maf'
    if not maf_path.exists():
        print(f"Missing file: {maf_path}", file=sys.stderr)
        return 1
    print('Searching for overlapping hg38 chr1 alignment blocks...')
    res = find_overlap_example(maf_path)
    if res is None:
        print('No overlap found within scanning window.', file=sys.stderr)
        return 2
    block1, block2 = res
    print('--- Overlap example: Block 1 ---')
    for l in block1:
        print(l)
    print('\n--- Overlap example: Block 2 ---')
    for l in block2:
        print(l)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
