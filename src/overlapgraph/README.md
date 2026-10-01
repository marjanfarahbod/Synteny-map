# Overlap graph - hg38 / mm39 1 Mb bin overlaps

## Overview

This folder contains code for building human (hg38) x mouse (mm39) overlap matrices in 1 Mb bins, treating them as bipartite graphs, and exporting/plotting their connected components.

The main pipeline processes the synteny net file from UCSC to generate an edge-line figure that links each human bin to the mouse bins it overlaps, colored by connected component.

This code was developed as part of an exploratory analysis and reflects the structure used during that process.

---

## Input Data

[UCSC for alignments](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm39/)

The synteny net file: [hg38.mm39.syn.net.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm39/hg38.mm39.syn.net.gz)

Chromosome lengths from UCSC: [mm39.chromsizes](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes), [hg38.chrom.sizes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes)

---

## Workflow

The pipeline for the edge-line figure (human bins 100-489, mouse bins 30-329) is:

```
data/external/hg38.mm39.syn.net
data/external/hg38.chrom.sizes
data/external/mm39.chrom.sizes
        │
        ▼  [1] src/overlapgraph/compute_paired_overlap_matrix_from_net.py
        │      --hg-start 100 --hg-regions 390 --mm-start 30 --mm-regions 300
        ▼
data/processed/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329.pkl
data/processed/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329.tsv   (same matrix as TSV, not used downstream)
        │
        ▼  [2] src/overlapgraph/export_cytoscape_components.py
        │      --matrix <the .pkl above>   (default --edge-threshold 10000)
        ▼
results/tables/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_cytoscape_edges.tsv
results/tables/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_cytoscape_nodes.tsv   (also produced; not used by the plot)
        │
        ▼  [3] src/overlapgraph/plot_cytoscape_edge_lines.py
        │      --edges <the _cytoscape_edges.tsv above>
        ▼
results/figures/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_edge_lines.pdf
```

Commands (from the repo root):

```
python3 src/overlapgraph/compute_paired_overlap_matrix_from_net.py --hg-start 100 --hg-regions 390 --mm-start 30 --mm-regions 300 --output-tsv data/processed/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329.tsv --output-pkl data/processed/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329.pkl
python3 src/overlapgraph/export_cytoscape_components.py --matrix data/processed/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329.pkl
python3 src/overlapgraph/plot_cytoscape_edge_lines.py --edges results/tables/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_cytoscape_edges.tsv --output results/figures/hg38_mm39_paired_overlap_1Mb_matrix_hg100_489_mm30_329_edge_lines.pdf
```

What each step does:

1. `compute_paired_overlap_matrix_from_net.py` parses the net file, keeps only paired (aligned) subsegments of each fill (gaps excluded), and builds a 390 x 300 bp-overlap matrix for human bins 100-489 x mouse bins 30-329.
2. `export_cytoscape_components.py` thresholds the matrix (edge if overlap > 10,000 bp), finds connected components, drops singleton vertices, and writes edge/node tables using absolute bin indices (offsets 100 and 30 are parsed from the `hg100_..._mm30_...` filename).
3. `plot_cytoscape_edge_lines.py` draws each edge as a line from human bin x at y=0 to mouse bin x at y=100, coloring lines by component with the `Paired` colormap.

Detailed descriptions of inputs, outputs, and functionality are provided in the header of each file.

---

## Notes

- `analyze_overlap_graph.py` produced the sibling `..._graph.pdf` (heatmap), but is not in the edge-line figure's dependency chain.
- `hg38_mm39_overlap_1Mb_matrix_490bins.pkl` (from `compute_overlap_matrix_pickle.py`, built from `hg38_mm39.synNet.human.bed12`, not the net file) is unrelated to this figure; it is the input to `analyze_bipartite_components.py`.
- The edge threshold was not passed explicitly in step 2, so the script's default (10,000 bp) applies.

---

## Output

The edge-line figure - see `results/figures/`.

---

## Requirements
- Python 3.9+ (tested with 3.12)
- numpy
- matplotlib

---
## Code Generation

This repository was developed with the assistance of OpenAI Codex.

Code was generated through iterative prompting and subsequently reviewed, modified, and validated by the author.

This README was created with the assistance of Claude Code (Anthropic).
