# Synteny map

## Overview

This repository contains code for exploring human (hg38) - mouse (mm39) synteny and comparing tissue expression across the two species.

This code was developed as part of an exploratory analysis and reflects the structure used during that process.

---

## Repository structure

Each folder has its own README with inputs, workflow, and outputs.

- [`src/syntenyMap/`](src/syntenyMap/README.md): synteny map figure of mm39 chromosomes onto hg38 chromosomes from the UCSC synteny net file, plus chain/MAF coverage counts.
- [`src/overlapgraph/`](src/overlapgraph/README.md): human x mouse 1 Mb bin overlap matrices, treated as bipartite graphs; connected components, Cytoscape export, and edge-line plots.
- [`src/tissueExpression/`](src/tissueExpression/README.md): GTEx and mouse expression preparation, one-to-one ortholog filtering, per-tissue statistics, and cross-species tissue comparisons.

---

## Input Data

Input data are not included in this repository. Scripts expect them under `data/external/` and write intermediate files to `data/interim/` and `data/processed/`, and figures/tables to `results/`. See each folder's README for sources.

---

## Usage pattern

Run scripts with `python3 src/<folder>/<script>.py` from the repository root. Default paths are relative to the repository root; see each script's header and `--help` for arguments.

---

## Requirements
- Python 3.9+
- numpy
- matplotlib

---
## Code Generation

This repository was developed with the assistance of OpenAI Codex.

Code was generated through iterative prompting and subsequently reviewed, modified, and validated by the author.

This README was created with the assistance of Claude Code (Anthropic).
