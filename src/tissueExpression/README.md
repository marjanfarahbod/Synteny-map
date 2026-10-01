# tissueExpression scripts

## Overview
This folder contains small, task-focused scripts for preparing GTEx/mouse expression data, filtering orthologs, computing statistics, and comparing tissues across species.

This code was developed as part of an exploratory analysis and reflects the structure used during that process. The goal in section 3, 2 and 1 was to extract the required samples and the ortholog genes present in datasets from both species and have them in the same order for down stream analyses. There was some back and forth and redundancies in the pipeline since I learned about duplicates and inconsistencies as I processed the data. 

## Data files and references


## Workflow

The overall pipeline is illustrated below:

![Workflow diagram](../../docs/codeDiagram_expressions02.png)

## Requirements
- Python 3.9+
- numpy, matplotlib

## Scripts and inputs/outputs
- `gtex_extract_4Tissues.py`: extracts GTEx samples for Liver/Lung/Muscle/Kidney from annotations (`data/external/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt`) and TPMs (`data/external/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct`); outputs `data/processed/GTEx_4Tissues.txt` and sample ID pickle `data/interim/GTEx_4Tissues_present_samples.pkl`.
- `filter_gtex_mouse_orthologs.py`: filters GTEx matrix to genes in `data/processed/mouse_human_ortholog_one2one.tsv`; outputs `data/processed/GTEx_4Tissues_mouseOrthologGenes.tsv`.
- `filter_mouse_one2one.py`: filters ortholog export to unique one-to-one pairs; outputs `data/processed/mouse_human_ortholog_one2one.tsv`.
- `filter_mouse_4tissue_fpkm.py`: filters mouse FPKM file (`data/external/ST6_Lietal.tsv`) to ortholog genes and four tissue columns; outputs `data/processed/mouse_orthologs_4Tissue_FPKM.tsv`.
- `map_shared_ortholog_expression.py` / `_V2.py`: keep shared ortholog genes between GTEx and mouse matrices using `mouse_human_ortholog_one2one.tsv`; write mapped TSVs for GTEx/mouse. V2 drops duplicate gene rows before mapping.
- `reorder_mapped_expression.py`: reorder mapped GTEx/mouse TSVs to match the order in `data/interim/one2oneGeneNameMapping.tsv`; outputs ordered TSVs.
- `build_one2one_gene_mapping.py`: build human→mouse gene name pairs from ortholog map and optionally a human mapped file; outputs `data/interim/one2oneGeneNameMapping.tsv`.
- `extract_human_mouse_genes.py`: extracts Ensembl/gene symbols from GTEx mapped file; outputs `data/processed/human15370genes_mouseOrtholog.pkl`.
- `column_sum_normalize.py`: column-normalize expression tables to sum to 1e6; appends `_columnSumNorm` by default.
- `compute_tissue_stats_GTEx.py` / `compute_tissue_stats_mouse.py`: compute per-tissue mean/variance vectors and generate histograms/heatmaps; output pickles (`GTEx4TissueMeanAndVar.pkl`, `Mouse4TissueMeanAndVar.pkl`) and figures in `results/figures`.
- `mouseGTEx_comparison.py`: compare GTEx vs mouse mean vectors (sorted heatmaps, correlations, Jaccard similarities); saves figures in `results/figures`.
- `GTEx_mouse_TF_stats.py`: TF-focused stats using Lambert TF list, builds TF mean table, and produces per-tissue heatmaps/correlations/Jaccard plots; outputs `data/processed/LambertTFsIDsNames.pkl`, `data/processed/LambertTFs_human_mouse_means.tsv`, and figures in `results/figures`.

## Usage pattern
Run scripts with `python3 path/to/script.py` from the repo root. Paths in defaults are relative to the repository root. Adjust arguments as needed (see each script’s `--help`).


## Code Generation

This repository was developed with the assistance of OpenAI Codex.

Code was generated through interative prompting and subsequently reviewed, modified, and validated by the author. 