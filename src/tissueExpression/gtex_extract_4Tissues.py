from __future__ import annotations

# Purpose: Extract GTEx samples for four tissues (Liver, Lung, Muscle, Kidney) and subset the RNA TPM
# matrix to those samples.
# Inputs:
#   - Annotations: data/external/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
#   - RNA TPM: data/external/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct
# Outputs:
#   - data/processed/GTEx_4Tissues.txt (subset TPM matrix)
#   - data/interim/GTEx_4Tissues_present_samples.pkl (dict of tissue -> sample IDs)

import csv
import pickle
from pathlib import Path


ANNOTATIONS_PATH = Path("data/external/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
RNA_PATH = Path("data/external/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct")
OUTPUT_PATH = Path("data/processed/GTEx_4Tissues.txt")
PRESENT_SAMPLES_PATH = Path("data/interim/GTEx_4Tissues_present_samples.pkl")

TISSUE_KEYWORDS = {
    "liver": "Liver",
    "lung": "Lung",
    "muscle": "Muscle",
    "kidney": "Kidney",
    "kideny": "Kidney",  # handle typo in request
}


def collect_sample_ids() -> tuple[dict[str, int], set[str], dict[str, set[str]]]:
    counts: dict[str, int] = {label: 0 for label in {"Liver", "Lung", "Muscle", "Kidney"}}
    sample_ids: set[str] = set()
    tissue_samples: dict[str, set[str]] = {label: set() for label in counts}

    with ANNOTATIONS_PATH.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Annotation file is missing a header row")
        if "SAMPID" not in reader.fieldnames:
            raise KeyError("Column 'SAMPID' not found in annotation file")
        if "SMTS" not in reader.fieldnames:
            raise KeyError("Column 'SMTS' not found in annotation file")

        for row in reader:
            tissue_note = row.get("SMTS", "")
            sampid = row.get("SAMPID", "")
            if not tissue_note or not sampid:
                continue
            normalized = tissue_note.lower()
            matched_labels = {
                TISSUE_KEYWORDS[keyword]
                for keyword in TISSUE_KEYWORDS
                if keyword in normalized
            }
            for label in matched_labels:
                counts[label] += 1
                sample_ids.add(sampid)
                tissue_samples[label].add(sampid)

    return counts, sample_ids, tissue_samples


def subset_rnaseq(
    sample_ids: set[str], tissue_samples: dict[str, set[str]]
) -> tuple[int, list[str], dict[str, list[str]]]:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    missing_ids: list[str] = []

    with RNA_PATH.open("r", newline="", encoding="utf-8") as infile, OUTPUT_PATH.open(
        "w", newline="", encoding="utf-8"
    ) as outfile:
        version_line = infile.readline().rstrip("\n")
        if not version_line:
            raise ValueError("RNA file is missing version line")
        dims_line = infile.readline().rstrip("\n")
        dims_parts = dims_line.split("\t")
        if len(dims_parts) < 2:
            raise ValueError("RNA file dims line is malformed")
        try:
            row_count = int(dims_parts[0])
        except ValueError as exc:
            raise ValueError("First value on dims line is not an integer") from exc

        header_line = infile.readline().rstrip("\n")
        header_fields = header_line.split("\t")
        if len(header_fields) < 3:
            raise ValueError("RNA file header must include Name, Description, and samples")
        if header_fields[0] != "Name" or header_fields[1] != "Description":
            raise ValueError("Expected first two RNA header columns to be Name and Description")

        available_samples = header_fields[2:]
        available_set = set(available_samples)
        missing_ids = sorted(sample_ids - available_set)
        present_samples = [sid for sid in available_samples if sid in sample_ids]
        present_indices = [idx for idx, sid in enumerate(available_samples) if sid in sample_ids]
        present_by_tissue = {
            tissue: sorted(sids & available_set) for tissue, sids in tissue_samples.items()
        }

        new_dims_line = f"{row_count}\t{len(present_samples)}"

        outfile.write(version_line + "\n")
        outfile.write(new_dims_line + "\n")
        header_out = ["Name", "Description"] + present_samples
        outfile.write("\t".join(header_out) + "\n")

        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")
        for row in reader:
            if len(row) < 2:
                continue
            core = row[:2]
            values = [row[2 + idx] for idx in present_indices]
            writer.writerow(core + values)

    return len(present_samples), missing_ids, present_by_tissue


def main() -> None:
    counts, sample_ids, tissue_samples = collect_sample_ids()
    print("Annotation counts per tissue:")
    for tissue in ("Liver", "Lung", "Muscle", "Kidney"):
        print(f"  {tissue}: {counts[tissue]}")

    sample_count, missing_ids, present_by_tissue = subset_rnaseq(sample_ids, tissue_samples)
    PRESENT_SAMPLES_PATH.parent.mkdir(parents=True, exist_ok=True)
    with PRESENT_SAMPLES_PATH.open("wb") as handle:
        pickle.dump(present_by_tissue, handle)

    print(f"Wrote samples: {sample_count}")
    print(f"Output: {OUTPUT_PATH}")
    print(f"Present sample lists saved to: {PRESENT_SAMPLES_PATH}")
    if missing_ids:
        print("Sample IDs missing from RNAseq header:")
        for sample_id in missing_ids:
            print(f"  {sample_id}")
    else:
        print("All selected sample IDs found in RNAseq file")


if __name__ == "__main__":
    main()
