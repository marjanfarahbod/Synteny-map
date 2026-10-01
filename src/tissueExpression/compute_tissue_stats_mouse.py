from __future__ import annotations

# Purpose: Compute per-tissue mean/variance vectors from mouse column-sum normalized FPKM data and
# mirror the GTEx stats workflow with comparable plots.
# Inputs:
#   - --input: data/processed/mouse_orthologs_4Tissue_FPKM_mapped_ordered_columnSumNorm.tsv (default)
# Outputs:
#   - --output: data/processed/Mouse4TissueMeanAndVar.pkl (means/variances)
#   - Figures in results/figures (histograms, mean/variance heatmaps)
# Notes:
#   - Tissue columns identified by prefixes mLi_, mLu_, mKi_, mMu_.

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pickle


TISSUES = ("Liver", "Lung", "Muscle", "Kidney")
PREFIXES = {"Liver": "mLi_", "Lung": "mLu_", "Muscle": "mMu_", "Kidney": "mKi_"}


def read_table(expr_path: Path) -> tuple[list[str], list[list[str]]]:
    """Read a TSV table with header and rows."""
    with expr_path.open("r", newline="", encoding="utf-8") as infile:
        header_line = infile.readline()
        if not header_line:
            raise ValueError(f"{expr_path} is empty or missing header")
        header = header_line.rstrip("\n").split("\t")
        reader = csv.reader(infile, delimiter="\t")
        rows = [row for row in reader if row]
    return header, rows


def compute_stats(header: list[str], rows: list[list[str]]) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Compute mean and variance per tissue across matching columns."""
    tissue_indices = {
        tissue: [idx for idx, name in enumerate(header) if name.startswith(prefix)]
        for tissue, prefix in PREFIXES.items()
    }

    means = {t: [] for t in TISSUES}
    variances = {t: [] for t in TISSUES}

    for row in rows:
        for tissue in TISSUES:
            idxs = tissue_indices[tissue]
            if not idxs:
                means[tissue].append(np.nan)
                variances[tissue].append(np.nan)
                continue
            vals = []
            for idx in idxs:
                try:
                    vals.append(float(row[idx]))
                except (ValueError, IndexError):
                    vals.append(0.0)
            arr = np.array(vals, dtype=float)
            means[tissue].append(float(np.mean(arr)))
            variances[tissue].append(float(np.var(arr)))

    means = {k: np.array(v, dtype=float) for k, v in means.items()}
    variances = {k: np.array(v, dtype=float) for k, v in variances.items()}
    return means, variances


def save_stats(output_path: Path, means: dict[str, np.ndarray], variances: dict[str, np.ndarray]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as handle:
        pickle.dump({"mean": means, "variance": variances}, handle)


def plot_heatmap(means: dict[str, np.ndarray], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    data = np.vstack([np.clip(means[t][:1000], None, 100) for t in TISSUES])
    plt.figure(figsize=(10, 4))
    plt.imshow(data, aspect="auto", cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Mean (normalized)")
    plt.yticks(range(len(TISSUES)), TISSUES)
    plt.xlabel("Gene index (first 1000)")
    plt.title("Mouse mean expression per tissue (first 1000)")
    plt.tight_layout()
    plt.savefig(output_dir / "Mouse4Tissue_mean_heatmap.pdf")
    plt.close()


def plot_mean_var_pdf(
    means: dict[str, np.ndarray], variances: dict[str, np.ndarray], tissue: str, output_dir: Path
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    mean_vals = means.get(tissue, np.array([]))
    var_vals = variances.get(tissue, np.array([]))
    if mean_vals.size:
        mean_vals = np.clip(mean_vals, None, 1000)
    if var_vals.size:
        var_vals = np.clip(var_vals, None, 1000)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(mean_vals[~np.isnan(mean_vals)], bins=50, color="tab:green", alpha=0.8)
    axes[0].set_title(f"{tissue} mean")
    axes[0].set_xlabel("Mean (clipped at 1000)")
    axes[0].set_ylabel("Count")
    axes[0].set_ylim(0, 12500)

    axes[1].hist(var_vals[~np.isnan(var_vals)], bins=50, color="tab:red", alpha=0.8)
    axes[1].set_title(f"{tissue} variance")
    axes[1].set_xlabel("Variance (clipped at 1000)")
    axes[1].set_ylabel("Count")
    axes[1].set_ylim(0, 12500)

    plt.tight_layout()
    outfile = output_dir / f"Mouse_{tissue}_mean_variance_hist.pdf"
    plt.savefig(outfile)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute per-gene mean and variance for mouse tissues and plot summaries."
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/processed/mouse_orthologs_4Tissue_FPKM_mapped_ordered_columnSumNorm.tsv"),
        help="Mouse expression TSV (column-sum normalized).",
    )
    parser.add_argument(
        "--output-pkl",
        type=Path,
        default=Path("data/processed/Mouse4TissueMeanAndVar.pkl"),
        help="Output pickle for means/variances.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory to save plots.",
    )
    args = parser.parse_args()

    header, rows = read_table(args.expression)
    means, variances = compute_stats(header, rows)
    save_stats(args.output_pkl, means, variances)
    plot_heatmap(means, args.fig_dir)
    for tissue in TISSUES:
        plot_mean_var_pdf(means, variances, tissue, args.fig_dir)
    print(f"Saved stats to {args.output_pkl}")
    print(f"Saved figures to {args.fig_dir}")


if __name__ == "__main__":
    main()
