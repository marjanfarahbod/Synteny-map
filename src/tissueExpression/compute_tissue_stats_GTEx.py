from __future__ import annotations

# Purpose: Compute per-tissue mean/variance vectors from GTEx column-sum normalized data and produce
# summary plots (histograms, heatmaps) plus pickled stats.
# Inputs:
#   - --input: data/processed/GTEx_4Tissues_mouseOrthologGenes_columnSumNorm.tsv (default)
# Outputs:
#   - --output: data/processed/GTEx4TissueMeanAndVar.pkl (means/variances)
#   - Figures in results/figures (per-tissue histograms, mean/variance heatmaps)
# Notes:
#   - Assumes GCT-like input with GTEx sample IDs in header; tissues inferred via provided sample lists.

import argparse
import csv
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TISSUES = ("Liver", "Lung", "Muscle", "Kidney")


def read_header(expr_path: Path) -> tuple[list[str], list[str]]:
    """Return header fields and any prefix lines (GCT meta)."""
    with expr_path.open("r", newline="", encoding="utf-8") as infile:
        first = infile.readline()
        if not first:
            raise ValueError(f"{expr_path} is empty")
        if first.startswith("#1."):
            dims = infile.readline()
            header_line = infile.readline()
            if not header_line:
                raise ValueError(f"{expr_path} missing header line")
            header = header_line.rstrip("\n").split("\t")
            prefix = [first, dims]
        else:
            header = first.rstrip("\n").split("\t")
            prefix = []
    return header, prefix


def load_present_samples(path: Path) -> dict[str, list[str]]:
    with path.open("rb") as handle:
        data = pickle.load(handle)
    # Ensure ordering is deterministic
    return {tissue: sorted(samples) for tissue, samples in data.items()}


def compute_stats(
    expr_path: Path,
    header: list[str],
    prefix: list[str],
    tissue_samples: dict[str, list[str]],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    sample_to_idx = {name: idx for idx, name in enumerate(header)}
    tissue_indices = {
        tissue: [sample_to_idx[s] for s in samples if s in sample_to_idx]
        for tissue, samples in tissue_samples.items()
    }

    means = {t: [] for t in TISSUES}
    variances = {t: [] for t in TISSUES}

    with expr_path.open("r", newline="", encoding="utf-8") as infile:
        # skip prefix lines already read
        if prefix:
            infile.readline()
            infile.readline()
            infile.readline()
        else:
            infile.readline()

        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            for tissue in TISSUES:
                idxs = tissue_indices.get(tissue, [])
                if not idxs:
                    means[tissue].append(np.nan)
                    variances[tissue].append(np.nan)
                    continue
                values = []
                for idx in idxs:
                    try:
                        values.append(float(row[idx]))
                    except (ValueError, IndexError):
                        values.append(0.0)
                arr = np.array(values, dtype=float)
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
    plt.title("Mean expression per tissue (first 1000)")
    plt.tight_layout()
    plt.savefig(output_dir / "GTEx4Tissue_mean_heatmap.pdf")
    plt.close()


def plot_variance_hist(variances: dict[str, np.ndarray], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    for ax, tissue in zip(axes.flat, TISSUES):
        ax.hist(variances[tissue][~np.isnan(variances[tissue])], bins=50, color="tab:blue", alpha=0.7)
        ax.set_title(f"{tissue} variance")
        ax.set_xlabel("Variance")
        ax.set_ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_dir / "GTEx4Tissue_variance_hist.png", dpi=300)
    plt.close()


def plot_liver_variance_pdf(variances: dict[str, np.ndarray], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    liver_var = variances.get("Liver", np.array([]))
    if liver_var.size:
        liver_var = np.clip(liver_var, None, 1000)
    plt.figure(figsize=(6, 4))
    plt.hist(liver_var[~np.isnan(liver_var)], bins=50, color="tab:red", alpha=0.8)
    plt.title("Liver variance")
    plt.xlabel("Variance")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_dir / "GTEx_Liver_variance_hist.pdf")
    plt.close()


def plot_liver_mean_pdf(means: dict[str, np.ndarray], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    liver_mean = means.get("Liver", np.array([]))
    if liver_mean.size:
        liver_mean = np.clip(liver_mean, None, 1000)
    plt.figure(figsize=(6, 4))
    plt.hist(liver_mean[~np.isnan(liver_mean)], bins=50, color="tab:green", alpha=0.8)
    plt.title("Liver mean")
    plt.xlabel("Mean")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_dir / "GTEx_Liver_mean_hist.pdf")
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
    outfile = output_dir / f"GTEx_{tissue}_mean_variance_hist.pdf"
    plt.savefig(outfile)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute per-gene mean and variance per tissue using column-sum-normalized GTEx file."
    )
    parser.add_argument(
        "--samples-pkl",
        type=Path,
        default=Path("data/interim/GTEx_4Tissues_present_samples.pkl"),
        help="Pickle with tissue->sample list.",
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("data/processed/GTEx_4Tissues_mouseOrthologGenes_mapped_ordered_columnSumNorm.tsv"),
        help="Column-sum normalized GTEx TSV (GCT-style).",
    )
    parser.add_argument(
        "--output-pkl",
        type=Path,
        default=Path("data/processed/GTEx4TissueMeanAndVar.pkl"),
        help="Output pickle for means/variances.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory to save plots.",
    )
    args = parser.parse_args()

    tissue_samples = load_present_samples(args.samples_pkl)
    header, prefix = read_header(args.expression)
    means, variances = compute_stats(args.expression, header, prefix, tissue_samples)
    save_stats(args.output_pkl, means, variances)
    plot_heatmap(means, args.fig_dir)
    for tissue in TISSUES:
        plot_mean_var_pdf(means, variances, tissue, args.fig_dir)
    print(f"Saved stats to {args.output_pkl}")
    print(f"Saved figures to {args.fig_dir}")


if __name__ == "__main__":
    main()
