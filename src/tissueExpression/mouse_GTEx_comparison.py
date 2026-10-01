from __future__ import annotations

# Purpose: Compare GTEx and mouse mean expression vectors per tissue by sorting GTEx values,
# reordering mouse accordingly, and plotting paired heatmaps plus correlation/Jaccard summaries.
# Inputs:
#   - --gtex-pkl: data/processed/GTEx4TissueMeanAndVar.pkl (default)
#   - --mouse-pkl: data/processed/Mouse4TissueMeanAndVar.pkl (default)
# Outputs:
#   - Heatmaps in results/figures (paired GTEx/mouse, Spearman/Pearson correlations, Jaccard at thresholds)
# Notes:
#   - Tissues assumed: Liver, Lung, Muscle.
#   - Correlations are visualized with fixed color limits (0.2–1).

import argparse
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TISSUES = ("Liver", "Lung", "Muscle")


def load_means(path: Path) -> dict[str, np.ndarray]:
    with path.open("rb") as handle:
        data = pickle.load(handle)
    return {t: np.array(data["mean"][t]) for t in TISSUES}


def plot_pair(tissue: str, gtex: np.ndarray, mouse: np.ndarray, out_path: Path) -> None:
    # Sort GTEx ascending and reorder mouse the same way
    order = np.argsort(gtex)
    gtex_sorted = np.clip(gtex[order], None, 100)
    mouse_reordered = np.clip(mouse[order], None, 100)
    data = np.vstack([gtex_sorted, mouse_reordered])

    plt.figure(figsize=(12, 2.5))
    plt.imshow(data, aspect="auto", cmap="pink", interpolation="nearest", rasterized=False)
    plt.colorbar(label="Mean (normalized)")
    plt.yticks([0, 1], [f"{tissue} GTEx (sorted)", f"{tissue} Mouse (reordered)"])
    plt.xlabel("Gene (sorted by GTEx mean)")
    plt.title(f"{tissue}: GTEx sorted vs Mouse reordered")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def plot_correlations(gtex_means: dict[str, np.ndarray], mouse_means: dict[str, np.ndarray], out_path: Path) -> None:
    # Group by species: first all GTEx tissues, then all Mouse tissues
    labels = [f"{t}_GTEx" for t in TISSUES] + [f"{t}_Mouse" for t in TISSUES]
    vectors = [gtex_means[t] for t in TISSUES] + [mouse_means[t] for t in TISSUES]
    n = len(vectors)
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            corr[i, j] = spearman_corr(vectors[i], vectors[j])

    vmin = 0.2
    vmax = 1
    plt.figure(figsize=(8, 6))
    x, y = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
    im = plt.pcolormesh(x, y, corr, cmap="YlGn", vmin=vmin, vmax=vmax, shading="auto")
    plt.colorbar(im, label="Spearman correlation")
    tick_pos = np.arange(n) + 0.5
    plt.xticks(tick_pos, labels, rotation=45, ha="right")
    plt.yticks(tick_pos, labels)
    for i in range(n):
        for j in range(n):
            plt.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=6, color="black")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def plot_correlations_pearson(
    gtex_means: dict[str, np.ndarray], mouse_means: dict[str, np.ndarray], cap: float, out_path: Path
) -> None:
    labels = [f"{t}_GTEx" for t in TISSUES] + [f"{t}_Mouse" for t in TISSUES]
    vectors = [np.clip(gtex_means[t], None, cap) for t in TISSUES] + [
        np.clip(mouse_means[t], None, cap) for t in TISSUES
    ]
    n = len(vectors)
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            corr[i, j] = pearson_corr(vectors[i], vectors[j])

    vmin = 0.2
    vmax = 1
    plt.figure(figsize=(8, 6))
    x, y = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
    im = plt.pcolormesh(x, y, corr, cmap="YlGn", vmin=vmin, vmax=vmax, shading="auto")
    plt.colorbar(im, label="Pearson correlation")
    tick_pos = np.arange(n) + 0.5
    plt.xticks(tick_pos, labels, rotation=45, ha="right")
    plt.yticks(tick_pos, labels)
    for i in range(n):
        for j in range(n):
            plt.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=6, color="black")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def plot_jaccard_heatmap(
    gtex_means: dict[str, np.ndarray], mouse_means: dict[str, np.ndarray], threshold: float, out_path: Path
) -> None:
    labels = [f"{t}_GTEx" for t in TISSUES] + [f"{t}_Mouse" for t in TISSUES]
    vectors = [(gtex_means[t] > threshold).astype(int) for t in TISSUES] + [
        (mouse_means[t] > threshold).astype(int) for t in TISSUES
    ]
    n = len(vectors)
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            inter = np.sum((vectors[i] == 1) & (vectors[j] == 1))
            union = np.sum((vectors[i] == 1) | (vectors[j] == 1))
            corr[i, j] = inter / union if union else np.nan
    vmin, vmax = np.nanmin(corr), np.nanmax(corr)
    plt.figure(figsize=(8, 6))
    x, y = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
    im = plt.pcolormesh(x, y, corr, cmap="YlGn", vmin=vmin, vmax=vmax, shading="auto")
    plt.colorbar(im, label="Jaccard similarity")
    tick_pos = np.arange(n) + 0.5
    plt.xticks(tick_pos, labels, rotation=45, ha="right")
    plt.yticks(tick_pos, labels)
    for i in range(n):
        for j in range(n):
            plt.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=6, color="black")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)
    plt.close()


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Spearman correlation without scipy."""
    if x.shape != y.shape:
        raise ValueError("Vectors must have the same shape")
    rx = _rank(x)
    ry = _rank(y)
    return np.corrcoef(rx, ry)[0, 1]


def _rank(a: np.ndarray) -> np.ndarray:
    """Return ranks with average handling for ties."""
    temp = a.argsort()
    ranks = np.empty_like(temp, dtype=float)
    ranks[temp] = np.arange(len(a))
    # average ranks for ties
    _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
    cumulative = np.cumsum(counts)
    start = cumulative - counts
    avg_rank = (start + cumulative - 1) / 2
    ranks = avg_rank[inv]
    return ranks


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.shape != y.shape:
        raise ValueError("Vectors must have the same shape")
    return float(np.corrcoef(x, y)[0, 1])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare GTEx and mouse mean vectors by sorting GTEx and reordering mouse accordingly."
    )
    parser.add_argument(
        "--gtex-pkl",
        type=Path,
        default=Path("data/processed/GTEx4TissueMeanAndVar.pkl"),
        help="Pickle with GTEx means/variances.",
    )
    parser.add_argument(
        "--mouse-pkl",
        type=Path,
        default=Path("data/processed/Mouse4TissueMeanAndVar.pkl"),
        help="Pickle with mouse means/variances.",
    )
    parser.add_argument(
        "--fig-dir",
        type=Path,
        default=Path("results/figures"),
        help="Directory to save comparison heatmaps.",
    )
    args = parser.parse_args()

    gtex_means = load_means(args.gtex_pkl)
    mouse_means = load_means(args.mouse_pkl)

    # Counts of mean values > 1
    for label, means in (("GTEx", gtex_means), ("Mouse", mouse_means)):
        for tissue in TISSUES:
            vec = means.get(tissue, np.array([]))
            count = int(np.sum(vec > 1))
            print(f"{label} {tissue}: {count} mean values > 1")

    for tissue in TISSUES:
        gtex_vec = gtex_means.get(tissue)
        mouse_vec = mouse_means.get(tissue)
        if gtex_vec is None or mouse_vec is None:
            continue
        out_path = args.fig_dir / f"{tissue}_GTExSorted_BW.pdf"
        plot_pair(tissue, gtex_vec, mouse_vec, out_path)
        print(f"Saved {out_path}")

    corr_path = args.fig_dir / "Mouse_GTEx_mean_spearman.svg"
    plot_correlations(gtex_means, mouse_means, corr_path)
    print(f"Saved {corr_path}")

    corr_pearson_path = args.fig_dir / "Mouse_GTEx_mean_pearson.svg"
    plot_correlations_pearson(gtex_means, mouse_means, 100, corr_pearson_path)
    print(f"Saved {corr_pearson_path}")

    corr_pearson_cap1000_path = args.fig_dir / "Mouse_GTEx_mean_pearson_cap1000.svg"
    plot_correlations_pearson(gtex_means, mouse_means, 1000, corr_pearson_cap1000_path)
    print(f"Saved {corr_pearson_cap1000_path}")

    for thr in (1, 2, 3):
        jac_path = args.fig_dir / f"Mouse_GTEx_jaccard_thr{thr}.pdf"
        plot_jaccard_heatmap(gtex_means, mouse_means, thr, jac_path)
        print(f"Saved {jac_path}")


if __name__ == "__main__":
    main()
