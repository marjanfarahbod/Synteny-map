from __future__ import annotations

# Purpose: Build TF-focused GTEx vs Mouse mean expression tables and visualizations.
# Inputs:
#   - --table: Lambert TF table (default: data/external/HumanTFsLambert2018_mmc2/Table S1. Related to Figure 1B-Table 1.tsv)
#   - --human-genes: data/processed/human15370genes_mouseOrtholog.pkl
#   - --mouse-genes: data/processed/mouse_orthologs_4Tissue_FPKM_mapped_ordered.tsv
#   - --gtex-stats: data/processed/GTEx4TissueMeanAndVar.pkl
#   - --mouse-stats: data/processed/Mouse4TissueMeanAndVar.pkl
# Outputs:
#   - TF ID/name pickle: data/processed/LambertTFsIDsNames.pkl
#   - TF mean table: data/processed/LambertTFs_human_mouse_means.tsv
#   - Figures in results/figures (per-tissue GTEx-sorted heatmaps, Pearson/Spearman correlations, Jaccard)
# Notes:
#   - Tissues: Liver, Lung, Muscle (Kidney excluded).
#   - Correlation heatmaps use color limits 0.2–1.

import argparse
import csv
import pickle
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

TISSUES = ("Liver", "Lung", "Muscle")


def extract_tf_ids(table_path: Path) -> list[tuple[str, str]]:
    tf_pairs: list[tuple[str, str]] = []
    with table_path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.reader(infile, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 4:
                continue
            if row[3].strip().lower() != "yes":
                continue
            gene_id = row[0].strip()
            gene_name = row[1].strip()
            if gene_id and gene_name:
                tf_pairs.append((gene_id, gene_name))
    return tf_pairs


def main() -> None:
    parser = argparse.ArgumentParser(description="GTEx vs Mouse TF stats helper.")
    parser.add_argument(
        "--table",
        type=Path,
        default=Path("data/external/HumanTFsLambert2018_mmc2/Table S1. Related to Figure 1B-Table 1.tsv"),
        help="Path to Lambert TF table (Table S1).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/LambertTFsIDsNames.pkl"),
        help="Output pickle for TF IDs/names.",
    )
    parser.add_argument(
        "--human-genes",
        type=Path,
        default=Path("data/processed/human15370genes_mouseOrtholog.pkl"),
        help="Path to human gene list pickle.",
    )
    parser.add_argument(
        "--mouse-genes",
        type=Path,
        default=Path("data/processed/mouse_orthologs_4Tissue_FPKM_mapped_ordered.tsv"),
        help="Mouse mapped ordered TSV (to read gene symbols).",
    )
    parser.add_argument(
        "--mouse-stats",
        type=Path,
        default=Path("data/processed/Mouse4TissueMeanAndVar.pkl"),
        help="Path to mouse stats pickle.",
    )
    parser.add_argument(
        "--gtex-stats",
        type=Path,
        default=Path("data/processed/GTEx4TissueMeanAndVar.pkl"),
        help="Path to GTEx stats pickle.",
    )
    parser.add_argument(
        "--output-table",
        type=Path,
        default=Path("data/processed/LambertTFs_human_mouse_means.tsv"),
        help="Output TSV for TF means/orthologs.",
    )
    args = parser.parse_args()

    tf_pairs = extract_tf_ids(args.table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wb") as handle:
        pickle.dump(tf_pairs, handle)
    print(f"Wrote {len(tf_pairs)} TF entries to {args.output}")

    tf_names = {name for _, name in tf_pairs}

    # Load gene lists
    with args.human_genes.open("rb") as f:
        human_genes = pickle.load(f)  # list of (ensembl, symbol)

    mouse_symbols: list[str] = []
    with args.mouse_genes.open("r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) >= 2:
                mouse_symbols.append(row[1])

    # Load stats
    with args.gtex_stats.open("rb") as f:
        gtex_stats = pickle.load(f)["mean"]
    with args.mouse_stats.open("rb") as f:
        mouse_stats = pickle.load(f)["mean"]

    rows = []
    for idx, (ens, sym) in enumerate(human_genes):
        if sym not in tf_names:
            continue
        mouse_sym = mouse_symbols[idx] if idx < len(mouse_symbols) else ""
        human_means = [gtex_stats[t][idx] if idx < len(gtex_stats[t]) else np.nan for t in TISSUES]
        mouse_means = [mouse_stats[t][idx] if idx < len(mouse_stats[t]) else np.nan for t in TISSUES]
        rows.append([ens, sym, mouse_sym] + human_means + mouse_means)

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    header = ["Ensembl_ID", "Human_symbol", "Mouse_symbol"] + [f"GTEx_{t}" for t in TISSUES] + [
        f"Mouse_{t}" for t in TISSUES
    ]
    with args.output_table.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"Wrote TF mean table: {len(rows)} rows -> {args.output_table}")

    # Heatmaps for each tissue: sort GTEx means, reorder mouse accordingly
    fig_dir = Path("results/figures")
    fig_dir.mkdir(parents=True, exist_ok=True)
    # Build arrays for TFs in rows list
    human_means_by_tissue = {t: np.array([r[3 + i] for r in rows], dtype=float) for i, t in enumerate(TISSUES)}
    mouse_means_by_tissue = {t: np.array([r[3 + len(TISSUES) + i] for r in rows], dtype=float) for i, t in enumerate(TISSUES)}
    for t in TISSUES:
        order = np.argsort(human_means_by_tissue[t])
        h_sorted = np.clip(human_means_by_tissue[t][order], None, 100)
        m_reordered = np.clip(mouse_means_by_tissue[t][order], None, 100)
        data_mat = np.vstack([h_sorted, m_reordered])
        plt.figure(figsize=(8, 2.5))
        plt.imshow(data_mat, aspect="auto", cmap="pink", interpolation="nearest", rasterized=False)
        plt.colorbar(label="Mean (normalized)")
        plt.yticks([0, 1], [f"{t} GTEx (sorted)", f"{t} Mouse (reordered)"])
        plt.xlabel("TFs (sorted by GTEx mean)")
        plt.title(f"{t}: GTEx sorted vs Mouse reordered (TFs)")
        plt.tight_layout()
        out_pdf = fig_dir / f"{t}_TF_GTExSorted_BW.pdf"
        plt.savefig(out_pdf)
        plt.close()
        print(f"Saved {out_pdf}")

    # Pearson correlation heatmap for TF mean vectors (capped)
    labels = [f"{t}_GTEx" for t in TISSUES] + [f"{t}_Mouse" for t in TISSUES]
    vectors = [np.clip(np.array([r[3 + i] for r in rows], dtype=float), None, 1000) for i, _ in enumerate(TISSUES)] + [
        np.clip(np.array([r[3 + len(TISSUES) + i] for r in rows], dtype=float), None, 1000)
        for i, _ in enumerate(TISSUES)
    ]
    n = len(vectors)
    corr = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            corr[i, j] = np.corrcoef(vectors[i], vectors[j])[0, 1]
    vmin, vmax = 0.2, 1
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
    corr_pdf = fig_dir / "Mouse_GTEx_mean_pearson_TF_noKidney_cap1000.svg"
    plt.savefig(corr_pdf)
    plt.close()
    print(f"Saved {corr_pdf}")

    # Spearman correlation heatmap for TF mean vectors (no cap)
    vectors_s = [np.array([r[3 + i] for r in rows], dtype=float) for i, _ in enumerate(TISSUES)] + [
        np.array([r[3 + len(TISSUES) + i] for r in rows], dtype=float) for i, _ in enumerate(TISSUES)
    ]
    corr_s = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            corr_s[i, j] = np.corrcoef(_rank(vectors_s[i]), _rank(vectors_s[j]))[0, 1]
    vmin_s, vmax_s = 0.2, 1
    plt.figure(figsize=(8, 6))
    x_s, y_s = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
    im = plt.pcolormesh(x_s, y_s, corr_s, cmap="YlGn", vmin=vmin_s, vmax=vmax_s, shading="auto")
    plt.colorbar(im, label="Spearman correlation")
    tick_pos_s = np.arange(n) + 0.5
    plt.xticks(tick_pos_s, labels, rotation=45, ha="right")
    plt.yticks(tick_pos_s, labels)
    for i in range(n):
        for j in range(n):
            plt.text(j, i, f"{corr_s[i,j]:.2f}", ha="center", va="center", fontsize=6, color="black")
    plt.tight_layout()
    corr_s_pdf = fig_dir / "Mouse_GTEx_mean_spearman_TF_noKidney.svg"
    plt.savefig(corr_s_pdf)
    plt.close()
    print(f"Saved {corr_s_pdf}")

    # Jaccard heatmaps at thresholds 1, 2, 3 using TF means
    def jaccard_heatmap(thresh: float) -> None:
        labels_j = [f"{t}_GTEx" for t in TISSUES] + [f"{t}_Mouse" for t in TISSUES]
        vectors_j = [(human_means_by_tissue[t] > thresh).astype(int) for t in TISSUES] + [
            (mouse_means_by_tissue[t] > thresh).astype(int) for t in TISSUES
        ]
        n_j = len(vectors_j)
        jac = np.zeros((n_j, n_j))
        for i in range(n_j):
            for j in range(n_j):
                inter = np.sum((vectors_j[i] == 1) & (vectors_j[j] == 1))
                union = np.sum((vectors_j[i] == 1) | (vectors_j[j] == 1))
                jac[i, j] = inter / union if union else np.nan
        vmin_j, vmax_j = np.nanmin(jac), np.nanmax(jac)
        plt.figure(figsize=(8, 6))
        im = plt.imshow(jac, cmap="YlGn", vmin=vmin_j, vmax=vmax_j)
        plt.colorbar(im, label="Jaccard similarity")
        plt.xticks(range(n_j), labels_j, rotation=45, ha="right")
        plt.yticks(range(n_j), labels_j)
        for i in range(n_j):
            for j in range(n_j):
                plt.text(j, i, f"{jac[i,j]:.2f}", ha="center", va="center", fontsize=6, color="black")
        plt.tight_layout()
        out_pdf = fig_dir / f"Mouse_GTEx_TF_jaccard_thr{thresh}.pdf"
        plt.savefig(out_pdf)
        plt.close()
        print(f"Saved {out_pdf}")

    for thr in (1, 2, 3):
        jaccard_heatmap(thr)

def _rank(a: np.ndarray) -> np.ndarray:
    temp = a.argsort()
    ranks = np.empty_like(temp, dtype=float)
    ranks[temp] = np.arange(len(a))
    _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
    cumulative = np.cumsum(counts)
    start = cumulative - counts
    avg_rank = (start + cumulative - 1) / 2
    ranks = avg_rank[inv]
    return ranks

if __name__ == "__main__":
    main()
