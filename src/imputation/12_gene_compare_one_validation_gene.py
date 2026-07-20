#!/usr/bin/env python
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def dense_vec(x):
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x)
    if x.ndim == 2:
        x = x[:, 0]
    return x.astype(float)

def zscore(x):
    x = np.asarray(x, dtype=float)
    s = np.nanstd(x)
    if not np.isfinite(s) or s <= 1e-12:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / s

def corr(fun, x, y):
    x = np.asarray(x); y = np.asarray(y)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan
    r = fun(x, y)
    return float(r.statistic if hasattr(r, "statistic") else r[0])

def find_gene(adata, gene):
    names = adata.var_names.astype(str).tolist()
    if gene in names:
        return gene
    lower = {g.lower(): g for g in names}
    return lower.get(gene.lower(), None)

def first(paths):
    for p in paths:
        if p.exists():
            return p
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default="/users/mjabin/projects/GeneBridge")
    ap.add_argument("--gene", default="ACTA2")
    ap.add_argument("--run-mode", choices=["full", "smoke"], default="full")
    args = ap.parse_args()

    root = Path(args.project_root)
    base = root / "data/processed/imputation_beta/Br8667"
    pred_dir = base / "benchmark/predictions"
    out_root = root / "outputs/imputation_beta/Br8667"
    fig_dir = out_root / "validation/method_comparison/figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    report_dir = out_root / "validation/method_comparison"
    report_dir.mkdir(parents=True, exist_ok=True)

    truth_path = base / "benchmark/truth/validation_truth_45genes.h5ad"
    truth = ad.read_h5ad(truth_path)

    priority = [
        args.gene, "ACTA2", "CLDN5", "PECAM1", "VWF", "RGS5", "PDGFRB",
        "MBP", "PLP1", "MOG", "MOBP", "MAG",
        "RORB", "TLE4", "FOXP2", "BCL11B", "CUX2", "CUX1", "SLC17A7",
        "AQP4", "GFAP", "ALDH1L1", "SLC1A3",
        "P2RY12", "CX3CR1", "C3",
    ]

    target = None
    truth_gene = None
    for g in priority:
        gg = find_gene(truth, g)
        if gg is not None:
            target = g.upper()
            truth_gene = gg
            break
    if truth_gene is None:
        raise ValueError("No priority target gene found in validation truth.")

    suffix = "_smoke" if args.run_mode == "smoke" else ""

    pred_files = {
        "VISTA": [
            pred_dir / f"vista_validation_predictions_45genes{suffix}.h5ad",
            pred_dir / "vista_validation_predictions_45genes.h5ad",
            pred_dir / "vista_validation_predictions_45genes_smoke.h5ad",
        ],
        "gimVI": [
            pred_dir / f"gimvi_validation_predictions_45genes{suffix}.h5ad",
            pred_dir / "gimvi_validation_predictions_45genes.h5ad",
            pred_dir / "gimvi_validation_predictions_45genes_smoke.h5ad",
        ],
        "SpaGE": [
            pred_dir / f"spage_validation_predictions_45genes{suffix}.h5ad",
            pred_dir / "spage_validation_predictions_45genes.h5ad",
            pred_dir / "spage_validation_predictions_45genes_smoke.h5ad",
        ],
        "Tangram": [
            pred_dir / f"tangram_validation_predictions_45genes_clusters{suffix}.h5ad",
            pred_dir / "tangram_validation_predictions_45genes_clusters.h5ad",
            pred_dir / "tangram_validation_predictions_45genes_clusters_smoke.h5ad",
            pred_dir / "tangram_validation_predictions_45genes.h5ad",
        ],
    }

    if "spatial" not in truth.obsm:
        raise ValueError("validation truth is missing obsm['spatial'].")

    xy = np.asarray(truth.obsm["spatial"])
    obs_raw = dense_vec(truth[:, truth_gene].X)
    obs_log = np.log1p(np.clip(obs_raw, 0, None))

    panels = [("Observed Xenium", obs_log, np.nan, np.nan)]
    rows = []

    for method in ["VISTA", "gimVI", "SpaGE", "Tangram"]:
        path = first(pred_files[method])
        if path is None:
            print(f"Missing {method} prediction; skipping.")
            continue

        pred = ad.read_h5ad(path)
        pred_gene = find_gene(pred, target)
        if pred_gene is None:
            print(f"{method} file exists but missing {target}; skipping.")
            continue

        common = truth.obs_names.intersection(pred.obs_names)
        obs = dense_vec(truth[common, truth_gene].X)
        pr = dense_vec(pred[common, pred_gene].X)

        obs_log_common = np.log1p(np.clip(obs, 0, None))
        pr_log_common = np.log1p(np.clip(pr, 0, None))

        full = np.full(truth.n_obs, np.nan)
        idx = truth.obs_names.get_indexer(common)
        full[idx] = pr_log_common

        sp = corr(spearmanr, obs_log_common, pr_log_common)
        pe = corr(pearsonr, obs_log_common, pr_log_common)

        panels.append((method, full, sp, pe))
        rows.append({
            "gene": target,
            "method": method,
            "prediction_file": str(path),
            "n_cells": len(common),
            "spearman_log1p": sp,
            "pearson_log1p": pe,
        })

    if len(panels) == 1:
        raise RuntimeError("No method prediction files found. Run methods first.")

    metrics = pd.DataFrame(rows)
    metrics_path = report_dir / f"one_gene_{target}_method_metrics.csv"
    metrics.to_csv(metrics_path, index=False)

    fig, axes = plt.subplots(1, len(panels), figsize=(4.0 * len(panels), 4.5), squeeze=False)
    last_scatter = None

    for i, (title, values, sp, pe) in enumerate(panels):
        ax = axes[0, i]
        zz = np.clip(zscore(values), -2.5, 2.5)
        m = np.isfinite(zz)
        last_scatter = ax.scatter(
            xy[m, 0], xy[m, 1], c=zz[m], s=1, linewidths=0,
            rasterized=True, vmin=-2.5, vmax=2.5
        )
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        if title != "Observed Xenium" and np.isfinite(sp):
            title = f"{title}\nSpearman={sp:.2f}"
        ax.set_title(title, fontsize=11)

    fig.colorbar(last_scatter, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02, label="z-score log1p")
    fig.suptitle(f"Held-out validation gene: {target}", y=1.03, fontsize=15)

    png = fig_dir / f"one_gene_{target}_observed_vs_all_methods.png"
    pdf = fig_dir / f"one_gene_{target}_observed_vs_all_methods.pdf"
    plt.savefig(png, dpi=250, bbox_inches="tight")
    plt.savefig(pdf, bbox_inches="tight")
    plt.close()

    print("Selected gene:", target)
    print("Saved metrics:", metrics_path)
    print("Saved figure:", png)
    print("Saved figure:", pdf)
    print(metrics)

if __name__ == "__main__":
    main()
