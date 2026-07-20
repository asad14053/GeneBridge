#!/usr/bin/env python
"""
07_compare_vista_Br8667.py

Immediate validation of VISTA-imputed Br8667 expression against the
300 genes measured by Xenium.

Important:
- This is reconstruction validation because these genes were available
  to the model during training.
- A held-out-gene experiment is required for unbiased imputation validation.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")
DATA_DIR = PROJECT_ROOT / "data/processed/imputation_beta/Br8667"
OUT_DIR = PROJECT_ROOT / "outputs/imputation_beta/Br8667/comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)

IMPUTED_PATH = DATA_DIR / "vista_Br8667_imputed_raw.h5ad"
XENIUM_PATH = DATA_DIR / "spatial_data_xenium_Br8667_vista.h5ad"

GENE_METRICS_PATH = OUT_DIR / "07_measured_gene_metrics.csv"
SUMMARY_PATH = OUT_DIR / "07_comparison_summary.csv"
CELLTYPE_PATH = OUT_DIR / "07_scClassify_pseudobulk_metrics.csv"
LAYER_PATH = OUT_DIR / "07_layer_pseudobulk_metrics.csv"
SPEARMAN_FIGURE = OUT_DIR / "07_gene_spearman_histogram.png"
MEAN_FIGURE = OUT_DIR / "07_observed_vs_imputed_gene_means.png"


def to_dense_float32(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if x.size < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan

    result = pearsonr(x, y)
    return float(result.statistic if hasattr(result, "statistic") else result[0])


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if x.size < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan

    result = spearmanr(x, y)
    return float(result.statistic if hasattr(result, "statistic") else result[0])


def pseudobulk_metrics(
    metadata: pd.DataFrame,
    observed: np.ndarray,
    predicted: np.ndarray,
    group_column: str,
    minimum_cells: int = 20,
) -> pd.DataFrame:
    if group_column not in metadata.columns:
        print(f"Skipping {group_column}: column not found.")
        return pd.DataFrame()

    groups = metadata[group_column].astype(str)
    records = []

    for group in sorted(groups.unique()):
        mask = groups.to_numpy() == group
        n_cells = int(mask.sum())

        if n_cells < minimum_cells:
            continue

        observed_mean = observed[mask].mean(axis=0)
        predicted_mean = predicted[mask].mean(axis=0)

        observed_log = np.log1p(observed_mean)
        predicted_log = np.log1p(np.clip(predicted_mean, 0, None))

        records.append(
            {
                "group_column": group_column,
                "group": group,
                "n_cells": n_cells,
                "pearson_across_measured_genes_log1p": safe_pearson(
                    observed_log,
                    predicted_log,
                ),
                "spearman_across_measured_genes": safe_spearman(
                    observed_mean,
                    predicted_mean,
                ),
                "mae_across_measured_genes_log1p": float(
                    np.mean(np.abs(predicted_log - observed_log))
                ),
                "rmse_across_measured_genes_log1p": float(
                    np.sqrt(np.mean((predicted_log - observed_log) ** 2))
                ),
            }
        )

    return pd.DataFrame(records)


def main() -> None:
    for path in [IMPUTED_PATH, XENIUM_PATH]:
        if not path.exists():
            raise FileNotFoundError(path)

    print("Loading original Xenium...")
    xenium = ad.read_h5ad(XENIUM_PATH)

    print("Opening imputed AnnData in backed mode...")
    imputed = ad.read_h5ad(IMPUTED_PATH, backed="r")

    print("\nOriginal Xenium:", xenium.shape)
    print("Imputed matrix:", imputed.shape)

    missing_cells = xenium.obs_names.difference(imputed.obs_names)
    if len(missing_cells) > 0:
        raise ValueError(
            f"{len(missing_cells)} Xenium cells are missing from the imputed file."
        )

    measured_genes = xenium.var_names[
        xenium.var_names.isin(imputed.var_names)
    ].tolist()

    if len(measured_genes) != xenium.n_vars:
        missing_genes = xenium.var_names.difference(imputed.var_names).tolist()
        raise ValueError(
            f"Expected all {xenium.n_vars} Xenium genes in the imputed file, "
            f"but {len(missing_genes)} are missing: {missing_genes[:20]}"
        )

    # Preserve original Xenium cell order and original Xenium gene order.
    cell_order = xenium.obs_names.tolist()
    gene_order = xenium.var_names.tolist()

    print(f"\nAligned cells: {len(cell_order):,}")
    print(f"Measured genes used for comparison: {len(gene_order):,}")

    observed = to_dense_float32(xenium[:, gene_order].X)

    # Load only the 300 measured genes rather than the full 28,907-gene matrix.




    # Load only the 300 measured genes rather than the full 28,907-gene matrix.
    #
    # The imputed AnnData is opened in backed mode. h5py does not permit
    # simultaneous fancy indexing on rows and columns. We therefore:
    #   1. load the selected genes for all cells;
    #   2. restore Xenium gene order in memory;
    #   3. restore Xenium cell order in memory.

    gene_positions = imputed.var_names.get_indexer(gene_order)

    if np.any(gene_positions < 0):
        missing = [
            gene_order[i]
            for i, position in enumerate(gene_positions)
            if position < 0
        ]
        raise ValueError(
            f"Measured genes missing from imputed matrix: {missing[:20]}"
        )

    # h5py requires fancy-index positions to be sorted.
    gene_sort_order = np.argsort(gene_positions)
    sorted_gene_positions = gene_positions[gene_sort_order]

    predicted_sorted = to_dense_float32(
        imputed[:, sorted_gene_positions].X
    )

    # Restore the original Xenium gene order.
    inverse_gene_sort = np.argsort(gene_sort_order)
    predicted_imputed_cell_order = predicted_sorted[:, inverse_gene_sort]

    cell_positions = imputed.obs_names.get_indexer(cell_order)

    if np.any(cell_positions < 0):
        missing = [
            cell_order[i]
            for i, position in enumerate(cell_positions)
            if position < 0
        ]
        raise ValueError(
            f"Xenium cells missing from imputed matrix: {missing[:20]}"
        )

    # Cell reordering now happens in memory.
    predicted = predicted_imputed_cell_order[cell_positions, :]
    predicted = np.clip(predicted, 0, None)



    if observed.shape != predicted.shape:
        raise ValueError(
            f"Observed and predicted shapes differ: "
            f"{observed.shape} versus {predicted.shape}"
        )

    observed_log = np.log1p(observed)
    predicted_log = np.log1p(predicted)

    records = []

    for j, gene in enumerate(gene_order):
        obs_gene = observed[:, j]
        pred_gene = predicted[:, j]
        obs_log_gene = observed_log[:, j]
        pred_log_gene = predicted_log[:, j]

        records.append(
            {
                "gene": gene,
                "n_cells": observed.shape[0],
                "observed_mean": float(obs_gene.mean()),
                "imputed_mean": float(pred_gene.mean()),
                "observed_variance": float(obs_gene.var()),
                "imputed_variance": float(pred_gene.var()),
                "observed_nonzero_fraction": float(np.mean(obs_gene > 0)),
                "imputed_positive_fraction": float(np.mean(pred_gene > 0)),
                "pearson_log1p_across_cells": safe_pearson(
                    obs_log_gene,
                    pred_log_gene,
                ),
                "spearman_across_cells": safe_spearman(
                    obs_gene,
                    pred_gene,
                ),
                "mae_log1p": float(
                    np.mean(np.abs(pred_log_gene - obs_log_gene))
                ),
                "rmse_log1p": float(
                    np.sqrt(np.mean((pred_log_gene - obs_log_gene) ** 2))
                ),
                "mean_bias_log1p": float(
                    np.mean(pred_log_gene - obs_log_gene)
                ),
            }
        )

    gene_metrics = pd.DataFrame(records).sort_values(
        "spearman_across_cells",
        ascending=False,
        na_position="last",
    )
    gene_metrics.to_csv(GENE_METRICS_PATH, index=False)

    observed_gene_means = observed.mean(axis=0)
    predicted_gene_means = predicted.mean(axis=0)

    global_pseudobulk_pearson = safe_pearson(
        np.log1p(observed_gene_means),
        np.log1p(predicted_gene_means),
    )
    global_pseudobulk_spearman = safe_spearman(
        observed_gene_means,
        predicted_gene_means,
    )

    summary = pd.DataFrame(
        [
            {"metric": "n_cells", "value": observed.shape[0]},
            {"metric": "n_measured_genes", "value": observed.shape[1]},
            {
                "metric": "median_gene_pearson_log1p_across_cells",
                "value": gene_metrics[
                    "pearson_log1p_across_cells"
                ].median(skipna=True),
            },
            {
                "metric": "median_gene_spearman_across_cells",
                "value": gene_metrics[
                    "spearman_across_cells"
                ].median(skipna=True),
            },
            {
                "metric": "median_gene_mae_log1p",
                "value": gene_metrics["mae_log1p"].median(skipna=True),
            },
            {
                "metric": "median_gene_rmse_log1p",
                "value": gene_metrics["rmse_log1p"].median(skipna=True),
            },
            {
                "metric": "global_pseudobulk_pearson_log1p_across_genes",
                "value": global_pseudobulk_pearson,
            },
            {
                "metric": "global_pseudobulk_spearman_across_genes",
                "value": global_pseudobulk_spearman,
            },
        ]
    )
    summary.to_csv(SUMMARY_PATH, index=False)

    celltype_metrics = pseudobulk_metrics(
        xenium.obs,
        observed,
        predicted,
        group_column="scClassify",
    )
    if not celltype_metrics.empty:
        celltype_metrics.to_csv(CELLTYPE_PATH, index=False)

    layer_metrics = pseudobulk_metrics(
        xenium.obs,
        observed,
        predicted,
        group_column="layer_annotation",
    )
    if not layer_metrics.empty:
        layer_metrics.to_csv(LAYER_PATH, index=False)

    # Figure 1: distribution of per-gene Spearman correlations.
    plt.figure(figsize=(7, 5))
    values = gene_metrics["spearman_across_cells"].dropna()
    plt.hist(values, bins=30)
    plt.xlabel("Spearman correlation across Xenium cells")
    plt.ylabel("Number of measured genes")
    plt.title("VISTA reconstruction of measured Xenium genes")
    plt.tight_layout()
    plt.savefig(SPEARMAN_FIGURE, dpi=200)
    plt.close()

    # Figure 2: measured-gene pseudobulk means.
    plt.figure(figsize=(6, 6))
    plt.scatter(
        np.log1p(observed_gene_means),
        np.log1p(predicted_gene_means),
        s=18,
        alpha=0.7,
    )
    low = min(
        float(np.log1p(observed_gene_means).min()),
        float(np.log1p(predicted_gene_means).min()),
    )
    high = max(
        float(np.log1p(observed_gene_means).max()),
        float(np.log1p(predicted_gene_means).max()),
    )
    plt.plot([low, high], [low, high], linestyle="--")
    plt.xlabel("Observed Xenium mean, log1p")
    plt.ylabel("VISTA-imputed mean, log1p")
    plt.title("Observed versus imputed measured-gene means")
    plt.tight_layout()
    plt.savefig(MEAN_FIGURE, dpi=200)
    plt.close()

    print("\nComparison complete.")
    print("\nSummary:")
    print(summary.to_string(index=False))

    print("\nTop 10 reconstructed genes:")
    print(
        gene_metrics[
            [
                "gene",
                "spearman_across_cells",
                "pearson_log1p_across_cells",
                "rmse_log1p",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )

    print("\nBottom 10 reconstructed genes:")
    print(
        gene_metrics[
            [
                "gene",
                "spearman_across_cells",
                "pearson_log1p_across_cells",
                "rmse_log1p",
            ]
        ]
        .tail(10)
        .to_string(index=False)
    )

    print("\nCreated:")
    for path in [
        GENE_METRICS_PATH,
        SUMMARY_PATH,
        CELLTYPE_PATH,
        LAYER_PATH,
        SPEARMAN_FIGURE,
        MEAN_FIGURE,
    ]:
        if path.exists():
            print(path)

    imputed.file.close()


if __name__ == "__main__":
    main()
