#!/usr/bin/env python
"""
07_check_imputed_output_Br8667.py

Day-1 QC for:
    vista_Br8667_imputed_raw.h5ad

Checks:
1. File dimensions, unique IDs, measured-gene coverage, and spatial coordinates.
2. NaN, infinite, negative, all-zero, and near-zero expression.
3. Separate statistics for the 300 measured Xenium genes and 28,607 unmeasured genes.
4. Reconstruction sanity check for the 300 measured genes.

The measured-gene comparison is not a held-out-gene test because these genes
were available to the full-panel VISTA model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")
DATA_DIR = PROJECT_ROOT / "data/processed/imputation_beta/Br8667"

DEFAULT_IMPUTED = DATA_DIR / "vista_Br8667_imputed_raw.h5ad"
DEFAULT_XENIUM = DATA_DIR / "spatial_data_xenium_Br8667_vista.h5ad"
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs/imputation_beta/Br8667/qc"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--imputed", type=Path, default=DEFAULT_IMPUTED)
    parser.add_argument("--xenium", type=Path, default=DEFAULT_XENIUM)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--chunk-size", type=int, default=128)
    parser.add_argument("--expected-cells", type=int, default=66164)
    parser.add_argument("--expected-genes", type=int, default=28907)
    parser.add_argument("--near-zero-threshold", type=float, default=1e-8)
    return parser.parse_args()


def to_numpy(matrix, dtype=np.float32):
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=dtype)


def safe_pearson(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x)[mask]
    y = np.asarray(y)[mask]
    if x.size < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan
    result = pearsonr(x, y)
    return float(result.statistic if hasattr(result, "statistic") else result[0])


def safe_spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x)[mask]
    y = np.asarray(y)[mask]
    if x.size < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan
    result = spearmanr(x, y)
    return float(result.statistic if hasattr(result, "statistic") else result[0])


def extract_backed_genes(adata, genes):
    """
    Read only selected genes from a backed h5ad.

    h5py cannot fancy-index rows and columns simultaneously, so this reads all
    rows and one sorted gene-position vector, then restores the requested order.
    """
    positions = adata.var_names.get_indexer(genes)

    if np.any(positions < 0):
        missing = [
            genes[i]
            for i, position in enumerate(positions)
            if position < 0
        ]
        raise ValueError(f"Genes missing from imputed file: {missing[:20]}")

    order = np.argsort(positions)
    sorted_positions = positions[order]

    matrix = to_numpy(
        adata[:, sorted_positions].X,
        dtype=np.float32,
    )

    return matrix[:, np.argsort(order)]


def scan_full_matrix(imputed, measured_genes, chunk_size, near_zero_threshold):
    """Calculate full-matrix QC without loading the complete matrix into RAM."""
    n_cells, n_genes = imputed.shape

    sums = np.zeros(n_genes, dtype=np.float64)
    sums_sq = np.zeros(n_genes, dtype=np.float64)
    finite_counts = np.zeros(n_genes, dtype=np.int64)
    positive_counts = np.zeros(n_genes, dtype=np.int64)
    zero_counts = np.zeros(n_genes, dtype=np.int64)
    negative_counts = np.zeros(n_genes, dtype=np.int64)
    invalid_counts = np.zeros(n_genes, dtype=np.int64)
    minima = np.full(n_genes, np.inf)
    maxima = np.full(n_genes, -np.inf)

    bad_cells = []

    total_invalid = 0
    total_negative = 0
    all_zero_cells = 0
    cells_with_invalid = 0
    cells_with_negative = 0

    print("\nScanning the full imputed matrix in chunks...")

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = to_numpy(imputed.X[start:end, :])

        finite = np.isfinite(chunk)
        invalid = ~finite
        negative = finite & (chunk < 0)
        positive = finite & (chunk > 0)
        zero = finite & (chunk == 0)
        safe = np.where(finite, chunk, 0.0)

        sums += safe.sum(axis=0, dtype=np.float64)
        sums_sq += (safe * safe).sum(axis=0, dtype=np.float64)
        finite_counts += finite.sum(axis=0)
        positive_counts += positive.sum(axis=0)
        zero_counts += zero.sum(axis=0)
        negative_counts += negative.sum(axis=0)
        invalid_counts += invalid.sum(axis=0)

        minima = np.minimum(
            minima,
            np.where(finite, chunk, np.inf).min(axis=0),
        )
        maxima = np.maximum(
            maxima,
            np.where(finite, chunk, -np.inf).max(axis=0),
        )

        cell_invalid = invalid.sum(axis=1)
        cell_negative = negative.sum(axis=1)
        cell_positive = positive.sum(axis=1)
        cell_finite = finite.sum(axis=1)

        cell_all_zero = (
            (cell_finite == n_genes)
            & (cell_positive == 0)
            & (cell_negative == 0)
        )

        total_invalid += int(invalid.sum())
        total_negative += int(negative.sum())
        all_zero_cells += int(cell_all_zero.sum())
        cells_with_invalid += int(np.sum(cell_invalid > 0))
        cells_with_negative += int(np.sum(cell_negative > 0))

        problem = cell_all_zero | (cell_invalid > 0) | (cell_negative > 0)

        for local_i in np.flatnonzero(problem):
            global_i = start + int(local_i)
            bad_cells.append(
                {
                    "cell_id": str(imputed.obs_names[global_i]),
                    "all_zero": bool(cell_all_zero[local_i]),
                    "invalid_entry_count": int(cell_invalid[local_i]),
                    "negative_entry_count": int(cell_negative[local_i]),
                    "positive_entry_count": int(cell_positive[local_i]),
                }
            )

        if end == n_cells or (start // chunk_size + 1) % 25 == 0:
            print(f"  Processed {end:,}/{n_cells:,} cells")

    denominator = np.maximum(finite_counts, 1)
    means = sums / denominator
    variances = np.maximum((sums_sq / denominator) - means**2, 0)

    minima[finite_counts == 0] = np.nan
    maxima[finite_counts == 0] = np.nan

    gene_stats = pd.DataFrame(
        {
            "gene": imputed.var_names.astype(str),
            "gene_group": [
                "measured_xenium" if g in measured_genes else "unmeasured_imputed"
                for g in imputed.var_names.astype(str)
            ],
            "finite_value_count": finite_counts,
            "invalid_value_count": invalid_counts,
            "negative_value_count": negative_counts,
            "zero_value_count": zero_counts,
            "positive_value_count": positive_counts,
            "mean_expression": means,
            "variance_expression": variances,
            "minimum_expression": minima,
            "maximum_expression": maxima,
            "positive_fraction": positive_counts / denominator,
            "zero_fraction": zero_counts / denominator,
        }
    )

    gene_stats["is_all_zero"] = (
        gene_stats["finite_value_count"].eq(n_cells)
        & gene_stats["positive_value_count"].eq(0)
        & gene_stats["negative_value_count"].eq(0)
    )
    gene_stats["has_invalid_values"] = gene_stats["invalid_value_count"] > 0
    gene_stats["has_negative_values"] = gene_stats["negative_value_count"] > 0
    gene_stats["is_near_zero_mean"] = (
        gene_stats["mean_expression"] <= near_zero_threshold
    )

    matrix_summary = {
        "invalid_matrix_entries": total_invalid,
        "negative_matrix_entries": total_negative,
        "all_zero_cells": all_zero_cells,
        "cells_with_invalid_values": cells_with_invalid,
        "cells_with_negative_values": cells_with_negative,
        "all_zero_genes": int(gene_stats["is_all_zero"].sum()),
        "near_zero_mean_genes": int(gene_stats["is_near_zero_mean"].sum()),
        "genes_with_invalid_values": int(
            gene_stats["has_invalid_values"].sum()
        ),
        "genes_with_negative_values": int(
            gene_stats["has_negative_values"].sum()
        ),
    }

    bad_cells = pd.DataFrame(
        bad_cells,
        columns=[
            "cell_id",
            "all_zero",
            "invalid_entry_count",
            "negative_entry_count",
            "positive_entry_count",
        ],
    )

    return gene_stats, bad_cells, matrix_summary


def reconstruction_metrics(xenium, imputed):
    """Compare the 300 measured Xenium genes with VISTA predictions."""
    genes = xenium.var_names.astype(str).tolist()

    if not set(genes).issubset(set(imputed.var_names)):
        missing = sorted(set(genes) - set(imputed.var_names))
        raise ValueError(f"Measured genes missing: {missing[:20]}")

    if not set(xenium.obs_names).issubset(set(imputed.obs_names)):
        raise ValueError("Some original Xenium cells are missing.")

    observed = to_numpy(xenium[:, genes].X)
    predicted_in_imputed_order = extract_backed_genes(imputed, genes)

    cell_positions = imputed.obs_names.get_indexer(xenium.obs_names)
    predicted = predicted_in_imputed_order[cell_positions, :]
    predicted = np.clip(predicted, 0, None)

    observed_log = np.log1p(observed)
    predicted_log = np.log1p(predicted)

    rows = []

    for j, gene in enumerate(genes):
        obs = observed[:, j]
        pred = predicted[:, j]
        obs_log = observed_log[:, j]
        pred_log = predicted_log[:, j]

        rows.append(
            {
                "gene": gene,
                "observed_mean": float(obs.mean()),
                "imputed_mean": float(pred.mean()),
                "observed_nonzero_fraction": float(np.mean(obs > 0)),
                "imputed_positive_fraction": float(np.mean(pred > 0)),
                "pearson_log1p_across_cells": safe_pearson(obs_log, pred_log),
                "spearman_across_cells": safe_spearman(obs, pred),
                "mae_log1p": float(np.mean(np.abs(pred_log - obs_log))),
                "rmse_log1p": float(
                    np.sqrt(np.mean((pred_log - obs_log) ** 2))
                ),
            }
        )

    metrics = pd.DataFrame(rows).sort_values(
        "spearman_across_cells",
        ascending=False,
        na_position="last",
    )

    summary = {
        "median_measured_gene_pearson_log1p": metrics[
            "pearson_log1p_across_cells"
        ].median(skipna=True),
        "median_measured_gene_spearman": metrics[
            "spearman_across_cells"
        ].median(skipna=True),
        "median_measured_gene_mae_log1p": metrics[
            "mae_log1p"
        ].median(skipna=True),
        "median_measured_gene_rmse_log1p": metrics[
            "rmse_log1p"
        ].median(skipna=True),
    }

    return metrics, summary


def save_figures(gene_stats, reconstruction, output_dir):
    measured = gene_stats[gene_stats["gene_group"] == "measured_xenium"]
    unmeasured = gene_stats[
        gene_stats["gene_group"] == "unmeasured_imputed"
    ]

    plt.figure(figsize=(8, 5))
    plt.hist(
        np.log1p(unmeasured["mean_expression"].clip(lower=0)),
        bins=60,
        alpha=0.65,
        label="Unmeasured/imputed genes",
    )
    plt.hist(
        np.log1p(measured["mean_expression"].clip(lower=0)),
        bins=40,
        alpha=0.65,
        label="Measured Xenium genes",
    )
    plt.xlabel("log1p mean predicted expression")
    plt.ylabel("Number of genes")
    plt.title("Br8667 predicted-expression distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "07_expression_distribution.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.hist(reconstruction["spearman_across_cells"].dropna(), bins=30)
    plt.xlabel("Spearman correlation across Xenium cells")
    plt.ylabel("Number of measured genes")
    plt.title("Measured-gene reconstruction correlations")
    plt.tight_layout()
    plt.savefig(output_dir / "07_reconstruction_correlations.png", dpi=200)
    plt.close()

    plt.figure(figsize=(6, 6))
    x = np.log1p(reconstruction["observed_mean"].clip(lower=0))
    y = np.log1p(reconstruction["imputed_mean"].clip(lower=0))
    plt.scatter(x, y, s=18, alpha=0.7)
    low = min(float(x.min()), float(y.min()))
    high = max(float(x.max()), float(y.max()))
    plt.plot([low, high], [low, high], linestyle="--")
    plt.xlabel("Observed Xenium mean, log1p")
    plt.ylabel("VISTA-predicted mean, log1p")
    plt.title("Observed versus reconstructed gene means")
    plt.tight_layout()
    plt.savefig(output_dir / "07_observed_vs_imputed_means.png", dpi=200)
    plt.close()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100)
    print("Day 1: Validate Br8667 VISTA-imputed output")
    print("=" * 100)

    for path in [args.imputed, args.xenium]:
        if not path.exists():
            raise FileNotFoundError(path)

    xenium = ad.read_h5ad(args.xenium)
    imputed = ad.read_h5ad(args.imputed, backed="r")

    measured_genes = set(xenium.var_names.astype(str))
    measured_missing = measured_genes - set(imputed.var_names.astype(str))

    has_spatial = "spatial" in imputed.obsm
    spatial_shape = tuple(imputed.obsm["spatial"].shape) if has_spatial else None
    spatial_finite = (
        bool(np.isfinite(np.asarray(imputed.obsm["spatial"])).all())
        if has_spatial
        else False
    )

    structural_checks = {
        "expected_cell_count": imputed.n_obs == args.expected_cells,
        "expected_gene_count": imputed.n_vars == args.expected_genes,
        "unique_cell_ids": imputed.obs_names.is_unique,
        "unique_gene_names": imputed.var_names.is_unique,
        "all_300_measured_genes_present": len(measured_missing) == 0,
        "all_xenium_cells_present": len(
            xenium.obs_names.difference(imputed.obs_names)
        ) == 0,
        "no_extra_imputed_cells": len(
            imputed.obs_names.difference(xenium.obs_names)
        ) == 0,
        "has_spatial_coordinates": has_spatial,
        "spatial_shape_is_n_by_2": spatial_shape == (imputed.n_obs, 2),
        "spatial_coordinates_are_finite": spatial_finite,
    }

    print("\nStructural checks:")
    for name, value in structural_checks.items():
        print(f"{name}: {value}")

    gene_stats, bad_cells, matrix_summary = scan_full_matrix(
        imputed,
        measured_genes,
        args.chunk_size,
        args.near_zero_threshold,
    )

    problem_mask = (
        gene_stats["is_all_zero"]
        | gene_stats["has_invalid_values"]
        | gene_stats["has_negative_values"]
        | gene_stats["is_near_zero_mean"]
    )
    bad_genes = gene_stats.loc[problem_mask].copy()

    reconstruction, reconstruction_summary = reconstruction_metrics(
        xenium,
        imputed,
    )

    measured_stats = gene_stats[
        gene_stats["gene_group"] == "measured_xenium"
    ]
    unmeasured_stats = gene_stats[
        gene_stats["gene_group"] == "unmeasured_imputed"
    ]

    rows = [
        {"section": "structure", "metric": "imputed_cells", "value": imputed.n_obs},
        {"section": "structure", "metric": "imputed_genes", "value": imputed.n_vars},
        {"section": "structure", "metric": "measured_genes", "value": xenium.n_vars},
        {
            "section": "structure",
            "metric": "unmeasured_genes",
            "value": imputed.n_vars - xenium.n_vars,
        },
        {"section": "structure", "metric": "spatial_shape", "value": str(spatial_shape)},
        {
            "section": "metadata",
            "metric": "training_history_verifiable_from_h5ad",
            "value": "No; confirm separately that this file came from the intended completed model.",
        },
    ]

    for name, value in structural_checks.items():
        rows.append(
            {"section": "structural_check", "metric": name, "value": value}
        )

    for name, value in matrix_summary.items():
        rows.append({"section": "matrix_qc", "metric": name, "value": value})

    rows.extend(
        [
            {
                "section": "gene_group",
                "metric": "measured_all_zero_genes",
                "value": int(measured_stats["is_all_zero"].sum()),
            },
            {
                "section": "gene_group",
                "metric": "unmeasured_all_zero_genes",
                "value": int(unmeasured_stats["is_all_zero"].sum()),
            },
            {
                "section": "gene_group",
                "metric": "measured_near_zero_genes",
                "value": int(measured_stats["is_near_zero_mean"].sum()),
            },
            {
                "section": "gene_group",
                "metric": "unmeasured_near_zero_genes",
                "value": int(unmeasured_stats["is_near_zero_mean"].sum()),
            },
        ]
    )

    for name, value in reconstruction_summary.items():
        rows.append(
            {"section": "reconstruction", "metric": name, "value": value}
        )

    summary = pd.DataFrame(rows)

    summary.to_csv(
        args.output_dir / "07_imputed_file_summary.csv",
        index=False,
    )
    gene_stats.to_csv(
        args.output_dir / "07_imputed_gene_statistics.csv",
        index=False,
    )
    bad_genes.to_csv(
        args.output_dir / "07_empty_or_invalid_genes.csv",
        index=False,
    )
    bad_cells.to_csv(
        args.output_dir / "07_invalid_or_empty_cells.csv",
        index=False,
    )
    reconstruction.to_csv(
        args.output_dir / "07_measured_gene_reconstruction.csv",
        index=False,
    )

    save_figures(gene_stats, reconstruction, args.output_dir)

    critical_pass = (
        all(structural_checks.values())
        and matrix_summary["invalid_matrix_entries"] == 0
        and matrix_summary["negative_matrix_entries"] == 0
        and matrix_summary["all_zero_cells"] == 0
    )

    print("\n" + "=" * 100)
    print(
        "PASS: all critical Day-1 checks passed"
        if critical_pass
        else "REVIEW REQUIRED: one or more critical checks failed"
    )
    print("=" * 100)

    print("\nMatrix summary:")
    for name, value in matrix_summary.items():
        print(f"{name}: {value}")

    print("\nReconstruction summary:")
    for name, value in reconstruction_summary.items():
        print(f"{name}: {value}")

    print("\nReports saved to:")
    print(args.output_dir)

    print(
        "\nReminder: reconstruction of the 300 measured genes is not "
        "the final held-out-gene imputation test."
    )

    imputed.file.close()


if __name__ == "__main__":
    main()
