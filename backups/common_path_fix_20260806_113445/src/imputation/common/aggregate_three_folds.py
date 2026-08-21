
#!/usr/bin/env python3
"""Combine three fold outputs into one 300-gene out-of-fold prediction."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge"),
    )
    parser.add_argument("--experiment", default="ex5")
    parser.add_argument("--experiment-label", default="Experiment 5")
    parser.add_argument("--model", default="vista")
    parser.add_argument("--model-label", default="VISTA")
    parser.add_argument("--fold-dir", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--n-folds", type=int, default=3)
    parser.add_argument("--n-clusters", type=int, default=10)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--seed", type=int, default=8667)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    fold_dir = (
        args.fold_dir.resolve()
        if args.fold_dir is not None
        else project_root
        / "data/processed/imputation_beta/Br8667/gene_folds_200_100"
    )
    output_root = (
        args.output_root.resolve()
        if args.output_root is not None
        else project_root / "outputs/imputation_beta/Br8667"
    )

    common_dir = project_root / "src/imputation_beta/common"
    sys.path.insert(0, str(common_dir))
    from benchmark_evaluation import (
        FOLD_METRIC_COLUMNS,
        calculate_cluster_metrics,
        plot_fold_metric_summary,
        plot_three_fold_metrics,
        to_dense_float32,
    )

    model_root = output_root / args.experiment / args.model
    combined_dir = model_root / "combined"
    figure_dir = combined_dir / "figures"
    combined_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    master_path = fold_dir / "gene_metrics_and_fold_assignment.h5ad"
    if not master_path.exists():
        raise FileNotFoundError(master_path)

    master = ad.read_h5ad(master_path, backed="r")
    master_gene_order = master.var_names.astype(str).tolist()
    master.file.close()

    prediction_parts = []
    native_parts = []
    truth_parts = []
    gene_metric_tables = []
    fold_metric_tables = []
    cell_order = None
    obs_metadata = None
    spatial_coordinates = None

    for fold in range(1, args.n_folds + 1):
        fold_output = model_root / f"fold_{fold}"
        prediction_path = fold_output / "predicted_heldout_genes.h5ad"
        truth_path = fold_dir / f"fold_{fold}_heldout_genes.h5ad"
        gene_metric_path = fold_output / "gene_level_metrics.csv"
        fold_metric_path = fold_output / "fold_level_metrics.csv"
        complete_path = fold_output / "run_complete.flag"

        for path in [
            prediction_path,
            truth_path,
            gene_metric_path,
            fold_metric_path,
            complete_path,
        ]:
            if not path.exists():
                raise FileNotFoundError(path)

        prediction = ad.read_h5ad(prediction_path)
        truth = ad.read_h5ad(truth_path)

        if not prediction.obs_names.equals(truth.obs_names):
            raise ValueError(f"Fold {fold}: prediction/truth cell order differs.")
        if not prediction.var_names.equals(truth.var_names):
            prediction = prediction[:, truth.var_names].copy()
        if prediction.shape != truth.shape:
            raise ValueError(
                f"Fold {fold}: prediction/truth shape differs: "
                f"{prediction.shape} vs {truth.shape}"
            )

        if cell_order is None:
            cell_order = prediction.obs_names.astype(str).tolist()
            obs_metadata = prediction.obs.copy()
            spatial_coordinates = np.asarray(prediction.obsm["spatial"]).copy()
        elif cell_order != prediction.obs_names.astype(str).tolist():
            raise ValueError(f"Fold {fold}: cell order differs from previous folds.")

        prediction_parts.append(
            pd.DataFrame(
                to_dense_float32(prediction.layers["count_scale"]),
                index=prediction.obs_names.astype(str),
                columns=prediction.var_names.astype(str),
            )
        )
        native_parts.append(
            pd.DataFrame(
                to_dense_float32(prediction.layers["native"]),
                index=prediction.obs_names.astype(str),
                columns=prediction.var_names.astype(str),
            )
        )
        truth_parts.append(
            pd.DataFrame(
                to_dense_float32(truth.X),
                index=truth.obs_names.astype(str),
                columns=truth.var_names.astype(str),
            )
        )

        gene_metric_tables.append(pd.read_csv(gene_metric_path))
        fold_metric_tables.append(pd.read_csv(fold_metric_path))

    predicted_table = pd.concat(prediction_parts, axis=1)
    native_table = pd.concat(native_parts, axis=1)
    truth_table = pd.concat(truth_parts, axis=1)

    for table_name, table in [
        ("prediction", predicted_table),
        ("native prediction", native_table),
        ("truth", truth_table),
    ]:
        if table.columns.duplicated().any():
            duplicated = table.columns[table.columns.duplicated()].tolist()
            raise ValueError(f"Duplicated genes in {table_name}: {duplicated[:20]}")
        if set(table.columns) != set(master_gene_order):
            missing = sorted(set(master_gene_order) - set(table.columns))
            extra = sorted(set(table.columns) - set(master_gene_order))
            raise ValueError(
                f"{table_name} does not cover the master 300 genes. "
                f"Missing={missing[:10]}, extra={extra[:10]}"
            )

    predicted_table = predicted_table.loc[cell_order, master_gene_order]
    native_table = native_table.loc[cell_order, master_gene_order]
    truth_table = truth_table.loc[cell_order, master_gene_order]

    oof_prediction = ad.AnnData(
        X=predicted_table.to_numpy(dtype=np.float32),
        obs=obs_metadata.copy(),
        var=pd.DataFrame(index=pd.Index(master_gene_order, name="gene")),
    )
    oof_prediction.layers["count_scale"] = predicted_table.to_numpy(dtype=np.float32)
    oof_prediction.layers["log1p"] = np.log1p(
        predicted_table.to_numpy(dtype=np.float32)
    ).astype(np.float32)
    oof_prediction.layers["native"] = native_table.to_numpy(dtype=np.float32)
    oof_prediction.obsm["spatial"] = spatial_coordinates.copy()
    oof_prediction.uns["benchmark"] = {
        "experiment": args.experiment,
        "experiment_label": args.experiment_label,
        "model": args.model,
        "model_label": args.model_label,
        "type": "300-gene out-of-fold prediction",
        "folds": args.n_folds,
        "output_scale": "nonnegative floating-point expected counts",
    }

    oof_truth = ad.AnnData(
        X=truth_table.to_numpy(dtype=np.float32),
        obs=obs_metadata.copy(),
        var=pd.DataFrame(index=pd.Index(master_gene_order, name="gene")),
    )
    oof_truth.obsm["spatial"] = spatial_coordinates.copy()
    oof_truth.uns["benchmark_role"] = "300-gene out-of-fold ground truth"

    prediction_path = combined_dir / "oof_predictions_300genes.h5ad"
    truth_path = combined_dir / "oof_ground_truth_300genes.h5ad"
    oof_prediction.write_h5ad(prediction_path, compression="gzip")
    oof_truth.write_h5ad(truth_path, compression="gzip")

    gene_metrics = pd.concat(gene_metric_tables, ignore_index=True)
    fold_metrics = pd.concat(fold_metric_tables, ignore_index=True).sort_values("fold")
    gene_metrics.to_csv(combined_dir / "gene_level_metrics_300genes.csv", index=False)
    fold_metrics.to_csv(combined_dir / "fold_level_metrics_3folds.csv", index=False)

    oof_nmi, oof_ari = calculate_cluster_metrics(
        truth_table.to_numpy(dtype=np.float32),
        predicted_table.to_numpy(dtype=np.float32),
        n_clusters=args.n_clusters,
        n_pcs=args.n_pcs,
        seed=args.seed,
    )

    summary = pd.DataFrame(
        [
            {
                "experiment": args.experiment,
                "experiment_label": args.experiment_label,
                "model": args.model,
                "model_label": args.model_label,
                "n_folds": args.n_folds,
                "n_cells": len(cell_order),
                "n_oof_genes": len(master_gene_order),
                "scc_mean": float(gene_metrics["scc"].mean(skipna=True)),
                "scc_sd_across_genes": float(gene_metrics["scc"].std(skipna=True)),
                "ssim_mean": float(gene_metrics["ssim"].mean(skipna=True)),
                "ssim_sd_across_genes": float(gene_metrics["ssim"].std(skipna=True)),
                "rmse_mean": float(gene_metrics["rmse"].mean(skipna=True)),
                "rmse_sd_across_genes": float(gene_metrics["rmse"].std(skipna=True)),
                "mae_mean": float(gene_metrics["mae"].mean(skipna=True)),
                "mae_sd_across_genes": float(gene_metrics["mae"].std(skipna=True)),
                "jsd_mean": float(gene_metrics["jsd"].mean(skipna=True)),
                "jsd_sd_across_genes": float(gene_metrics["jsd"].std(skipna=True)),
                "moran_abs_error_mean": float(
                    gene_metrics["moran_abs_error"].mean(skipna=True)
                ),
                "moran_abs_error_sd_across_genes": float(
                    gene_metrics["moran_abs_error"].std(skipna=True)
                ),
                "nmi": oof_nmi,
                "ari": oof_ari,
                "fold_nmi_mean": float(fold_metrics["nmi"].mean()),
                "fold_nmi_sd": float(fold_metrics["nmi"].std()),
                "fold_ari_mean": float(fold_metrics["ari"].mean()),
                "fold_ari_sd": float(fold_metrics["ari"].std()),
                "training_seconds_total": float(fold_metrics["training_seconds"].sum()),
                "prediction_seconds_total": float(fold_metrics["prediction_seconds"].sum()),
            }
        ]
    )
    summary.to_csv(combined_dir / "model_experiment_summary.csv", index=False)

    # Reuse the same standard labels and plotting functions that future models use.
    plot_three_fold_metrics(
        fold_metrics,
        figure_dir / "three_fold_evaluation_metrics.png",
        title=f"{args.experiment_label} | {args.model_label} | Three-fold evaluation",
    )

    summary_for_plot = summary.copy()
    plot_fold_metric_summary(
        summary_for_plot,
        figure_dir / "oof_300gene_evaluation_metrics.png",
        title=f"{args.experiment_label} | {args.model_label} | 300-gene OOF summary",
    )

    manifest = {
        "experiment": args.experiment,
        "model": args.model,
        "prediction": str(prediction_path.resolve()),
        "truth": str(truth_path.resolve()),
        "gene_metrics": str(
            (combined_dir / "gene_level_metrics_300genes.csv").resolve()
        ),
        "fold_metrics": str(
            (combined_dir / "fold_level_metrics_3folds.csv").resolve()
        ),
        "summary": str(
            (combined_dir / "model_experiment_summary.csv").resolve()
        ),
        "oof_nmi": oof_nmi,
        "oof_ari": oof_ari,
    }
    (combined_dir / "combined_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    (combined_dir / "aggregate_complete.flag").write_text("SUCCESS\n", encoding="utf-8")

    print("=" * 100)
    print("SUCCESS: three VISTA folds aggregated")
    print("OOF prediction:", prediction_path)
    print("Model summary:", combined_dir / "model_experiment_summary.csv")
    print("Three-fold plot:", figure_dir / "three_fold_evaluation_metrics.png")
    print("=" * 100)


if __name__ == "__main__":
    main()
