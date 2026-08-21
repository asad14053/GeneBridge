
#!/usr/bin/env python3
"""Run one VISTA fold for the GeneBridge 3-fold benchmark on HGCC.

Standard output contract:
- X and layers["count_scale"]: nonnegative floating-point expected counts
- layers["log1p"]: log1p(count_scale)
- layers["native"]: untouched VISTA held-out prediction before count calibration
- exactly the same Xenium cells and 100 held-out genes as the fold truth file
"""

from __future__ import annotations

import argparse
import fcntl
import gc
import json
import os
import random
import re
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import torch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge"),
    )
    parser.add_argument("--experiment", default="ex5")
    parser.add_argument("--experiment-label", default="Experiment 5")
    parser.add_argument("--fold", type=int, required=True, choices=[1, 2, 3])
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--fold-dir", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)

    parser.add_argument(
        "--reference-mode",
        choices=["benchmark_300", "full_reference"],
        default="benchmark_300",
        help=(
            "benchmark_300 keeps the complete 300-gene benchmark panel in the "
            "snRNA reference. full_reference retains every reference gene and "
            "is substantially more expensive."
        ),
    )
    parser.add_argument("--run-mode", choices=["smoke", "full"], default="full")
    parser.add_argument("--smoke-epochs", type=int, default=2)
    parser.add_argument("--full-epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--n-latent", type=int, default=32)
    parser.add_argument("--neighbor-size", type=int, default=20)
    parser.add_argument("--seed", type=int, default=8667)
    parser.add_argument("--cpus", type=int, default=8)

    parser.add_argument("--use-gpu", action="store_true")
    parser.add_argument("--skip-model-save", action="store_true")
    parser.add_argument("--use-pseudocorrelation-filter", action="store_true")
    parser.add_argument("--pseudo-correlation-threshold", type=float, default=0.5)
    parser.add_argument("--pseudo-pvalue-threshold", type=float, default=0.05)
    parser.add_argument("--minimum-selected-genes", type=int, default=10)

    parser.add_argument("--k-spatial", type=int, default=15)
    parser.add_argument("--ssim-grid-size", type=int, default=128)
    parser.add_argument("--n-clusters", type=int, default=10)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--plot-genes", type=int, default=10)

    return parser.parse_args()


def ensure_vista_patch(project_root: Path) -> None:
    """Safely apply the no-validation-loader patch used by the demo notebook."""
    model_file = project_root / "src/imputation/VISTA/vista/_model.py"
    if not model_file.exists():
        print(f"VISTA source patch skipped; file not found: {model_file}")
        return

    lock_path = model_file.with_suffix(".patch.lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        text = model_file.read_text(encoding="utf-8")

        old_block = (
            "            val = ds.val_dataloader()\n"
            "            val_dls.append(val)\n"
            "            val.mode = i"
        )
        new_block = (
            "            val = ds.val_dataloader()\n"
            "            if val is not None:\n"
            "                val_dls.append(val)\n"
            "                val.mode = i"
        )

        if old_block in text:
            backup = model_file.with_suffix(".py.before_vista_no_val_patch")
            if not backup.exists():
                backup.write_text(text, encoding="utf-8")
            model_file.write_text(text.replace(old_block, new_block), encoding="utf-8")
            print(f"Patched VISTA validation loader: {model_file}")
        elif new_block in text:
            print(f"VISTA validation-loader patch already present: {model_file}")
        else:
            print(
                "WARNING: VISTA validation-loader block was not recognized. "
                "Training may fail when validation_size=0."
            )


def mean_expression(adata: ad.AnnData) -> np.ndarray:
    if sparse.issparse(adata.X):
        return np.asarray(adata.X.mean(axis=0)).ravel()
    return np.asarray(adata.X).mean(axis=0)


def calculate_pseudo_correlation(
    adata_sc: ad.AnnData,
    adata_st: ad.AnnData,
    *,
    celltype: str = "scClassify",
    p_value_threshold: float = 0.05,
    correlation_threshold: float = 0.5,
    exclude_types=("Other",),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    import scipy.stats as stats

    if celltype not in adata_sc.obs.columns or celltype not in adata_st.obs.columns:
        raise KeyError(
            f"{celltype!r} is required in both snRNA and Xenium metadata "
            "when pseudocorrelation filtering is enabled."
        )

    overlap = sorted(set(adata_sc.var_names.astype(str)) & set(adata_st.var_names.astype(str)))
    if not overlap:
        raise ValueError("No overlapping genes for pseudocorrelation filtering.")

    sc_overlap = adata_sc[:, overlap]
    st_overlap = adata_st[:, overlap]

    common_types = sorted(
        (
            set(sc_overlap.obs[celltype].astype(str))
            & set(st_overlap.obs[celltype].astype(str))
        )
        - set(exclude_types)
    )
    if len(common_types) < 3:
        raise ValueError(f"At least three common cell types are required; found {common_types}")

    pseudo_sc = []
    pseudo_st = []
    for label in common_types:
        sc_subset = sc_overlap[sc_overlap.obs[celltype].astype(str) == label]
        st_subset = st_overlap[st_overlap.obs[celltype].astype(str) == label]
        pseudo_sc.append(mean_expression(sc_subset))
        pseudo_st.append(mean_expression(st_subset))

    pseudo_sc = np.asarray(pseudo_sc)
    pseudo_st = np.asarray(pseudo_st)

    records = []
    for index, gene in enumerate(overlap):
        x = pseudo_sc[:, index]
        y = pseudo_st[:, index]
        if np.allclose(x, x[0]) or np.allclose(y, y[0]):
            correlation, pvalue = np.nan, np.nan
        else:
            correlation, pvalue = stats.pearsonr(x, y)
        records.append((gene, correlation, pvalue))

    table = pd.DataFrame(records, columns=["gene", "pearson", "pvalue"]).set_index("gene")
    selected = table.loc[
        (table["pvalue"] < p_value_threshold)
        & (table["pearson"] > correlation_threshold)
    ]
    return selected, table


def history_to_table(model) -> pd.DataFrame:
    history = getattr(model, "history", None)
    if history is None or len(history) == 0:
        return pd.DataFrame()

    tables = []
    for metric_name, values in history.items():
        if isinstance(values, pd.DataFrame):
            raw = values.to_numpy(dtype=object).reshape(-1)
        elif isinstance(values, pd.Series):
            raw = values.to_numpy(dtype=object).reshape(-1)
        else:
            raw = np.asarray(values, dtype=object).reshape(-1)

        numeric = []
        for value in raw:
            try:
                if torch.is_tensor(value):
                    if value.numel() != 1:
                        continue
                    value = value.detach().cpu().item()
                value = float(value)
                if np.isfinite(value):
                    numeric.append(value)
            except (TypeError, ValueError):
                continue

        if numeric:
            tables.append(
                pd.DataFrame(
                    {
                        "epoch": np.arange(1, len(numeric) + 1),
                        "metric": str(metric_name),
                        "value": numeric,
                    }
                )
            )

    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def save_training_history(history: pd.DataFrame, output_dir: Path) -> None:
    if history.empty:
        return
    history.to_csv(output_dir / "training_history.csv", index=False)

    for metric, frame in history.groupby("metric", sort=True):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(frame["epoch"], frame["value"], linewidth=2)
        ax.set_xlabel("Recorded epoch")
        ax.set_ylabel(metric.replace("_", " "))
        ax.set_title(f"VISTA {metric.replace('_', ' ')}")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", metric)
        fig.savefig(output_dir / f"training_{safe}.png", dpi=220, bbox_inches="tight")
        plt.close(fig)


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

    output_dir = output_root / args.experiment / "vista" / f"fold_{args.fold}"
    model_dir = output_dir / "model"
    figure_dir = output_dir / "figures"
    for directory in [output_dir, model_dir, figure_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    common_dir = project_root / "src/imputation/common"
    vista_source_root = project_root / "src/imputation/VISTA"
    sys.path.insert(0, str(common_dir))
    if vista_source_root.exists():
        sys.path.insert(0, str(vista_source_root))

    from benchmark_evaluation import (
        assert_nonnegative_count_like,
        build_fold_summary,
        calculate_cluster_metrics,
        calculate_gene_metrics,
        plot_fold_metric_summary,
        plot_gene_metric_distributions,
        plot_ten_gene_maps,
        row_sums,
        select_or_load_plot_genes,
        to_dense_float32,
    )

    ensure_vista_patch(project_root)
    from vista import GIMVI_GCN

    seed = int(args.seed + args.fold)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(max(1, int(args.cpus)))
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.set_float32_matmul_precision("high")

    if args.use_gpu and not torch.cuda.is_available():
        raise RuntimeError("--use-gpu was requested, but CUDA is unavailable.")

    observed_path = fold_dir / f"fold_{args.fold}_observed_genes.h5ad"
    heldout_path = fold_dir / f"fold_{args.fold}_heldout_genes.h5ad"
    master_path = fold_dir / "gene_metrics_and_fold_assignment.h5ad"

    for path in [observed_path, heldout_path, args.reference]:
        if not path.exists():
            raise FileNotFoundError(path)

    spatial_observed = sc.read_h5ad(observed_path)
    spatial_truth = sc.read_h5ad(heldout_path)
    seq_data = sc.read_h5ad(args.reference)

    spatial_observed.var_names = spatial_observed.var_names.astype(str)
    spatial_truth.var_names = spatial_truth.var_names.astype(str)
    seq_data.var_names = seq_data.var_names.astype(str)
    spatial_observed.obs_names = spatial_observed.obs_names.astype(str)
    spatial_truth.obs_names = spatial_truth.obs_names.astype(str)
    seq_data.obs_names = seq_data.obs_names.astype(str)

    if spatial_observed.n_vars != 200:
        raise ValueError(f"Fold {args.fold} observed file has {spatial_observed.n_vars} genes, expected 200.")
    if spatial_truth.n_vars != 100:
        raise ValueError(f"Fold {args.fold} truth file has {spatial_truth.n_vars} genes, expected 100.")
    if not spatial_observed.obs_names.equals(spatial_truth.obs_names):
        raise ValueError("Observed and held-out Xenium cell IDs/order differ.")
    if not set(spatial_observed.var_names).isdisjoint(set(spatial_truth.var_names)):
        raise ValueError("Observed and held-out gene sets overlap.")
    if "spatial" not in spatial_observed.obsm:
        raise KeyError("Observed Xenium input lacks obsm['spatial'].")

    assert_nonnegative_count_like(spatial_observed.X, "Observed Xenium fold input")
    assert_nonnegative_count_like(spatial_truth.X, "Held-out Xenium ground truth")

    observed_genes = spatial_observed.var_names.tolist()
    heldout_genes = spatial_truth.var_names.tolist()
    benchmark_gene_set = set(observed_genes) | set(heldout_genes)

    if len(benchmark_gene_set) != 300:
        raise ValueError(f"Fold union has {len(benchmark_gene_set)} genes, expected 300.")

    if master_path.exists():
        master = ad.read_h5ad(master_path, backed="r")
        master_order = master.var_names.astype(str).tolist()
        master.file.close()
        if set(master_order) != benchmark_gene_set:
            raise ValueError("Master 300-gene file and fold union contain different genes.")
        benchmark_genes = master_order
    else:
        benchmark_genes = observed_genes + heldout_genes

    missing_reference = pd.Index(benchmark_genes).difference(seq_data.var_names)
    if len(missing_reference):
        raise ValueError(
            f"{len(missing_reference)} benchmark genes are absent from the snRNA reference: "
            f"{missing_reference[:20].tolist()}"
        )

    # Keep raw counts in X. The benchmark_300 mode substantially reduces the
    # VISTA decode matrix while retaining all 200 observed and 100 held-out genes.
    if args.reference_mode == "benchmark_300":
        seq_data_model = seq_data[:, benchmark_genes].copy()
    else:
        seq_data_model = seq_data.copy()

    # VISTA setup metadata.
    if "names" not in seq_data_model.obs.columns:
        seq_data_model.obs["names"] = seq_data_model.obs_names.astype(str)
    if "names" not in spatial_observed.obs.columns:
        spatial_observed.obs["names"] = spatial_observed.obs_names.astype(str)
    seq_data_model.obs["ind_x"] = seq_data_model.obs_names.astype(str)
    spatial_observed.obs["ind_x"] = spatial_observed.obs_names.astype(str)
    if "batch" not in spatial_observed.obs.columns:
        spatial_observed.obs["batch"] = "Br8667"

    coordinates = np.asarray(spatial_observed.obsm["spatial"], dtype=np.float32)
    if coordinates.shape[0] != spatial_observed.n_obs or coordinates.shape[1] < 2:
        raise ValueError(f"Invalid spatial coordinate shape: {coordinates.shape}")
    coordinates = coordinates[:, :2]
    spatial_observed.obsm["spatial"] = coordinates
    spatial_truth.obsm["spatial"] = coordinates.copy()

    # For fair model comparison, the default is all 200 observed genes.
    info_genes = observed_genes.copy()
    if args.use_pseudocorrelation_filter:
        selected, all_statistics = calculate_pseudo_correlation(
            seq_data_model,
            spatial_observed,
            p_value_threshold=args.pseudo_pvalue_threshold,
            correlation_threshold=args.pseudo_correlation_threshold,
        )
        all_statistics.to_csv(output_dir / "pseudocorrelation_all_observed_genes.csv")
        selected.to_csv(output_dir / "pseudocorrelation_selected_genes.csv")
        info_genes = selected.index.astype(str).tolist()
        if len(info_genes) < args.minimum_selected_genes:
            raise ValueError(
                f"Only {len(info_genes)} genes passed pseudocorrelation filtering."
            )

    if set(info_genes) & set(heldout_genes):
        raise ValueError("Held-out genes leaked into VISTA spatial input.")
    if not set(info_genes).issubset(set(observed_genes)):
        raise ValueError("VISTA information genes are not a subset of the 200 observed genes.")

    spatial_model = spatial_observed[:, info_genes].copy()

    GIMVI_GCN.setup_anndata(
        spatial_model,
        batch_key="batch",
        obs_names="names",
    )
    GIMVI_GCN.setup_anndata(seq_data_model)

    max_epochs = args.smoke_epochs if args.run_mode == "smoke" else args.full_epochs
    use_gpu = "cuda:0" if args.use_gpu else False

    print("=" * 100)
    print("VISTA 3-fold benchmark")
    print("=" * 100)
    print("Experiment:", args.experiment)
    print("Fold:", args.fold)
    print("Reference:", args.reference)
    print("Reference mode:", args.reference_mode)
    print("Observed Xenium:", spatial_observed.shape)
    print("Held-out truth:", spatial_truth.shape)
    print("VISTA spatial input:", spatial_model.shape)
    print("VISTA snRNA input:", seq_data_model.shape)
    print("Epochs:", max_epochs)
    print("GPU:", use_gpu)
    print("Output:", output_dir)

    model = GIMVI_GCN(
        seq_data_model,
        spatial_model,
        n_latent=args.n_latent,
        neighbor_size=args.neighbor_size,
    )

    training_start = time.time()
    model.train(
        max_epochs=max_epochs,
        train_size=1.0,
        validation_size=0.0,
        use_gpu=use_gpu,
        batch_size=args.batch_size,
    )
    training_seconds = time.time() - training_start

    history = history_to_table(model)
    save_training_history(history, figure_dir)

    if not args.skip_model_save:
        model.save(
            str(model_dir),
            overwrite=True,
            save_anndata=False,
        )

    prediction_start = time.time()
    result = model.get_imputed_values(
        normalized=False,
        batch_size=args.batch_size,
    )
    full_native = result[0] if isinstance(result, (tuple, list)) else result
    full_native = np.asarray(full_native, dtype=np.float32)
    prediction_seconds = time.time() - prediction_start

    expected_shape = (spatial_model.n_obs, seq_data_model.n_vars)
    if full_native.shape != expected_shape:
        raise ValueError(
            f"Unexpected VISTA prediction shape {full_native.shape}; expected {expected_shape}."
        )
    if not np.isfinite(full_native).all() or np.any(full_native < 0):
        raise ValueError("VISTA native prediction contains invalid or negative values.")

    reference_index = pd.Index(seq_data_model.var_names.astype(str))
    observed_positions = reference_index.get_indexer(observed_genes)
    heldout_positions = reference_index.get_indexer(heldout_genes)
    if np.any(observed_positions < 0) or np.any(heldout_positions < 0):
        raise ValueError("Observed or held-out benchmark genes are missing from VISTA output.")

    native_observed = full_native[:, observed_positions]
    native_heldout = full_native[:, heldout_positions]

    # Convert the native nonnegative prediction to Xenium count scale using only
    # the 200 observed genes. Held-out ground truth is never used for calibration.
    observed_count_matrix = to_dense_float32(spatial_observed[:, observed_genes].X)
    observed_library = observed_count_matrix.sum(axis=1)
    native_observed_library = native_observed.sum(axis=1)

    scale_factor = np.divide(
        observed_library,
        native_observed_library,
        out=np.zeros_like(observed_library, dtype=np.float32),
        where=native_observed_library > 1e-12,
    )
    predicted_counts = native_heldout * scale_factor[:, None]
    predicted_counts = np.clip(predicted_counts, 0, None).astype(np.float32, copy=False)

    del full_native, native_observed
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    prediction_adata = ad.AnnData(
        X=predicted_counts,
        obs=spatial_truth.obs.copy(),
        var=spatial_truth.var.copy(),
    )
    prediction_adata.layers["count_scale"] = predicted_counts.copy()
    prediction_adata.layers["log1p"] = np.log1p(predicted_counts).astype(np.float32)
    prediction_adata.layers["native"] = native_heldout.astype(np.float32, copy=False)
    prediction_adata.obsm["spatial"] = coordinates.copy()
    prediction_adata.obs["vista_count_scale_factor"] = scale_factor

    prediction_adata.uns["benchmark"] = {
        "experiment": args.experiment,
        "experiment_label": args.experiment_label,
        "model": "vista",
        "model_label": "VISTA",
        "fold": int(args.fold),
        "observed_gene_count": 200,
        "heldout_gene_count": 100,
        "output_scale": "nonnegative floating-point expected counts",
        "native_output": "VISTA get_imputed_values(normalized=False)",
        "count_scale_calibration": (
            "Per-cell ratio of observed Xenium counts across the 200 observed genes "
            "to VISTA native predictions across the same 200 genes."
        ),
        "reference": str(args.reference.resolve()),
        "reference_mode": args.reference_mode,
        "seed": seed,
        "epochs": int(max_epochs),
        "batch_size": int(args.batch_size),
        "n_latent": int(args.n_latent),
        "neighbor_size": int(args.neighbor_size),
        "pseudocorrelation_filter": bool(args.use_pseudocorrelation_filter),
    }

    prediction_path = output_dir / "predicted_heldout_genes.h5ad"
    prediction_adata.write_h5ad(prediction_path, compression="gzip")

    observed_truth = to_dense_float32(spatial_truth[:, heldout_genes].X)
    gene_metrics = calculate_gene_metrics(
        observed_truth,
        predicted_counts,
        heldout_genes,
        coordinates,
        k_neighbors=args.k_spatial,
        ssim_grid_size=args.ssim_grid_size,
        n_jobs=args.cpus,
    )

    fold_metric_columns = [
        "log_mean_expression",
        "dispersion",
        "detection_fraction",
        "morans_I",
        "heldout_fold",
        "metric_stratum",
    ]
    available_fold_columns = [
        column for column in fold_metric_columns if column in spatial_truth.var.columns
    ]
    if available_fold_columns:
        fold_gene_metadata = spatial_truth.var[available_fold_columns].copy()
        fold_gene_metadata["gene"] = spatial_truth.var_names.astype(str)
        gene_metrics = gene_metrics.merge(fold_gene_metadata, on="gene", how="left")

    gene_metrics.insert(0, "fold", args.fold)
    gene_metrics.insert(0, "model", "vista")
    gene_metrics.insert(0, "experiment", args.experiment)
    gene_metrics.to_csv(output_dir / "gene_level_metrics.csv", index=False)

    nmi, ari = calculate_cluster_metrics(
        observed_truth,
        predicted_counts,
        n_clusters=args.n_clusters,
        n_pcs=args.n_pcs,
        seed=args.seed,
    )

    fold_summary = build_fold_summary(
        gene_metrics,
        experiment=args.experiment,
        model_key="vista",
        model_label="VISTA",
        fold=args.fold,
        n_cells=spatial_truth.n_obs,
        n_observed_genes=spatial_observed.n_vars,
        n_heldout_genes=spatial_truth.n_vars,
        nmi=nmi,
        ari=ari,
        training_seconds=training_seconds,
        prediction_seconds=prediction_seconds,
    )
    fold_summary.to_csv(output_dir / "fold_level_metrics.csv", index=False)

    plot_fold_metric_summary(
        fold_summary,
        figure_dir / "evaluation_metrics.png",
        title=f"{args.experiment_label} | VISTA | Fold {args.fold}",
    )
    plot_gene_metric_distributions(
        gene_metrics,
        figure_dir / "gene_metric_distributions.png",
        title=f"{args.experiment_label} | VISTA | Fold {args.fold} | 100 held-out genes",
    )

    shared_plot_gene_path = (
        fold_dir / "plot_genes" / f"fold_{args.fold}_plot_genes_10.csv"
    )
    selected_plot_genes = select_or_load_plot_genes(
        spatial_truth,
        shared_plot_gene_path,
        n_genes=args.plot_genes,
        seed=args.seed,
    )
    if len(selected_plot_genes) != 10:
        raise ValueError("The fair-comparison spatial figure requires exactly 10 genes.")

    pd.DataFrame(
        {
            "plot_order": np.arange(1, len(selected_plot_genes) + 1),
            "gene": selected_plot_genes,
        }
    ).to_csv(output_dir / "selected_plot_genes.csv", index=False)

    plot_ten_gene_maps(
        observed_truth,
        predicted_counts,
        heldout_genes,
        selected_plot_genes,
        coordinates,
        gene_metrics,
        figure_dir / "ten_gene_observed_vs_prediction.png",
        figure_dir / "ten_gene_observed_vs_prediction.pdf",
        experiment_label=args.experiment_label,
        model_label="VISTA",
        fold=args.fold,
    )

    run_config = {
        "experiment": args.experiment,
        "experiment_label": args.experiment_label,
        "model": "vista",
        "fold": args.fold,
        "reference": str(args.reference.resolve()),
        "observed_input": str(observed_path.resolve()),
        "heldout_truth": str(heldout_path.resolve()),
        "prediction": str(prediction_path.resolve()),
        "output_dir": str(output_dir.resolve()),
        "run_mode": args.run_mode,
        "epochs": max_epochs,
        "batch_size": args.batch_size,
        "n_latent": args.n_latent,
        "neighbor_size": args.neighbor_size,
        "reference_mode": args.reference_mode,
        "use_gpu": bool(args.use_gpu),
        "seed": seed,
        "metric_labels": {
            "scc": "SCC",
            "ssim": "SSIM",
            "rmse": "RMSE",
            "mae": "MAE",
            "jsd": "Jensen–Shannon divergence",
            "moran_abs_error": "Moran's I absolute error",
            "nmi": "NMI",
            "ari": "ARI",
        },
    }
    (output_dir / "run_config.json").write_text(
        json.dumps(run_config, indent=2),
        encoding="utf-8",
    )

    runtime = {
        "training_seconds": training_seconds,
        "prediction_seconds": prediction_seconds,
        "all_zero_predicted_cells": int(np.sum(predicted_counts.sum(axis=1) == 0)),
        "all_zero_predicted_genes": int(np.sum(predicted_counts.sum(axis=0) == 0)),
        "scale_factor_min": float(np.min(scale_factor)),
        "scale_factor_median": float(np.median(scale_factor)),
        "scale_factor_max": float(np.max(scale_factor)),
    }
    (output_dir / "runtime_summary.json").write_text(
        json.dumps(runtime, indent=2),
        encoding="utf-8",
    )

    (output_dir / "run_complete.flag").write_text(
        f"SUCCESS {time.strftime('%Y-%m-%d %H:%M:%S')}\n",
        encoding="utf-8",
    )

    print("=" * 100)
    print("SUCCESS")
    print("Prediction:", prediction_path)
    print("Gene metrics:", output_dir / "gene_level_metrics.csv")
    print("Fold metrics:", output_dir / "fold_level_metrics.csv")
    print("Ten-gene plot:", figure_dir / "ten_gene_observed_vs_prediction.png")
    print("=" * 100)


if __name__ == "__main__":
    main()
