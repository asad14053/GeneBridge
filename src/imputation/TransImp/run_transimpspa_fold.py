#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gc
import json
import os
import random
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch

from scipy import sparse
from sklearn.neighbors import NearestNeighbors


MODEL_KEY = "transimpspa"
MODEL_LABEL = "TransImpSpa"


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--project-root",
        type=Path,
        default=Path(
            "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
        ),
    )

    p.add_argument("--experiment", required=True)
    p.add_argument("--experiment-label", required=True)

    p.add_argument(
        "--fold",
        type=int,
        required=True,
        choices=[1, 2, 3],
    )

    p.add_argument(
        "--reference",
        type=Path,
        required=True,
    )

    p.add_argument(
        "--fold-dir",
        type=Path,
        default=None,
    )

    p.add_argument(
        "--output-root",
        type=Path,
        default=None,
    )

    # Paper/default TransImpSpa settings.
    p.add_argument("--mapping-lowdim", type=int, default=256)
    p.add_argument("--epochs", type=int, default=2000)
    p.add_argument("--lr", type=float, default=0.01)
    p.add_argument("--weight-decay", type=float, default=0.01)
    p.add_argument("--clip-max", type=float, default=10.0)
    p.add_argument("--wt-spa", type=float, default=1.0)

    # Squidpy generic spatial_neighbors default.
    p.add_argument("--spatial-neighbors", type=int, default=6)

    p.add_argument("--seed", type=int, default=8667)
    p.add_argument("--cpus", type=int, default=16)

    # Shared benchmark evaluator.
    p.add_argument("--k-spatial", type=int, default=15)
    p.add_argument("--ssim-grid-size", type=int, default=128)
    p.add_argument("--n-clusters", type=int, default=10)
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--plot-genes", type=int, default=10)

    return p.parse_args()


def dense_float32(X):
    if sparse.issparse(X):
        X = X.toarray()

    return np.asarray(
        X,
        dtype=np.float32,
    )


def assert_finite_nonnegative(X, name):
    if sparse.issparse(X):
        values = X.data
    else:
        values = np.asarray(X).ravel()

    if values.size == 0:
        return

    if not np.isfinite(values).all():
        raise FloatingPointError(
            f"{name} contains NaN/Inf."
        )

    minimum = float(values.min())

    if minimum < -1e-7:
        raise ValueError(
            f"{name} contains negative values; min={minimum}"
        )


def normalized_reference(
    reference_counts,
    observed_indices,
    target_sum=10000.0,
):
    """
    Normalize all 300 reference genes using only the
    200 observed genes to compute each reference-cell library size.
    """

    X = np.asarray(
        reference_counts,
        dtype=np.float64,
    )

    library = X[
        :,
        observed_indices,
    ].sum(axis=1)

    keep = (
        np.isfinite(library)
        & (library > 0)
    )

    X = X[keep]
    library = library[keep]

    X_norm = np.log1p(
        X / library[:, None] * target_sum
    )

    if not np.isfinite(X_norm).all():
        raise FloatingPointError(
            "Reference normalization produced NaN/Inf."
        )

    return (
        X_norm.astype(np.float32),
        keep,
    )


def normalized_spatial(
    observed_counts,
    target_sum=10000.0,
):
    X = np.asarray(
        observed_counts,
        dtype=np.float64,
    )

    library = X.sum(axis=1)

    if (
        not np.isfinite(library).all()
        or np.any(library <= 0)
    ):
        raise ValueError(
            "Spatial observed-gene library contains "
            "zero/non-finite cells."
        )

    X_norm = np.log1p(
        X / library[:, None] * target_sum
    )

    if not np.isfinite(X_norm).all():
        raise FloatingPointError(
            "Spatial normalization produced NaN/Inf."
        )

    return (
        X_norm.astype(np.float32),
        library.astype(np.float64),
    )


def make_spatial_knn(
    coordinates,
    n_neighbors=6,
):
    """
    Squidpy-equivalent generic kNN connectivity:
    physical Euclidean coordinates, no self edge.

    Symmetrize so an edge is retained if either observation
    considers the other a nearest neighbor.
    """

    nn = NearestNeighbors(
        n_neighbors=n_neighbors + 1,
        metric="euclidean",
        n_jobs=-1,
    )

    nn.fit(coordinates)

    indices = nn.kneighbors(
        coordinates,
        return_distance=False,
    )[:, 1:]

    n = coordinates.shape[0]

    rows = np.repeat(
        np.arange(n),
        n_neighbors,
    )

    cols = indices.reshape(-1)

    data = np.ones(
        rows.shape[0],
        dtype=np.float32,
    )

    adjacency = sparse.coo_matrix(
        (data, (rows, cols)),
        shape=(n, n),
        dtype=np.float32,
    ).tocsr()

    adjacency = adjacency.maximum(
        adjacency.T
    )

    adjacency.setdiag(0)
    adjacency.eliminate_zeros()

    return adjacency.tocoo()


def extract_prediction(result, n_cells, n_genes):
    """
    Current documentation returns prediction directly when
    uncertainty simulation is disabled.

    This also tolerates package versions returning the prediction
    as element 0 of a tuple/list.
    """

    candidates = []

    if isinstance(result, (tuple, list)):
        candidates.extend(result)
    else:
        candidates.append(result)

    for candidate in candidates:
        try:
            array = np.asarray(candidate)
        except Exception:
            continue

        if array.shape == (n_cells, n_genes):
            return array.astype(
                np.float32,
                copy=False,
            )

    shapes = []

    for candidate in candidates:
        try:
            shapes.append(
                np.asarray(candidate).shape
            )
        except Exception:
            shapes.append(
                str(type(candidate))
            )

    raise ValueError(
        "Could not identify TransImpSpa prediction. "
        f"Expected {(n_cells, n_genes)}, got {shapes}"
    )


def main():
    args = parse_args()

    os.environ["OMP_NUM_THREADS"] = str(args.cpus)
    os.environ["MKL_NUM_THREADS"] = str(args.cpus)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.cpus)
    os.environ["NUMEXPR_NUM_THREADS"] = str(args.cpus)

    torch.set_num_threads(
        args.cpus
    )

    project_root = (
        args.project_root.resolve()
    )

    fold_dir = (
        args.fold_dir.resolve()
        if args.fold_dir
        else (
            project_root
            / "data/processed/imputation_beta/Br8667"
            / "gene_folds_200_100"
        )
    )

    output_root = (
        args.output_root.resolve()
        if args.output_root
        else (
            project_root
            / "outputs/imputation_beta/Br8667"
        )
    )

    output_dir = (
        output_root
        / args.experiment
        / MODEL_KEY
        / f"fold_{args.fold}"
    )

    figure_dir = output_dir / "figures"
    diagnostics_dir = output_dir / "diagnostics"

    for directory in [
        output_dir,
        figure_dir,
        diagnostics_dir,
    ]:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    common_dir = (
        project_root
        / "src/imputation/common"
    )

    sys.path.insert(
        0,
        str(common_dir),
    )

    from benchmark_evaluation import (
        assert_nonnegative_count_like,
        build_fold_summary,
        calculate_cluster_metrics,
        calculate_gene_metrics,
        plot_fold_metric_summary,
        plot_gene_metric_distributions,
        plot_ten_gene_maps,
        select_or_load_plot_genes,
        to_dense_float32,
    )

    transpa_root = (
        project_root
        / "src/imputation/TransImp/tranSpa"
    )

    if not (
        transpa_root
        / "transpa/util.py"
    ).is_file():
        raise FileNotFoundError(
            f"TranSpa repository missing: {transpa_root}"
        )

    sys.path.insert(
        0,
        str(transpa_root),
    )

    from transpa.util import expTransImp

    seed = args.seed + args.fold

    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)

    observed_path = (
        fold_dir
        / f"fold_{args.fold}_observed_genes.h5ad"
    )

    heldout_path = (
        fold_dir
        / f"fold_{args.fold}_heldout_genes.h5ad"
    )

    spatial_obs = sc.read_h5ad(
        observed_path
    )

    spatial_truth = sc.read_h5ad(
        heldout_path
    )

    reference = sc.read_h5ad(
        args.reference
    )

    for x in [
        spatial_obs,
        spatial_truth,
        reference,
    ]:
        x.obs_names = (
            x.obs_names.astype(str)
        )

        x.var_names = (
            x.var_names.astype(str)
        )

    if spatial_obs.n_vars != 200:
        raise ValueError(
            f"Expected 200 observed genes; got {spatial_obs.n_vars}"
        )

    if spatial_truth.n_vars != 100:
        raise ValueError(
            f"Expected 100 held-out genes; got {spatial_truth.n_vars}"
        )

    if not spatial_obs.obs_names.equals(
        spatial_truth.obs_names
    ):
        raise ValueError(
            "Observed and held-out Xenium cell order differs."
        )

    observed_genes = (
        spatial_obs.var_names.tolist()
    )

    heldout_genes = (
        spatial_truth.var_names.tolist()
    )

    if not set(observed_genes).isdisjoint(
        heldout_genes
    ):
        raise ValueError(
            "Observed and held-out genes overlap."
        )

    benchmark_genes = (
        observed_genes
        + heldout_genes
    )

    if len(set(benchmark_genes)) != 300:
        raise ValueError(
            "Fold union is not exactly 300 unique genes."
        )

    missing_reference = (
        pd.Index(benchmark_genes)
        .difference(
            reference.var_names
        )
    )

    if len(missing_reference):
        raise ValueError(
            "Benchmark genes missing from reference: "
            f"{missing_reference[:20].tolist()}"
        )

    if "spatial" not in spatial_obs.obsm:
        raise KeyError(
            "Xenium input lacks obsm['spatial']."
        )

    assert_nonnegative_count_like(
        spatial_obs.X,
        "TransImpSpa observed Xenium counts",
    )

    assert_nonnegative_count_like(
        spatial_truth.X,
        "TransImpSpa held-out Xenium truth",
    )

    reference_300 = dense_float32(
        reference[
            :,
            benchmark_genes,
        ].X
    )

    spatial_200 = dense_float32(
        spatial_obs[
            :,
            observed_genes,
        ].X
    )

    assert_finite_nonnegative(
        reference_300,
        "Reference benchmark counts",
    )

    assert_finite_nonnegative(
        spatial_200,
        "Spatial observed counts",
    )

    observed_indices = np.arange(
        0,
        200,
        dtype=int,
    )

    (
        reference_norm,
        reference_keep,
    ) = normalized_reference(
        reference_300,
        observed_indices,
        target_sum=10000.0,
    )

    (
        spatial_norm,
        xenium_observed_library,
    ) = normalized_spatial(
        spatial_200,
        target_sum=10000.0,
    )

    n_removed_reference = int(
        (~reference_keep).sum()
    )

    if reference_norm.shape[0] < 10:
        raise ValueError(
            "Too few reference cells after observed-gene filtering."
        )

    # Strict fold preservation: all 200 observed genes must be usable.
    ref_train_var = np.var(
        reference_norm[:, :200],
        axis=0,
    )

    spa_train_var = np.var(
        spatial_norm,
        axis=0,
    )

    bad = (
        ~np.isfinite(ref_train_var)
        | ~np.isfinite(spa_train_var)
        | (ref_train_var <= 1e-12)
        | (spa_train_var <= 1e-12)
    )

    if np.any(bad):
        bad_genes = [
            observed_genes[i]
            for i in np.where(bad)[0]
        ]

        raise ValueError(
            "TransImpSpa requires all 200 benchmark training genes; "
            f"zero/non-finite variance genes: {bad_genes}"
        )

    coordinates = np.asarray(
        spatial_obs.obsm["spatial"],
        dtype=np.float32,
    )[:, :2]

    if (
        coordinates.shape
        != (spatial_obs.n_obs, 2)
    ):
        raise ValueError(
            f"Unexpected spatial shape: {coordinates.shape}"
        )

    if not np.isfinite(
        coordinates
    ).all():
        raise ValueError(
            "Spatial coordinates contain NaN/Inf."
        )

    adjacency = make_spatial_knn(
        coordinates,
        n_neighbors=args.spatial_neighbors,
    )

    sparse.save_npz(
        diagnostics_dir
        / "spatial_connectivities_6nn.npz",
        adjacency.tocsr(),
    )

    pd.DataFrame(
        {
            "cell_id":
                spatial_obs.obs_names,

            "observed_200_library":
                xenium_observed_library,
        }
    ).to_csv(
        diagnostics_dir
        / "observed_library_sizes.csv",
        index=False,
    )

    reference_df = pd.DataFrame(
        reference_norm,
        columns=benchmark_genes,
    )

    spatial_df = pd.DataFrame(
        spatial_norm,
        index=spatial_obs.obs_names,
        columns=observed_genes,
    )

    device = torch.device(
        "cuda:0"
        if torch.cuda.is_available()
        else "cpu"
    )

    print("=" * 100)
    print("STANDARD TransImpSpa 3-FOLD BENCHMARK")
    print("=" * 100)

    print("Experiment:", args.experiment)
    print("Fold:", args.fold)
    print("Reference raw shape:", reference.shape)
    print("Reference model shape:", reference_df.shape)
    print("Removed zero-observed reference cells:", n_removed_reference)
    print("Spatial model shape:", spatial_df.shape)
    print("Held-out truth:", spatial_truth.shape)

    print("Train genes:", len(observed_genes))
    print("Test genes:", len(heldout_genes))

    print("signature_mode: cell")
    print("mapping_mode: lowrank")
    print("mapping_lowdim:", args.mapping_lowdim)

    print("Spatial neighbors:", args.spatial_neighbors)
    print("Adjacency nnz:", adjacency.nnz)

    print("wt_spa:", args.wt_spa)
    print("epochs:", args.epochs)
    print("lr:", args.lr)
    print("weight_decay:", args.weight_decay)
    print("clip_max:", args.clip_max)

    print("Device:", device)
    print("Seed:", seed)

    print("=" * 100)

    start = time.time()

    result = expTransImp(
        df_ref=reference_df,
        df_tgt=spatial_df,
        train_gene=observed_genes,
        test_gene=heldout_genes,

        signature_mode="cell",
        mapping_mode="lowrank",
        mapping_lowdim=args.mapping_lowdim,

        spa_adj=adjacency,

        lr=args.lr,
        weight_decay=args.weight_decay,
        n_epochs=args.epochs,
        clip_max=args.clip_max,
        wt_spa=args.wt_spa,

        n_simulation=None,

        device=device,
        seed=seed,
    )

    prediction_seconds = (
        time.time() - start
    )

    native_prediction = extract_prediction(
        result,
        spatial_obs.n_obs,
        len(heldout_genes),
    )

    assert_finite_nonnegative(
        native_prediction,
        "TransImpSpa native prediction",
    )

    native_prediction = np.maximum(
        native_prediction,
        0,
    )

    # Native output is on the same log-normalized scale used
    # for the fitted translation.
    normalized_abundance = np.expm1(
        native_prediction.astype(
            np.float64
        )
    )

    if (
        not np.isfinite(
            normalized_abundance
        ).all()
    ):
        raise FloatingPointError(
            "expm1(native TransImpSpa) produced NaN/Inf."
        )

    normalized_abundance = np.maximum(
        normalized_abundance,
        0,
    )

    # Reverse our leakage-safe normalization.
    predicted_counts = (
        normalized_abundance
        * (
            xenium_observed_library[:, None]
            / 10000.0
        )
    ).astype(
        np.float32
    )

    assert_finite_nonnegative(
        predicted_counts,
        "TransImpSpa count-scale predictions",
    )

    truth_counts = to_dense_float32(
        spatial_truth[
            :,
            heldout_genes,
        ].X
    )

    prediction_adata = ad.AnnData(
        X=predicted_counts,
        obs=spatial_truth.obs.copy(),
        var=spatial_truth.var.copy(),
    )

    prediction_adata.layers[
        "count_scale"
    ] = predicted_counts.copy()

    prediction_adata.layers[
        "log1p"
    ] = np.log1p(
        predicted_counts
    ).astype(np.float32)

    prediction_adata.layers[
        "native"
    ] = native_prediction.copy()

    prediction_adata.obsm[
        "spatial"
    ] = coordinates.copy()

    prediction_adata.uns[
        "benchmark"
    ] = {
        "model":
            MODEL_KEY,

        "model_label":
            MODEL_LABEL,

        "variant":
            "TransImpSpa",

        "signature_mode":
            "cell",

        "mapping_mode":
            "lowrank",

        "mapping_lowdim":
            int(args.mapping_lowdim),

        "wt_spa":
            float(args.wt_spa),

        "spatial_neighbors":
            int(args.spatial_neighbors),

        "epochs":
            int(args.epochs),

        "lr":
            float(args.lr),

        "weight_decay":
            float(args.weight_decay),

        "clip_max":
            float(args.clip_max),

        "normalization":
            (
                "Both modalities log1p(count / "
                "200-observed-gene library * 10000); "
                "reference held-out genes use the same "
                "reference observed-gene size factor."
            ),

        "heldout_truth_used_for_training":
            False,

        "heldout_truth_used_for_normalization":
            False,

        "device":
            str(device),

        "seed":
            int(seed),
    }

    prediction_path = (
        output_dir
        / "predicted_heldout_genes.h5ad"
    )

    prediction_adata.write_h5ad(
        prediction_path,
        compression="gzip",
    )

    gene_metrics = calculate_gene_metrics(
        truth_counts,
        predicted_counts,
        heldout_genes,
        coordinates,

        k_neighbors=args.k_spatial,
        ssim_grid_size=args.ssim_grid_size,
        n_jobs=args.cpus,
    )

    gene_metrics.insert(
        0,
        "fold",
        args.fold,
    )

    gene_metrics.insert(
        0,
        "model",
        MODEL_KEY,
    )

    gene_metrics.insert(
        0,
        "experiment",
        args.experiment,
    )

    gene_metrics.to_csv(
        output_dir
        / "gene_level_metrics.csv",
        index=False,
    )

    nmi, ari = calculate_cluster_metrics(
        truth_counts,
        predicted_counts,

        n_clusters=args.n_clusters,
        n_pcs=args.n_pcs,
        seed=args.seed,
    )

    fold_summary = build_fold_summary(
        gene_metrics,

        experiment=args.experiment,
        model_key=MODEL_KEY,
        model_label=MODEL_LABEL,
        fold=args.fold,

        n_cells=spatial_truth.n_obs,
        n_observed_genes=200,
        n_heldout_genes=100,

        nmi=nmi,
        ari=ari,

        training_seconds=prediction_seconds,
        prediction_seconds=0.0,
    )

    fold_summary.to_csv(
        output_dir
        / "fold_level_metrics.csv",
        index=False,
    )

    plot_fold_metric_summary(
        fold_summary,

        figure_dir
        / "evaluation_metrics.png",

        title=(
            f"{args.experiment_label} "
            f"| TransImpSpa | Fold {args.fold}"
        ),
    )

    plot_gene_metric_distributions(
        gene_metrics,

        figure_dir
        / "gene_metric_distributions.png",

        title=(
            f"{args.experiment_label} "
            f"| TransImpSpa | Fold {args.fold}"
        ),
    )

    plot_gene_path = (
        fold_dir
        / "plot_genes"
        / f"fold_{args.fold}_plot_genes_10.csv"
    )

    selected_plot_genes = (
        select_or_load_plot_genes(
            spatial_truth,
            plot_gene_path,
            n_genes=args.plot_genes,
            seed=args.seed,
        )
    )

    if len(selected_plot_genes) != 10:
        raise ValueError(
            "Expected exactly 10 shared visualization genes."
        )

    pd.DataFrame(
        {
            "plot_order":
                np.arange(1, 11),

            "gene":
                selected_plot_genes,
        }
    ).to_csv(
        output_dir
        / "selected_plot_genes.csv",
        index=False,
    )

    plot_ten_gene_maps(
        truth_counts,
        predicted_counts,
        heldout_genes,
        selected_plot_genes,
        coordinates,
        gene_metrics,

        figure_dir
        / "ten_gene_observed_vs_prediction.png",

        figure_dir
        / "ten_gene_observed_vs_prediction.pdf",

        experiment_label=args.experiment_label,
        model_label=MODEL_LABEL,
        fold=args.fold,
    )

    config = {
        "model":
            MODEL_KEY,

        "model_label":
            MODEL_LABEL,

        "experiment":
            args.experiment,

        "fold":
            args.fold,

        "reference":
            str(
                args.reference.resolve()
            ),

        "n_cells":
            int(spatial_obs.n_obs),

        "n_reference_cells_used":
            int(reference_df.shape[0]),

        "n_reference_cells_removed":
            n_removed_reference,

        "train_genes":
            200,

        "test_genes":
            100,

        "signature_mode":
            "cell",

        "mapping_mode":
            "lowrank",

        "mapping_lowdim":
            args.mapping_lowdim,

        "wt_spa":
            args.wt_spa,

        "spatial_neighbors":
            args.spatial_neighbors,

        "epochs":
            args.epochs,

        "lr":
            args.lr,

        "weight_decay":
            args.weight_decay,

        "clip_max":
            args.clip_max,

        "prediction_seconds":
            prediction_seconds,

        "device":
            str(device),

        "seed":
            seed,
    }

    with (
        output_dir
        / "run_config.json"
    ).open(
        "w"
    ) as f:
        json.dump(
            config,
            f,
            indent=2,
        )

    success = (
        "SUCCESS\n"
        f"model={MODEL_KEY}\n"
        f"experiment={args.experiment}\n"
        f"fold={args.fold}\n"
    )

    (
        output_dir
        / "complete.flag"
    ).write_text(success)

    (
        output_dir
        / "run_complete.flag"
    ).write_text(success)

    print()
    print("=" * 100)
    print("SUCCESS")
    print("Model:", MODEL_LABEL)
    print("Experiment:", args.experiment)
    print("Fold:", args.fold)
    print("Prediction:", prediction_path)
    print(
        "Runtime minutes:",
        f"{prediction_seconds / 60:.2f}",
    )
    print(
        "run_complete.flag:",
        output_dir / "run_complete.flag",
    )
    print("=" * 100)

    gc.collect()


if __name__ == "__main__":
    main()
