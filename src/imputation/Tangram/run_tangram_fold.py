#!/usr/bin/env python3
"""Run one Tangram fold for the GeneBridge 3-fold imputation benchmark on HGCC.

Output contract:
- X and layers["count_scale"]: nonnegative floating-point expected Xenium counts
- layers["log1p"]: log1p(count_scale)
- layers["native"]: native Tangram projection before count calibration
- same Xenium cells and 100 held-out genes as the fold truth

Tangram configuration follows the Br8667 validation template:
- cluster mode
- uniform density prior
- log-normalized input, target sum 1e4
- 500 epochs
"""

from __future__ import annotations

import argparse
import gc
import inspect
import json
import random
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


MODEL_KEY = "tangram"
MODEL_LABEL = "Tangram"

CLUSTER_LABEL_CANDIDATES = [
    "cellType_broad_k",
    "cellType_broad_hc",
    "cellType_hc",
    "cell_type",
    "cell_type_annotation",
    "scClassify",
    "labels",
]


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(
            "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
        ),
    )

    parser.add_argument("--experiment", required=True)
    parser.add_argument("--experiment-label", required=True)

    parser.add_argument(
        "--fold",
        type=int,
        required=True,
        choices=[1, 2, 3],
    )

    parser.add_argument(
        "--reference",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--fold-dir",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
    )

    parser.add_argument(
        "--tangram-mode",
        choices=["clusters", "cells"],
        default="clusters",
    )

    parser.add_argument(
        "--cluster-label",
        default="auto",
    )

    parser.add_argument(
        "--min-reference-cells-per-cluster",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--density-prior",
        default="uniform",
    )

    parser.add_argument(
        "--normalization-target-sum",
        type=float,
        default=1e4,
    )

    parser.add_argument(
        "--run-mode",
        choices=["smoke", "full"],
        default="full",
    )

    parser.add_argument(
        "--smoke-epochs",
        type=int,
        default=50,
    )

    parser.add_argument(
        "--full-epochs",
        type=int,
        default=500,
    )

    parser.add_argument(
        "--device",
        default="cpu",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=8667,
    )

    parser.add_argument(
        "--cpus",
        type=int,
        default=16,
    )

    parser.add_argument(
        "--k-spatial",
        type=int,
        default=15,
    )

    parser.add_argument(
        "--ssim-grid-size",
        type=int,
        default=128,
    )

    parser.add_argument(
        "--n-clusters",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
    )

    parser.add_argument(
        "--plot-genes",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--skip-map-save",
        action="store_true",
    )

    return parser.parse_args()


def first_existing_column(frame, candidates):
    return next(
        (
            column
            for column in candidates
            if column in frame.columns
        ),
        None,
    )


def dense_float32(matrix):
    if sparse.issparse(matrix):
        matrix = matrix.toarray()

    return np.asarray(
        matrix,
        dtype=np.float32,
    )


def matrix_values(matrix):
    if sparse.issparse(matrix):
        return np.asarray(
            matrix.data
        )

    return np.asarray(
        matrix
    ).ravel()


def assert_finite_nonnegative(
    matrix,
    name,
    atol=1e-7,
):
    values = matrix_values(
        matrix
    )

    if values.size == 0:
        return

    if not np.isfinite(
        values
    ).all():
        raise ValueError(
            f"{name} contains NaN or infinity."
        )

    minimum = float(
        np.min(values)
    )

    if minimum < -atol:
        raise ValueError(
            f"{name} contains negative values; "
            f"minimum={minimum}."
        )


def normalized_log_matrix_from_adata(
    adata,
    genes,
    library_sizes,
    target_sum,
):
    genes = [
        str(gene)
        for gene in genes
    ]

    positions = (
        adata.var_names
        .astype(str)
        .get_indexer(genes)
    )

    if np.any(
        positions < 0
    ):
        missing = [
            genes[index]
            for index, position
            in enumerate(positions)
            if position < 0
        ]

        raise ValueError(
            "Missing genes during Tangram "
            f"normalization: {missing[:20]}"
        )

    counts = dense_float32(
        adata[
            :,
            genes,
        ].X
    )

    library_sizes = np.asarray(
        library_sizes,
        dtype=np.float32,
    ).reshape(-1)

    if (
        library_sizes.shape[0]
        != adata.n_obs
    ):
        raise ValueError(
            "Library-size vector does not "
            "match AnnData observations."
        )

    if not np.isfinite(
        library_sizes
    ).all():
        raise ValueError(
            "Library sizes contain NaN or infinity."
        )

    if np.any(
        library_sizes <= 0
    ):
        raise ValueError(
            "Library sizes contain "
            f"{int(np.sum(library_sizes <= 0))} "
            "zero/negative cells."
        )

    normalized = np.log1p(
        counts
        / library_sizes[:, None]
        * float(target_sum)
    ).astype(
        np.float32,
        copy=False,
    )

    assert_finite_nonnegative(
        normalized,
        "Tangram normalized matrix",
    )

    return normalized


def make_adata(
    matrix,
    obs,
    genes,
):
    out = ad.AnnData(
        X=np.asarray(
            matrix,
            dtype=np.float32,
        ),
        obs=obs.copy(),
        var=pd.DataFrame(
            index=pd.Index(
                genes,
                dtype=str,
            )
        ),
    )

    out.obs_names = (
        obs.index
        .astype(str)
    )

    out.var_names = pd.Index(
        genes,
        dtype=str,
    )

    return out


def parse_training_history(history):
    if history is None:
        return pd.DataFrame()

    if isinstance(
        history,
        pd.DataFrame,
    ):
        return history.copy()

    if isinstance(
        history,
        dict,
    ):
        return pd.DataFrame(
            {
                key: pd.Series(value)
                for key, value
                in history.items()
            }
        )

    if isinstance(
        history,
        (list, tuple),
    ):
        if (
            len(history)
            and isinstance(
                history[0],
                dict,
            )
        ):
            return pd.DataFrame(
                history
            )

        return pd.DataFrame(
            {
                "record": np.arange(
                    1,
                    len(history) + 1,
                ),
                "value": history,
            }
        )

    return pd.DataFrame(
        {
            "training_history_repr": [
                repr(history)
            ]
        }
    )


def save_training_diagnostics(
    ad_map,
    output_dir,
):
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if (
        "train_genes_df"
        in ad_map.uns
    ):
        scores = (
            ad_map.uns[
                "train_genes_df"
            ]
        )

        if isinstance(
            scores,
            pd.DataFrame,
        ):
            scores = scores.copy()
            scores.index.name = "gene"

            scores.to_csv(
                output_dir
                / "training_gene_scores.csv"
            )
        else:
            pd.DataFrame(
                scores
            ).to_csv(
                output_dir
                / "training_gene_scores.csv",
                index=False,
            )

    history = parse_training_history(
        ad_map.uns.get(
            "training_history"
        )
    )

    history.to_csv(
        output_dir
        / "training_history.csv",
        index=False,
    )

    numeric_columns = (
        history
        .select_dtypes(
            include=[np.number]
        )
        .columns
        .tolist()
    )

    for column in numeric_columns:

        values = pd.to_numeric(
            history[column],
            errors="coerce",
        ).to_numpy(
            dtype=float
        )

        finite_mask = np.isfinite(
            values
        )

        if not finite_mask.any():
            continue

        fig, ax = plt.subplots(
            figsize=(8, 5)
        )

        ax.plot(
            np.arange(
                1,
                len(values) + 1,
            ),
            values,
            linewidth=2,
        )

        ax.set_xlabel(
            "Recorded step"
        )

        ax.set_ylabel(
            column
        )

        ax.set_title(
            "Tangram training history: "
            f"{column}"
        )

        ax.grid(
            alpha=0.25
        )

        fig.tight_layout()

        safe = "".join(
            character
            if (
                character.isalnum()
                or character in "._-"
            )
            else "_"
            for character in column
        )

        fig.savefig(
            output_dir
            / f"training_{safe}.png",
            dpi=220,
            bbox_inches="tight",
        )

        plt.close(fig)


def main():
    args = parse_args()

    project_root = (
        args.project_root
        .resolve()
    )

    fold_dir = (
        args.fold_dir.resolve()
        if args.fold_dir
        is not None
        else (
            project_root
            / "data"
            / "processed"
            / "imputation_beta"
            / "Br8667"
            / "gene_folds_200_100"
        )
    )

    output_root = (
        args.output_root.resolve()
        if args.output_root
        is not None
        else (
            project_root
            / "outputs"
            / "imputation_beta"
            / "Br8667"
        )
    )

    output_dir = (
        output_root
        / args.experiment
        / MODEL_KEY
        / f"fold_{args.fold}"
    )

    figure_dir = (
        output_dir
        / "figures"
    )

    diagnostics_dir = (
        output_dir
        / "diagnostics"
    )

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
        / "src"
        / "imputation"
        / "common"
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
        row_sums,
        select_or_load_plot_genes,
        to_dense_float32,
    )

    try:
        import tangram as tg
    except ImportError as error:
        raise ImportError(
            "Tangram is not importable. "
            "Activate/install tangram-sc before "
            "submitting the arrays."
        ) from error

    seed = int(
        args.seed
        + args.fold
    )

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    torch.set_num_threads(
        max(
            1,
            int(args.cpus),
        )
    )

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(
            seed
        )

    requested_device = str(
        args.device
    )

    if (
        requested_device.startswith(
            "cuda"
        )
        and not torch.cuda.is_available()
    ):
        raise RuntimeError(
            f"Device {requested_device!r} "
            "requested but CUDA is unavailable."
        )

    observed_path = (
        fold_dir
        / (
            f"fold_{args.fold}"
            "_observed_genes.h5ad"
        )
    )

    heldout_path = (
        fold_dir
        / (
            f"fold_{args.fold}"
            "_heldout_genes.h5ad"
        )
    )

    master_path = (
        fold_dir
        / "gene_metrics_and_fold_assignment.h5ad"
    )

    for path in [
        observed_path,
        heldout_path,
        args.reference,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                path
            )

    spatial_observed = (
        sc.read_h5ad(
            observed_path
        )
    )

    spatial_truth = (
        sc.read_h5ad(
            heldout_path
        )
    )

    seq_data = (
        sc.read_h5ad(
            args.reference
        )
    )

    for current in [
        spatial_observed,
        spatial_truth,
        seq_data,
    ]:
        current.obs_names = (
            current.obs_names
            .astype(str)
        )

        current.var_names = (
            current.var_names
            .astype(str)
        )

    if (
        spatial_observed.n_vars
        != 200
    ):
        raise ValueError(
            f"Observed input has "
            f"{spatial_observed.n_vars} genes; "
            "expected 200."
        )

    if (
        spatial_truth.n_vars
        != 100
    ):
        raise ValueError(
            f"Held-out truth has "
            f"{spatial_truth.n_vars} genes; "
            "expected 100."
        )

    if not (
        spatial_observed
        .obs_names
        .equals(
            spatial_truth
            .obs_names
        )
    ):
        raise ValueError(
            "Observed and held-out "
            "Xenium cell IDs/order differ."
        )

    if not set(
        spatial_observed.var_names
    ).isdisjoint(
        set(
            spatial_truth.var_names
        )
    ):
        raise ValueError(
            "Observed and held-out "
            "gene sets overlap."
        )

    if (
        "spatial"
        not in spatial_observed.obsm
    ):
        raise KeyError(
            "Observed Xenium lacks "
            "obsm['spatial']."
        )

    assert_nonnegative_count_like(
        spatial_observed.X,
        "Observed Xenium fold input",
    )

    assert_nonnegative_count_like(
        spatial_truth.X,
        "Held-out Xenium ground truth",
    )

    assert_nonnegative_count_like(
        seq_data.X,
        "Tangram reference input",
    )

    observed_genes = (
        spatial_observed
        .var_names
        .astype(str)
        .tolist()
    )

    heldout_genes = (
        spatial_truth
        .var_names
        .astype(str)
        .tolist()
    )

    benchmark_gene_set = (
        set(observed_genes)
        | set(heldout_genes)
    )

    if (
        len(benchmark_gene_set)
        != 300
    ):
        raise ValueError(
            f"Fold union has "
            f"{len(benchmark_gene_set)} genes; "
            "expected 300."
        )

    if master_path.exists():

        master = ad.read_h5ad(
            master_path,
            backed="r",
        )

        benchmark_genes = (
            master
            .var_names
            .astype(str)
            .tolist()
        )

        master.file.close()

        if (
            set(benchmark_genes)
            != benchmark_gene_set
        ):
            raise ValueError(
                "Master 300-gene file and "
                "fold union contain different genes."
            )

    else:
        benchmark_genes = (
            observed_genes
            + heldout_genes
        )

    missing_reference = (
        pd.Index(
            benchmark_genes
        )
        .difference(
            seq_data.var_names
        )
    )

    if len(
        missing_reference
    ):
        raise ValueError(
            f"{len(missing_reference)} benchmark "
            "genes are absent from reference: "
            f"{missing_reference[:20].tolist()}"
        )

    coordinates = np.asarray(
        spatial_observed
        .obsm["spatial"],
        dtype=np.float32,
    )

    if (
        coordinates.ndim != 2
        or coordinates.shape[0]
        != spatial_observed.n_obs
        or coordinates.shape[1] < 2
    ):
        raise ValueError(
            "Invalid spatial coordinate shape: "
            f"{coordinates.shape}"
        )

    coordinates = (
        coordinates[:, :2]
    )

    spatial_truth.obsm[
        "spatial"
    ] = coordinates.copy()

    observed_library = (
        row_sums(
            spatial_observed.X
        )
        .astype(
            np.float32,
            copy=False,
        )
    )

    if np.any(
        observed_library <= 0
    ):
        raise ValueError(
            f"Fold {args.fold} still has "
            f"{int(np.sum(observed_library <= 0))} "
            "zero-count cells across the "
            "200 observed genes."
        )

    reference_library = (
        row_sums(
            seq_data.X
        )
        .astype(
            np.float32,
            copy=False,
        )
    )

    bad_reference_cells = (
        ~np.isfinite(
            reference_library
        )
        | (
            reference_library
            <= 0
        )
    )

    if np.any(
        bad_reference_cells
    ):
        raise ValueError(
            "Reference contains "
            f"{int(bad_reference_cells.sum())} "
            "nonfinite/zero-total cells."
        )

    cluster_label = None

    if (
        args.tangram_mode
        == "clusters"
    ):
        if (
            args.cluster_label
            == "auto"
        ):
            cluster_label = (
                first_existing_column(
                    seq_data.obs,
                    CLUSTER_LABEL_CANDIDATES,
                )
            )
        elif (
            args.cluster_label
            in seq_data.obs.columns
        ):
            cluster_label = (
                args.cluster_label
            )

        if cluster_label is None:
            raise KeyError(
                "Tangram cluster mode requires "
                "a reference cluster label. "
                f"Tried: {CLUSTER_LABEL_CANDIDATES}"
            )

        labels = (
            seq_data.obs[
                cluster_label
            ]
            .astype(str)
        )

        label_counts = (
            labels
            .value_counts()
        )

        keep_reference = (
            labels
            .map(label_counts)
            .to_numpy()
            >= int(
                args.min_reference_cells_per_cluster
            )
        )

        n_removed = int(
            (~keep_reference).sum()
        )

        if n_removed:
            seq_data = (
                seq_data[
                    keep_reference
                ]
                .copy()
            )

            reference_library = (
                reference_library[
                    keep_reference
                ]
            )

        print(
            "Reference cluster label:",
            cluster_label,
        )

        print(
            "Reference clusters:",
            seq_data.obs[
                cluster_label
            ]
            .astype(str)
            .nunique(),
        )

        print(
            "Removed small-cluster cells:",
            n_removed,
        )

    reference_normalized = (
        normalized_log_matrix_from_adata(
            seq_data,
            benchmark_genes,
            reference_library,
            target_sum=(
                args.normalization_target_sum
            ),
        )
    )

    spatial_normalized = (
        normalized_log_matrix_from_adata(
            spatial_observed,
            observed_genes,
            observed_library,
            target_sum=(
                args.normalization_target_sum
            ),
        )
    )

    ad_sc = make_adata(
        reference_normalized,
        seq_data.obs,
        benchmark_genes,
    )

    ad_sp = make_adata(
        spatial_normalized,
        spatial_observed.obs,
        observed_genes,
    )

    ad_sp.obsm[
        "spatial"
    ] = coordinates.copy()

    if (
        cluster_label
        is not None
    ):
        ad_sc.obs[
            cluster_label
        ] = (
            seq_data.obs[
                cluster_label
            ]
            .astype(str)
            .to_numpy()
        )

    assert_finite_nonnegative(
        ad_sc.X,
        "Tangram reference normalized input",
    )

    assert_finite_nonnegative(
        ad_sp.X,
        "Tangram spatial normalized input",
    )

    reference_training_sum = (
        np.asarray(
            ad_sc[
                :,
                observed_genes,
            ].X.sum(
                axis=0
            )
        )
        .reshape(-1)
    )

    spatial_training_sum = (
        np.asarray(
            ad_sp.X.sum(
                axis=0
            )
        )
        .reshape(-1)
    )

    if np.any(
        reference_training_sum
        <= 0
    ):
        bad = np.asarray(
            observed_genes
        )[
            reference_training_sum
            <= 0
        ].tolist()

        raise ValueError(
            "Observed Tangram genes are "
            "all-zero in reference: "
            f"{bad[:20]}"
        )

    if np.any(
        spatial_training_sum
        <= 0
    ):
        bad = np.asarray(
            observed_genes
        )[
            spatial_training_sum
            <= 0
        ].tolist()

        raise ValueError(
            "Observed Tangram genes are "
            "all-zero in Xenium: "
            f"{bad[:20]}"
        )

    tg.pp_adatas(
        ad_sc,
        ad_sp,
        genes=observed_genes,
        gene_to_lowercase=False,
    )

    actual_training_genes = [
        str(gene)
        for gene
        in ad_sc.uns[
            "training_genes"
        ]
    ]

    removed_training_genes = sorted(
        set(observed_genes)
        - set(actual_training_genes)
    )

    if removed_training_genes:

        pd.DataFrame(
            {
                "removed_gene":
                    removed_training_genes
            }
        ).to_csv(
            diagnostics_dir
            / "removed_training_genes.csv",
            index=False,
        )

        raise ValueError(
            "Tangram preprocessing removed "
            f"{len(removed_training_genes)} "
            "of the 200 observed genes: "
            f"{removed_training_genes[:20]}"
        )

    if (
        len(actual_training_genes)
        != 200
    ):
        raise ValueError(
            "Tangram retained "
            f"{len(actual_training_genes)} "
            "training genes; expected 200."
        )

    missing_after_preprocess = (
        pd.Index(
            heldout_genes
        )
        .difference(
            ad_sc.var_names
            .astype(str)
        )
    )

    if len(
        missing_after_preprocess
    ):
        raise ValueError(
            "Held-out genes disappeared "
            "after Tangram preprocessing: "
            f"{missing_after_preprocess[:20].tolist()}"
        )

    max_epochs = (
        args.smoke_epochs
        if args.run_mode
        == "smoke"
        else args.full_epochs
    )

    map_kwargs = {
        "adata_sc": ad_sc,
        "adata_sp": ad_sp,
        "mode": args.tangram_mode,
        "density_prior": args.density_prior,
        "num_epochs": int(
            max_epochs
        ),
        "device": requested_device,
        "random_state": seed,
    }

    if (
        args.tangram_mode
        == "clusters"
    ):
        map_kwargs[
            "cluster_label"
        ] = cluster_label

    print(
        "=" * 100
    )

    print(
        "STANDARD Tangram 3-FOLD BENCHMARK"
    )

    print(
        "=" * 100
    )

    print(
        "Experiment:",
        args.experiment,
    )

    print(
        "Experiment label:",
        args.experiment_label,
    )

    print(
        "Fold:",
        args.fold,
    )

    print(
        "Reference:",
        args.reference,
    )

    print(
        "Tangram reference input:",
        ad_sc.shape,
    )

    print(
        "Tangram spatial input:",
        ad_sp.shape,
    )

    print(
        "Held-out truth:",
        spatial_truth.shape,
    )

    print(
        "Tangram mode:",
        args.tangram_mode,
    )

    print(
        "Cluster label:",
        cluster_label,
    )

    print(
        "Density prior:",
        args.density_prior,
    )

    print(
        "Normalization target:",
        args.normalization_target_sum,
    )

    print(
        "Epochs:",
        max_epochs,
    )

    print(
        "Device:",
        requested_device,
    )

    print(
        "Output:",
        output_dir,
    )

    print(
        "=" * 100
    )

    training_start = (
        time.time()
    )

    ad_map = (
        tg.map_cells_to_space(
            **map_kwargs
        )
    )

    training_seconds = (
        time.time()
        - training_start
    )

    assert_finite_nonnegative(
        ad_map.X,
        "Tangram mapping matrix",
    )

    save_training_diagnostics(
        ad_map,
        diagnostics_dir,
    )

    if not (
        args.skip_map_save
    ):
        ad_map.write_h5ad(
            output_dir
            / "tangram_map.h5ad",
            compression="gzip",
        )

    project_kwargs = {
        "adata_map": ad_map,
        "adata_sc": ad_sc,
    }

    project_signature = (
        inspect.signature(
            tg.project_genes
        )
    )

    if (
        args.tangram_mode
        == "clusters"
        and "cluster_label"
        in project_signature.parameters
    ):
        project_kwargs[
            "cluster_label"
        ] = cluster_label

    prediction_start = (
        time.time()
    )

    ad_ge = tg.project_genes(
        **project_kwargs
    )

    prediction_seconds = (
        time.time()
        - prediction_start
    )

    ad_ge.obs_names = (
        ad_ge.obs_names
        .astype(str)
    )

    ad_ge.var_names = (
        ad_ge.var_names
        .astype(str)
    )

    if not (
        ad_ge.obs_names
        .equals(
            ad_sp.obs_names
        )
    ):
        missing_cells = (
            ad_sp.obs_names
            .difference(
                ad_ge.obs_names
            )
        )

        if len(
            missing_cells
        ):
            raise ValueError(
                "Tangram projection is missing "
                f"{len(missing_cells)} "
                "Xenium cells."
            )

        ad_ge = (
            ad_ge[
                ad_sp.obs_names
            ]
            .copy()
        )

    lower_to_actual = {}

    for gene in (
        ad_ge.var_names
        .astype(str)
    ):
        lower = (
            gene.lower()
        )

        if (
            lower
            not in lower_to_actual
        ):
            lower_to_actual[
                lower
            ] = gene

    actual_projected_genes = []
    missing_projected = []

    for gene in benchmark_genes:

        if (
            gene
            in ad_ge.var_names
        ):
            actual_projected_genes.append(
                gene
            )

        elif (
            gene.lower()
            in lower_to_actual
        ):
            actual_projected_genes.append(
                lower_to_actual[
                    gene.lower()
                ]
            )

        else:
            missing_projected.append(
                gene
            )

    if missing_projected:
        raise ValueError(
            "Tangram projection is missing "
            "benchmark genes: "
            f"{missing_projected[:20]}"
        )

    projected = (
        ad_ge[
            :,
            actual_projected_genes,
        ]
        .copy()
    )

    projected.var_names = (
        pd.Index(
            benchmark_genes,
            dtype=str,
        )
    )

    native_full = (
        dense_float32(
            projected.X
        )
    )

    expected_shape = (
        spatial_observed.n_obs,
        300,
    )

    if (
        native_full.shape
        != expected_shape
    ):
        raise ValueError(
            "Unexpected Tangram projection "
            f"shape {native_full.shape}; "
            f"expected {expected_shape}."
        )

    assert_finite_nonnegative(
        native_full,
        "Tangram native projection",
    )

    benchmark_index = (
        pd.Index(
            benchmark_genes
        )
    )

    observed_positions = (
        benchmark_index
        .get_indexer(
            observed_genes
        )
    )

    heldout_positions = (
        benchmark_index
        .get_indexer(
            heldout_genes
        )
    )

    if (
        np.any(
            observed_positions < 0
        )
        or np.any(
            heldout_positions < 0
        )
    ):
        raise ValueError(
            "Observed/held-out genes "
            "cannot be indexed in projection."
        )

    native_observed = (
        native_full[
            :,
            observed_positions,
        ]
    )

    native_heldout = (
        native_full[
            :,
            heldout_positions,
        ]
    )

    native_observed_library = (
        native_observed
        .sum(
            axis=1
        )
        .astype(
            np.float32,
            copy=False,
        )
    )

    zero_native_cells = (
        native_observed_library
        <= 1e-12
    )

    if np.any(
        zero_native_cells
    ):
        pd.DataFrame(
            {
                "cell_id":
                    spatial_observed
                    .obs_names[
                        zero_native_cells
                    ]
                    .astype(str),

                "observed_xenium_library":
                    observed_library[
                        zero_native_cells
                    ],

                "native_tangram_library":
                    native_observed_library[
                        zero_native_cells
                    ],
            }
        ).to_csv(
            diagnostics_dir
            / "zero_native_observed_library_cells.csv",
            index=False,
        )

        raise ValueError(
            "Tangram produced zero projected "
            "expression across the 200 observed "
            f"genes for {int(zero_native_cells.sum())} "
            "Xenium cells."
        )

    scale_factor = (
        observed_library
        / native_observed_library
    )

    if (
        not np.isfinite(
            scale_factor
        ).all()
        or np.any(
            scale_factor < 0
        )
    ):
        raise ValueError(
            "Tangram count-scale calibration "
            "factors contain invalid values."
        )

    predicted_counts = (
        native_heldout
        * scale_factor[:, None]
    )

    predicted_counts = (
        np.clip(
            predicted_counts,
            0,
            None,
        )
        .astype(
            np.float32,
            copy=False,
        )
    )

    if not np.isfinite(
        predicted_counts
    ).all():
        raise ValueError(
            "Tangram count-scale predictions "
            "contain NaN or infinity."
        )

    quantiles = np.quantile(
        scale_factor,
        [
            0,
            0.01,
            0.25,
            0.50,
            0.75,
            0.99,
            1.0,
        ],
    )

    pd.DataFrame(
        {
            "quantile": [
                "0",
                "0.01",
                "0.25",
                "0.50",
                "0.75",
                "0.99",
                "1.0",
            ],
            "scale_factor":
                quantiles,
        }
    ).to_csv(
        diagnostics_dir
        / "count_scale_factor_quantiles.csv",
        index=False,
    )

    prediction_adata = (
        ad.AnnData(
            X=predicted_counts,
            obs=spatial_truth.obs.copy(),
            var=spatial_truth.var.copy(),
        )
    )

    prediction_adata.layers[
        "count_scale"
    ] = predicted_counts.copy()

    prediction_adata.layers[
        "log1p"
    ] = np.log1p(
        predicted_counts
    ).astype(
        np.float32
    )

    prediction_adata.layers[
        "native"
    ] = native_heldout.astype(
        np.float32,
        copy=False,
    )

    prediction_adata.obsm[
        "spatial"
    ] = coordinates.copy()

    prediction_adata.obs[
        "tangram_count_scale_factor"
    ] = scale_factor.astype(
        np.float32
    )

    prediction_adata.uns[
        "benchmark"
    ] = {
        "experiment":
            args.experiment,

        "experiment_label":
            args.experiment_label,

        "model":
            MODEL_KEY,

        "model_label":
            MODEL_LABEL,

        "fold":
            int(args.fold),

        "observed_gene_count":
            200,

        "heldout_gene_count":
            100,

        "output_scale":
            "nonnegative floating-point expected counts",

        "native_output":
            (
                "Tangram project_genes output "
                "from log-normalized reference expression"
            ),

        "count_scale_calibration":
            (
                "Per-cell ratio of observed Xenium "
                "counts across the 200 observed genes "
                "to Tangram native projection summed "
                "across the same 200 genes."
            ),

        "reference":
            str(
                args.reference.resolve()
            ),

        "tangram_mode":
            args.tangram_mode,

        "cluster_label":
            cluster_label,

        "density_prior":
            args.density_prior,

        "normalization_target_sum":
            float(
                args.normalization_target_sum
            ),

        "epochs":
            int(max_epochs),

        "device":
            requested_device,

        "seed":
            seed,
    }

    prediction_path = (
        output_dir
        / "predicted_heldout_genes.h5ad"
    )

    prediction_adata.write_h5ad(
        prediction_path,
        compression="gzip",
    )

    observed_truth = (
        to_dense_float32(
            spatial_truth[
                :,
                heldout_genes,
            ].X
        )
    )

    gene_metrics = (
        calculate_gene_metrics(
            observed_truth,
            predicted_counts,
            heldout_genes,
            coordinates,
            k_neighbors=(
                args.k_spatial
            ),
            ssim_grid_size=(
                args.ssim_grid_size
            ),
            n_jobs=args.cpus,
        )
    )

    fold_metadata_columns = [
        "log_mean_expression",
        "dispersion",
        "detection_fraction",
        "morans_I",
        "heldout_fold",
        "metric_stratum",
    ]

    available_columns = [
        column
        for column
        in fold_metadata_columns
        if (
            column
            in spatial_truth.var.columns
        )
    ]

    if available_columns:

        fold_gene_metadata = (
            spatial_truth.var[
                available_columns
            ]
            .copy()
        )

        fold_gene_metadata[
            "gene"
        ] = (
            spatial_truth
            .var_names
            .astype(str)
        )

        gene_metrics = (
            gene_metrics.merge(
                fold_gene_metadata,
                on="gene",
                how="left",
            )
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

    nmi, ari = (
        calculate_cluster_metrics(
            observed_truth,
            predicted_counts,
            n_clusters=(
                args.n_clusters
            ),
            n_pcs=(
                args.n_pcs
            ),
            seed=args.seed,
        )
    )

    fold_summary = (
        build_fold_summary(
            gene_metrics,
            experiment=(
                args.experiment
            ),
            model_key=(
                MODEL_KEY
            ),
            model_label=(
                MODEL_LABEL
            ),
            fold=args.fold,
            n_cells=(
                spatial_truth.n_obs
            ),
            n_observed_genes=(
                spatial_observed.n_vars
            ),
            n_heldout_genes=(
                spatial_truth.n_vars
            ),
            nmi=nmi,
            ari=ari,
            training_seconds=(
                training_seconds
            ),
            prediction_seconds=(
                prediction_seconds
            ),
        )
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
            f"| Tangram | Fold {args.fold}"
        ),
    )

    plot_gene_metric_distributions(
        gene_metrics,
        figure_dir
        / "gene_metric_distributions.png",
        title=(
            f"{args.experiment_label} "
            f"| Tangram | Fold {args.fold} "
            "| 100 held-out genes"
        ),
    )

    shared_plot_gene_path = (
        fold_dir
        / "plot_genes"
        / (
            f"fold_{args.fold}"
            "_plot_genes_10.csv"
        )
    )

    selected_plot_genes = (
        select_or_load_plot_genes(
            spatial_truth,
            shared_plot_gene_path,
            n_genes=args.plot_genes,
            seed=args.seed,
        )
    )

    if (
        len(selected_plot_genes)
        != 10
    ):
        raise ValueError(
            "Fair-comparison spatial figure "
            "requires exactly 10 genes."
        )

    pd.DataFrame(
        {
            "plot_order":
                np.arange(
                    1,
                    len(
                        selected_plot_genes
                    ) + 1,
                ),

            "gene":
                selected_plot_genes,
        }
    ).to_csv(
        output_dir
        / "selected_plot_genes.csv",
        index=False,
    )

    plot_ten_gene_maps(
        observed_truth,
        predicted_counts,
        heldout_genes,
        selected_plot_genes,
        coordinates,
        gene_metrics,
        figure_dir
        / "ten_gene_observed_vs_prediction.png",
        figure_dir
        / "ten_gene_observed_vs_prediction.pdf",
        experiment_label=(
            args.experiment_label
        ),
        model_label=(
            MODEL_LABEL
        ),
        fold=args.fold,
    )

    run_config = {
        "experiment":
            args.experiment,

        "experiment_label":
            args.experiment_label,

        "model":
            MODEL_KEY,

        "model_label":
            MODEL_LABEL,

        "fold":
            int(args.fold),

        "reference":
            str(
                args.reference.resolve()
            ),

        "observed_input":
            str(
                observed_path.resolve()
            ),

        "heldout_truth":
            str(
                heldout_path.resolve()
            ),

        "prediction":
            str(
                prediction_path.resolve()
            ),

        "output_dir":
            str(
                output_dir.resolve()
            ),

        "run_mode":
            args.run_mode,

        "tangram_mode":
            args.tangram_mode,

        "cluster_label":
            cluster_label,

        "density_prior":
            args.density_prior,

        "normalization_target_sum":
            float(
                args.normalization_target_sum
            ),

        "epochs":
            int(max_epochs),

        "device":
            requested_device,

        "seed":
            seed,

        "cpus":
            int(args.cpus),

        "training_seconds":
            float(
                training_seconds
            ),

        "prediction_seconds":
            float(
                prediction_seconds
            ),

        "scale_factor_quantiles": {
            "min":
                float(quantiles[0]),

            "p01":
                float(quantiles[1]),

            "p25":
                float(quantiles[2]),

            "median":
                float(quantiles[3]),

            "p75":
                float(quantiles[4]),

            "p99":
                float(quantiles[5]),

            "max":
                float(quantiles[6]),
        },
    }

    with (
        output_dir
        / "run_config.json"
    ).open(
        "w",
        encoding="utf-8",
    ) as handle:

        json.dump(
            run_config,
            handle,
            indent=2,
        )

    (
        output_dir
        / "complete.flag"
    ).write_text(
        (
            "SUCCESS\n"
            f"experiment={args.experiment}\n"
            f"fold={args.fold}\n"
        ),
        encoding="utf-8",
    )

    print(
        "=" * 100
    )

    print(
        "SUCCESS"
    )

    print(
        "Prediction:",
        prediction_path,
    )

    print(
        "Gene metrics:",
        output_dir
        / "gene_level_metrics.csv",
    )

    print(
        "Fold metrics:",
        output_dir
        / "fold_level_metrics.csv",
    )

    print(
        "Training minutes:",
        f"{training_seconds / 60:.2f}",
    )

    print(
        "Projection minutes:",
        f"{prediction_seconds / 60:.2f}",
    )

    print(
        "Scale factor min/median/max:",
        f"{quantiles[0]:.6g}",
        f"{quantiles[3]:.6g}",
        f"{quantiles[6]:.6g}",
    )

    print(
        "=" * 100
    )

    del (
        native_full,
        native_observed,
        native_heldout,
        predicted_counts,
        ad_ge,
        projected,
    )

    gc.collect()

    if torch.cuda.is_available():
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
