#!/usr/bin/env python3

"""
SpaGE three-fold GeneBridge benchmark.

Spatial model input
-------------------
Exactly 200 observed Xenium genes before SpaGE-specific alignment filtering.

Evaluation
----------
Exactly 100 held-out Xenium genes.

SpaGE normalization
-------------------
Reference:
    log1p(count / cell_required_300_total * 1,000,000)

Spatial:
    log1p(count / observed_200_total * median_observed_200_total)

Alignment genes
---------------
Observed genes only, retained when:
    reference detected cells >= 10
    reference normalized variance > threshold
    spatial normalized variance > threshold

Prediction
----------
SpaGE internally predicts all 300 benchmark genes.

The predictions for the 200 observed genes are used ONLY for
count-scale calibration.

The 100 held-out Xenium truth values never enter:
    normalization
    alignment
    SpaGE
    calibration

Output
------
X:
    calibrated expected Xenium counts for 100 held-out genes

layers["count_scale"]:
    same expected counts

layers["log1p"]:
    log1p(expected counts)

layers["native"]:
    raw SpaGE prediction on its native log1p-reference scale
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import random
import sys
import time
import warnings
from pathlib import Path

import anndata as ad

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc

from scipy import sparse


MODEL_KEY = "spage"
MODEL_LABEL = "SpaGE"


# =============================================================================
# Arguments
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(
            "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
        ),
    )

    parser.add_argument(
        "--experiment",
        required=True,
    )

    parser.add_argument(
        "--experiment-label",
        required=True,
    )

    parser.add_argument(
        "--fold",
        required=True,
        type=int,
        choices=[1, 2, 3],
    )

    parser.add_argument(
        "--reference",
        required=True,
        type=Path,
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
        "--n-pv",
        type=int,
        default=50,
    )

    parser.add_argument(
        "--min-reference-detected-cells",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--variance-threshold",
        type=float,
        default=1e-8,
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

    # Shared benchmark settings.
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

    return parser.parse_args()


# =============================================================================
# Helpers
# =============================================================================

def row_sums(matrix):
    if sparse.issparse(matrix):
        return np.asarray(
            matrix.sum(axis=1)
        ).ravel()

    return np.asarray(
        matrix
    ).sum(axis=1)


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
    tolerance=1e-7,
):
    values = matrix_values(
        matrix
    )

    if values.size == 0:
        return

    if not np.isfinite(
        values
    ).all():
        raise FloatingPointError(
            f"{name} contains NaN/Inf."
        )

    minimum = float(
        np.min(values)
    )

    if minimum < -tolerance:
        raise ValueError(
            f"{name} contains negative values. "
            f"Minimum={minimum}"
        )


def matrix_to_dataframe(
    adata,
    genes,
):
    genes = [
        str(gene)
        for gene in genes
    ]

    positions = (
        adata.var_names
        .get_indexer(
            genes
        )
    )

    if np.any(
        positions < 0
    ):
        missing = [
            genes[i]
            for i, position
            in enumerate(positions)
            if position < 0
        ]

        raise ValueError(
            "Genes missing from AnnData: "
            f"{missing[:20]}"
        )

    matrix = (
        adata[
            :,
            genes,
        ].X
    )

    if sparse.issparse(
        matrix
    ):
        matrix = (
            matrix
            .toarray()
        )

    matrix = np.asarray(
        matrix,
        dtype=np.float32,
    )

    return pd.DataFrame(
        matrix,
        index=(
            adata.obs_names
            .astype(str)
        ),
        columns=genes,
    )


def normalize_reference_log_cpm(
    raw_counts,
):
    """
    SpaGE tutorial reference normalization.

        log1p(
            count /
            required-gene cell total *
            1,000,000
        )
    """

    values = raw_counts.to_numpy(
        dtype=np.float64,
        copy=True,
    )

    library_sizes = (
        values.sum(
            axis=1
        )
    )

    if (
        not np.isfinite(
            library_sizes
        ).all()
    ):
        raise FloatingPointError(
            "Reference library sizes contain NaN/Inf."
        )

    if np.any(
        library_sizes <= 0
    ):
        raise ValueError(
            "Reference normalization received "
            "zero-library cells."
        )

    normalized = np.log1p(
        (
            values
            / library_sizes[:, None]
        )
        * 1_000_000.0
    )

    if not np.isfinite(
        normalized
    ).all():
        raise FloatingPointError(
            "Reference SpaGE normalization "
            "generated NaN/Inf."
        )

    return pd.DataFrame(
        normalized.astype(
            np.float32
        ),
        index=raw_counts.index.copy(),
        columns=raw_counts.columns.copy(),
    )


def normalize_spatial_input(
    raw_counts,
):
    """
    SpaGE tutorial spatial normalization.

        log1p(
            count /
            200-observed-gene library *
            median observed library
        )
    """

    values = raw_counts.to_numpy(
        dtype=np.float64,
        copy=True,
    )

    library_sizes = (
        values.sum(
            axis=1
        )
    )

    if (
        not np.isfinite(
            library_sizes
        ).all()
    ):
        raise FloatingPointError(
            "Spatial library sizes contain NaN/Inf."
        )

    if np.any(
        library_sizes <= 0
    ):
        raise ValueError(
            "Spatial input contains zero-library cells."
        )

    scale_factor = float(
        np.median(
            library_sizes
        )
    )

    if (
        not np.isfinite(
            scale_factor
        )
        or scale_factor <= 0
    ):
        raise ValueError(
            "Invalid median spatial library size."
        )

    normalized = np.log1p(
        (
            values
            / library_sizes[:, None]
        )
        * scale_factor
    )

    if not np.isfinite(
        normalized
    ).all():
        raise FloatingPointError(
            "Spatial SpaGE normalization "
            "generated NaN/Inf."
        )

    return (
        pd.DataFrame(
            normalized.astype(
                np.float32
            ),
            index=raw_counts.index.copy(),
            columns=raw_counts.columns.copy(),
        ),
        library_sizes.astype(
            np.float64
        ),
        scale_factor,
    )


# =============================================================================
# Main
# =============================================================================

def main():

    args = parse_args()

    # -------------------------------------------------------------------------
    # Thread control
    # -------------------------------------------------------------------------

    os.environ[
        "OMP_NUM_THREADS"
    ] = str(
        args.cpus
    )

    os.environ[
        "MKL_NUM_THREADS"
    ] = str(
        args.cpus
    )

    os.environ[
        "OPENBLAS_NUM_THREADS"
    ] = str(
        args.cpus
    )

    os.environ[
        "NUMEXPR_NUM_THREADS"
    ] = str(
        args.cpus
    )


    project_root = (
        args.project_root
        .resolve()
    )

    fold_dir = (
        args.fold_dir.resolve()
        if args.fold_dir is not None
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
        if args.output_root is not None
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


    # -------------------------------------------------------------------------
    # Shared benchmark functions
    # -------------------------------------------------------------------------

    common_dir = (
        project_root
        / "src"
        / "imputation"
        / "common"
    )

    sys.path.insert(
        0,
        str(
            common_dir
        ),
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


    # -------------------------------------------------------------------------
    # SpaGE
    # -------------------------------------------------------------------------

    spage_repo = (
        project_root
        / "src"
        / "imputation"
        / "SpaGE"
    )

    required_spage_file = (
        spage_repo
        / "SpaGE"
        / "main.py"
    )

    if not required_spage_file.is_file():
        raise FileNotFoundError(
            "SpaGE source is missing:\n"
            f"{required_spage_file}"
        )

    sys.path.insert(
        0,
        str(
            spage_repo
        ),
    )

    from SpaGE.main import SpaGE


    # -------------------------------------------------------------------------
    # Reproducibility
    # -------------------------------------------------------------------------

    seed = (
        int(
            args.seed
        )
        + int(
            args.fold
        )
    )

    random.seed(
        seed
    )

    np.random.seed(
        seed
    )


    # -------------------------------------------------------------------------
    # Paths
    # -------------------------------------------------------------------------

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
        if not path.is_file():
            raise FileNotFoundError(
                path
            )


    # -------------------------------------------------------------------------
    # Load
    # -------------------------------------------------------------------------

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

    for adata in [
        spatial_observed,
        spatial_truth,
        seq_data,
    ]:
        adata.obs_names = (
            adata.obs_names
            .astype(str)
        )

        adata.var_names = (
            adata.var_names
            .astype(str)
        )


    # -------------------------------------------------------------------------
    # Fold integrity
    # -------------------------------------------------------------------------

    if (
        spatial_observed.n_vars
        != 200
    ):
        raise ValueError(
            "Observed fold does not contain "
            "exactly 200 genes."
        )

    if (
        spatial_truth.n_vars
        != 100
    ):
        raise ValueError(
            "Held-out fold does not contain "
            "exactly 100 genes."
        )

    if not (
        spatial_observed
        .obs_names
        .equals(
            spatial_truth.obs_names
        )
    ):
        raise ValueError(
            "Observed/truth Xenium cells "
            "or cell order differ."
        )

    observed_genes = (
        spatial_observed
        .var_names
        .tolist()
    )

    heldout_genes = (
        spatial_truth
        .var_names
        .tolist()
    )

    if not set(
        observed_genes
    ).isdisjoint(
        heldout_genes
    ):
        raise ValueError(
            "Observed and held-out genes overlap."
        )

    if (
        len(
            set(
                observed_genes
            )
            | set(
                heldout_genes
            )
        )
        != 300
    ):
        raise ValueError(
            "Fold union is not exactly 300 genes."
        )

    if (
        "spatial"
        not in spatial_observed.obsm
    ):
        raise KeyError(
            "Observed Xenium lacks "
            "obsm['spatial']."
        )


    # -------------------------------------------------------------------------
    # Count / NaN checks
    # -------------------------------------------------------------------------

    assert_nonnegative_count_like(
        spatial_observed.X,
        "SpaGE observed Xenium input",
    )

    assert_nonnegative_count_like(
        spatial_truth.X,
        "SpaGE held-out Xenium truth",
    )


    # -------------------------------------------------------------------------
    # Canonical 300-gene order
    # -------------------------------------------------------------------------

    benchmark_gene_set = (
        set(
            observed_genes
        )
        | set(
            heldout_genes
        )
    )

    if master_path.is_file():

        master = (
            ad.read_h5ad(
                master_path,
                backed="r",
            )
        )

        benchmark_genes = (
            master.var_names
            .astype(str)
            .tolist()
        )

        master.file.close()

        if (
            set(
                benchmark_genes
            )
            != benchmark_gene_set
        ):
            raise ValueError(
                "Master 300-gene panel differs "
                "from fold union."
            )

    else:

        benchmark_genes = (
            observed_genes
            + heldout_genes
        )


    # -------------------------------------------------------------------------
    # All benchmark genes in reference
    # -------------------------------------------------------------------------

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
            "genes absent from reference: "
            f"{missing_reference[:20].tolist()}"
        )


    reference_benchmark = (
        seq_data[
            :,
            benchmark_genes,
        ]
    )

    assert_nonnegative_count_like(
        reference_benchmark.X,
        "SpaGE 300-gene reference counts",
    )


    # -------------------------------------------------------------------------
    # Coordinates
    # -------------------------------------------------------------------------

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
            "Invalid spatial coordinates: "
            f"{coordinates.shape}"
        )

    coordinates = (
        coordinates[
            :,
            :2,
        ]
    )

    if not np.isfinite(
        coordinates
    ).all():
        raise ValueError(
            "Spatial coordinates contain NaN/Inf."
        )

    spatial_truth.obsm[
        "spatial"
    ] = coordinates.copy()


    # -------------------------------------------------------------------------
    # Raw DataFrames
    # -------------------------------------------------------------------------

    reference_raw = (
        matrix_to_dataframe(
            seq_data,
            benchmark_genes,
        )
    )

    spatial_raw = (
        matrix_to_dataframe(
            spatial_observed,
            observed_genes,
        )
    )

    assert_finite_nonnegative(
        reference_raw.to_numpy(),
        "SpaGE reference raw matrix",
    )

    assert_finite_nonnegative(
        spatial_raw.to_numpy(),
        "SpaGE spatial raw matrix",
    )


    # -------------------------------------------------------------------------
    # Remove reference cells empty across required 300 genes
    # -------------------------------------------------------------------------

    reference_library = (
        reference_raw.sum(
            axis=1
        )
        .to_numpy(
            dtype=np.float64
        )
    )

    reference_keep = (
        np.isfinite(
            reference_library
        )
        & (
            reference_library > 0
        )
    )

    removed_reference_cells = int(
        (
            ~reference_keep
        ).sum()
    )

    reference_raw = (
        reference_raw.loc[
            reference_keep
        ]
        .copy()
    )

    if (
        reference_raw.shape[0]
        < 50
    ):
        raise ValueError(
            "Too few nonempty reference cells "
            "remain for SpaGE."
        )


    # -------------------------------------------------------------------------
    # Normalize
    # -------------------------------------------------------------------------

    reference_normalized = (
        normalize_reference_log_cpm(
            reference_raw
        )
    )

    (
        spatial_normalized,
        observed_library,
        spatial_scale_factor,
    ) = (
        normalize_spatial_input(
            spatial_raw
        )
    )


    # -------------------------------------------------------------------------
    # Alignment-gene filtering
    # -------------------------------------------------------------------------

    reference_detection_counts = (
        (
            reference_raw[
                observed_genes
            ]
            > 0
        )
        .sum(
            axis=0
        )
    )

    reference_std = (
        reference_normalized[
            observed_genes
        ]
        .std(
            axis=0,
            ddof=0,
        )
    )

    spatial_std = (
        spatial_normalized[
            observed_genes
        ]
        .std(
            axis=0,
            ddof=0,
        )
    )

    alignment_table = pd.DataFrame(
        {
            "gene":
                observed_genes,

            "reference_detected_cells": [
                int(
                    reference_detection_counts[
                        gene
                    ]
                )
                for gene in observed_genes
            ],

            "reference_std": [
                float(
                    reference_std[
                        gene
                    ]
                )
                for gene in observed_genes
            ],

            "spatial_std": [
                float(
                    spatial_std[
                        gene
                    ]
                )
                for gene in observed_genes
            ],
        }
    )

    alignment_table[
        "passes_reference_detection"
    ] = (
        alignment_table[
            "reference_detected_cells"
        ]
        >= args.min_reference_detected_cells
    )

    alignment_table[
        "passes_variance"
    ] = (
        np.isfinite(
            alignment_table[
                "reference_std"
            ]
        )
        & np.isfinite(
            alignment_table[
                "spatial_std"
            ]
        )
        & (
            alignment_table[
                "reference_std"
            ]
            > args.variance_threshold
        )
        & (
            alignment_table[
                "spatial_std"
            ]
            > args.variance_threshold
        )
    )

    alignment_table[
        "used_for_alignment"
    ] = (
        alignment_table[
            "passes_reference_detection"
        ]
        & alignment_table[
            "passes_variance"
        ]
    )

    alignment_genes = (
        alignment_table.loc[
            alignment_table[
                "used_for_alignment"
            ],
            "gene",
        ]
        .astype(str)
        .tolist()
    )

    alignment_table.to_csv(
        diagnostics_dir
        / "spage_alignment_genes.csv",
        index=False,
    )

    if (
        len(
            alignment_genes
        )
        < 2
    ):
        raise ValueError(
            "Too few valid SpaGE alignment genes."
        )


    # -------------------------------------------------------------------------
    # Build SpaGE matrices
    # -------------------------------------------------------------------------

    spage_spatial_input = (
        spatial_normalized[
            alignment_genes
        ]
        .astype(
            np.float32
        )
    )

    # All 300 are retained in RNA_data so SpaGE can predict:
    #   200 observed genes -> calibration only
    #   100 held-out genes -> benchmark prediction
    spage_reference_input = (
        reference_normalized[
            benchmark_genes
        ]
        .astype(
            np.float32
        )
    )

    requested_n_pv = int(
        args.n_pv
    )

    n_pv = min(
        requested_n_pv,
        len(
            alignment_genes
        ),
        spage_spatial_input.shape[0],
        spage_reference_input.shape[0],
    )

    if n_pv < 2:
        raise ValueError(
            f"Effective SpaGE n_pv={n_pv}."
        )

    assert_finite_nonnegative(
        spage_spatial_input.to_numpy(),
        "SpaGE normalized spatial alignment matrix",
    )


    # -------------------------------------------------------------------------
    # Log
    # -------------------------------------------------------------------------

    print(
        "=" * 100
    )

    print(
        "STANDARD SpaGE 3-FOLD BENCHMARK"
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
        "Reference raw shape:",
        seq_data.shape,
    )

    print(
        "Reference 300-gene cells retained:",
        reference_raw.shape,
    )

    print(
        "Reference zero-library cells removed:",
        removed_reference_cells,
    )

    print(
        "Spatial observed input:",
        spatial_observed.shape,
    )

    print(
        "Held-out truth:",
        spatial_truth.shape,
    )

    print(
        "Alignment genes retained:",
        len(
            alignment_genes
        ),
        "/ 200",
    )

    print(
        "Requested n_pv:",
        requested_n_pv,
    )

    print(
        "n_pv passed to SpaGE:",
        n_pv,
    )

    print(
        "SpaGE genes_to_predict:",
        len(
            benchmark_genes
        ),
    )

    print(
        "Spatial normalization median library:",
        spatial_scale_factor,
    )

    print(
        "Seed:",
        seed,
    )

    print(
        "Output:",
        output_dir,
    )

    print(
        "=" * 100
    )


    # -------------------------------------------------------------------------
    # Run SpaGE
    # -------------------------------------------------------------------------

    prediction_start = (
        time.time()
    )

    with warnings.catch_warnings():
        warnings.simplefilter(
            "default"
        )

        native_prediction_all = (
            SpaGE(
                Spatial_data=(
                    spage_spatial_input
                ),
                RNA_data=(
                    spage_reference_input
                ),
                n_pv=n_pv,
                genes_to_predict=(
                    benchmark_genes
                ),
            )
        )

    prediction_seconds = (
        time.time()
        - prediction_start
    )


    # -------------------------------------------------------------------------
    # Align output
    # -------------------------------------------------------------------------

    if not isinstance(
        native_prediction_all,
        pd.DataFrame,
    ):
        native_prediction_all = (
            pd.DataFrame(
                native_prediction_all,
                columns=benchmark_genes,
            )
        )

    if (
        native_prediction_all.shape[0]
        != spatial_observed.n_obs
    ):
        raise ValueError(
            "SpaGE output cell count is wrong: "
            f"{native_prediction_all.shape}"
        )

    native_prediction_all.index = (
        spatial_observed
        .obs_names
        .copy()
    )

    missing_prediction_genes = (
        pd.Index(
            benchmark_genes
        )
        .difference(
            native_prediction_all.columns
        )
    )

    if len(
        missing_prediction_genes
    ):
        raise ValueError(
            "SpaGE prediction missing genes: "
            f"{missing_prediction_genes[:20].tolist()}"
        )

    native_prediction_all = (
        native_prediction_all[
            benchmark_genes
        ]
        .astype(
            np.float32
        )
    )

    native_values = (
        native_prediction_all
        .to_numpy(
            dtype=np.float32
        )
    )

    assert_finite_nonnegative(
        native_values,
        "SpaGE native 300-gene prediction",
    )

    native_values = np.maximum(
        native_values,
        0,
    )


    # -------------------------------------------------------------------------
    # Back-transform SpaGE native log1p-reference abundance
    # -------------------------------------------------------------------------

    abundance_all = np.expm1(
        native_values.astype(
            np.float64
        )
    )

    if (
        not np.isfinite(
            abundance_all
        ).all()
    ):
        raise FloatingPointError(
            "expm1(SpaGE native prediction) "
            "generated NaN/Inf."
        )

    if np.any(
        abundance_all < -1e-8
    ):
        raise ValueError(
            "Back-transformed SpaGE abundance "
            "contains negative values."
        )

    abundance_all = np.maximum(
        abundance_all,
        0,
    )


    gene_to_index = {
        gene: index
        for index, gene
        in enumerate(
            benchmark_genes
        )
    }

    observed_indices = np.asarray(
        [
            gene_to_index[
                gene
            ]
            for gene in observed_genes
        ],
        dtype=int,
    )

    heldout_indices = np.asarray(
        [
            gene_to_index[
                gene
            ]
            for gene in heldout_genes
        ],
        dtype=int,
    )


    # -------------------------------------------------------------------------
    # Leakage-free count calibration
    # -------------------------------------------------------------------------

    predicted_observed_abundance = (
        abundance_all[
            :,
            observed_indices,
        ]
    )

    predicted_observed_library = (
        predicted_observed_abundance.sum(
            axis=1
        )
    )

    bad_native_library = (
        ~np.isfinite(
            predicted_observed_library
        )
        | (
            predicted_observed_library
            <= 1e-12
        )
    )

    if np.any(
        bad_native_library
    ):

        pd.DataFrame(
            {
                "cell_id":
                    spatial_observed
                    .obs_names[
                        bad_native_library
                    ]
                    .astype(str),

                "observed_xenium_library":
                    observed_library[
                        bad_native_library
                    ],

                "spage_predicted_observed_library":
                    predicted_observed_library[
                        bad_native_library
                    ],
            }
        ).to_csv(
            diagnostics_dir
            / "invalid_spage_calibration_cells.csv",
            index=False,
        )

        raise FloatingPointError(
            "SpaGE predicted zero/invalid "
            "observed-gene abundance for "
            f"{int(bad_native_library.sum())} cells."
        )

    count_scale_factor = (
        observed_library
        / predicted_observed_library
    )

    if (
        not np.isfinite(
            count_scale_factor
        ).all()
        or np.any(
            count_scale_factor <= 0
        )
    ):
        raise FloatingPointError(
            "SpaGE count-scale factors contain "
            "NaN/Inf/nonpositive values."
        )


    # -------------------------------------------------------------------------
    # Held-out prediction
    # -------------------------------------------------------------------------

    heldout_abundance = (
        abundance_all[
            :,
            heldout_indices,
        ]
    )

    predicted_counts = (
        heldout_abundance
        * count_scale_factor[
            :,
            None,
        ]
    ).astype(
        np.float32
    )

    assert_finite_nonnegative(
        predicted_counts,
        "SpaGE held-out count-scale prediction",
    )

    predicted_counts = np.maximum(
        predicted_counts,
        0,
    ).astype(
        np.float32,
        copy=False,
    )

    native_heldout = (
        native_values[
            :,
            heldout_indices,
        ]
        .astype(
            np.float32,
            copy=False,
        )
    )


    # -------------------------------------------------------------------------
    # Calibration diagnostics
    # -------------------------------------------------------------------------

    probabilities = np.asarray(
        [
            0,
            0.01,
            0.25,
            0.50,
            0.75,
            0.99,
            1.0,
        ],
        dtype=float,
    )

    quantiles = np.quantile(
        count_scale_factor,
        probabilities,
    )

    pd.DataFrame(
        {
            "quantile":
                probabilities,

            "count_scale_factor":
                quantiles,
        }
    ).to_csv(
        diagnostics_dir
        / "count_scale_factor_quantiles.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Prediction H5AD
    # -------------------------------------------------------------------------

    prediction_adata = ad.AnnData(
        X=predicted_counts,

        obs=(
            spatial_truth
            .obs
            .copy()
        ),

        var=(
            spatial_truth
            .var
            .copy()
        ),
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
    ] = native_heldout.copy()

    prediction_adata.obsm[
        "spatial"
    ] = coordinates.copy()

    prediction_adata.obs[
        "spage_count_scale_factor"
    ] = (
        count_scale_factor
        .astype(
            np.float32
        )
    )

    prediction_adata.uns[
        "benchmark"
    ] = {
        "model":
            MODEL_KEY,

        "model_label":
            MODEL_LABEL,

        "experiment":
            args.experiment,

        "experiment_label":
            args.experiment_label,

        "fold":
            int(
                args.fold
            ),

        "observed_gene_count":
            200,

        "heldout_gene_count":
            100,

        "alignment_gene_count":
            int(
                len(
                    alignment_genes
                )
            ),

        "requested_n_pv":
            int(
                requested_n_pv
            ),

        "n_pv_passed_to_spage":
            int(
                n_pv
            ),

        "minimum_reference_detected_cells":
            int(
                args.min_reference_detected_cells
            ),

        "reference_normalization":
            (
                "log1p(count / required_300_gene_cell_total * 1e6)"
            ),

        "spatial_normalization":
            (
                "log1p(count / observed_200_gene_cell_total * "
                "median_observed_200_gene_cell_total)"
            ),

        "native_scale":
            "SpaGE weighted reference log1p-CPM-like expression",

        "count_scale_calibration":
            (
                "Per-cell actual Xenium 200-observed-gene library "
                "divided by back-transformed SpaGE prediction across "
                "the same 200 observed genes."
            ),

        "heldout_truth_used_for_model":
            False,

        "heldout_truth_used_for_calibration":
            False,

        "reference":
            str(
                args.reference.resolve()
            ),

        "seed":
            int(
                seed
            ),
    }


    prediction_path = (
        output_dir
        / "predicted_heldout_genes.h5ad"
    )

    prediction_adata.write_h5ad(
        prediction_path,
        compression="gzip",
    )


    # -------------------------------------------------------------------------
    # Shared count-scale evaluation
    # -------------------------------------------------------------------------

    truth_counts = (
        to_dense_float32(
            spatial_truth[
                :,
                heldout_genes,
            ].X
        )
    )

    gene_metrics = (
        calculate_gene_metrics(
            truth_counts,
            predicted_counts,
            heldout_genes,
            coordinates,
            k_neighbors=(
                args.k_spatial
            ),
            ssim_grid_size=(
                args.ssim_grid_size
            ),
            n_jobs=(
                args.cpus
            ),
        )
    )


    # -------------------------------------------------------------------------
    # Fold metadata
    # -------------------------------------------------------------------------

    metadata_columns = [
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
        in metadata_columns
        if column
        in spatial_truth.var.columns
    ]

    if available_columns:

        metadata = (
            spatial_truth.var[
                available_columns
            ]
            .copy()
        )

        metadata[
            "gene"
        ] = (
            spatial_truth
            .var_names
            .astype(str)
        )

        gene_metrics = (
            gene_metrics.merge(
                metadata,
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


    # -------------------------------------------------------------------------
    # NMI / ARI
    # -------------------------------------------------------------------------

    nmi, ari = (
        calculate_cluster_metrics(
            truth_counts,
            predicted_counts,
            n_clusters=(
                args.n_clusters
            ),
            n_pcs=(
                args.n_pcs
            ),
            seed=(
                args.seed
            ),
        )
    )


    # -------------------------------------------------------------------------
    # Fold summary
    # -------------------------------------------------------------------------

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

            fold=(
                args.fold
            ),

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

            training_seconds=0.0,

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


    # -------------------------------------------------------------------------
    # Same benchmark plots
    # -------------------------------------------------------------------------

    plot_fold_metric_summary(
        fold_summary,

        figure_dir
        / "evaluation_metrics.png",

        title=(
            f"{args.experiment_label} "
            f"| SpaGE | Fold {args.fold}"
        ),
    )

    plot_gene_metric_distributions(
        gene_metrics,

        figure_dir
        / "gene_metric_distributions.png",

        title=(
            f"{args.experiment_label} "
            f"| SpaGE | Fold {args.fold} "
            "| 100 held-out genes"
        ),
    )


    # -------------------------------------------------------------------------
    # Same frozen 10-gene fold visualization
    # -------------------------------------------------------------------------

    plot_gene_path = (
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
            plot_gene_path,
            n_genes=(
                args.plot_genes
            ),
            seed=(
                args.seed
            ),
        )
    )

    if len(
        selected_plot_genes
    ) != 10:
        raise ValueError(
            "Expected exactly 10 shared "
            "visualization genes."
        )

    pd.DataFrame(
        {
            "plot_order":
                np.arange(
                    1,
                    11,
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

        experiment_label=(
            args.experiment_label
        ),

        model_label=(
            MODEL_LABEL
        ),

        fold=(
            args.fold
        ),
    )


    # -------------------------------------------------------------------------
    # Configuration manifest
    # -------------------------------------------------------------------------

    run_config = {
        "model":
            MODEL_KEY,

        "model_label":
            MODEL_LABEL,

        "experiment":
            args.experiment,

        "experiment_label":
            args.experiment_label,

        "fold":
            int(
                args.fold
            ),

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

        "reference_cells_original":
            int(
                seq_data.n_obs
            ),

        "reference_cells_removed_zero_300":
            int(
                removed_reference_cells
            ),

        "reference_cells_used":
            int(
                reference_raw.shape[0]
            ),

        "alignment_genes":
            int(
                len(
                    alignment_genes
                )
            ),

        "requested_n_pv":
            int(
                requested_n_pv
            ),

        "n_pv_passed_to_spage":
            int(
                n_pv
            ),

        "min_reference_detected_cells":
            int(
                args.min_reference_detected_cells
            ),

        "variance_threshold":
            float(
                args.variance_threshold
            ),

        "spatial_scale_factor":
            float(
                spatial_scale_factor
            ),

        "prediction_seconds":
            float(
                prediction_seconds
            ),

        "seed":
            int(
                seed
            ),

        "count_scale_factor_quantiles": {
            "min":
                float(
                    quantiles[0]
                ),

            "p01":
                float(
                    quantiles[1]
                ),

            "p25":
                float(
                    quantiles[2]
                ),

            "median":
                float(
                    quantiles[3]
                ),

            "p75":
                float(
                    quantiles[4]
                ),

            "p99":
                float(
                    quantiles[5]
                ),

            "max":
                float(
                    quantiles[6]
                ),
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


    # -------------------------------------------------------------------------
    # BOTH completion flags
    # -------------------------------------------------------------------------

    success_text = (
        "SUCCESS\n"
        f"model={MODEL_KEY}\n"
        f"experiment={args.experiment}\n"
        f"fold={args.fold}\n"
    )

    (
        output_dir
        / "complete.flag"
    ).write_text(
        success_text,
        encoding="utf-8",
    )

    (
        output_dir
        / "run_complete.flag"
    ).write_text(
        success_text,
        encoding="utf-8",
    )


    # -------------------------------------------------------------------------
    # Final log
    # -------------------------------------------------------------------------

    print()
    print(
        "=" * 100
    )

    print(
        "SUCCESS"
    )

    print(
        "Model:",
        MODEL_LABEL,
    )

    print(
        "Experiment:",
        args.experiment,
    )

    print(
        "Fold:",
        args.fold,
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
        "Alignment genes:",
        len(
            alignment_genes
        ),
    )

    print(
        "SpaGE prediction minutes:",
        f"{prediction_seconds / 60:.2f}",
    )

    print(
        "Count scale min/median/max:",
        f"{quantiles[0]:.6g}",
        f"{quantiles[3]:.6g}",
        f"{quantiles[6]:.6g}",
    )

    print(
        "run_complete.flag:",
        output_dir
        / "run_complete.flag",
    )

    print(
        "=" * 100
    )

    gc.collect()


if __name__ == "__main__":
    main()
