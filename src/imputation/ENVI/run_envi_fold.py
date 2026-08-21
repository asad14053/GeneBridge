#!/usr/bin/env python3

"""
Run one ENVI fold for the GeneBridge three-fold imputation benchmark.

Output contract
---------------
X:
    calibrated nonnegative expected Xenium counts

layers["count_scale"]:
    same calibrated expected counts

layers["log1p"]:
    log1p(count_scale)

layers["native"]:
    untouched ENVI held-out prediction

Spatial input:
    exactly 200 observed Xenium genes

Evaluation:
    exactly 100 held-out Xenium genes

Held-out Xenium truth is NEVER supplied to ENVI and is NEVER used
for count-scale calibration.
"""

from __future__ import annotations

import argparse
import gc
import inspect
import json
import os
import random
import sys
import time
from pathlib import Path


# ------------------------------------------------------------
# Force JAX CPU mode before importing scenvi/jax
# ------------------------------------------------------------

os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("JAX_PLATFORM_NAME", "cpu")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")


import anndata as ad

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc

from scipy import sparse


MODEL_KEY = "envi"
MODEL_LABEL = "ENVI"


# ============================================================
# Arguments
# ============================================================

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

    # ENVI parameters from validated notebook
    parser.add_argument(
        "--run-mode",
        choices=["smoke", "full"],
        default="full",
    )

    parser.add_argument(
        "--smoke-training-steps",
        type=int,
        default=500,
    )

    parser.add_argument(
        "--full-training-steps",
        type=int,
        default=16000,
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=128,
    )

    parser.add_argument(
        "--num-hvg",
        type=int,
        default=2048,
    )

    parser.add_argument(
        "--num-cov-genes",
        type=int,
        default=64,
    )

    parser.add_argument(
        "--k-nearest",
        type=int,
        default=8,
    )

    parser.add_argument(
        "--covet-batch-size",
        type=int,
        default=256,
    )

    parser.add_argument(
        "--spatial-dist",
        default="pois",
    )

    parser.add_argument(
        "--sc-dist",
        default="nb",
    )

    parser.add_argument(
        "--stable-eps",
        type=float,
        default=1e-6,
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

    # Shared benchmark parameters
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


# ============================================================
# Matrix helpers
# ============================================================

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


def row_sums(matrix):
    if sparse.issparse(matrix):
        return np.asarray(
            matrix.sum(axis=1)
        ).ravel()

    return np.asarray(
        matrix
    ).sum(axis=1)


def gene_sums(matrix):
    if sparse.issparse(matrix):
        return np.asarray(
            matrix.sum(axis=0)
        ).ravel()

    return np.asarray(
        matrix
    ).sum(axis=0)


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

    finite = np.isfinite(
        values
    )

    if not finite.all():
        raise ValueError(
            f"{name} contains "
            f"{int((~finite).sum())} NaN/Inf values."
        )

    minimum = float(
        np.min(values)
    )

    if minimum < -atol:
        raise ValueError(
            f"{name} contains negative values. "
            f"Minimum={minimum}"
        )


# ============================================================
# Reference-gene selection
# ============================================================

def choose_reference_genes(
    seq_data,
    benchmark_genes,
    num_hvg,
):
    """
    Preserve uploaded ENVI template design:

        all required benchmark genes
        +
        additional reference HVGs
    """

    benchmark_set = set(
        benchmark_genes
    )

    if (
        "highly_variable"
        in seq_data.var.columns
    ):
        hvg_mask = (
            seq_data.var[
                "highly_variable"
            ]
            .fillna(False)
            .astype(bool)
            .to_numpy()
        )

        hvg_candidates = (
            seq_data.var_names[
                hvg_mask
            ]
            .astype(str)
            .tolist()
        )

        hvg_source = (
            "reference.var['highly_variable']"
        )

    else:

        # Same deterministic fallback principle
        # as uploaded ENVI notebook.
        hvg_candidates = (
            seq_data.var_names
            .astype(str)
            .tolist()
        )

        hvg_source = (
            "reference var-order fallback"
        )

    extra_hvg = []

    for gene in hvg_candidates:

        if gene in benchmark_set:
            continue

        extra_hvg.append(
            gene
        )

        if len(extra_hvg) >= num_hvg:
            break

    selected = list(
        dict.fromkeys(
            benchmark_genes
            + extra_hvg
        )
    )

    return (
        selected,
        extra_hvg,
        hvg_source,
    )


# ============================================================
# JAX parameter check
# ============================================================

def check_jax_params_finite(
    params,
    jax,
):
    leaves = (
        jax.tree_util
        .tree_leaves(
            params
        )
    )

    bad = 0

    for leaf in leaves:

        values = np.asarray(
            leaf
        )

        if not np.isfinite(
            values
        ).all():

            bad += 1

    if bad:
        raise FloatingPointError(
            "ENVI training generated "
            f"non-finite model parameters "
            f"in {bad} JAX parameter leaves."
        )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    os.environ[
        "OMP_NUM_THREADS"
    ] = str(args.cpus)

    os.environ[
        "MKL_NUM_THREADS"
    ] = str(args.cpus)

    os.environ[
        "OPENBLAS_NUM_THREADS"
    ] = str(args.cpus)

    os.environ[
        "NUMEXPR_NUM_THREADS"
    ] = str(args.cpus)

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


    # --------------------------------------------------------
    # Shared frozen evaluator
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # ENVI imports
    # --------------------------------------------------------

    import jax
    import scenvi

    print(
        "scenvi:",
        getattr(
            scenvi,
            "__version__",
            "version-not-exposed",
        ),
    )

    print(
        "JAX:",
        getattr(
            jax,
            "__version__",
            "unknown",
        ),
    )

    print(
        "JAX backend:",
        jax.default_backend(),
    )

    print(
        "JAX devices:",
        jax.devices(),
    )

    if (
        jax.default_backend()
        .lower()
        != "cpu"
    ):
        raise RuntimeError(
            "Expected ENVI/JAX CPU backend, "
            f"found {jax.default_backend()}."
        )


    # --------------------------------------------------------
    # Reproducibility
    # --------------------------------------------------------

    seed = int(
        args.seed
        + args.fold
    )

    random.seed(
        seed
    )

    np.random.seed(
        seed
    )


    # --------------------------------------------------------
    # Input paths
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Fold integrity
    # --------------------------------------------------------

    if (
        spatial_observed.n_vars
        != 200
    ):
        raise ValueError(
            f"Observed fold has "
            f"{spatial_observed.n_vars} genes. "
            "Expected 200."
        )

    if (
        spatial_truth.n_vars
        != 100
    ):
        raise ValueError(
            f"Held-out fold has "
            f"{spatial_truth.n_vars} genes. "
            "Expected 100."
        )

    if not (
        spatial_observed
        .obs_names
        .equals(
            spatial_truth.obs_names
        )
    ):
        raise ValueError(
            "Observed and held-out "
            "Xenium cells/order differ."
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


    # --------------------------------------------------------
    # Strict count / NaN validation
    # --------------------------------------------------------

    assert_nonnegative_count_like(
        spatial_observed.X,
        "Observed Xenium fold input",
    )

    assert_nonnegative_count_like(
        spatial_truth.X,
        "Held-out Xenium truth",
    )

    assert_nonnegative_count_like(
        seq_data.X,
        "ENVI reference input",
    )


    # --------------------------------------------------------
    # Gene identities
    # --------------------------------------------------------

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
            f"{len(benchmark_gene_set)} genes. "
            "Expected 300."
        )


    # --------------------------------------------------------
    # Canonical 300-gene order
    # --------------------------------------------------------

    if master_path.exists():

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
            set(benchmark_genes)
            != benchmark_gene_set
        ):
            raise ValueError(
                "Master 300-gene panel "
                "does not match fold union."
            )

    else:

        benchmark_genes = (
            observed_genes
            + heldout_genes
        )


    # --------------------------------------------------------
    # All 300 benchmark genes MUST exist in reference
    # --------------------------------------------------------

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
            "genes absent from ENVI reference: "
            f"{missing_reference[:20].tolist()}"
        )


    # --------------------------------------------------------
    # Spatial coordinates
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Common nonempty target-cell audit
    # --------------------------------------------------------

    observed_library = (
        row_sums(
            spatial_observed.X
        )
        .astype(
            np.float32
        )
    )

    bad_spatial_cells = (
        ~np.isfinite(
            observed_library
        )
        | (
            observed_library
            <= 0
        )
    )

    if np.any(
        bad_spatial_cells
    ):
        raise ValueError(
            f"Fold {args.fold} still has "
            f"{int(bad_spatial_cells.sum())} "
            "invalid/zero Xenium cells "
            "across the 200 observed genes."
        )


    # --------------------------------------------------------
    # ENVI reference gene set:
    # 300 benchmark genes + up to 2048 HVG genes
    # --------------------------------------------------------

    (
        selected_reference_genes,
        extra_hvg,
        hvg_source,
    ) = choose_reference_genes(
        seq_data,
        benchmark_genes,
        args.num_hvg,
    )

    sc_model = (
        seq_data[
            :,
            selected_reference_genes,
        ]
        .copy()
    )

    spatial_model = (
        spatial_observed[
            :,
            observed_genes,
        ]
        .copy()
    )


    # --------------------------------------------------------
    # Remove reference cells empty on ENVI-selected genes only
    # --------------------------------------------------------

    selected_reference_library = (
        row_sums(
            sc_model.X
        )
    )

    keep_reference = (
        np.isfinite(
            selected_reference_library
        )
        & (
            selected_reference_library
            > 0
        )
    )

    removed_reference_cells = int(
        (~keep_reference).sum()
    )

    if removed_reference_cells:

        sc_model = (
            sc_model[
                keep_reference
            ]
            .copy()
        )


    # --------------------------------------------------------
    # Every benchmark gene must be informative in reference
    # --------------------------------------------------------

    required_reference_sums = (
        gene_sums(
            sc_model[
                :,
                benchmark_genes,
            ].X
        )
    )

    zero_reference_genes = (
        np.asarray(
            benchmark_genes
        )[
            required_reference_sums
            <= 0
        ]
        .tolist()
    )

    if zero_reference_genes:
        raise ValueError(
            "Benchmark genes all-zero in "
            "ENVI reference: "
            f"{zero_reference_genes[:20]}"
        )


    # --------------------------------------------------------
    # All 200 training genes must have spatial signal
    # --------------------------------------------------------

    spatial_gene_sums = (
        gene_sums(
            spatial_model.X
        )
    )

    zero_spatial_genes = (
        np.asarray(
            observed_genes
        )[
            spatial_gene_sums
            <= 0
        ]
        .tolist()
    )

    if zero_spatial_genes:
        raise ValueError(
            "Observed Xenium genes all-zero "
            "across target cells: "
            f"{zero_spatial_genes[:20]}"
        )


    # --------------------------------------------------------
    # Required ENVI metadata
    # --------------------------------------------------------

    if (
        "batch"
        not in spatial_model.obs.columns
    ):
        spatial_model.obs[
            "batch"
        ] = "Br8667"

    else:
        spatial_model.obs[
            "batch"
        ] = (
            spatial_model.obs[
                "batch"
            ]
            .astype(str)
        )

    if (
        "cell_type"
        not in spatial_model.obs.columns
    ):
        spatial_model.obs[
            "cell_type"
        ] = "unknown"


    # --------------------------------------------------------
    # Dense float32 input as in uploaded ENVI template
    # --------------------------------------------------------

    spatial_model.X = (
        dense_float32(
            spatial_model.X
        )
    )

    sc_model.X = (
        dense_float32(
            sc_model.X
        )
    )

    assert_finite_nonnegative(
        spatial_model.X,
        "ENVI spatial model input",
    )

    assert_finite_nonnegative(
        sc_model.X,
        "ENVI reference model input",
    )


    # --------------------------------------------------------
    # Input audit
    # --------------------------------------------------------

    pd.DataFrame(
        {
            "gene":
                selected_reference_genes,

            "role": [
                (
                    "benchmark"
                    if gene
                    in benchmark_gene_set
                    else "additional_reference_gene"
                )
                for gene
                in selected_reference_genes
            ],
        }
    ).to_csv(
        diagnostics_dir
        / "envi_reference_gene_manifest.csv",
        index=False,
    )

    pd.DataFrame(
        {
            "item": [
                "HVG source",
                "requested additional HVGs",
                "actual additional genes",
                "selected reference genes",
                "reference cells original",
                "reference cells removed zero selected genes",
                "reference cells retained",
            ],

            "value": [
                hvg_source,
                args.num_hvg,
                len(extra_hvg),
                len(selected_reference_genes),
                seq_data.n_obs,
                removed_reference_cells,
                sc_model.n_obs,
            ],
        }
    ).to_csv(
        diagnostics_dir
        / "envi_input_audit.csv",
        index=False,
    )


    # --------------------------------------------------------
    # Full run parameters
    # --------------------------------------------------------

    training_steps = (
        args.smoke_training_steps
        if args.run_mode
        == "smoke"
        else args.full_training_steps
    )


    # --------------------------------------------------------
    # Version-safe ENVI constructor
    # --------------------------------------------------------

    requested_kwargs = {
        "spatial_data":
            spatial_model,

        "sc_data":
            sc_model,

        "spatial_key":
            "spatial",

        "batch_key":
            "batch",

        "k_nearest":
            int(
                args.k_nearest
            ),

        "covet_batch_size":
            int(
                args.covet_batch_size
            ),

        "num_cov_genes":
            int(
                args.num_cov_genes
            ),

        "num_HVG":
            int(
                args.num_hvg
            ),

        # The 200 observed genes overlap spatial and are
        # automatically retained by ENVI.
        # Explicitly protect the 100 held-out genes.
        "sc_genes":
            heldout_genes,

        "spatial_dist":
            args.spatial_dist,

        "sc_dist":
            args.sc_dist,

        # Current ENVI exposes this specifically to
        # improve numerical stability.
        "stable_eps":
            float(
                args.stable_eps
            ),
    }

    constructor_signature = (
        inspect.signature(
            scenvi.ENVI
        )
    )

    accepted_constructor_args = set(
        constructor_signature.parameters
    )

    essential_args = {
        "spatial_data",
        "sc_data",
        "spatial_key",
        "batch_key",
        "k_nearest",
        "num_cov_genes",
        "num_HVG",
        "sc_genes",
        "spatial_dist",
        "sc_dist",
    }

    missing_essential = sorted(
        essential_args
        - accepted_constructor_args
    )

    if missing_essential:
        raise RuntimeError(
            "Installed ENVI version is "
            "incompatible with template. "
            f"Missing arguments: "
            f"{missing_essential}"
        )

    envi_kwargs = {
        key: value
        for key, value
        in requested_kwargs.items()
        if key
        in accepted_constructor_args
    }

    dropped_optional = sorted(
        set(
            requested_kwargs
        )
        - set(
            envi_kwargs
        )
    )


    # --------------------------------------------------------
    # Log configuration
    # --------------------------------------------------------

    print(
        "=" * 100
    )

    print(
        "STANDARD ENVI 3-FOLD BENCHMARK"
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
        "ENVI reference model input:",
        sc_model.shape,
    )

    print(
        "ENVI spatial model input:",
        spatial_model.shape,
    )

    print(
        "Held-out truth:",
        spatial_truth.shape,
    )

    print(
        "Training steps:",
        training_steps,
    )

    print(
        "Batch size:",
        args.batch_size,
    )

    print(
        "num_HVG:",
        args.num_hvg,
    )

    print(
        "num_cov_genes:",
        args.num_cov_genes,
    )

    print(
        "k_nearest:",
        args.k_nearest,
    )

    print(
        "covet_batch_size:",
        args.covet_batch_size,
    )

    print(
        "spatial_dist:",
        args.spatial_dist,
    )

    print(
        "sc_dist:",
        args.sc_dist,
    )

    print(
        "stable_eps:",
        (
            args.stable_eps
            if "stable_eps"
            in envi_kwargs
            else "unsupported/default"
        ),
    )

    print(
        "Seed:",
        seed,
    )

    print(
        "Dropped optional kwargs:",
        dropped_optional,
    )

    print(
        "Output:",
        output_dir,
    )

    print(
        "=" * 100
    )


    # --------------------------------------------------------
    # Initialize ENVI
    # --------------------------------------------------------

    init_start = (
        time.time()
    )

    envi_model = (
        scenvi.ENVI(
            **envi_kwargs
        )
    )

    initialization_seconds = (
        time.time()
        - init_start
    )


    # --------------------------------------------------------
    # Verify ENVI did not change our training-gene set
    # --------------------------------------------------------

    if (
        envi_model.spatial_data.n_vars
        != 200
    ):
        raise ValueError(
            "ENVI retained "
            f"{envi_model.spatial_data.n_vars} "
            "spatial genes. Expected 200."
        )

    retained_spatial_genes = (
        envi_model
        .spatial_data
        .var_names
        .astype(str)
        .tolist()
    )

    if (
        set(
            retained_spatial_genes
        )
        != set(
            observed_genes
        )
    ):
        raise ValueError(
            "ENVI changed the "
            "200 observed-gene set."
        )


    # --------------------------------------------------------
    # Verify 100 held-out genes survived preprocessing
    # --------------------------------------------------------

    missing_heldout_after_init = (
        pd.Index(
            heldout_genes
        )
        .difference(
            envi_model
            .sc_data
            .var_names
            .astype(str)
        )
    )

    if len(
        missing_heldout_after_init
    ):
        raise ValueError(
            "ENVI preprocessing dropped "
            "held-out genes: "
            f"{missing_heldout_after_init[:20].tolist()}"
        )


    # --------------------------------------------------------
    # Seed ENVI/JAX
    # --------------------------------------------------------

    try:
        jax_key = (
            jax.random.key(
                seed
            )
        )

    except AttributeError:
        jax_key = (
            jax.random.PRNGKey(
                seed
            )
        )


    # --------------------------------------------------------
    # Train
    # --------------------------------------------------------

    train_signature = (
        inspect.signature(
            envi_model.train
        )
    )

    train_kwargs = {
        "training_steps":
            int(
                training_steps
            ),

        "batch_size":
            int(
                args.batch_size
            ),
    }

    if (
        "key"
        in train_signature.parameters
    ):
        train_kwargs[
            "key"
        ] = jax_key

    training_start = (
        time.time()
    )

    envi_model.train(
        **train_kwargs
    )

    training_seconds = (
        time.time()
        - training_start
    )


    # --------------------------------------------------------
    # NaN check immediately after training
    # --------------------------------------------------------

    check_jax_params_finite(
        envi_model.params,
        jax,
    )

    if (
        "envi_latent"
        not in envi_model
        .spatial_data
        .obsm
    ):
        raise KeyError(
            "ENVI training did not create "
            "spatial_data.obsm['envi_latent']."
        )

    latent = np.asarray(
        envi_model
        .spatial_data
        .obsm[
            "envi_latent"
        ],
        dtype=np.float32,
    )

    if not np.isfinite(
        latent
    ).all():
        raise FloatingPointError(
            "ENVI spatial latent contains "
            "NaN/Inf after training."
        )


    # --------------------------------------------------------
    # Impute
    # --------------------------------------------------------

    prediction_start = (
        time.time()
    )

    envi_model.impute_genes()

    prediction_seconds = (
        time.time()
        - prediction_start
    )

    if (
        "imputation"
        not in envi_model
        .spatial_data
        .obsm
    ):
        raise KeyError(
            "ENVI did not create "
            "spatial_data.obsm['imputation']."
        )


    # --------------------------------------------------------
    # Read native ENVI imputation
    # --------------------------------------------------------

    imputation = (
        envi_model
        .spatial_data
        .obsm[
            "imputation"
        ]
    )

    if isinstance(
        imputation,
        pd.DataFrame,
    ):

        imputation_df = (
            imputation.copy()
        )

    else:

        imputation_df = (
            pd.DataFrame(
                np.asarray(
                    imputation
                ),

                index=(
                    envi_model
                    .spatial_data
                    .obs_names
                    .astype(str)
                ),

                columns=(
                    envi_model
                    .sc_data
                    .var_names
                    .astype(str)
                ),
            )
        )

    imputation_df.index = (
        imputation_df
        .index
        .astype(str)
    )

    imputation_df.columns = (
        imputation_df
        .columns
        .astype(str)
    )


    # --------------------------------------------------------
    # Cell alignment
    # --------------------------------------------------------

    if not (
        spatial_observed
        .obs_names
        .equals(
            pd.Index(
                imputation_df.index
            )
        )
    ):

        missing_cells = (
            spatial_observed
            .obs_names
            .difference(
                imputation_df.index
            )
        )

        if len(
            missing_cells
        ):
            raise ValueError(
                "ENVI imputation missing "
                f"{len(missing_cells)} "
                "Xenium cells."
            )

        imputation_df = (
            imputation_df.loc[
                spatial_observed.obs_names
            ]
        )


    # --------------------------------------------------------
    # Require complete benchmark output
    # --------------------------------------------------------

    missing_imputed = (
        pd.Index(
            benchmark_genes
        )
        .difference(
            imputation_df.columns
        )
    )

    if len(
        missing_imputed
    ):
        raise ValueError(
            "ENVI imputation missing "
            "benchmark genes: "
            f"{missing_imputed[:20].tolist()}"
        )


    # --------------------------------------------------------
    # Native predictions for observed + heldout
    # --------------------------------------------------------

    native_observed = (
        imputation_df.loc[
            :,
            observed_genes,
        ]
        .to_numpy(
            dtype=np.float32
        )
    )

    native_heldout = (
        imputation_df.loc[
            :,
            heldout_genes,
        ]
        .to_numpy(
            dtype=np.float32
        )
    )

    assert_finite_nonnegative(
        native_observed,
        "ENVI native observed-gene prediction",
    )

    assert_finite_nonnegative(
        native_heldout,
        "ENVI native held-out prediction",
    )


    # --------------------------------------------------------
    # Count-scale calibration using ONLY observed genes
    # --------------------------------------------------------

    native_observed_library = (
        native_observed
        .sum(
            axis=1,
            dtype=np.float64,
        )
    )

    invalid_native_library = (
        ~np.isfinite(
            native_observed_library
        )
        | (
            native_observed_library
            <= 1e-12
        )
    )

    if np.any(
        invalid_native_library
    ):

        pd.DataFrame(
            {
                "cell_id":
                    spatial_observed
                    .obs_names[
                        invalid_native_library
                    ]
                    .astype(str),

                "observed_xenium_library":
                    observed_library[
                        invalid_native_library
                    ],

                "native_envi_observed_library":
                    native_observed_library[
                        invalid_native_library
                    ],
            }
        ).to_csv(
            diagnostics_dir
            / "invalid_native_observed_library_cells.csv",
            index=False,
        )

        raise FloatingPointError(
            "ENVI produced invalid/zero "
            "native prediction across "
            "the 200 observed genes for "
            f"{int(invalid_native_library.sum())} cells."
        )


    scale_factor = (
        observed_library
        .astype(np.float64)
        / native_observed_library
    )

    if (
        not np.isfinite(
            scale_factor
        ).all()
        or np.any(
            scale_factor <= 0
        )
    ):
        raise FloatingPointError(
            "ENVI count-scale factors "
            "contain invalid values."
        )


    predicted_counts = (
        native_heldout
        .astype(np.float64)
        * scale_factor[
            :,
            None,
        ]
    ).astype(
        np.float32
    )

    assert_finite_nonnegative(
        predicted_counts,
        "ENVI calibrated held-out prediction",
    )

    # Only numerical tolerance after the strict
    # nonnegative check above.
    predicted_counts = (
        np.maximum(
            predicted_counts,
            0,
        )
        .astype(
            np.float32,
            copy=False,
        )
    )


    # --------------------------------------------------------
    # Scale-factor audit
    # --------------------------------------------------------

    probabilities = np.asarray(
        [
            0,
            0.01,
            0.25,
            0.50,
            0.75,
            0.99,
            1.0,
        ]
    )

    quantiles = (
        np.quantile(
            scale_factor,
            probabilities,
        )
    )

    pd.DataFrame(
        {
            "quantile":
                probabilities,

            "scale_factor":
                quantiles,
        }
    ).to_csv(
        diagnostics_dir
        / "count_scale_factor_quantiles.csv",
        index=False,
    )


    # --------------------------------------------------------
    # Prediction AnnData
    # --------------------------------------------------------

    prediction_adata = (
        ad.AnnData(
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
    )

    prediction_adata.layers[
        "count_scale"
    ] = predicted_counts.copy()

    prediction_adata.layers[
        "log1p"
    ] = (
        np.log1p(
            predicted_counts
        )
        .astype(
            np.float32
        )
    )

    prediction_adata.layers[
        "native"
    ] = (
        native_heldout
        .astype(
            np.float32,
            copy=False,
        )
    )

    prediction_adata.obsm[
        "spatial"
    ] = coordinates.copy()

    prediction_adata.obs[
        "envi_count_scale_factor"
    ] = (
        scale_factor
        .astype(
            np.float32
        )
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
            int(
                args.fold
            ),

        "observed_gene_count":
            200,

        "heldout_gene_count":
            100,

        "output_scale":
            "nonnegative floating-point expected Xenium counts",

        "native_output":
            "ENVI spatial_data.obsm['imputation']",

        "count_scale_calibration":
            (
                "Per-cell observed Xenium 200-gene library "
                "divided by ENVI native prediction across "
                "the same 200 observed genes."
            ),

        "heldout_truth_used_for_calibration":
            False,

        "reference":
            str(
                args.reference.resolve()
            ),

        "training_steps":
            int(
                training_steps
            ),

        "batch_size":
            int(
                args.batch_size
            ),

        "num_hvg":
            int(
                args.num_hvg
            ),

        "num_cov_genes":
            int(
                args.num_cov_genes
            ),

        "k_nearest":
            int(
                args.k_nearest
            ),

        "covet_batch_size":
            int(
                args.covet_batch_size
            ),

        "spatial_dist":
            args.spatial_dist,

        "sc_dist":
            args.sc_dist,

        "stable_eps":
            (
                float(
                    args.stable_eps
                )
                if "stable_eps"
                in envi_kwargs
                else None
            ),

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


    # --------------------------------------------------------
    # Same benchmark evaluation as VISTA/gimVI/Tangram
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Preserve fold metadata
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # NMI / ARI
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Fold metric summary
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Same metrics plots
    # --------------------------------------------------------

    plot_fold_metric_summary(
        fold_summary,

        figure_dir
        / "evaluation_metrics.png",

        title=(
            f"{args.experiment_label} "
            f"| ENVI | Fold {args.fold}"
        ),
    )


    plot_gene_metric_distributions(
        gene_metrics,

        figure_dir
        / "gene_metric_distributions.png",

        title=(
            f"{args.experiment_label} "
            f"| ENVI | Fold {args.fold} "
            "| 100 held-out genes"
        ),
    )


    # --------------------------------------------------------
    # Same frozen 10-gene fold visualization
    # --------------------------------------------------------

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

    if (
        len(
            selected_plot_genes
        )
        != 10
    ):
        raise ValueError(
            "Expected exactly 10 "
            "shared visualization genes."
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


    # --------------------------------------------------------
    # Run manifest
    # --------------------------------------------------------

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

        "run_mode":
            args.run_mode,

        "training_steps":
            int(
                training_steps
            ),

        "batch_size":
            int(
                args.batch_size
            ),

        "num_hvg":
            int(
                args.num_hvg
            ),

        "num_cov_genes":
            int(
                args.num_cov_genes
            ),

        "k_nearest":
            int(
                args.k_nearest
            ),

        "covet_batch_size":
            int(
                args.covet_batch_size
            ),

        "spatial_dist":
            args.spatial_dist,

        "sc_dist":
            args.sc_dist,

        "seed":
            seed,

        "cpus":
            int(
                args.cpus
            ),

        "reference_cells_removed_zero_selected":
            removed_reference_cells,

        "initialization_seconds":
            float(
                initialization_seconds
            ),

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


    # --------------------------------------------------------
    # BOTH completion flags
    # --------------------------------------------------------

    success_text = (
        "SUCCESS\n"
        f"experiment={args.experiment}\n"
        f"fold={args.fold}\n"
        f"model={MODEL_KEY}\n"
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


    # --------------------------------------------------------
    # Final log
    # --------------------------------------------------------

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
        "Initialization minutes:",
        f"{initialization_seconds / 60:.2f}",
    )

    print(
        "Training minutes:",
        f"{training_seconds / 60:.2f}",
    )

    print(
        "Imputation minutes:",
        f"{prediction_seconds / 60:.2f}",
    )

    print(
        "Scale factor min/median/max:",
        f"{quantiles[0]:.6g}",
        f"{quantiles[3]:.6g}",
        f"{quantiles[6]:.6g}",
    )

    print(
        "complete.flag:",
        output_dir
        / "complete.flag",
    )

    print(
        "run_complete.flag:",
        output_dir
        / "run_complete.flag",
    )

    print(
        "=" * 100
    )

    del (
        native_observed,
        native_heldout,
        predicted_counts,
        truth_counts,
    )

    gc.collect()


if __name__ == "__main__":
    main()
