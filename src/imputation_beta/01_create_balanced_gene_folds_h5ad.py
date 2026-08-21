#!/usr/bin/env python3
"""Create three balanced 200-observed/100-held-out Xenium gene folds."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler


FEATURES = [
    "log_mean_expression",
    "dispersion",
    "detection_fraction",
    "morans_I",
]

LOG = logging.getLogger("gene_folds")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)

    parser.add_argument(
        "--counts-layer",
        default="auto",
        help="auto, X, raw, or the name of a layer",
    )

    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--expected-genes", type=int, default=300)
    parser.add_argument("--n-folds", type=int, default=3)
    parser.add_argument("--heldout-per-fold", type=int, default=100)

    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--k-neighbors", type=int, default=15)
    parser.add_argument("--n-strata", type=int, default=30)
    parser.add_argument("--seed", type=int, default=14053)
    parser.add_argument("--n-jobs", type=int, default=-1)

    return parser.parse_args()


def as_float_matrix(matrix):
    if sparse.issparse(matrix):
        matrix = matrix.tocsr().astype(np.float64, copy=True)
        matrix.eliminate_zeros()
        return matrix

    return np.asarray(matrix, dtype=np.float64)


def get_counts(adata: ad.AnnData, source: str):
    if source == "auto":
        if "counts" in adata.layers:
            return (
                as_float_matrix(adata.layers["counts"]),
                "layers['counts']",
            )

        if adata.raw is not None:
            raw_names = pd.Index(adata.raw.var_names.astype(str))
            positions = raw_names.get_indexer(
                adata.var_names.astype(str)
            )

            if np.all(positions >= 0):
                return (
                    as_float_matrix(adata.raw.X[:, positions]),
                    "raw.X",
                )

        return as_float_matrix(adata.X), "X"

    if source == "X":
        return as_float_matrix(adata.X), "X"

    if source == "raw":
        if adata.raw is None:
            raise ValueError("adata.raw is absent")

        raw_names = pd.Index(adata.raw.var_names.astype(str))
        positions = raw_names.get_indexer(
            adata.var_names.astype(str)
        )

        if np.any(positions < 0):
            raise ValueError(
                "Some current genes are absent from adata.raw"
            )

        return (
            as_float_matrix(adata.raw.X[:, positions]),
            "raw.X",
        )

    if source not in adata.layers:
        raise KeyError(
            f"Layer {source!r} not found. "
            f"Available layers: {list(adata.layers.keys())}"
        )

    return (
        as_float_matrix(adata.layers[source]),
        f"layers['{source}']",
    )


def validate_input(
    adata: ad.AnnData,
    counts,
    args: argparse.Namespace,
) -> np.ndarray:
    if adata.n_vars != args.expected_genes:
        raise ValueError(
            f"Expected {args.expected_genes} genes, "
            f"but found {adata.n_vars}"
        )

    if (
        args.expected_genes
        != args.n_folds * args.heldout_per_fold
    ):
        raise ValueError(
            "expected_genes must equal "
            "n_folds × heldout_per_fold"
        )

    if not adata.var_names.is_unique:
        raise ValueError("Gene names are not unique")

    if not adata.obs_names.is_unique:
        raise ValueError("Cell IDs are not unique")

    if counts.shape != adata.shape:
        raise ValueError(
            f"Counts shape {counts.shape} does not match "
            f"AnnData shape {adata.shape}"
        )

    if args.spatial_key not in adata.obsm:
        raise KeyError(
            f"adata.obsm[{args.spatial_key!r}] is missing"
        )

    coordinates = np.asarray(
        adata.obsm[args.spatial_key],
        dtype=np.float64,
    )

    if (
        coordinates.shape[0] != adata.n_obs
        or coordinates.shape[1] < 2
    ):
        raise ValueError(
            f"Invalid spatial-coordinate shape: "
            f"{coordinates.shape}"
        )

    coordinates = coordinates[:, :2]

    if not np.isfinite(coordinates).all():
        raise ValueError(
            "Spatial coordinates contain NaN or Inf"
        )

    values = (
        counts.data
        if sparse.issparse(counts)
        else counts
    )

    if not np.isfinite(values).all():
        raise ValueError(
            "Count matrix contains NaN or Inf"
        )

    if np.min(values, initial=0) < 0:
        raise ValueError(
            "Count matrix contains negative values"
        )

    return coordinates


def log_normalize(counts, target_sum: float):
    if sparse.issparse(counts):
        matrix = counts.tocsr().astype(
            np.float64,
            copy=True,
        )

        totals = np.asarray(
            matrix.sum(axis=1)
        ).ravel()

        scale = np.divide(
            target_sum,
            totals,
            out=np.zeros_like(
                totals,
                dtype=np.float64,
            ),
            where=totals > 0,
        )

        matrix = sparse.diags(scale).dot(matrix).tocsr()
        matrix.data = np.log1p(matrix.data)

        return matrix

    matrix = np.asarray(
        counts,
        dtype=np.float64,
    ).copy()

    totals = matrix.sum(axis=1)

    scale = np.divide(
        target_sum,
        totals,
        out=np.zeros_like(
            totals,
            dtype=np.float64,
        ),
        where=totals > 0,
    )

    matrix *= scale[:, None]

    return np.log1p(matrix)


def calculate_mean_variance(matrix):
    if sparse.issparse(matrix):
        mean = np.asarray(
            matrix.mean(axis=0)
        ).ravel()

        squared_mean = np.asarray(
            matrix.power(2).mean(axis=0)
        ).ravel()

        variance = np.maximum(
            squared_mean - mean**2,
            0.0,
        )

        return mean, variance

    return (
        np.asarray(matrix.mean(axis=0)),
        np.asarray(matrix.var(axis=0)),
    )


def calculate_detection_fraction(counts):
    if sparse.issparse(counts):
        matrix = counts.tocsc(copy=True)
        matrix.eliminate_zeros()

        n_detected = np.diff(matrix.indptr)
    else:
        n_detected = np.count_nonzero(
            counts > 0,
            axis=0,
        )

    detection_fraction = (
        n_detected / counts.shape[0]
    )

    return n_detected, detection_fraction


def construct_spatial_weights(
    coordinates: np.ndarray,
    k_neighbors: int,
    n_jobs: int,
):
    if coordinates.shape[0] < 2:
        raise ValueError(
            "At least two cells are required"
        )

    k_neighbors = min(
        k_neighbors,
        coordinates.shape[0] - 1,
    )

    graph = kneighbors_graph(
        coordinates,
        n_neighbors=k_neighbors,
        mode="connectivity",
        include_self=False,
        n_jobs=n_jobs,
    ).tocsr()

    # Make graph symmetric.
    graph = graph.maximum(graph.T).tocsr()

    graph.setdiag(0)
    graph.eliminate_zeros()

    row_sums = np.asarray(
        graph.sum(axis=1)
    ).ravel()

    if np.any(row_sums == 0):
        raise ValueError(
            "Spatial graph contains isolated cells"
        )

    # Row-standardized spatial weights.
    weights = sparse.diags(
        1.0 / row_sums
    ).dot(graph).tocsr()

    return weights


def get_gene_vector(
    matrix,
    gene_index: int,
) -> np.ndarray:
    if sparse.issparse(matrix):
        return (
            matrix[:, gene_index]
            .toarray()
            .ravel()
        )

    return np.asarray(
        matrix[:, gene_index]
    ).ravel()


def calculate_morans_i(
    expression,
    weights,
) -> np.ndarray:
    n_cells, n_genes = expression.shape

    weight_sum = float(weights.sum())
    scale = n_cells / weight_sum

    results = np.zeros(
        n_genes,
        dtype=np.float64,
    )

    for gene_index in range(n_genes):
        expression_vector = get_gene_vector(
            expression,
            gene_index,
        )

        centered = (
            expression_vector
            - expression_vector.mean()
        )

        denominator = float(
            centered @ centered
        )

        if denominator == 0:
            moran = 0.0
        else:
            numerator = float(
                centered @ (
                    weights @ centered
                )
            )

            moran = (
                scale
                * numerator
                / denominator
            )

        results[gene_index] = moran

        if (
            (gene_index + 1) % 25 == 0
            or gene_index + 1 == n_genes
        ):
            LOG.info(
                "Moran's I completed: %d/%d genes",
                gene_index + 1,
                n_genes,
            )

    return results


def calculate_gene_metrics(
    adata: ad.AnnData,
    counts,
    normalized_expression,
    spatial_weights,
) -> pd.DataFrame:
    mean_expression, expression_variance = (
        calculate_mean_variance(
            normalized_expression
        )
    )

    (
        n_detected_cells,
        detection_fraction,
    ) = calculate_detection_fraction(counts)

    metrics = pd.DataFrame(
        index=adata.var_names.astype(str)
    )

    metrics["gene"] = metrics.index

    metrics["log_mean_expression"] = (
        mean_expression
    )

    metrics["dispersion"] = (
        expression_variance
        / (mean_expression + 1e-8)
    )

    metrics["detection_fraction"] = (
        detection_fraction
    )

    metrics["morans_I"] = calculate_morans_i(
        normalized_expression,
        spatial_weights,
    )

    metrics["n_detected_cells"] = (
        n_detected_cells
    )

    metrics["log_expression_variance"] = (
        expression_variance
    )

    if not np.isfinite(
        metrics[FEATURES].to_numpy()
    ).all():
        raise ValueError(
            "Calculated metrics contain NaN or Inf"
        )

    return metrics


def assign_balanced_folds(
    metrics: pd.DataFrame,
    n_folds: int,
    capacity: int,
    n_strata: int,
    seed: int,
) -> pd.DataFrame:
    # Convert each feature to percentile rank.
    # This gives the four metrics comparable scales.
    ranked_features = metrics[FEATURES].rank(
        method="average",
        pct=True,
    )

    standardized_features = (
        StandardScaler().fit_transform(
            ranked_features
        )
    )

    # Create temporary strata of genes with similar
    # four-dimensional metric profiles.
    n_strata = min(
        n_strata,
        len(metrics) // n_folds,
    )

    strata = KMeans(
        n_clusters=n_strata,
        random_state=seed,
        n_init=50,
    ).fit_predict(
        standardized_features
    )

    random_generator = np.random.default_rng(seed)

    assignment = np.full(
        len(metrics),
        -1,
        dtype=int,
    )

    fold_sizes = np.zeros(
        n_folds,
        dtype=int,
    )

    fold_feature_sums = np.zeros(
        (
            n_folds,
            standardized_features.shape[1],
        ),
        dtype=np.float64,
    )

    labels, label_sizes = np.unique(
        strata,
        return_counts=True,
    )

    # Process the largest strata first.
    for label in labels[
        np.argsort(-label_sizes)
    ]:
        gene_indices = np.where(
            strata == label
        )[0]

        random_generator.shuffle(
            gene_indices
        )

        local_stratum_counts = np.zeros(
            n_folds,
            dtype=int,
        )

        for gene_index in gene_indices:
            available_folds = np.where(
                fold_sizes < capacity
            )[0]

            # First priority:
            # distribute genes from the same stratum evenly.
            minimum_local_count = (
                local_stratum_counts[
                    available_folds
                ].min()
            )

            available_folds = available_folds[
                local_stratum_counts[
                    available_folds
                ]
                == minimum_local_count
            ]

            # Second priority:
            # keep each fold close to the global
            # four-feature average.
            scores = []

            for fold_index in available_folds:
                projected_mean = (
                    fold_feature_sums[fold_index]
                    + standardized_features[
                        gene_index
                    ]
                ) / (
                    fold_sizes[fold_index] + 1
                )

                scores.append(
                    float(
                        projected_mean
                        @ projected_mean
                    )
                )

            scores = np.asarray(scores)

            best_folds = available_folds[
                np.isclose(
                    scores,
                    scores.min(),
                )
            ]

            selected_fold = int(
                random_generator.choice(
                    best_folds
                )
            )

            assignment[gene_index] = (
                selected_fold
            )

            fold_sizes[selected_fold] += 1

            fold_feature_sums[
                selected_fold
            ] += standardized_features[
                gene_index
            ]

            local_stratum_counts[
                selected_fold
            ] += 1

    if np.any(assignment < 0):
        raise RuntimeError(
            "Some genes were not assigned to a fold"
        )

    if not np.all(fold_sizes == capacity):
        raise RuntimeError(
            f"Invalid fold sizes: "
            f"{fold_sizes.tolist()}"
        )

    metrics = metrics.copy()

    metrics["metric_stratum"] = (
        strata + 1
    )

    metrics["heldout_fold"] = (
        assignment + 1
    )

    for feature in FEATURES:
        metrics[
            f"percentile_{feature}"
        ] = ranked_features[
            feature
        ].to_numpy()

    return metrics


def add_metrics_to_var(
    adata: ad.AnnData,
    metrics: pd.DataFrame,
) -> None:
    aligned_metrics = metrics.loc[
        adata.var_names.astype(str)
    ]

    for column in metrics.columns:
        if column == "gene":
            continue

        adata.var[column] = (
            aligned_metrics[column].to_numpy()
        )


def make_leakage_safe(
    adata: ad.AnnData,
    spatial_key: str,
) -> None:
    # Variable slicing does not subset .raw,
    # so remove .raw explicitly.
    adata.raw = None

    # Preserve physical spatial coordinates only.
    spatial_coordinates = np.asarray(
        adata.obsm[spatial_key]
    ).copy()

    for key in list(adata.obsm.keys()):
        del adata.obsm[key]

    adata.obsm[spatial_key] = (
        spatial_coordinates
    )

    adata.obsp.clear()
    adata.varm.clear()
    adata.varp.clear()
    adata.uns.clear()


def create_balance_summary(
    metrics: pd.DataFrame,
    n_folds: int,
) -> pd.DataFrame:
    rows = []

    for fold in range(1, n_folds + 1):
        fold_metrics = metrics[
            metrics["heldout_fold"] == fold
        ]

        for metric in FEATURES:
            rows.append(
                {
                    "fold": fold,
                    "metric": metric,
                    "n": len(fold_metrics),
                    "mean": fold_metrics[
                        metric
                    ].mean(),
                    "median": fold_metrics[
                        metric
                    ].median(),
                    "std": fold_metrics[
                        metric
                    ].std(ddof=0),
                }
            )

    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format=(
            "%(asctime)s | "
            "%(levelname)s | "
            "%(message)s"
        ),
    )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    LOG.info(
        "Reading input: %s",
        args.input,
    )

    adata = ad.read_h5ad(args.input)

    adata.var_names = (
        adata.var_names.astype(str)
    )

    adata.obs_names = (
        adata.obs_names.astype(str)
    )

    counts, counts_source = get_counts(
        adata,
        args.counts_layer,
    )

    coordinates = validate_input(
        adata,
        counts,
        args,
    )

    LOG.info(
        "Input shape: %d cells × %d genes",
        adata.n_obs,
        adata.n_vars,
    )

    LOG.info(
        "Counts source: %s",
        counts_source,
    )

    normalized_expression = log_normalize(
        counts,
        args.target_sum,
    )

    spatial_weights = construct_spatial_weights(
        coordinates,
        args.k_neighbors,
        args.n_jobs,
    )

    metrics = calculate_gene_metrics(
        adata,
        counts,
        normalized_expression,
        spatial_weights,
    )

    metrics = assign_balanced_folds(
        metrics,
        n_folds=args.n_folds,
        capacity=args.heldout_per_fold,
        n_strata=args.n_strata,
        seed=args.seed,
    )

    balance_summary = create_balance_summary(
        metrics,
        args.n_folds,
    )

    LOG.info(
        "Fold balance summary:\n%s",
        balance_summary.to_string(index=False),
    )

    # Master H5AD containing all 300 genes,
    # metrics, strata and fold assignments.
    master = adata.copy()
    master.raw = None

    add_metrics_to_var(
        master,
        metrics,
    )

    master.uns["fold_design"] = {
        "features": FEATURES,
        "counts_source": counts_source,
        "spatial_key": args.spatial_key,
        "k_neighbors": args.k_neighbors,
        "n_strata": args.n_strata,
        "n_folds": args.n_folds,
        "heldout_per_fold": (
            args.heldout_per_fold
        ),
        "seed": args.seed,
        "description": (
            "The same folds must be reused for "
            "Experiment 5, 5.1, and 5.3."
        ),
    }

    master.uns["fold_balance_summary"] = {
        column: balance_summary[column].tolist()
        for column in balance_summary.columns
    }

    master_path = (
        args.output_dir
        / "gene_metrics_and_fold_assignment.h5ad"
    )

    master.write_h5ad(
        master_path,
        compression="gzip",
    )

    LOG.info(
        "Wrote master file: %s",
        master_path,
    )

    manifest = {
        "input": str(args.input.resolve()),
        "counts_source": counts_source,
        "seed": args.seed,
        "features": FEATURES,
        "folds": {},
    }

    for fold in range(
        1,
        args.n_folds + 1,
    ):
        heldout_genes = metrics.index[
            metrics["heldout_fold"] == fold
        ].tolist()

        observed_genes = metrics.index[
            metrics["heldout_fold"] != fold
        ].tolist()

        observed = adata[
            :,
            observed_genes,
        ].copy()

        heldout = adata[
            :,
            heldout_genes,
        ].copy()

        make_leakage_safe(
            observed,
            args.spatial_key,
        )

        make_leakage_safe(
            heldout,
            args.spatial_key,
        )

        add_metrics_to_var(
            observed,
            metrics,
        )

        add_metrics_to_var(
            heldout,
            metrics,
        )

        observed.var["split_role"] = (
            "observed_model_input"
        )

        heldout.var["split_role"] = (
            "heldout_ground_truth"
        )

        shared_metadata = {
            "fold": fold,
            "seed": args.seed,
            "features": FEATURES,
            "counts_source_used_for_folding": (
                counts_source
            ),
            "observed_gene_count": (
                len(observed_genes)
            ),
            "heldout_gene_count": (
                len(heldout_genes)
            ),
            "observed_genes": observed_genes,
            "heldout_genes": heldout_genes,
        }

        observed.uns["fold_design"] = {
            **shared_metadata,
            "role": "observed_model_input",
        }

        heldout.uns["fold_design"] = {
            **shared_metadata,
            "role": "heldout_ground_truth",
            "warning": (
                "Never provide this H5AD "
                "to the imputation model."
            ),
        }

        observed_path = (
            args.output_dir
            / f"fold_{fold}_observed_genes.h5ad"
        )

        heldout_path = (
            args.output_dir
            / f"fold_{fold}_heldout_genes.h5ad"
        )

        observed.write_h5ad(
            observed_path,
            compression="gzip",
        )

        heldout.write_h5ad(
            heldout_path,
            compression="gzip",
        )

        manifest["folds"][str(fold)] = {
            "observed": str(
                observed_path.resolve()
            ),
            "heldout": str(
                heldout_path.resolve()
            ),
            "n_observed": len(
                observed_genes
            ),
            "n_heldout": len(
                heldout_genes
            ),
        }

        LOG.info(
            "Fold %d: observed=%s, heldout=%s",
            fold,
            observed.shape,
            heldout.shape,
        )

    # Final held-out coverage and overlap checks.
    heldout_sets = [
        set(
            metrics.index[
                metrics["heldout_fold"] == fold
            ]
        )
        for fold in range(
            1,
            args.n_folds + 1,
        )
    ]

    if (
        len(set.union(*heldout_sets))
        != args.expected_genes
    ):
        raise RuntimeError(
            "Held-out folds do not cover "
            "all 300 genes"
        )

    if any(
        heldout_sets[first]
        & heldout_sets[second]
        for first in range(args.n_folds)
        for second in range(
            first + 1,
            args.n_folds,
        )
    ):
        raise RuntimeError(
            "Held-out folds overlap"
        )

    manifest_path = (
        args.output_dir
        / "fold_manifest.json"
    )

    with manifest_path.open("w") as handle:
        json.dump(
            manifest,
            handle,
            indent=2,
        )

    LOG.info(
        "SUCCESS: fold H5AD files created in %s",
        args.output_dir,
    )


if __name__ == "__main__":
    main()
