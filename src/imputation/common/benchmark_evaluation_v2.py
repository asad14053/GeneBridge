
#!/usr/bin/env python3
"""Shared evaluation utilities for the GeneBridge 3-fold imputation benchmark.

The same functions and labels should be reused by VISTA, gimVI, Tangram,
ENVI, TransImp, and SpaGE so their results are directly comparable.
"""

from __future__ import annotations

import json
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable, Sequence

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import uniform_filter
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr
from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler


GENE_METRIC_COLUMNS = [
    "scc",
    "ssim",
    "rmse",
    "mae",
    "jsd",
    "moran_abs_error",
]

FOLD_METRIC_COLUMNS = [
    "scc_mean",
    "ssim_mean",
    "rmse_mean",
    "mae_mean",
    "jsd_mean",
    "moran_abs_error_mean",
    "nmi",
    "ari",
]

METRIC_LABELS = {
    "scc_mean": "SCC",
    "ssim_mean": "SSIM",
    "rmse_mean": "RMSE",
    "mae_mean": "MAE",
    "jsd_mean": "Jensen–Shannon divergence",
    "moran_abs_error_mean": "Moran's I absolute error",
    "nmi": "NMI",
    "ari": "ARI",
}

HIGHER_IS_BETTER = ["scc_mean", "ssim_mean", "nmi", "ari"]
LOWER_IS_BETTER = [
    "rmse_mean",
    "mae_mean",
    "jsd_mean",
    "moran_abs_error_mean",
]


def to_dense_float32(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def matrix_values(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    return np.asarray(matrix).ravel()


def assert_nonnegative_count_like(matrix, name: str, atol: float = 1e-6) -> None:
    """Accept integer dtype or floating dtype whose values are integer-like."""
    values = matrix_values(matrix)
    if values.size == 0:
        return
    if not np.isfinite(values).all():
        raise ValueError(f"{name} contains NaN or infinity.")
    if np.min(values) < -atol:
        raise ValueError(f"{name} contains negative values.")
    if not np.allclose(values, np.rint(values), atol=atol, rtol=0):
        sample = values[np.abs(values - np.rint(values)) > atol][:10]
        raise ValueError(
            f"{name} is not raw count-like. Example non-integer values: {sample.tolist()}"
        )


def row_sums(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel()
    return np.asarray(matrix).sum(axis=1)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]
    if x.size < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan
    result = spearmanr(x, y)
    return float(result.statistic if hasattr(result, "statistic") else result[0])


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    sd = float(np.std(values))
    if not np.isfinite(sd) or sd <= 1e-12:
        return np.zeros_like(values, dtype=np.float64)
    return (values - float(np.mean(values))) / sd


class SpatialRasterizer:
    """Rasterize irregular Xenium cell coordinates onto one fixed grid."""

    def __init__(self, coordinates: np.ndarray, grid_size: int = 128):
        coordinates = np.asarray(coordinates, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError(f"Invalid spatial coordinate shape: {coordinates.shape}")
        self.coordinates = coordinates[:, :2]
        self.grid_size = int(grid_size)
        if self.grid_size < 16:
            raise ValueError("grid_size must be at least 16.")

        x = self.coordinates[:, 0]
        y = self.coordinates[:, 1]
        x_min, x_max = float(np.min(x)), float(np.max(x))
        y_min, y_max = float(np.min(y)), float(np.max(y))

        x_scaled = (x - x_min) / max(x_max - x_min, 1e-12)
        y_scaled = (y - y_min) / max(y_max - y_min, 1e-12)

        x_bin = np.clip((x_scaled * (self.grid_size - 1)).astype(int), 0, self.grid_size - 1)
        y_bin = np.clip((y_scaled * (self.grid_size - 1)).astype(int), 0, self.grid_size - 1)

        self.flat_index = y_bin * self.grid_size + x_bin
        self.bin_counts = np.bincount(
            self.flat_index,
            minlength=self.grid_size * self.grid_size,
        ).astype(np.float64)
        self.occupied = self.bin_counts.reshape(self.grid_size, self.grid_size) > 0

    def rasterize(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=np.float64)
        if values.shape[0] != self.flat_index.shape[0]:
            raise ValueError("Expression vector and spatial coordinates have different cell counts.")
        sums = np.bincount(
            self.flat_index,
            weights=values,
            minlength=self.grid_size * self.grid_size,
        ).astype(np.float64)
        means = np.divide(
            sums,
            self.bin_counts,
            out=np.zeros_like(sums),
            where=self.bin_counts > 0,
        )
        return means.reshape(self.grid_size, self.grid_size)


def local_ssim(
    observed_image: np.ndarray,
    predicted_image: np.ndarray,
    occupied_mask: np.ndarray,
    window_size: int = 7,
) -> float:
    """Standard local-window SSIM, averaged over occupied spatial bins."""
    x = np.asarray(observed_image, dtype=np.float64)
    y = np.asarray(predicted_image, dtype=np.float64)

    combined_min = float(min(np.min(x), np.min(y)))
    combined_max = float(max(np.max(x), np.max(y)))
    data_range = combined_max - combined_min

    if data_range <= 1e-12:
        return 1.0 if np.allclose(x, y) else 0.0

    x = (x - combined_min) / data_range
    y = (y - combined_min) / data_range

    size = int(window_size)
    mu_x = uniform_filter(x, size=size, mode="reflect")
    mu_y = uniform_filter(y, size=size, mode="reflect")

    sigma_x = np.maximum(uniform_filter(x * x, size=size, mode="reflect") - mu_x * mu_x, 0)
    sigma_y = np.maximum(uniform_filter(y * y, size=size, mode="reflect") - mu_y * mu_y, 0)
    sigma_xy = uniform_filter(x * y, size=size, mode="reflect") - mu_x * mu_y

    c1 = 0.01 ** 2
    c2 = 0.03 ** 2

    numerator = (2 * mu_x * mu_y + c1) * (2 * sigma_xy + c2)
    denominator = (mu_x * mu_x + mu_y * mu_y + c1) * (sigma_x + sigma_y + c2)
    ssim_map = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=np.abs(denominator) > 1e-15,
    )

    mask = np.asarray(occupied_mask, dtype=bool)
    if not np.any(mask):
        return np.nan
    return float(np.mean(ssim_map[mask]))


def build_spatial_weights(
    coordinates: np.ndarray,
    k_neighbors: int = 15,
    n_jobs: int = 1,
):
    coordinates = np.asarray(coordinates, dtype=np.float64)[:, :2]
    if coordinates.shape[0] < 2:
        raise ValueError("At least two cells are required for Moran's I.")
    k = min(int(k_neighbors), coordinates.shape[0] - 1)
    graph = kneighbors_graph(
        coordinates,
        n_neighbors=k,
        mode="connectivity",
        include_self=False,
        n_jobs=n_jobs,
    ).tocsr()
    graph = graph.maximum(graph.T).tocsr()
    graph.setdiag(0)
    graph.eliminate_zeros()
    row_sum = np.asarray(graph.sum(axis=1)).ravel()
    if np.any(row_sum == 0):
        raise ValueError("The spatial graph contains isolated cells.")
    return sparse.diags(1.0 / row_sum).dot(graph).tocsr()


def morans_i(values: np.ndarray, weights) -> float:
    values = np.asarray(values, dtype=np.float64)
    centered = values - float(np.mean(values))
    denominator = float(centered @ centered)
    if denominator <= 1e-15:
        return 0.0
    numerator = float(centered @ (weights @ centered))
    return float((values.size / float(weights.sum())) * numerator / denominator)


def js_divergence(observed: np.ndarray, predicted: np.ndarray) -> float:
    """Jensen–Shannon divergence between nonnegative cell-wise gene profiles."""
    p = np.clip(np.asarray(observed, dtype=np.float64), 0, None)
    q = np.clip(np.asarray(predicted, dtype=np.float64), 0, None)

    p_sum = float(p.sum())
    q_sum = float(q.sum())

    if p_sum <= 1e-15 and q_sum <= 1e-15:
        return 0.0
    if p_sum <= 1e-15 or q_sum <= 1e-15:
        return 1.0

    eps = 1e-12
    p = (p + eps) / float((p + eps).sum())
    q = (q + eps) / float((q + eps).sum())
    return float(jensenshannon(p, q, base=2.0) ** 2)


def calculate_gene_metrics(
    observed_counts: np.ndarray,
    predicted_counts: np.ndarray,
    genes: Sequence[str],
    coordinates: np.ndarray,
    *,
    k_neighbors: int = 15,
    ssim_grid_size: int = 128,
    n_jobs: int = 1,
) -> pd.DataFrame:
    observed_counts = np.asarray(observed_counts, dtype=np.float32)
    predicted_counts = np.asarray(predicted_counts, dtype=np.float32)

    if observed_counts.shape != predicted_counts.shape:
        raise ValueError(
            f"Observed/predicted shape mismatch: {observed_counts.shape} vs {predicted_counts.shape}"
        )
    if observed_counts.shape[1] != len(genes):
        raise ValueError("Gene list length does not match matrix columns.")
    if np.any(predicted_counts < 0) or not np.isfinite(predicted_counts).all():
        raise ValueError("Predictions must be finite and nonnegative.")

    observed_log = np.log1p(observed_counts)
    predicted_log = np.log1p(predicted_counts)

    rasterizer = SpatialRasterizer(coordinates, grid_size=ssim_grid_size)
    weights = build_spatial_weights(coordinates, k_neighbors=k_neighbors, n_jobs=n_jobs)

    rows = []
    for index, gene in enumerate(genes):
        obs = observed_counts[:, index]
        pred = predicted_counts[:, index]
        obs_log = observed_log[:, index]
        pred_log = predicted_log[:, index]

        observed_image = rasterizer.rasterize(obs_log)
        predicted_image = rasterizer.rasterize(pred_log)

        observed_moran = morans_i(obs_log, weights)
        predicted_moran = morans_i(pred_log, weights)

        rows.append(
            {
                "gene": str(gene),
                "scc": safe_spearman(obs, pred),
                "ssim": local_ssim(
                    observed_image,
                    predicted_image,
                    rasterizer.occupied,
                ),
                "rmse": float(np.sqrt(np.mean((pred - obs) ** 2))),
                "mae": float(np.mean(np.abs(pred - obs))),
                "jsd": js_divergence(obs, pred),
                "moran_observed": observed_moran,
                "moran_predicted": predicted_moran,
                "moran_abs_error": abs(predicted_moran - observed_moran),
                "observed_mean": float(np.mean(obs)),
                "predicted_mean": float(np.mean(pred)),
                "observed_variance": float(np.var(obs)),
                "predicted_variance": float(np.var(pred)),
                "observed_detection_fraction": float(np.mean(obs > 0)),
                "predicted_positive_fraction": float(np.mean(pred > 0)),
            }
        )
    return pd.DataFrame(rows)


def cluster_expression(
    count_matrix: np.ndarray,
    *,
    n_clusters: int = 10,
    n_pcs: int = 30,
    seed: int = 8667,
) -> np.ndarray:
    """One frozen clustering pipeline shared by truth and every model."""
    x = np.log1p(np.clip(np.asarray(count_matrix, dtype=np.float32), 0, None))
    x = StandardScaler().fit_transform(x)

    max_pcs = min(int(n_pcs), x.shape[1], x.shape[0] - 1)
    if max_pcs < 2:
        raise ValueError("Not enough cells/genes for clustering.")
    latent = PCA(
        n_components=max_pcs,
        svd_solver="randomized",
        random_state=seed,
    ).fit_transform(x)

    clusters = min(int(n_clusters), max(2, latent.shape[0] - 1))
    return MiniBatchKMeans(
        n_clusters=clusters,
        random_state=seed,
        n_init=20,
        batch_size=min(4096, max(256, latent.shape[0])),
    ).fit_predict(latent)


def calculate_cluster_metrics(
    observed_counts: np.ndarray,
    predicted_counts: np.ndarray,
    *,
    n_clusters: int = 10,
    n_pcs: int = 30,
    seed: int = 8667,
) -> tuple[float, float]:
    truth_labels = cluster_expression(
        observed_counts,
        n_clusters=n_clusters,
        n_pcs=n_pcs,
        seed=seed,
    )
    predicted_labels = cluster_expression(
        predicted_counts,
        n_clusters=n_clusters,
        n_pcs=n_pcs,
        seed=seed,
    )
    nmi = normalized_mutual_info_score(truth_labels, predicted_labels)
    ari = adjusted_rand_score(truth_labels, predicted_labels)
    return float(nmi), float(ari)


def build_fold_summary(
    gene_metrics: pd.DataFrame,
    *,
    experiment: str,
    model_key: str,
    model_label: str,
    fold: int,
    n_cells: int,
    n_observed_genes: int,
    n_heldout_genes: int,
    nmi: float,
    ari: float,
    training_seconds: float,
    prediction_seconds: float,
) -> pd.DataFrame:
    row = {
        "experiment": experiment,
        "model": model_key,
        "model_label": model_label,
        "fold": int(fold),
        "n_cells": int(n_cells),
        "n_observed_genes": int(n_observed_genes),
        "n_heldout_genes": int(n_heldout_genes),
        "scc_mean": float(gene_metrics["scc"].mean(skipna=True)),
        "scc_median": float(gene_metrics["scc"].median(skipna=True)),
        "ssim_mean": float(gene_metrics["ssim"].mean(skipna=True)),
        "ssim_median": float(gene_metrics["ssim"].median(skipna=True)),
        "rmse_mean": float(gene_metrics["rmse"].mean(skipna=True)),
        "mae_mean": float(gene_metrics["mae"].mean(skipna=True)),
        "jsd_mean": float(gene_metrics["jsd"].mean(skipna=True)),
        "moran_abs_error_mean": float(gene_metrics["moran_abs_error"].mean(skipna=True)),
        "nmi": float(nmi),
        "ari": float(ari),
        "training_seconds": float(training_seconds),
        "prediction_seconds": float(prediction_seconds),
    }
    return pd.DataFrame([row])


def _atomic_write_csv(dataframe: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".csv",
        prefix=path.stem + ".",
        dir=path.parent,
        delete=False,
    ) as handle:
        temp_path = Path(handle.name)
        dataframe.to_csv(handle, index=False)
    os.replace(temp_path, path)


def select_or_load_plot_genes(
    heldout_adata: ad.AnnData,
    path: Path,
    *,
    n_genes: int = 10,
    seed: int = 8667,
) -> list[str]:
    """Create one model-independent 10-gene list per fold.

    Genes are selected as medoids of 10 clusters in the four fold-construction
    metrics. The resulting CSV is shared by every imputation model.
    """
    path = Path(path)
    if path.exists():
        table = pd.read_csv(path)
        genes = table["gene"].astype(str).tolist()
        missing = pd.Index(genes).difference(heldout_adata.var_names.astype(str))
        if len(missing):
            raise ValueError(f"Stored plot-gene list has missing genes: {missing.tolist()}")
        return genes[:n_genes]

    features = [
        "log_mean_expression",
        "dispersion",
        "detection_fraction",
        "morans_I",
    ]
    available = [feature for feature in features if feature in heldout_adata.var.columns]

    var = heldout_adata.var.copy()
    var.index = heldout_adata.var_names.astype(str)

    if len(available) == 4 and heldout_adata.n_vars >= n_genes:
        ranked = var[available].rank(method="average", pct=True)
        scaled = StandardScaler().fit_transform(ranked)
        labels = KMeans(
            n_clusters=n_genes,
            random_state=seed,
            n_init=50,
        ).fit_predict(scaled)

        selected = []
        for cluster in range(n_genes):
            positions = np.where(labels == cluster)[0]
            center = scaled[positions].mean(axis=0)
            distance = np.sum((scaled[positions] - center) ** 2, axis=1)
            selected.append(var.index[positions[int(np.argmin(distance))]])
    else:
        matrix = to_dense_float32(heldout_adata.X)
        order = np.argsort(-np.var(matrix, axis=0))
        selected = var.index[order[:n_genes]].tolist()

    selected_table = pd.DataFrame(
        {
            "plot_order": np.arange(1, len(selected) + 1),
            "gene": selected,
        }
    )
    for feature in available:
        selected_table[feature] = var.loc[selected, feature].to_numpy()

    _atomic_write_csv(selected_table, path)
    return selected


def plot_fold_metric_summary(
    summary,
    output_path,
    title="Evaluation metrics",
):
    """Plot benchmark metrics without mixing incompatible units."""

    from pathlib import Path

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    output_path = Path(output_path)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if isinstance(summary, pd.Series):
        frame = summary.to_frame().T
    else:
        frame = pd.DataFrame(summary).copy()

    if frame.empty:
        raise ValueError(
            "Cannot plot an empty metric summary."
        )

    row = frame.iloc[0]

    metric_aliases = {
        "SCC": [
            "scc_mean",
            "scc",
        ],
        "SSIM": [
            "ssim_mean",
            "ssim",
        ],
        "NMI": [
            "nmi",
            "oof_nmi",
            "fold_nmi_mean",
        ],
        "ARI": [
            "ari",
            "oof_ari",
            "fold_ari_mean",
        ],
        "RMSE": [
            "rmse_mean",
            "rmse",
        ],
        "MAE": [
            "mae_mean",
            "mae",
        ],
        "Jensen–Shannon divergence": [
            "jsd_mean",
            "jsd",
        ],
        "Moran's I absolute error": [
            "moran_abs_error_mean",
            "moran_abs_error",
        ],
    }

    def get_metric(label):
        for column in metric_aliases[label]:
            if (
                column in row.index
                and pd.notna(row[column])
            ):
                return float(row[column])

        return np.nan

    groups = [
        {
            "title": "Higher is better",
            "subtitle": "Dimensionless similarity scores",
            "ylabel": "Score",
            "metrics": [
                "SCC",
                "SSIM",
                "NMI",
                "ARI",
            ],
        },
        {
            "title": "Lower is better",
            "subtitle": "Expected Xenium count scale",
            "ylabel": "Absolute count error",
            "metrics": [
                "RMSE",
                "MAE",
            ],
        },
        {
            "title": "Lower is better",
            "subtitle": "Dimensionless spatial/distribution error",
            "ylabel": "Divergence / absolute error",
            "metrics": [
                "Jensen–Shannon divergence",
                "Moran's I absolute error",
            ],
        },
    ]

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(20, 6),
    )

    for axis, group in zip(
        axes,
        groups,
    ):
        labels = []
        values = []

        for label in group["metrics"]:
            value = get_metric(label)

            if np.isfinite(value):
                labels.append(label)
                values.append(value)

        if not values:
            axis.axis("off")
            continue

        positions = np.arange(
            len(values)
        )

        bars = axis.bar(
            positions,
            values,
        )

        axis.set_xticks(
            positions
        )

        axis.set_xticklabels(
            labels,
            rotation=22,
            ha="right",
        )

        axis.set_title(
            group["title"]
            + "\n"
            + group["subtitle"]
        )

        axis.set_ylabel(
            group["ylabel"]
        )

        axis.grid(
            axis="y",
            alpha=0.25,
        )

        finite_maximum = max(values)

        upper_limit = (
            1.0
            if (
                group["ylabel"] == "Score"
                and finite_maximum <= 1.0
            )
            else finite_maximum * 1.22
        )

        if upper_limit <= 0:
            upper_limit = 1.0

        axis.set_ylim(
            0,
            upper_limit,
        )

        for bar, value in zip(
            bars,
            values,
        ):
            axis.text(
                bar.get_x()
                + bar.get_width() / 2,
                bar.get_height(),
                f"{value:.3g}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    figure.suptitle(
        title,
        fontsize=13,
    )

    figure.tight_layout(
        rect=[0, 0, 1, 0.94]
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )



def plot_gene_metric_distributions(
    gene_metrics: pd.DataFrame,
    output_path: Path,
    *,
    title: str,
) -> None:
    labels = {
        "scc": "SCC",
        "ssim": "SSIM",
        "rmse": "RMSE",
        "mae": "MAE",
        "jsd": "Jensen–Shannon divergence",
        "moran_abs_error": "Moran's I absolute error",
    }
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, metric in zip(axes.ravel(), GENE_METRIC_COLUMNS):
        values = gene_metrics[metric].replace([np.inf, -np.inf], np.nan).dropna()
        ax.hist(values, bins=20)
        ax.set_xlabel(labels[metric])
        ax.set_ylabel("Held-out genes")
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle(title)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_ten_gene_maps(
    observed_counts,
    predicted_counts,
    gene_names,
    selected_genes,
    coordinates,
    gene_metrics,
    png_path,
    pdf_path,
    *,
    experiment_label,
    model_label,
    fold=None,
):
    """
    Plot observed and predicted spatial expression for a fixed gene panel.

    The historical function name is retained for backward compatibility,
    but the function now supports any number of selected genes.

    Visualization scale:
        log1p(expected Xenium counts)

    For each gene, observed and predicted panels use the same gene-specific
    color limits, permitting direct within-gene visual comparison.
    """

    from pathlib import Path

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import sparse


    def to_dense(matrix):
        if sparse.issparse(matrix):
            matrix = matrix.toarray()

        return np.asarray(
            matrix,
            dtype=np.float32,
        )


    png_path = Path(png_path)
    pdf_path = Path(pdf_path)

    png_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pdf_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    observed = to_dense(
        observed_counts
    )

    predicted = to_dense(
        predicted_counts
    )

    coordinates = np.asarray(
        coordinates,
        dtype=np.float32,
    )

    gene_names = [
        str(gene)
        for gene in gene_names
    ]

    selected_genes = list(
        dict.fromkeys(
            str(gene)
            for gene in selected_genes
        )
    )

    if observed.shape != predicted.shape:
        raise ValueError(
            "Observed and predicted matrices differ in shape: "
            f"{observed.shape} versus {predicted.shape}"
        )

    if observed.ndim != 2:
        raise ValueError(
            f"Expression matrices must be two-dimensional: {observed.shape}"
        )

    if observed.shape[1] != len(gene_names):
        raise ValueError(
            "Number of gene names does not match matrix columns: "
            f"{len(gene_names)} versus {observed.shape[1]}"
        )

    if coordinates.ndim != 2 or coordinates.shape[0] != observed.shape[0]:
        raise ValueError(
            "Spatial-coordinate rows must match expression rows: "
            f"{coordinates.shape} versus {observed.shape}"
        )

    if coordinates.shape[1] < 2:
        raise ValueError(
            f"At least two spatial-coordinate columns are required: {coordinates.shape}"
        )

    coordinates = coordinates[:, :2]

    gene_to_index = {
        gene: index
        for index, gene in enumerate(gene_names)
    }

    available_genes = [
        gene
        for gene in selected_genes
        if gene in gene_to_index
    ]

    missing_genes = [
        gene
        for gene in selected_genes
        if gene not in gene_to_index
    ]

    if missing_genes:
        print(
            "WARNING: selected highlight genes missing from this matrix:",
            missing_genes,
        )

    if not available_genes:
        raise ValueError(
            "None of the selected highlight genes are available."
        )

    metrics = pd.DataFrame(
        gene_metrics
    ).copy()

    if not metrics.empty and "gene" in metrics.columns:
        metrics["gene"] = metrics["gene"].astype(str)
        metrics = (
            metrics
            .drop_duplicates("gene")
            .set_index("gene")
        )
    else:
        metrics = pd.DataFrame()

    n_genes = len(
        available_genes
    )

    figure_height = max(
        4.0,
        3.15 * n_genes,
    )

    figure, axes = plt.subplots(
        n_genes,
        2,
        figsize=(13, figure_height),
        squeeze=False,
    )

    point_size = float(
        np.clip(
            60000.0 / max(observed.shape[0], 1),
            0.25,
            4.0,
        )
    )

    for row_index, gene in enumerate(
        available_genes
    ):
        matrix_index = gene_to_index[
            gene
        ]

        observed_gene = np.clip(
            observed[:, matrix_index],
            0,
            None,
        )

        predicted_gene = np.clip(
            predicted[:, matrix_index],
            0,
            None,
        )

        observed_log = np.log1p(
            observed_gene
        )

        predicted_log = np.log1p(
            predicted_gene
        )

        combined_values = np.concatenate(
            [
                observed_log,
                predicted_log,
            ]
        )

        finite_values = combined_values[
            np.isfinite(combined_values)
        ]

        if finite_values.size == 0:
            vmin = 0.0
            vmax = 1.0
        else:
            vmin = float(
                np.nanpercentile(
                    finite_values,
                    1,
                )
            )

            vmax = float(
                np.nanpercentile(
                    finite_values,
                    99,
                )
            )

            if not np.isfinite(vmin):
                vmin = 0.0

            if not np.isfinite(vmax):
                vmax = 1.0

            if vmax <= vmin:
                vmax = vmin + 1e-6

        observed_axis = axes[
            row_index,
            0,
        ]

        predicted_axis = axes[
            row_index,
            1,
        ]

        observed_scatter = observed_axis.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            c=observed_log,
            s=point_size,
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )

        predicted_scatter = predicted_axis.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            c=predicted_log,
            s=point_size,
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )

        expected_region = None
        metric_text = []

        if (
            not metrics.empty
            and gene in metrics.index
        ):
            metric_row = metrics.loc[
                gene
            ]

            if isinstance(
                metric_row,
                pd.DataFrame,
            ):
                metric_row = metric_row.iloc[0]

            for region_column in (
                "expected_region",
                "layer",
                "region",
            ):
                if (
                    region_column in metric_row.index
                    and pd.notna(
                        metric_row[
                            region_column
                        ]
                    )
                ):
                    expected_region = str(
                        metric_row[
                            region_column
                        ]
                    )

                    break

            metric_labels = [
                ("SCC", "scc"),
                ("SSIM", "ssim"),
                ("RMSE", "rmse"),
                ("MAE", "mae"),
                ("JSD", "jsd"),
                ("Moran error", "moran_abs_error"),
            ]

            for label, column in metric_labels:
                if (
                    column in metric_row.index
                    and pd.notna(
                        metric_row[column]
                    )
                ):
                    metric_text.append(
                        f"{label}={float(metric_row[column]):.3f}"
                    )

        gene_label = gene

        if expected_region:
            gene_label += (
                f" — {expected_region}"
            )

        observed_axis.set_title(
            f"{gene_label}\nObserved Xenium",
            fontsize=10,
        )

        predicted_axis.set_title(
            f"{gene_label}\nPredicted: {model_label}",
            fontsize=10,
        )

        if metric_text:
            predicted_axis.set_xlabel(
                " | ".join(
                    metric_text
                ),
                fontsize=7.5,
            )

        for axis in (
            observed_axis,
            predicted_axis,
        ):
            axis.set_aspect(
                "equal"
            )

            axis.set_xticks([])
            axis.set_yticks([])

            for spine in axis.spines.values():
                spine.set_visible(
                    False
                )

        colorbar = figure.colorbar(
            predicted_scatter,
            ax=[
                observed_axis,
                predicted_axis,
            ],
            fraction=0.018,
            pad=0.012,
        )

        colorbar.set_label(
            "log1p(expected counts)",
            fontsize=8,
        )

        colorbar.ax.tick_params(
            labelsize=7
        )

    fold_label = (
        f"Fold {fold}"
        if fold is not None
        else "300-gene OOF"
    )

    figure.suptitle(
        f"{experiment_label} | {model_label} | {fold_label}\n"
        f"Fixed advisor highlight panel: observed versus predicted",
        fontsize=14,
        y=0.998,
    )

    figure.subplots_adjust(
        top=0.975,
        bottom=0.02,
        left=0.03,
        right=0.93,
        hspace=0.42,
        wspace=0.08,
    )

    figure.savefig(
        png_path,
        dpi=220,
        bbox_inches="tight",
    )

    figure.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

    print(
        "Highlight-gene plot saved:",
        png_path,
    )

    print(
        "Genes plotted:",
        available_genes,
    )



def plot_three_fold_metrics(
    fold_metrics,
    output_path,
    title="Three-fold evaluation",
):
    """Plot fold metrics in groups with compatible units."""

    from pathlib import Path

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    output_path = Path(output_path)
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    frame = pd.DataFrame(
        fold_metrics
    ).copy()

    if frame.empty:
        raise ValueError(
            "Cannot plot an empty fold-metric table."
        )

    if "fold" not in frame.columns:
        frame["fold"] = np.arange(
            1,
            len(frame) + 1,
        )

    frame = frame.sort_values(
        "fold"
    )

    metric_aliases = {
        "SCC": [
            "scc_mean",
            "scc",
        ],
        "SSIM": [
            "ssim_mean",
            "ssim",
        ],
        "NMI": [
            "nmi",
        ],
        "ARI": [
            "ari",
        ],
        "RMSE": [
            "rmse_mean",
            "rmse",
        ],
        "MAE": [
            "mae_mean",
            "mae",
        ],
        "Jensen–Shannon divergence": [
            "jsd_mean",
            "jsd",
        ],
        "Moran's I absolute error": [
            "moran_abs_error_mean",
            "moran_abs_error",
        ],
    }

    def resolve_column(label):
        for column in metric_aliases[label]:
            if column in frame.columns:
                return column

        return None

    groups = [
        {
            "title": "Higher is better",
            "subtitle": "Dimensionless similarity scores",
            "ylabel": "Score",
            "metrics": [
                "SCC",
                "SSIM",
                "NMI",
                "ARI",
            ],
        },
        {
            "title": "Lower is better",
            "subtitle": "Expected Xenium count scale",
            "ylabel": "Absolute count error",
            "metrics": [
                "RMSE",
                "MAE",
            ],
        },
        {
            "title": "Lower is better",
            "subtitle": "Dimensionless spatial/distribution error",
            "ylabel": "Divergence / absolute error",
            "metrics": [
                "Jensen–Shannon divergence",
                "Moran's I absolute error",
            ],
        },
    ]

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(20, 6),
    )

    fold_values = frame[
        "fold"
    ].to_numpy()

    for axis, group in zip(
        axes,
        groups,
    ):
        plotted = 0

        for label in group["metrics"]:
            column = resolve_column(
                label
            )

            if column is None:
                continue

            values = pd.to_numeric(
                frame[column],
                errors="coerce",
            ).to_numpy(dtype=float)

            axis.plot(
                fold_values,
                values,
                marker="o",
                linewidth=2,
                label=label,
            )

            plotted += 1

        if plotted == 0:
            axis.axis("off")
            continue

        axis.set_title(
            group["title"]
            + "\n"
            + group["subtitle"]
        )

        axis.set_xlabel(
            "Fold"
        )

        axis.set_ylabel(
            group["ylabel"]
        )

        axis.set_xticks(
            fold_values
        )

        axis.grid(
            alpha=0.25,
        )

        axis.legend(
            fontsize=8,
        )

        if group["ylabel"] == "Score":
            axis.set_ylim(
                0,
                1,
            )

    figure.suptitle(
        title,
        fontsize=13,
    )

    figure.tight_layout(
        rect=[0, 0, 1, 0.94]
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

