#!/usr/bin/env python
"""
09_plot_imputed_genes_Br8667.py

Plot spatial expression patterns for a small set of VISTA-imputed genes
in Xenium Br8667.

Inputs
------
1. Full VISTA-imputed matrix:
   data/processed/imputation_beta/Br8667/
       vista_Br8667_imputed_raw.h5ad

2. Original measured Xenium matrix:
   data/processed/imputation_beta/Br8667/
       spatial_data_xenium_Br8667_vista.h5ad

3. Day-1 gene statistics:
   outputs/imputation_beta/Br8667/qc/
       07_imputed_gene_statistics.csv

Outputs
-------
outputs/imputation_beta/Br8667/qc/imputed_gene_plots/
    09_selected_imputed_genes.csv
    09_imputed_gene_layer_means.csv
    09_imputed_gene_celltype_means.csv
    09_<GENE>_spatial.png
    09_imputed_gene_spatial_maps.pdf

Purpose
-------
The script plots only genes that were NOT measured by the original 300-gene
Xenium panel. Therefore, these maps show VISTA-imputed spatial patterns rather
than reconstruction of measured Xenium genes.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")
DATA_DIR = PROJECT_ROOT / "data/processed/imputation_beta/Br8667"
QC_DIR = PROJECT_ROOT / "outputs/imputation_beta/Br8667/qc"

DEFAULT_IMPUTED = DATA_DIR / "vista_Br8667_imputed_raw.h5ad"
DEFAULT_XENIUM = DATA_DIR / "spatial_data_xenium_Br8667_vista.h5ad"
DEFAULT_GENE_STATS = QC_DIR / "07_imputed_gene_statistics.csv"
DEFAULT_OUTPUT_DIR = QC_DIR / "imputed_gene_plots"

# Candidate genes chosen to represent distinct cortical layers and cell classes.
# The script keeps a candidate only when it is present in the imputed matrix,
# absent from the measured Xenium panel, and has non-empty predicted expression.
PREFERRED_MARKERS = [
    ("CUX2", "Upper-layer excitatory neuron"),
    ("RORB", "Layer 4 excitatory neuron"),
    ("BCL11B", "Deep-layer excitatory neuron"),
    ("TLE4", "Layer 6 corticothalamic neuron"),
    ("FOXP2", "Deep-layer neuron"),
    ("SLC17A7", "Excitatory neuron"),
    ("GAD1", "Inhibitory neuron"),
    ("GAD2", "Inhibitory neuron"),
    ("MBP", "Mature oligodendrocyte / white matter"),
    ("PLP1", "Mature oligodendrocyte / white matter"),
    ("MOG", "Mature oligodendrocyte"),
    ("AQP4", "Astrocyte"),
    ("GFAP", "Astrocyte"),
    ("P2RY12", "Microglia"),
    ("CX3CR1", "Microglia"),
    ("CLDN5", "Endothelial cell"),
    ("FLT1", "Endothelial cell"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot spatial patterns of selected VISTA-imputed Br8667 genes."
    )

    parser.add_argument(
        "--imputed",
        type=Path,
        default=DEFAULT_IMPUTED,
        help="Full VISTA-imputed h5ad file.",
    )

    parser.add_argument(
        "--xenium",
        type=Path,
        default=DEFAULT_XENIUM,
        help="Original measured Xenium h5ad file.",
    )

    parser.add_argument(
        "--gene-stats",
        type=Path,
        default=DEFAULT_GENE_STATS,
        help="Day-1 gene-statistics CSV.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for spatial plots and summary tables.",
    )

    parser.add_argument(
        "--genes",
        nargs="*",
        default=None,
        help=(
            "Optional manual gene list, for example: "
            "--genes CUX2 RORB TLE4 MBP AQP4 P2RY12"
        ),
    )

    parser.add_argument(
        "--n-genes",
        type=int,
        default=6,
        help="Number of genes to plot when using automatic selection. Default: 6.",
    )

    parser.add_argument(
        "--point-size",
        type=float,
        default=1.0,
        help="Scatter-point size. Default: 1.0.",
    )

    parser.add_argument(
        "--upper-quantile",
        type=float,
        default=0.99,
        help=(
            "Upper quantile used to cap each gene's color scale. "
            "Default: 0.99."
        ),
    )

    parser.add_argument(
        "--invert-y",
        action="store_true",
        help="Invert the y-axis when needed to match image-style orientation.",
    )

    return parser.parse_args()


def to_numpy(matrix, dtype=np.float32) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()

    return np.asarray(matrix, dtype=dtype)


def extract_backed_genes(
    adata: ad.AnnData,
    genes: list[str],
) -> np.ndarray:
    """
    Read selected gene columns from a backed AnnData file.

    h5py requires fancy-index positions to be sorted. The requested gene order
    is restored after loading.
    """
    positions = adata.var_names.get_indexer(genes)

    if np.any(positions < 0):
        missing = [
            genes[i]
            for i, position in enumerate(positions)
            if position < 0
        ]
        raise ValueError(
            f"Selected genes missing from imputed file: {missing}"
        )

    sort_order = np.argsort(positions)
    sorted_positions = positions[sort_order]

    matrix_sorted = to_numpy(
        adata[:, sorted_positions].X,
        dtype=np.float32,
    )

    inverse_order = np.argsort(sort_order)
    return matrix_sorted[:, inverse_order]


def make_case_insensitive_name_map(names) -> dict[str, str]:
    name_map: dict[str, str] = {}

    for name in names:
        text = str(name)
        name_map.setdefault(text.upper(), text)

    return name_map


def choose_genes(
    imputed: ad.AnnData,
    xenium: ad.AnnData,
    gene_stats: pd.DataFrame,
    requested_genes: list[str] | None,
    n_genes: int,
) -> pd.DataFrame:
    """
    Select unmeasured genes for visualization.

    Selection priority:
    1. User-requested genes, when provided.
    2. Biologically interpretable marker candidates.
    3. High-variance unmeasured genes as fallbacks.
    """
    imputed_map = make_case_insensitive_name_map(
        imputed.var_names
    )

    measured_upper = {
        str(gene).upper()
        for gene in xenium.var_names
    }

    stats = gene_stats.copy()

    required_columns = {
        "gene",
        "gene_group",
        "mean_expression",
        "variance_expression",
        "positive_fraction",
        "is_all_zero",
        "has_invalid_values",
        "has_negative_values",
    }

    missing_columns = (
        required_columns
        - set(stats.columns)
    )

    if missing_columns:
        raise KeyError(
            "Gene-statistics file is missing columns: "
            f"{sorted(missing_columns)}"
        )

    stats["gene"] = stats["gene"].astype(str)
    stats["gene_upper"] = stats["gene"].str.upper()

    eligible = stats.loc[
        (stats["gene_group"] == "unmeasured_imputed")
        & (~stats["is_all_zero"].astype(bool))
        & (~stats["has_invalid_values"].astype(bool))
        & (~stats["has_negative_values"].astype(bool))
        & (stats["mean_expression"] > 0)
        & (stats["variance_expression"] > 0)
    ].copy()

    eligible_by_upper = eligible.set_index(
        "gene_upper",
        drop=False,
    )

    selected: list[dict] = []
    used_upper: set[str] = set()

    def add_gene(
        requested_name: str,
        source: str,
        description: str,
    ) -> None:
        upper = requested_name.upper()

        if upper in used_upper:
            return

        if upper not in imputed_map:
            print(
                f"Skipping {requested_name}: "
                "not present in the imputed matrix."
            )
            return

        if upper in measured_upper:
            print(
                f"Skipping {requested_name}: "
                "already measured by the Xenium panel."
            )
            return

        if upper not in eligible_by_upper.index:
            print(
                f"Skipping {requested_name}: "
                "predicted expression is empty, invalid, or constant."
            )
            return

        row = eligible_by_upper.loc[upper]

        # A duplicate uppercase name is unlikely, but handle it safely.
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]

        selected.append(
            {
                "gene": str(row["gene"]),
                "selection_source": source,
                "biological_description": description,
                "mean_expression": float(
                    row["mean_expression"]
                ),
                "variance_expression": float(
                    row["variance_expression"]
                ),
                "positive_fraction": float(
                    row["positive_fraction"]
                ),
            }
        )

        used_upper.add(upper)

    if requested_genes:
        for gene in requested_genes:
            add_gene(
                requested_name=gene,
                source="manual",
                description="User-selected gene",
            )

    else:
        for gene, description in PREFERRED_MARKERS:
            add_gene(
                requested_name=gene,
                source="biological_candidate",
                description=description,
            )

            if len(selected) >= n_genes:
                break

    if len(selected) < n_genes:
        fallback = eligible.copy()

        # Favor genes with high spatial/cell-level variation while requiring
        # nonzero average predicted expression.
        fallback["selection_score"] = (
            fallback["variance_expression"]
            / (
                fallback["mean_expression"]
                + 1e-8
            )
        )

        fallback = fallback.sort_values(
            ["selection_score", "variance_expression"],
            ascending=False,
        )

        for _, row in fallback.iterrows():
            upper = str(row["gene"]).upper()

            if upper in used_upper:
                continue

            selected.append(
                {
                    "gene": str(row["gene"]),
                    "selection_source": "high_variance_fallback",
                    "biological_description": (
                        "Automatically selected high-variation "
                        "unmeasured gene"
                    ),
                    "mean_expression": float(
                        row["mean_expression"]
                    ),
                    "variance_expression": float(
                        row["variance_expression"]
                    ),
                    "positive_fraction": float(
                        row["positive_fraction"]
                    ),
                }
            )

            used_upper.add(upper)

            if len(selected) >= n_genes:
                break

    selected_table = pd.DataFrame(selected)

    if selected_table.empty:
        raise ValueError(
            "No eligible unmeasured genes were selected."
        )

    if requested_genes is None:
        selected_table = selected_table.iloc[
            :n_genes
        ].copy()

    return selected_table


def safe_filename(text: str) -> str:
    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        text,
    )


def grouped_means(
    metadata: pd.DataFrame,
    expression: np.ndarray,
    genes: list[str],
    group_column: str,
) -> pd.DataFrame:
    if group_column not in metadata.columns:
        return pd.DataFrame()

    groups = metadata[group_column].astype(str)
    records: list[dict] = []

    for group in sorted(groups.unique()):
        mask = groups.to_numpy() == group

        if not np.any(mask):
            continue

        means = expression[mask].mean(axis=0)

        for gene_index, gene in enumerate(genes):
            records.append(
                {
                    "group_column": group_column,
                    "group": group,
                    "gene": gene,
                    "n_cells": int(mask.sum()),
                    "mean_imputed_expression": float(
                        means[gene_index]
                    ),
                }
            )

    return pd.DataFrame(records)


def plot_one_gene(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    gene: str,
    description: str,
    output_path: Path,
    point_size: float,
    upper_quantile: float,
    invert_y: bool,
    pdf: PdfPages,
) -> None:
    log_values = np.log1p(
        np.clip(values, 0, None)
    )

    positive = log_values[
        log_values > 0
    ]

    if positive.size > 0:
        vmax = float(
            np.quantile(
                positive,
                upper_quantile,
            )
        )
    else:
        vmax = 1.0

    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(
            np.nanmax(log_values)
        )

    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    figure, axis = plt.subplots(
        figsize=(7.5, 7.5)
    )

    scatter = axis.scatter(
        x,
        y,
        c=log_values,
        s=point_size,
        linewidths=0,
        rasterized=True,
        vmin=0,
        vmax=vmax,
    )

    axis.set_aspect("equal")
    axis.set_xlabel("X coordinate")
    axis.set_ylabel("Y coordinate")
    axis.set_title(
        f"{gene}: VISTA-imputed spatial expression\n"
        f"{description}"
    )

    if invert_y:
        axis.invert_yaxis()

    colorbar = figure.colorbar(
        scatter,
        ax=axis,
        fraction=0.046,
        pad=0.04,
    )
    colorbar.set_label(
        "log1p predicted count rate"
    )

    figure.tight_layout()
    figure.savefig(
        output_path,
        dpi=250,
    )
    pdf.savefig(
        figure,
        dpi=250,
    )
    plt.close(figure)


def main() -> None:
    args = parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    for path in [
        args.imputed,
        args.xenium,
        args.gene_stats,
    ]:
        if not path.exists():
            raise FileNotFoundError(path)

    print("=" * 100)
    print(
        "Day 1 extension: Plot VISTA-imputed "
        "spatial gene-expression patterns"
    )
    print("=" * 100)

    print("\nOpening imputed file in backed mode...")
    imputed = ad.read_h5ad(
        args.imputed,
        backed="r",
    )

    print("Opening original Xenium file...")
    xenium = ad.read_h5ad(
        args.xenium,
        backed="r",
    )

    print("Loading Day-1 gene statistics...")
    gene_stats = pd.read_csv(
        args.gene_stats
    )

    if "spatial" not in imputed.obsm:
        raise KeyError(
            "The imputed file does not contain "
            "obsm['spatial']."
        )

    coordinates = np.asarray(
        imputed.obsm["spatial"],
        dtype=np.float32,
    )

    if coordinates.shape != (
        imputed.n_obs,
        2,
    ):
        raise ValueError(
            "Expected spatial coordinates with shape "
            f"({imputed.n_obs}, 2), received "
            f"{coordinates.shape}."
        )

    selected = choose_genes(
        imputed=imputed,
        xenium=xenium,
        gene_stats=gene_stats,
        requested_genes=args.genes,
        n_genes=args.n_genes,
    )

    genes = selected["gene"].tolist()

    print("\nSelected unmeasured genes:")
    print(
        selected[
            [
                "gene",
                "selection_source",
                "biological_description",
                "mean_expression",
                "variance_expression",
            ]
        ].to_string(index=False)
    )

    print("\nLoading selected predicted-expression columns...")
    expression = extract_backed_genes(
        imputed,
        genes,
    )

    if expression.shape != (
        imputed.n_obs,
        len(genes),
    ):
        raise ValueError(
            "Unexpected selected-expression shape: "
            f"{expression.shape}"
        )

    selected["minimum_expression"] = (
        expression.min(axis=0)
    )
    selected["maximum_expression"] = (
        expression.max(axis=0)
    )
    selected["median_expression"] = (
        np.median(expression, axis=0)
    )
    selected["zero_cell_count"] = (
        (expression == 0).sum(axis=0)
    )
    selected["positive_cell_count"] = (
        (expression > 0).sum(axis=0)
    )

    selected_path = (
        args.output_dir
        / "09_selected_imputed_genes.csv"
    )

    selected.to_csv(
        selected_path,
        index=False,
    )

    layer_means = grouped_means(
        metadata=imputed.obs,
        expression=expression,
        genes=genes,
        group_column="layer_annotation",
    )

    layer_path = (
        args.output_dir
        / "09_imputed_gene_layer_means.csv"
    )

    layer_means.to_csv(
        layer_path,
        index=False,
    )

    celltype_means = grouped_means(
        metadata=imputed.obs,
        expression=expression,
        genes=genes,
        group_column="scClassify",
    )

    celltype_path = (
        args.output_dir
        / "09_imputed_gene_celltype_means.csv"
    )

    celltype_means.to_csv(
        celltype_path,
        index=False,
    )

    pdf_path = (
        args.output_dir
        / "09_imputed_gene_spatial_maps.pdf"
    )

    x = coordinates[:, 0]
    y = coordinates[:, 1]

    with PdfPages(pdf_path) as pdf:
        for gene_index, row in selected.iterrows():
            gene = str(row["gene"])
            description = str(
                row["biological_description"]
            )

            png_path = (
                args.output_dir
                / (
                    "09_"
                    + safe_filename(gene)
                    + "_spatial.png"
                )
            )

            print(f"Plotting {gene}...")

            plot_one_gene(
                x=x,
                y=y,
                values=expression[
                    :,
                    gene_index,
                ],
                gene=gene,
                description=description,
                output_path=png_path,
                point_size=args.point_size,
                upper_quantile=args.upper_quantile,
                invert_y=args.invert_y,
                pdf=pdf,
            )

    print("\nCreated:")
    print(selected_path)
    print(layer_path)
    print(celltype_path)
    print(pdf_path)

    for gene in genes:
        print(
            args.output_dir
            / (
                "09_"
                + safe_filename(gene)
                + "_spatial.png"
            )
        )

    print(
        "\nThese genes were verified as absent from "
        "the original measured Xenium panel."
    )

    imputed.file.close()
    xenium.file.close()


if __name__ == "__main__":
    main()
