#!/usr/bin/env python3

from __future__ import annotations

import ast
import shutil
import textwrap
from datetime import datetime
from pathlib import Path


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

COMMON_DIR = (
    PROJECT_ROOT
    / "src"
    / "imputation"
    / "common"
)

TARGET_FILES = [
    COMMON_DIR / "benchmark_evaluation.py",
    COMMON_DIR / "benchmark_evaluation_v2.py",
]


NEW_FUNCTION = r'''
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
'''


def replace_function(
    source: str,
    function_name: str,
    replacement: str,
) -> str:
    tree = ast.parse(
        source
    )

    candidates = [
        node
        for node in tree.body
        if isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef,
            ),
        )
        and node.name == function_name
    ]

    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected exactly one {function_name}() definition; "
            f"found {len(candidates)}."
        )

    node = candidates[0]

    lines = source.splitlines(
        keepends=True
    )

    start = node.lineno - 1
    end = node.end_lineno

    replacement_text = (
        textwrap.dedent(
            replacement
        ).strip()
        + "\n\n"
    )

    lines[start:end] = [
        replacement_text
    ]

    return "".join(
        lines
    )


def main() -> None:
    timestamp = datetime.now().strftime(
        "%Y%m%d_%H%M%S"
    )

    patched = []

    for target_file in TARGET_FILES:
        if not target_file.exists():
            print(
                "Skipping missing file:",
                target_file,
            )

            continue

        backup_file = target_file.with_name(
            target_file.name
            + f".before_highlight_plot_{timestamp}"
        )

        shutil.copy2(
            target_file,
            backup_file,
        )

        source = target_file.read_text(
            encoding="utf-8"
        )

        updated = replace_function(
            source=source,
            function_name="plot_ten_gene_maps",
            replacement=NEW_FUNCTION,
        )

        target_file.write_text(
            updated,
            encoding="utf-8",
        )

        patched.append(
            target_file
        )

        print(
            "Patched:",
            target_file,
        )

        print(
            "Backup:",
            backup_file,
        )

    if not patched:
        raise RuntimeError(
            "No benchmark evaluation files were patched."
        )


if __name__ == "__main__":
    main()
