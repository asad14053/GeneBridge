#!/usr/bin/env python3

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ADVISOR_GENES = [
    "MBP",
    "RELN",
    "CUX2",
    "ADCYAP1",
    "VAMP1",
    "RORB",
    "PCP4",
    "NPTX1",
    "CCK",
    "PVALB",
    "ENC1",
]

EXPERIMENTS = {
    "Ex-5": "ex5",
    "Ex-5.1": "ex5_1",
    "Ex-5.3": "ex5_3",
}

METRIC_ALIASES = {
    "SCC": ["SCC", "scc"],
    "SSIM": ["SSIM", "ssim"],
    "RMSE": ["RMSE", "rmse"],
    "MAE": ["MAE", "mae"],
    "JSD": ["JSD", "jsd"],
    "Moran_error": [
        "Moran_error",
        "moran_error",
        "moran_abs_error",
        "Moran_abs_error",
    ],
}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--project-root",
        type=Path,
        required=True,
    )

    parser.add_argument(
        "--model",
        required=True,
    )

    parser.add_argument(
        "--model-label",
        required=True,
    )

    parser.add_argument(
        "--aggregate-subdir",
        default="combined_v2",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )

    return parser.parse_args()


def find_column(df, aliases):
    mapping = {
        str(c).lower(): c
        for c in df.columns
    }

    for alias in aliases:
        if alias.lower() in mapping:
            return mapping[
                alias.lower()
            ]

    return None


def find_gene_column(df):
    col = find_column(
        df,
        [
            "gene",
            "gene_name",
            "gene_symbol",
            "feature",
            "var_name",
        ],
    )

    if col is None:
        raise KeyError(
            f"Gene column not found. "
            f"Columns: {df.columns.tolist()}"
        )

    return col


def matrix_from_adata(adata):
    preferred_layers = [
        "count_scale",
        "counts",
        "raw_counts",
    ]

    for layer in preferred_layers:
        if layer in adata.layers:
            print(
                f"Using layer '{layer}' "
                f"for {adata.shape}"
            )
            return adata.layers[layer]

    print(
        f"Using X for {adata.shape}"
    )

    return adata.X


def extract_gene_vector(
    adata,
    matrix,
    gene,
):
    idx = adata.var_names.get_loc(
        gene
    )

    values = matrix[:, idx]

    if sparse.issparse(values):
        values = values.toarray()

    values = np.asarray(
        values
    ).reshape(-1)

    return values.astype(
        np.float64,
        copy=False,
    )


def spatial_coordinates(adata):
    if "spatial" in adata.obsm:
        coords = np.asarray(
            adata.obsm["spatial"]
        )

        if (
            coords.ndim == 2
            and coords.shape[1] >= 2
        ):
            return coords[:, :2]

    candidates = [
        ("x_centroid", "y_centroid"),
        ("x", "y"),
        ("X", "Y"),
        ("spatial_x", "spatial_y"),
        ("center_x", "center_y"),
    ]

    for x_col, y_col in candidates:
        if (
            x_col in adata.obs.columns
            and y_col in adata.obs.columns
        ):
            return np.column_stack(
                [
                    adata.obs[
                        x_col
                    ].to_numpy(),
                    adata.obs[
                        y_col
                    ].to_numpy(),
                ]
            )

    raise KeyError(
        "Could not locate spatial coordinates."
    )


def load_metric_tables(
    root,
    model,
    aggregate_subdir,
):
    frames = []

    for label, exp in EXPERIMENTS.items():

        path = (
            root
            / exp
            / model
            / aggregate_subdir
            / "gene_level_metrics_300genes.csv"
        )

        df = pd.read_csv(path)

        gene_col = find_gene_column(
            df
        )

        df = df.copy()

        df["gene"] = (
            df[gene_col]
            .astype(str)
        )

        available = [
            gene
            for gene in ADVISOR_GENES
            if gene in set(df["gene"])
        ]

        sub = (
            df[
                df["gene"]
                .isin(available)
            ]
            .copy()
        )

        # Aggregated metric tables may already contain these
        # metadata columns. Overwrite them rather than inserting
        # duplicate columns.
        sub["experiment"] = label
        sub["model"] = model

        # Keep metadata columns first for consistent downstream files.
        front_columns = [
            "experiment",
            "model",
            "gene",
        ]

        remaining_columns = [
            column
            for column in sub.columns
            if column not in front_columns
        ]

        sub = sub[
            front_columns + remaining_columns
        ]

        frames.append(sub)

    return pd.concat(
        frames,
        ignore_index=True,
    )


def make_metric_heatmaps(
    metrics_df,
    figure_dir,
    model_label,
):
    for metric, aliases in METRIC_ALIASES.items():

        col = find_column(
            metrics_df,
            aliases,
        )

        if col is None:
            print(
                f"Skipping {metric}: "
                "column not found"
            )
            continue

        table = (
            metrics_df[
                [
                    "gene",
                    "experiment",
                    col,
                ]
            ]
            .pivot(
                index="gene",
                columns="experiment",
                values=col,
            )
            .reindex(
                ADVISOR_GENES
            )
        )

        table = table.reindex(
            columns=[
                "Ex-5",
                "Ex-5.1",
                "Ex-5.3",
            ]
        )

        fig, ax = plt.subplots(
            figsize=(6.5, 8.5),
        )

        values = table.to_numpy(
            dtype=float
        )

        im = ax.imshow(
            values,
            aspect="auto",
        )

        ax.set_xticks(
            np.arange(
                table.shape[1]
            )
        )

        ax.set_xticklabels(
            table.columns
        )

        ax.set_yticks(
            np.arange(
                table.shape[0]
            )
        )

        ax.set_yticklabels(
            table.index
        )

        for i in range(
            table.shape[0]
        ):
            for j in range(
                table.shape[1]
            ):
                value = values[i, j]

                if np.isfinite(value):
                    ax.text(
                        j,
                        i,
                        f"{value:.3f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                    )

        fig.colorbar(
            im,
            ax=ax,
            label=metric,
        )

        ax.set_title(
            f"{model_label}: advisor 11 genes — {metric}"
        )

        fig.tight_layout()

        fig.savefig(
            figure_dir
            / f"advisor_11_{metric}_heatmap.png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.close(fig)


def main():
    args = parse_args()

    output_dir = args.output_dir

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure_dir = (
        output_dir
        / "figures"
    )

    figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    root = (
        args.project_root
        / "outputs"
        / "imputation_beta"
        / "Br8667"
    )

    print("=" * 100)
    print("ADVISOR 11-GENE COMPARISON")
    print("=" * 100)
    print("Model:", args.model_label)

    metrics_df = load_metric_tables(
        root,
        args.model,
        args.aggregate_subdir,
    )

    metric_output = (
        output_dir
        / "advisor_11_gene_metrics_long.csv"
    )

    metrics_df.to_csv(
        metric_output,
        index=False,
    )

    summary_rows = []

    for exp in [
        "Ex-5",
        "Ex-5.1",
        "Ex-5.3",
    ]:
        sub = metrics_df[
            metrics_df[
                "experiment"
            ] == exp
        ]

        row = {
            "experiment": exp,
            "n_genes": sub[
                "gene"
            ].nunique(),
        }

        for metric, aliases in METRIC_ALIASES.items():
            col = find_column(
                sub,
                aliases,
            )

            if col is not None:
                row[
                    f"{metric}_mean"
                ] = pd.to_numeric(
                    sub[col],
                    errors="coerce",
                ).mean()

        summary_rows.append(row)

    summary_df = pd.DataFrame(
        summary_rows
    )

    summary_output = (
        output_dir
        / "advisor_11_gene_mean_metrics_by_experiment.csv"
    )

    summary_df.to_csv(
        summary_output,
        index=False,
    )

    make_metric_heatmaps(
        metrics_df,
        figure_dir,
        args.model_label,
    )

    prediction_adatas = {}
    prediction_matrices = {}

    ground_truth = None
    ground_truth_matrix = None

    for label, exp in EXPERIMENTS.items():

        prediction_path = (
            root
            / exp
            / args.model
            / args.aggregate_subdir
            / "oof_predictions_300genes.h5ad"
        )

        truth_path = (
            root
            / exp
            / args.model
            / args.aggregate_subdir
            / "oof_ground_truth_300genes.h5ad"
        )

        print(
            f"\nLoading prediction: {prediction_path}"
        )

        pred = ad.read_h5ad(
            prediction_path
        )

        prediction_adatas[
            label
        ] = pred

        prediction_matrices[
            label
        ] = matrix_from_adata(
            pred
        )

        if ground_truth is None:

            print(
                f"Loading ground truth: {truth_path}"
            )

            ground_truth = ad.read_h5ad(
                truth_path
            )

            ground_truth_matrix = (
                matrix_from_adata(
                    ground_truth
                )
            )

    common_cells = list(
        ground_truth.obs_names
    )

    for label, pred in prediction_adatas.items():

        pred_cells = set(
            pred.obs_names
        )

        common_cells = [
            cell
            for cell in common_cells
            if cell in pred_cells
        ]

    print(
        "\nCommon cells:",
        f"{len(common_cells):,}",
    )

    ground_truth = ground_truth[
        common_cells
    ].copy()

    ground_truth_matrix = (
        matrix_from_adata(
            ground_truth
        )
    )

    for label in prediction_adatas:

        prediction_adatas[
            label
        ] = prediction_adatas[
            label
        ][common_cells].copy()

        prediction_matrices[
            label
        ] = matrix_from_adata(
            prediction_adatas[
                label
            ]
        )

    coords = spatial_coordinates(
        ground_truth
    )

    available_genes = [
        gene
        for gene in ADVISOR_GENES
        if (
            gene in ground_truth.var_names
            and all(
                gene in prediction_adatas[label].var_names
                for label in prediction_adatas
            )
        )
    ]

    missing_genes = [
        gene
        for gene in ADVISOR_GENES
        if gene not in available_genes
    ]

    print(
        "Genes plotted:",
        available_genes,
    )

    if missing_genes:
        print(
            "Missing/skipped:",
            missing_genes,
        )

    for gene in available_genes:

        truth = extract_gene_vector(
            ground_truth,
            ground_truth_matrix,
            gene,
        )

        vectors = {
            "Ground truth": truth,
        }

        for label in [
            "Ex-5",
            "Ex-5.1",
            "Ex-5.3",
        ]:

            vectors[label] = (
                extract_gene_vector(
                    prediction_adatas[
                        label
                    ],
                    prediction_matrices[
                        label
                    ],
                    gene,
                )
            )

        all_values = np.concatenate(
            list(
                vectors.values()
            )
        )

        transformed = np.log1p(
            np.clip(
                all_values,
                a_min=0,
                a_max=None,
            )
        )

        vmax = np.nanpercentile(
            transformed,
            99.5,
        )

        if (
            not np.isfinite(vmax)
            or vmax <= 0
        ):
            vmax = 1.0

        fig, axes = plt.subplots(
            1,
            4,
            figsize=(18, 5),
        )

        for ax, title in zip(
            axes,
            [
                "Ground truth",
                "Ex-5",
                "Ex-5.1",
                "Ex-5.3",
            ],
        ):

            values = np.log1p(
                np.clip(
                    vectors[title],
                    a_min=0,
                    a_max=None,
                )
            )

            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                s=1,
                vmin=0,
                vmax=vmax,
                rasterized=True,
            )

            ax.set_title(title)
            ax.set_aspect("equal")
            ax.axis("off")

            # Spatial image convention.
            ax.invert_yaxis()

        fig.colorbar(
            scatter,
            ax=axes.ravel().tolist(),
            shrink=0.75,
            label="log1p(count-scale expression)",
        )

        fig.suptitle(
            f"{args.model_label} — {gene}",
            fontsize=14,
        )

        fig.savefig(
            figure_dir
            / f"{gene}_groundtruth_vs_experiments.png",
            dpi=250,
            bbox_inches="tight",
        )

        plt.close(fig)

        print(
            "Saved spatial comparison:",
            gene,
        )

    print("\nAdvisor metric summary:")
    print(
        summary_df.to_string(
            index=False
        )
    )

    print("\nSaved:")
    print(metric_output)
    print(summary_output)
    print(figure_dir)

    print("\nSUCCESS")


if __name__ == "__main__":
    main()
