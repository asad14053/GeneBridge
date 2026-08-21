#!/usr/bin/env python3

"""
Compare the fixed advisor highlight-gene panel across Experiments 5, 5.1,
and 5.3 for one imputation model.

The spatial figure contains:

    Observed Xenium | Experiment 5 | Experiment 5.1 | Experiment 5.3

Each gene uses a shared color scale across observed and all predictions.
The displayed expression scale is log1p(expected Xenium counts).

This script is model-agnostic and can be reused for:
    VISTA, gimVI, Tangram, ENVI, TransImp, and SpaGE.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

DEFAULT_OUTPUT_ROOT = (
    PROJECT_ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

DEFAULT_HIGHLIGHT_FILE = (
    PROJECT_ROOT
    / "data"
    / "metadata"
    / "imputation"
    / "advisor_highlight_genes.csv"
)

EXPERIMENTS = [
    {
        "key": "ex5",
        "label": "Experiment 5",
        "short_label": "Ex-5",
    },
    {
        "key": "ex5_1",
        "label": "Experiment 5.1",
        "short_label": "Ex-5.1",
    },
    {
        "key": "ex5_3",
        "label": "Experiment 5.3",
        "short_label": "Ex-5.3",
    },
]

METRIC_COLUMNS = [
    "scc",
    "ssim",
    "rmse",
    "mae",
    "jsd",
    "moran_abs_error",
]

HIGHER_IS_BETTER = {
    "scc": True,
    "ssim": True,
    "rmse": False,
    "mae": False,
    "jsd": False,
    "moran_abs_error": False,
}

METRIC_LABELS = {
    "scc": "SCC",
    "ssim": "SSIM",
    "rmse": "RMSE",
    "mae": "MAE",
    "jsd": "JSD",
    "moran_abs_error": "Moran's I absolute error",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--model",
        required=True,
        help="Model directory key, for example vista or gimvi.",
    )

    parser.add_argument(
        "--model-label",
        required=True,
        help="Display label, for example VISTA or gimVI.",
    )

    parser.add_argument(
        "--combined-dir",
        default="auto",
        help=(
            "Combined result directory name. "
            "Use auto, combined_v2, or combined. "
            "Auto prefers combined_v2."
        ),
    )

    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
    )

    parser.add_argument(
        "--highlight-file",
        type=Path,
        default=DEFAULT_HIGHLIGHT_FILE,
    )

    return parser.parse_args()


def to_dense_float32(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()

    return np.asarray(
        matrix,
        dtype=np.float32,
    )


def resolve_combined_dir(
    model_root: Path,
    requested: str,
) -> Path:
    if requested != "auto":
        candidate = model_root / requested

        if not candidate.is_dir():
            raise FileNotFoundError(candidate)

        return candidate

    candidates = [
        model_root / "combined_v2",
        model_root / "combined",
    ]

    for candidate in candidates:
        required_files = [
            candidate / "oof_predictions_300genes.h5ad",
            candidate / "oof_ground_truth_300genes.h5ad",
            candidate / "gene_level_metrics_300genes.csv",
        ]

        if all(path.is_file() for path in required_files):
            return candidate

    raise FileNotFoundError(
        "Could not find a complete combined_v2 or combined directory under:\n"
        f"{model_root}"
    )


def load_highlight_panel(
    highlight_file: Path,
) -> pd.DataFrame:
    if not highlight_file.is_file():
        raise FileNotFoundError(highlight_file)

    highlight = pd.read_csv(highlight_file)

    required_columns = {
        "plot_order",
        "gene",
        "expected_region",
        "source_category",
    }

    missing_columns = required_columns - set(highlight.columns)

    if missing_columns:
        raise KeyError(
            "Highlight-gene file is missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    highlight["gene"] = (
        highlight["gene"]
        .astype(str)
        .str.strip()
    )

    # NTNG2 is not in the 300-gene Xenium benchmark.
    highlight = highlight.loc[
        highlight["gene"].str.upper() != "NTNG2"
    ].copy()

    highlight = (
        highlight
        .sort_values("plot_order")
        .drop_duplicates("gene")
        .reset_index(drop=True)
    )

    highlight["plot_order"] = np.arange(
        1,
        len(highlight) + 1,
    )

    if len(highlight) != 11:
        raise ValueError(
            "Expected exactly 11 advisor-highlight genes after excluding "
            f"NTNG2; found {len(highlight)}:\n"
            f"{highlight['gene'].tolist()}"
        )

    return highlight


def load_experiment(
    output_root: Path,
    model: str,
    experiment: dict[str, str],
    requested_combined_dir: str,
) -> dict:
    model_root = (
        output_root
        / experiment["key"]
        / model
    )

    combined_dir = resolve_combined_dir(
        model_root=model_root,
        requested=requested_combined_dir,
    )

    prediction_path = (
        combined_dir
        / "oof_predictions_300genes.h5ad"
    )

    truth_path = (
        combined_dir
        / "oof_ground_truth_300genes.h5ad"
    )

    metrics_path = (
        combined_dir
        / "gene_level_metrics_300genes.csv"
    )

    for path in [
        prediction_path,
        truth_path,
        metrics_path,
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)

    prediction = ad.read_h5ad(prediction_path)
    truth = ad.read_h5ad(truth_path)

    prediction.obs_names = prediction.obs_names.astype(str)
    prediction.var_names = prediction.var_names.astype(str)
    truth.obs_names = truth.obs_names.astype(str)
    truth.var_names = truth.var_names.astype(str)

    if not prediction.obs_names.equals(truth.obs_names):
        raise ValueError(
            f"{experiment['label']}: prediction and truth cell order differs."
        )

    if not prediction.var_names.equals(truth.var_names):
        missing = truth.var_names.difference(
            prediction.var_names
        )

        if len(missing):
            raise ValueError(
                f"{experiment['label']}: prediction is missing genes: "
                f"{missing[:20].tolist()}"
            )

        prediction = prediction[
            :,
            truth.var_names,
        ].copy()

    if prediction.shape != truth.shape:
        raise ValueError(
            f"{experiment['label']}: prediction/truth shapes differ: "
            f"{prediction.shape} versus {truth.shape}"
        )

    if "count_scale" in prediction.layers:
        predicted_counts = to_dense_float32(
            prediction.layers["count_scale"]
        )

        prediction_source = "layers['count_scale']"
    else:
        predicted_counts = to_dense_float32(
            prediction.X
        )

        prediction_source = "X"

    truth_counts = to_dense_float32(
        truth.X
    )

    if "spatial" in truth.obsm:
        coordinates = np.asarray(
            truth.obsm["spatial"],
            dtype=np.float32,
        )
    elif "spatial" in prediction.obsm:
        coordinates = np.asarray(
            prediction.obsm["spatial"],
            dtype=np.float32,
        )
    else:
        raise KeyError(
            f"{experiment['label']}: obsm['spatial'] is missing."
        )

    metrics = pd.read_csv(metrics_path)

    if "gene" not in metrics.columns:
        raise KeyError(
            f"{metrics_path} does not contain a gene column."
        )

    metrics["gene"] = (
        metrics["gene"]
        .astype(str)
        .str.strip()
    )

    metrics = metrics.drop_duplicates(
        "gene"
    )

    missing_metric_columns = [
        metric
        for metric in METRIC_COLUMNS
        if metric not in metrics.columns
    ]

    if missing_metric_columns:
        raise KeyError(
            f"{experiment['label']} metric file is missing: "
            + ", ".join(missing_metric_columns)
        )

    return {
        **experiment,
        "combined_dir": combined_dir,
        "prediction_path": prediction_path,
        "truth_path": truth_path,
        "metrics_path": metrics_path,
        "prediction_source": prediction_source,
        "obs_names": truth.obs_names.astype(str).tolist(),
        "var_names": truth.var_names.astype(str).tolist(),
        "truth_counts": truth_counts,
        "predicted_counts": predicted_counts,
        "coordinates": coordinates[:, :2],
        "metrics": metrics,
    }


def validate_experiment_alignment(
    experiment_data: list[dict],
) -> None:
    reference = experiment_data[0]

    for current in experiment_data[1:]:
        if current["obs_names"] != reference["obs_names"]:
            raise ValueError(
                f"Cell order differs between {reference['label']} "
                f"and {current['label']}."
            )

        if current["var_names"] != reference["var_names"]:
            raise ValueError(
                f"Gene order differs between {reference['label']} "
                f"and {current['label']}."
            )

        if current["coordinates"].shape != reference["coordinates"].shape:
            raise ValueError(
                f"Coordinate shapes differ between {reference['label']} "
                f"and {current['label']}."
            )

        coordinate_difference = float(
            np.max(
                np.abs(
                    current["coordinates"]
                    - reference["coordinates"]
                )
            )
        )

        if coordinate_difference > 1e-6:
            raise ValueError(
                f"Spatial coordinates differ between experiments; "
                f"maximum difference={coordinate_difference}"
            )

        truth_difference = float(
            np.max(
                np.abs(
                    current["truth_counts"]
                    - reference["truth_counts"]
                )
            )
        )

        if truth_difference > 1e-6:
            raise ValueError(
                f"OOF ground truth differs between {reference['label']} "
                f"and {current['label']}; maximum difference={truth_difference}"
            )


def build_metric_tables(
    experiment_data: list[dict],
    highlight: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tables = []

    for current in experiment_data:
        table = current["metrics"].merge(
            highlight[
                [
                    "plot_order",
                    "gene",
                    "expected_region",
                    "source_category",
                ]
            ],
            on="gene",
            how="inner",
            validate="one_to_one",
        )

        # Aggregated gene-metric files may already contain experiment
        # and model metadata. Overwrite the experiment labels safely
        # instead of inserting duplicate columns.
        table["experiment"] = current["key"]
        table["experiment_label"] = current["label"]
        table["experiment_short_label"] = current["short_label"]

        front_columns = [
            "experiment",
            "experiment_label",
            "experiment_short_label",
        ]

        table = table[
            front_columns
            + [
                column
                for column in table.columns
                if column not in front_columns
            ]
        ]

        tables.append(table)

    metrics_long = pd.concat(
        tables,
        ignore_index=True,
    )

    metrics_long = metrics_long.sort_values(
        [
            "plot_order",
            "experiment",
        ]
    )

    summary_records = []

    for current in experiment_data:
        subset = metrics_long.loc[
            metrics_long["experiment"]
            == current["key"]
        ]

        record = {
            "experiment": current["key"],
            "experiment_label": current["label"],
            "experiment_short_label": current["short_label"],
            "n_highlight_genes": int(subset["gene"].nunique()),
        }

        for metric in METRIC_COLUMNS:
            record[f"{metric}_mean_11genes"] = float(
                subset[metric].mean()
            )

            record[f"{metric}_sd_11genes"] = float(
                subset[metric].std()
            )

        summary_records.append(record)

    summary = pd.DataFrame(
        summary_records
    )

    winner_records = []

    for gene in highlight["gene"]:
        gene_table = metrics_long.loc[
            metrics_long["gene"] == gene
        ]

        for metric in METRIC_COLUMNS:
            valid = gene_table.dropna(
                subset=[metric]
            )

            if valid.empty:
                continue

            if HIGHER_IS_BETTER[metric]:
                winner_index = valid[metric].idxmax()
            else:
                winner_index = valid[metric].idxmin()

            winner_row = valid.loc[
                winner_index
            ]

            winner_records.append(
                {
                    "gene": gene,
                    "expected_region": winner_row[
                        "expected_region"
                    ],
                    "metric": metric,
                    "metric_label": METRIC_LABELS[
                        metric
                    ],
                    "higher_is_better": HIGHER_IS_BETTER[
                        metric
                    ],
                    "best_experiment": winner_row[
                        "experiment"
                    ],
                    "best_experiment_label": winner_row[
                        "experiment_label"
                    ],
                    "best_value": float(
                        winner_row[metric]
                    ),
                }
            )

    winners = pd.DataFrame(
        winner_records
    )

    return metrics_long, summary, winners


def compact_metric_text(
    metrics: pd.DataFrame,
    gene: str,
) -> str:
    rows = metrics.loc[
        metrics["gene"] == gene
    ]

    if rows.empty:
        return ""

    row = rows.iloc[0]

    return (
        f"SCC={row['scc']:.3f} | "
        f"SSIM={row['ssim']:.3f}\n"
        f"RMSE={row['rmse']:.2f} | "
        f"MAE={row['mae']:.2f}"
    )


def create_spatial_comparison_plot(
    experiment_data: list[dict],
    highlight: pd.DataFrame,
    model_label: str,
    output_png: Path,
    output_pdf: Path,
) -> None:
    reference = experiment_data[0]

    gene_names = reference["var_names"]

    gene_to_index = {
        gene: index
        for index, gene in enumerate(gene_names)
    }

    selected_genes = highlight["gene"].tolist()

    missing_genes = [
        gene
        for gene in selected_genes
        if gene not in gene_to_index
    ]

    if missing_genes:
        raise ValueError(
            "The following advisor genes are absent from the "
            "300-gene OOF benchmark: "
            + ", ".join(missing_genes)
        )

    coordinates = reference["coordinates"]
    observed_counts = reference["truth_counts"]

    n_genes = len(selected_genes)
    n_columns = 1 + len(experiment_data)

    figure_width = 5.0 * n_columns
    figure_height = max(
        8.0,
        2.85 * n_genes,
    )

    figure, axes = plt.subplots(
        n_genes,
        n_columns,
        figsize=(
            figure_width,
            figure_height,
        ),
        squeeze=False,
    )

    point_size = float(
        np.clip(
            60000.0 / max(
                observed_counts.shape[0],
                1,
            ),
            0.25,
            4.0,
        )
    )

    column_titles = [
        "Observed Xenium",
        *[
            f"{current['short_label']}\n{model_label}"
            for current in experiment_data
        ],
    ]

    for column_index, title in enumerate(
        column_titles
    ):
        axes[0, column_index].set_title(
            title,
            fontsize=12,
            fontweight="bold",
            pad=10,
        )

    for row_index, highlight_row in highlight.iterrows():
        gene = str(
            highlight_row["gene"]
        )

        expected_region = str(
            highlight_row["expected_region"]
        )

        gene_index = gene_to_index[
            gene
        ]

        observed_gene = np.clip(
            observed_counts[:, gene_index],
            0,
            None,
        )

        predicted_gene_values = [
            np.clip(
                current["predicted_counts"][
                    :,
                    gene_index,
                ],
                0,
                None,
            )
            for current in experiment_data
        ]

        log_values = [
            np.log1p(observed_gene),
            *[
                np.log1p(values)
                for values in predicted_gene_values
            ],
        ]

        combined = np.concatenate(
            log_values
        )

        finite = combined[
            np.isfinite(combined)
        ]

        if finite.size == 0:
            vmin = 0.0
            vmax = 1.0
        else:
            vmin = float(
                np.percentile(
                    finite,
                    1,
                )
            )

            vmax = float(
                np.percentile(
                    finite,
                    99,
                )
            )

            if vmax <= vmin:
                vmax = vmin + 1e-6

        row_scatter = None

        for column_index, values in enumerate(
            log_values
        ):
            axis = axes[
                row_index,
                column_index,
            ]

            row_scatter = axis.scatter(
                coordinates[:, 0],
                coordinates[:, 1],
                c=values,
                s=point_size,
                vmin=vmin,
                vmax=vmax,
                linewidths=0,
                rasterized=True,
            )

            axis.set_aspect("equal")
            axis.set_xticks([])
            axis.set_yticks([])

            for spine in axis.spines.values():
                spine.set_visible(False)

            if column_index == 0:
                axis.set_ylabel(
                    f"{gene}\n{expected_region}",
                    rotation=0,
                    ha="right",
                    va="center",
                    fontsize=10,
                    labelpad=15,
                )
            else:
                current = experiment_data[
                    column_index - 1
                ]

                axis.set_xlabel(
                    compact_metric_text(
                        current["metrics"],
                        gene,
                    ),
                    fontsize=7.5,
                )

        colorbar = figure.colorbar(
            row_scatter,
            ax=axes[row_index, :].tolist(),
            fraction=0.012,
            pad=0.008,
        )

        colorbar.set_label(
            "log1p(expected counts)",
            fontsize=7,
        )

        colorbar.ax.tick_params(
            labelsize=6,
        )

    figure.suptitle(
        f"{model_label}: fixed 11-gene comparison across "
        "Experiments 5, 5.1, and 5.3\n"
        "Observed Xenium versus 300-gene out-of-fold predictions",
        fontsize=15,
        y=0.998,
    )

    figure.subplots_adjust(
        top=0.965,
        bottom=0.025,
        left=0.10,
        right=0.95,
        hspace=0.48,
        wspace=0.08,
    )

    figure.savefig(
        output_png,
        dpi=200,
        bbox_inches="tight",
    )

    figure.savefig(
        output_pdf,
        bbox_inches="tight",
    )

    plt.close(figure)


def write_text_report(
    output_file: Path,
    model_label: str,
    highlight: pd.DataFrame,
    summary: pd.DataFrame,
    experiment_data: list[dict],
) -> None:
    lines = [
        f"{model_label.upper()} ADVISOR 11-GENE COMPARISON REPORT",
        "=" * 110,
        "",
        "Experiments compared:",
        "  Experiment 5   : Br8667-only Huuki reference",
        "  Experiment 5.1 : all 10 Huuki brains, including Br8667",
        "  Experiment 5.3 : 9 Huuki brains, excluding Br8667",
        "",
        "Visualization scale:",
        "  log1p(expected Xenium counts)",
        "",
        "The observed and all three predicted panels use one shared",
        "gene-specific color scale within each row.",
        "",
        "Genes:",
    ]

    for _, row in highlight.iterrows():
        lines.append(
            f"  {int(row['plot_order']):2d}. "
            f"{row['gene']:<8} "
            f"{row['expected_region']}"
        )

    lines.extend(
        [
            "",
            "MEAN PERFORMANCE ACROSS THE 11 HIGHLIGHT GENES",
            "-" * 110,
        ]
    )

    display_columns = [
        "experiment_short_label",
        "scc_mean_11genes",
        "ssim_mean_11genes",
        "rmse_mean_11genes",
        "mae_mean_11genes",
        "jsd_mean_11genes",
        "moran_abs_error_mean_11genes",
    ]

    lines.append(
        summary[
            display_columns
        ].to_string(
            index=False
        )
    )

    lines.extend(
        [
            "",
            "INPUT DIRECTORIES",
            "-" * 110,
        ]
    )

    for current in experiment_data:
        lines.append(
            f"{current['short_label']}: "
            f"{current['combined_dir']}"
        )

    lines.append("")

    output_file.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_arguments()

    output_root = args.output_root.resolve()

    highlight = load_highlight_panel(
        args.highlight_file.resolve()
    )

    experiment_data = [
        load_experiment(
            output_root=output_root,
            model=args.model,
            experiment=experiment,
            requested_combined_dir=args.combined_dir,
        )
        for experiment in EXPERIMENTS
    ]

    validate_experiment_alignment(
        experiment_data
    )

    gene_set = set(
        experiment_data[0]["var_names"]
    )

    missing_highlight_genes = [
        gene
        for gene in highlight["gene"]
        if gene not in gene_set
    ]

    if missing_highlight_genes:
        raise ValueError(
            "Highlight genes missing from the 300-gene benchmark: "
            + ", ".join(missing_highlight_genes)
        )

    metrics_long, summary, winners = build_metric_tables(
        experiment_data=experiment_data,
        highlight=highlight,
    )

    report_dir = (
        output_root
        / "comparisons"
        / args.model
        / "ex5_ex5_1_ex5_3"
        / "advisor_11_genes"
    )

    report_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    png_path = (
        report_dir
        / "advisor_11_gene_ex5_ex5_1_ex5_3_comparison.png"
    )

    pdf_path = (
        report_dir
        / "advisor_11_gene_ex5_ex5_1_ex5_3_comparison.pdf"
    )

    metrics_long_path = (
        report_dir
        / "advisor_11_gene_metrics_by_experiment.csv"
    )

    summary_path = (
        report_dir
        / "advisor_11_gene_metric_summary.csv"
    )

    winners_path = (
        report_dir
        / "advisor_11_gene_best_experiment_by_metric.csv"
    )

    report_path = (
        report_dir
        / "advisor_11_gene_comparison_report.txt"
    )

    manifest_path = (
        report_dir
        / "advisor_11_gene_comparison_manifest.json"
    )

    create_spatial_comparison_plot(
        experiment_data=experiment_data,
        highlight=highlight,
        model_label=args.model_label,
        output_png=png_path,
        output_pdf=pdf_path,
    )

    metrics_long.to_csv(
        metrics_long_path,
        index=False,
    )

    summary.to_csv(
        summary_path,
        index=False,
    )

    winners.to_csv(
        winners_path,
        index=False,
    )

    write_text_report(
        output_file=report_path,
        model_label=args.model_label,
        highlight=highlight,
        summary=summary,
        experiment_data=experiment_data,
    )

    manifest = {
        "model": args.model,
        "model_label": args.model_label,
        "highlight_gene_count": len(highlight),
        "highlight_genes": highlight["gene"].tolist(),
        "experiments": {
            current["key"]: {
                "label": current["label"],
                "combined_dir": str(
                    current["combined_dir"].resolve()
                ),
                "prediction": str(
                    current["prediction_path"].resolve()
                ),
                "truth": str(
                    current["truth_path"].resolve()
                ),
                "gene_metrics": str(
                    current["metrics_path"].resolve()
                ),
                "prediction_source": current[
                    "prediction_source"
                ],
            }
            for current in experiment_data
        },
        "outputs": {
            "png": str(png_path.resolve()),
            "pdf": str(pdf_path.resolve()),
            "metrics_long": str(metrics_long_path.resolve()),
            "metric_summary": str(summary_path.resolve()),
            "metric_winners": str(winners_path.resolve()),
            "text_report": str(report_path.resolve()),
        },
    }

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        ),
        encoding="utf-8",
    )

    print("=" * 110)
    print("ADVISOR 11-GENE CROSS-EXPERIMENT REPORT COMPLETE")
    print("=" * 110)
    print("Model:", args.model_label)
    print("Genes:", highlight["gene"].tolist())
    print()
    print("Spatial comparison PNG:", png_path)
    print("Spatial comparison PDF:", pdf_path)
    print("Gene-level metrics:", metrics_long_path)
    print("11-gene summary:", summary_path)
    print("Best experiment table:", winners_path)
    print("Text report:", report_path)
    print("Manifest:", manifest_path)
    print()
    print("11-gene metric summary:")
    print(
        summary.to_string(
            index=False
        )
    )
    print("=" * 110)


if __name__ == "__main__":
    main()
