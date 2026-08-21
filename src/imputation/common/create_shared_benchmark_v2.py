#!/usr/bin/env python3

from __future__ import annotations

import ast
import re
import shutil
import textwrap
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

BENCHMARK_V1 = (
    COMMON_DIR
    / "benchmark_evaluation.py"
)

BENCHMARK_V2 = (
    COMMON_DIR
    / "benchmark_evaluation_v2.py"
)

AGGREGATOR_V1 = (
    COMMON_DIR
    / "aggregate_three_folds.py"
)

AGGREGATOR_V2 = (
    COMMON_DIR
    / "aggregate_three_folds_v2.py"
)


PLOT_FOLD_METRIC_SUMMARY = r'''
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
'''


PLOT_THREE_FOLD_METRICS = r'''
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
'''


def replace_functions(
    source: str,
    replacements: dict[str, str],
) -> str:
    tree = ast.parse(
        source
    )

    function_nodes = {
        node.name: node
        for node in tree.body
        if isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef,
            ),
        )
    }

    missing = sorted(
        set(replacements)
        - set(function_nodes)
    )

    if missing:
        raise RuntimeError(
            "Required plotting functions were not found: "
            + ", ".join(missing)
        )

    lines = source.splitlines(
        keepends=True
    )

    nodes_to_replace = sorted(
        (
            function_nodes[name]
            for name in replacements
        ),
        key=lambda node: node.lineno,
        reverse=True,
    )

    for node in nodes_to_replace:
        replacement = textwrap.dedent(
            replacements[node.name]
        ).strip() + "\n\n"

        start = node.lineno - 1
        end = node.end_lineno

        lines[start:end] = [
            replacement
        ]

    return "".join(
        lines
    )


def build_benchmark_v2() -> None:
    if not BENCHMARK_V1.is_file():
        raise FileNotFoundError(
            BENCHMARK_V1
        )

    shutil.copy2(
        BENCHMARK_V1,
        BENCHMARK_V2,
    )

    source = BENCHMARK_V2.read_text(
        encoding="utf-8"
    )

    source = replace_functions(
        source,
        {
            "plot_fold_metric_summary": (
                PLOT_FOLD_METRIC_SUMMARY
            ),
            "plot_three_fold_metrics": (
                PLOT_THREE_FOLD_METRICS
            ),
        },
    )

    BENCHMARK_V2.write_text(
        source,
        encoding="utf-8",
    )


def build_aggregator_v2() -> None:
    if not AGGREGATOR_V1.is_file():
        raise FileNotFoundError(
            AGGREGATOR_V1
        )

    source = AGGREGATOR_V1.read_text(
        encoding="utf-8"
    )

    source = source.replace(
        'project_root / "src/imputation_beta/common"',
        'project_root / "src/imputation/common"',
    )

    source = source.replace(
        'project_root / "src/imputation/common"',
        'project_root / "src/imputation/common"',
    )

    source = re.sub(
        r"from\s+benchmark_evaluation\s+import",
        "from benchmark_evaluation_v2 import",
        source,
    )

    seed_argument = (
        '    parser.add_argument('
        '"--seed", type=int, default=8667)\n'
    )

    combined_argument = (
        '    parser.add_argument(\n'
        '        "--combined-dir-name",\n'
        '        default="combined_v2",\n'
        '        help="Output subdirectory under the model directory.",\n'
        '    )\n'
    )

    if "--combined-dir-name" not in source:
        if seed_argument not in source:
            raise RuntimeError(
                "Could not locate the aggregator seed argument."
            )

        source = source.replace(
            seed_argument,
            seed_argument + combined_argument,
            1,
        )

    old_combined = (
        '    combined_dir = model_root / "combined"\n'
    )

    new_combined = (
        "    combined_dir = "
        "model_root / args.combined_dir_name\n"
    )

    if old_combined in source:
        source = source.replace(
            old_combined,
            new_combined,
            1,
        )
    elif new_combined not in source:
        raise RuntimeError(
            "Could not locate combined_dir assignment."
        )

    if "src/imputation_beta/common" in source:
        raise RuntimeError(
            "The legacy common-directory path remains "
            "in aggregate_three_folds_v2.py."
        )

    AGGREGATOR_V2.write_text(
        source,
        encoding="utf-8",
    )


def main() -> None:
    COMMON_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    build_benchmark_v2()
    build_aggregator_v2()

    print("Created:")
    print(BENCHMARK_V2)
    print(AGGREGATOR_V2)


if __name__ == "__main__":
    main()
