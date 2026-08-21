#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# =============================================================================
# PATHS
# =============================================================================

PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

INPUT_FILE = (
    BASE_DIR
    / "final_visualizations"
    / "panel_level"
    / "across_methods_by_expression"
    / "all_methods_all_experiments_gene_metrics.csv"
)

OUT_DIR = (
    BASE_DIR
    / "final_visualizations"
    / "panel_level"
    / "rmse_mae_logscale"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =============================================================================
# CONFIG
# =============================================================================

METHODS = [
    ("vista", "VISTA"),
    ("gimvi", "gimVI"),
    ("tangram", "Tangram"),
    ("envi", "ENVI"),
    ("spage", "SpaGE"),
    ("transimpspa", "TransImpSpa"),
]

EXPERIMENTS = [
    ("ex5", "Experiment 5"),
    ("ex5_1", "Experiment 5.1"),
    ("ex5_3", "Experiment 5.3"),
]

STRATA = [
    "Low",
    "Medium",
    "High",
]

METRICS = [
    "RMSE",
    "MAE",
]


# =============================================================================
# LOAD + VALIDATE
# =============================================================================

def load_data():

    if not INPUT_FILE.is_file():
        raise FileNotFoundError(
            f"Input not found:\n{INPUT_FILE}"
        )

    print("Input:")
    print(INPUT_FILE)

    df = pd.read_csv(
        INPUT_FILE
    )

    required = {
        "gene",
        "experiment",
        "method",
        "expression_stratum",
        "RMSE",
        "MAE",
    }

    missing = required - set(
        df.columns
    )

    if missing:
        raise KeyError(
            f"Missing columns: {sorted(missing)}"
        )

    for metric in METRICS:

        df[metric] = pd.to_numeric(
            df[metric],
            errors="raise",
        )

        values = df[
            metric
        ].to_numpy(
            dtype=float
        )

        if not np.isfinite(
            values
        ).all():
            raise ValueError(
                f"{metric}: non-finite values detected."
            )

        if np.min(values) < 0:
            raise ValueError(
                f"{metric}: negative values detected."
            )

        # plotting-only transformation
        df[
            f"log1p_{metric}"
        ] = np.log1p(
            values
        )

    return df


# =============================================================================
# RAW-ERROR TAIL SUMMARY
# =============================================================================

def make_tail_summary(
    df,
):

    rows = []

    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        for stratum in STRATA:

            for (
                method,
                method_label,
            ) in METHODS:

                block = df.loc[
                    (
                        df[
                            "experiment"
                        ]
                        == experiment
                    )
                    &
                    (
                        df[
                            "method"
                        ]
                        == method
                    )
                    &
                    (
                        df[
                            "expression_stratum"
                        ].astype(str)
                        == stratum
                    )
                ]

                if len(block) != 100:
                    raise ValueError(
                        f"{experiment}/{method}/{stratum}: "
                        f"expected 100 genes, found {len(block)}"
                    )

                for metric in METRICS:

                    values = (
                        block[
                            metric
                        ]
                        .to_numpy(
                            dtype=float
                        )
                    )

                    rows.append(
                        {
                            "experiment":
                                experiment,

                            "experiment_label":
                                experiment_label,

                            "expression_stratum":
                                stratum,

                            "method":
                                method,

                            "method_label":
                                method_label,

                            "metric":
                                metric,

                            "n_genes":
                                len(values),

                            "mean":
                                float(
                                    np.mean(values)
                                ),

                            "median":
                                float(
                                    np.median(values)
                                ),

                            "q75":
                                float(
                                    np.quantile(
                                        values,
                                        0.75,
                                    )
                                ),

                            "q90":
                                float(
                                    np.quantile(
                                        values,
                                        0.90,
                                    )
                                ),

                            "q95":
                                float(
                                    np.quantile(
                                        values,
                                        0.95,
                                    )
                                ),

                            "q99":
                                float(
                                    np.quantile(
                                        values,
                                        0.99,
                                    )
                                ),

                            "maximum":
                                float(
                                    np.max(values)
                                ),
                        }
                    )

    return pd.DataFrame(
        rows
    )


# =============================================================================
# PLOT
# =============================================================================

def make_plot(
    df,
    experiment,
    experiment_label,
):

    subset = df.loc[
        df[
            "experiment"
        ]
        == experiment
    ].copy()


    fig, axes = plt.subplots(
        3,
        2,
        figsize=(
            13,
            13,
        ),
        squeeze=False,
    )


    positions = np.arange(
        1,
        len(METHODS) + 1,
        dtype=float,
    )

    method_labels = [
        label
        for key, label in METHODS
    ]

    rng = np.random.default_rng(
        8667
    )


    for row_idx, stratum in enumerate(
        STRATA
    ):

        for col_idx, metric in enumerate(
            METRICS
        ):

            ax = axes[
                row_idx,
                col_idx
            ]

            transformed_metric = (
                f"log1p_{metric}"
            )

            distributions = []


            for (
                method,
                method_label,
            ) in METHODS:

                values = (
                    subset.loc[
                        (
                            subset[
                                "method"
                            ]
                            == method
                        )
                        &
                        (
                            subset[
                                "expression_stratum"
                            ].astype(str)
                            == stratum
                        ),
                        transformed_metric,
                    ]
                    .to_numpy(
                        dtype=float
                    )
                )

                if len(values) != 100:
                    raise ValueError(
                        f"{experiment}/{method}/{stratum}: "
                        f"expected 100 values, found {len(values)}"
                    )

                distributions.append(
                    values
                )


            # -----------------------------------------------------------------
            # Violin
            # -----------------------------------------------------------------

            violin = ax.violinplot(
                distributions,
                positions=positions,
                widths=0.80,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )

            for body in violin[
                "bodies"
            ]:
                body.set_alpha(
                    0.30
                )


            # -----------------------------------------------------------------
            # Boxplot
            # -----------------------------------------------------------------

            ax.boxplot(
                distributions,
                positions=positions,
                widths=0.24,
                showfliers=False,
                patch_artist=False,
                medianprops={
                    "linewidth": 1.6,
                },
            )


            # -----------------------------------------------------------------
            # Individual genes
            # -----------------------------------------------------------------

            for (
                position,
                values,
            ) in zip(
                positions,
                distributions,
            ):

                jitter = rng.uniform(
                    -0.10,
                    0.10,
                    size=len(values),
                )

                ax.scatter(
                    np.full(
                        len(values),
                        position,
                    )
                    + jitter,
                    values,
                    s=8,
                    alpha=0.22,
                    linewidths=0,
                )


            # -----------------------------------------------------------------
            # Labels
            # -----------------------------------------------------------------

            if row_idx == 0:

                ax.set_title(
                    (
                        f"log1p({metric})\n"
                        "lower is better"
                    ),
                    fontsize=12,
                    fontweight="bold",
                )


            if col_idx == 0:

                ax.set_ylabel(
                    (
                        f"{stratum} expression\n"
                        "Transformed error"
                    ),
                    fontsize=11,
                    fontweight="bold",
                )


            ax.set_xticks(
                positions
            )

            ax.set_xticklabels(
                method_labels,
                rotation=35,
                ha="right",
                fontsize=9,
            )

            ax.grid(
                axis="y",
                alpha=0.20,
            )

            ax.spines[
                "top"
            ].set_visible(False)

            ax.spines[
                "right"
            ].set_visible(False)


    fig.suptitle(
        (
            f"{experiment_label}: RMSE and MAE across six imputation methods\n"
            "log1p transformation reveals central distributions while "
            "raw-scale outliers are retained separately"
        ),
        fontsize=15,
        fontweight="bold",
        y=0.995,
    )


    fig.tight_layout(
        rect=[
            0,
            0,
            1,
            0.94,
        ]
    )


    png = (
        OUT_DIR
        / f"{experiment}_rmse_mae_log1p.png"
    )

    pdf = (
        OUT_DIR
        / f"{experiment}_rmse_mae_log1p.pdf"
    )


    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


    return (
        png,
        pdf,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 100)
    print(
        "RMSE / MAE LOG-SCALE CROSS-METHOD COMPARISON"
    )
    print("=" * 100)


    df = load_data()


    # -------------------------------------------------------------------------
    # Save raw tail statistics
    # -------------------------------------------------------------------------

    tail_summary = (
        make_tail_summary(
            df
        )
    )

    tail_file = (
        OUT_DIR
        / "rmse_mae_tail_summary.csv"
    )

    tail_summary.to_csv(
        tail_file,
        index=False,
    )


    figures = {}


    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        print()
        print(
            "Plotting:",
            experiment_label,
        )

        png, pdf = make_plot(
            df,
            experiment,
            experiment_label,
        )

        figures[
            experiment
        ] = {
            "png":
                str(png),

            "pdf":
                str(pdf),
        }


    # -------------------------------------------------------------------------
    # Compact catastrophic-tail table
    # -------------------------------------------------------------------------

    tail_rank = (
        tail_summary[
            [
                "experiment",
                "expression_stratum",
                "method_label",
                "metric",
                "median",
                "q95",
                "q99",
                "maximum",
            ]
        ]
        .copy()
    )


    tail_rank.to_csv(
        OUT_DIR
        / "rmse_mae_median_tail_comparison.csv",
        index=False,
    )


    print()
    print("=" * 100)
    print("TAIL SUMMARY")
    print("=" * 100)

    print(
        tail_rank.to_string(
            index=False
        )
    )


    manifest = {
        "input":
            str(INPUT_FILE),

        "metrics":
            [
                "RMSE",
                "MAE",
            ],

        "plot_transform":
            "log1p(raw metric)",

        "raw_scale_outputs_preserved":
            True,

        "tail_statistics":
            [
                "median",
                "q95",
                "q99",
                "maximum",
            ],

        "figures":
            figures,
    }


    with (
        OUT_DIR
        / "manifest.json"
    ).open(
        "w"
    ) as handle:

        json.dump(
            manifest,
            handle,
            indent=2,
        )


    (
        OUT_DIR
        / "plot_complete.flag"
    ).write_text(
        "PASS\n"
    )


    print()
    print("=" * 100)
    print("DONE")
    print("=" * 100)

    print(
        "Tail summary:",
        tail_file,
    )


if __name__ == "__main__":
    main()
