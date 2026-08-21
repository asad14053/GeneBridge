#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

OUT_DIR = (
    BASE_DIR
    / "final_visualizations"
    / "panel_level"
    / "across_methods_by_expression"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

COMBINED_DIR = "combined_v2"


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


METRICS = [
    "SCC",
    "SSIM",
    "RMSE",
    "MAE",
    "JSD",
    "Moran_error",
]


SOURCE_COLUMNS = {
    "SCC": "scc",
    "SSIM": "ssim",
    "RMSE": "rmse",
    "MAE": "mae",
    "JSD": "jsd",
    "Moran_error": "moran_abs_error",
}


HIGHER_IS_BETTER = {
    "SCC": True,
    "SSIM": True,
    "RMSE": False,
    "MAE": False,
    "JSD": False,
    "Moran_error": False,
}


STRATA = [
    "Low",
    "Medium",
    "High",
]


# =============================================================================
# FIXED EXPRESSION STRATA
# =============================================================================

def build_expression_strata():

    path = (
        BASE_DIR
        / "ex5"
        / "vista"
        / COMBINED_DIR
        / "gene_level_metrics_300genes.csv"
    )

    df = pd.read_csv(path)

    required = {
        "gene",
        "observed_mean",
    }

    missing = required - set(df.columns)

    if missing:
        raise KeyError(
            f"Missing columns: {sorted(missing)}"
        )


    genes = (
        df[
            [
                "gene",
                "observed_mean",
            ]
        ]
        .copy()
    )

    genes["gene"] = (
        genes["gene"]
        .astype(str)
    )

    genes["observed_mean"] = (
        pd.to_numeric(
            genes["observed_mean"],
            errors="raise",
        )
    )


    if len(genes) != 300:
        raise ValueError(
            f"Expected 300 genes; found {len(genes)}"
        )


    genes[
        "log1p_observed_mean"
    ] = np.log1p(
        genes[
            "observed_mean"
        ]
    )


    ranks = (
        genes[
            "log1p_observed_mean"
        ]
        .rank(
            method="first",
            ascending=True,
        )
    )


    genes[
        "expression_stratum"
    ] = pd.qcut(
        ranks,
        q=3,
        labels=STRATA,
    )


    counts = (
        genes[
            "expression_stratum"
        ]
        .value_counts()
        .reindex(
            STRATA
        )
    )


    if not (
        counts.to_numpy()
        == 100
    ).all():
        raise ValueError(
            "Expected exactly 100 genes per stratum."
        )


    return genes


# =============================================================================
# LOAD ALL METHODS
# =============================================================================

def load_all_metrics(
    strata,
):

    frames = []


    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        for (
            method,
            method_label,
        ) in METHODS:

            path = (
                BASE_DIR
                / experiment
                / method
                / COMBINED_DIR
                / "gene_level_metrics_300genes.csv"
            )


            if not path.is_file():
                raise FileNotFoundError(path)


            print(
                f"Loading: "
                f"{experiment_label} | {method_label}"
            )


            source = pd.read_csv(
                path
            )


            out = pd.DataFrame(
                {
                    "gene":
                        source[
                            "gene"
                        ].astype(str),

                    "experiment":
                        experiment,

                    "experiment_label":
                        experiment_label,

                    "method":
                        method,

                    "method_label":
                        method_label,
                }
            )


            for (
                metric,
                source_column,
            ) in SOURCE_COLUMNS.items():

                if source_column not in source.columns:

                    raise KeyError(
                        f"{source_column} missing from:\n"
                        f"{path}"
                    )


                out[
                    metric
                ] = pd.to_numeric(
                    source[
                        source_column
                    ],
                    errors="raise",
                )


            out = out.merge(
                strata[
                    [
                        "gene",
                        "observed_mean",
                        "log1p_observed_mean",
                        "expression_stratum",
                    ]
                ],
                on="gene",
                how="left",
                validate="one_to_one",
            )


            if (
                out[
                    "expression_stratum"
                ]
                .isna()
                .any()
            ):
                raise ValueError(
                    f"Missing expression stratum: "
                    f"{experiment}/{method}"
                )


            frames.append(
                out
            )


    combined = pd.concat(
        frames,
        ignore_index=True,
    )


    # -------------------------------------------------------------------------
    # Strict validation
    # -------------------------------------------------------------------------

    expected_rows = (
        3
        * 6
        * 300
    )


    if len(combined) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} rows; "
            f"found {len(combined)}"
        )


    for metric in METRICS:

        values = (
            combined[
                metric
            ]
            .to_numpy(
                dtype=float
            )
        )

        if not np.isfinite(
            values
        ).all():
            raise ValueError(
                f"{metric} contains non-finite values."
            )


    return combined


# =============================================================================
# SUMMARY + RANKS
# =============================================================================

def create_summary(
    data,
):

    rows = []


    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        for stratum in STRATA:

            for metric in METRICS:

                for (
                    method,
                    method_label,
                ) in METHODS:

                    values = (
                        data.loc[
                            (
                                data[
                                    "experiment"
                                ]
                                == experiment
                            )
                            &
                            (
                                data[
                                    "method"
                                ]
                                == method
                            )
                            &
                            (
                                data[
                                    "expression_stratum"
                                ]
                                .astype(str)
                                == stratum
                            ),
                            metric,
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

                            "metric":
                                metric,

                            "method":
                                method,

                            "method_label":
                                method_label,

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

                            "std":
                                float(
                                    np.std(
                                        values,
                                        ddof=1,
                                    )
                                ),

                            "q25":
                                float(
                                    np.quantile(
                                        values,
                                        0.25,
                                    )
                                ),

                            "q75":
                                float(
                                    np.quantile(
                                        values,
                                        0.75,
                                    )
                                ),

                            "higher_is_better":
                                HIGHER_IS_BETTER[
                                    metric
                                ],
                        }
                    )


    summary = pd.DataFrame(
        rows
    )


    # -------------------------------------------------------------------------
    # Rank by median within each:
    # experiment × expression stratum × metric
    # -------------------------------------------------------------------------

    rank_frames = []


    for (
        experiment,
        stratum,
        metric,
    ), block in summary.groupby(
        [
            "experiment",
            "expression_stratum",
            "metric",
        ],
        observed=True,
    ):

        block = block.copy()


        ascending = (
            not HIGHER_IS_BETTER[
                metric
            ]
        )


        block = block.sort_values(
            "median",
            ascending=ascending,
        )


        block[
            "median_rank"
        ] = np.arange(
            1,
            len(block) + 1,
        )


        rank_frames.append(
            block
        )


    ranked = pd.concat(
        rank_frames,
        ignore_index=True,
    )


    return (
        summary,
        ranked,
    )


# =============================================================================
# FIGURE
# =============================================================================

def make_experiment_plot(
    data,
    experiment,
    experiment_label,
):

    subset = (
        data.loc[
            data[
                "experiment"
            ]
            == experiment
        ]
        .copy()
    )


    fig, axes = plt.subplots(
        3,
        6,
        figsize=(
            25,
            12,
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
                            ]
                            .astype(str)
                            == stratum
                        ),
                        metric,
                    ]
                    .to_numpy(
                        dtype=float
                    )
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
                widths=0.22,
                showfliers=False,
                patch_artist=False,
                medianprops={
                    "linewidth": 1.5,
                },
            )


            # -----------------------------------------------------------------
            # Gene points
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
                    s=6,
                    alpha=0.20,
                    linewidths=0,
                )


            # -----------------------------------------------------------------
            # Labels
            # -----------------------------------------------------------------

            display_metric = (
                "Moran error"
                if metric
                == "Moran_error"
                else metric
            )


            direction = (
                "higher is better"
                if HIGHER_IS_BETTER[
                    metric
                ]
                else
                "lower is better"
            )


            if row_idx == 0:

                ax.set_title(
                    f"{display_metric}\n"
                    f"({direction})",
                    fontsize=11,
                    fontweight="bold",
                )


            if col_idx == 0:

                ax.set_ylabel(
                    f"{stratum} expression\n"
                    "Raw metric value",
                    fontsize=10,
                    fontweight="bold",
                )


            ax.set_xticks(
                positions
            )


            ax.set_xticklabels(
                method_labels,
                rotation=45,
                ha="right",
                fontsize=7,
            )


            ax.grid(
                axis="y",
                alpha=0.20,
            )


            # -----------------------------------------------------------------
            # Defined ranges
            # -----------------------------------------------------------------

            if metric == "SCC":

                ax.set_ylim(
                    -1,
                    1,
                )

                ax.axhline(
                    0,
                    linewidth=0.7,
                    alpha=0.5,
                )


            elif metric in {
                "SSIM",
                "JSD",
            }:

                all_values = np.concatenate(
                    distributions
                )


                if (
                    np.min(
                        all_values
                    )
                    >= 0
                    and
                    np.max(
                        all_values
                    )
                    <= 1
                ):

                    ax.set_ylim(
                        0,
                        1,
                    )


            ax.spines[
                "top"
            ].set_visible(False)

            ax.spines[
                "right"
            ].set_visible(False)


    fig.suptitle(
        (
            f"{experiment_label}: Direct comparison "
            "of six imputation methods\n"
            "300 Xenium genes stratified by fixed "
            "observed-expression tertiles"
        ),
        fontsize=16,
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
        / f"{experiment}_methods_by_expression.png"
    )


    pdf = (
        OUT_DIR
        / f"{experiment}_methods_by_expression.pdf"
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
# COMPACT RANK TABLE
# =============================================================================

def create_rank_matrix(
    ranked,
):

    tables = []


    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        block = (
            ranked.loc[
                ranked[
                    "experiment"
                ]
                == experiment
            ]
            .copy()
        )


        block[
            "metric_stratum"
        ] = (
            block[
                "metric"
            ]
            .astype(str)
            + "__"
            + block[
                "expression_stratum"
            ]
            .astype(str)
        )


        matrix = (
            block.pivot(
                index="method_label",
                columns="metric_stratum",
                values="median_rank",
            )
        )


        matrix[
            "mean_rank"
        ] = (
            matrix.mean(
                axis=1
            )
        )


        matrix = (
            matrix.sort_values(
                "mean_rank"
            )
        )


        matrix.insert(
            0,
            "experiment",
            experiment_label,
        )


        matrix.insert(
            1,
            "method",
            matrix.index,
        )


        matrix = (
            matrix.reset_index(
                drop=True
            )
        )


        matrix.to_csv(
            OUT_DIR
            / f"{experiment}_median_rank_matrix.csv",
            index=False,
        )


        tables.append(
            matrix[
                [
                    "experiment",
                    "method",
                    "mean_rank",
                ]
            ]
        )


    all_ranks = pd.concat(
        tables,
        ignore_index=True,
    )


    return all_ranks


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 110)
    print(
        "DIRECT CROSS-METHOD COMPARISON "
        "BY EXPRESSION STRATUM"
    )
    print("=" * 110)


    strata = (
        build_expression_strata()
    )


    strata.to_csv(
        OUT_DIR
        / "fixed_gene_expression_strata.csv",
        index=False,
    )


    data = load_all_metrics(
        strata
    )


    data.to_csv(
        OUT_DIR
        / "all_methods_all_experiments_gene_metrics.csv",
        index=False,
    )


    summary, ranked = (
        create_summary(
            data
        )
    )


    summary.to_csv(
        OUT_DIR
        / "cross_method_metric_summary.csv",
        index=False,
    )


    ranked.to_csv(
        OUT_DIR
        / "cross_method_metric_ranks.csv",
        index=False,
    )


    figures = {}


    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        png, pdf = (
            make_experiment_plot(
                data,
                experiment,
                experiment_label,
            )
        )


        figures[
            experiment
        ] = {
            "png":
                str(png),

            "pdf":
                str(pdf),
        }


    overall_rank = (
        create_rank_matrix(
            ranked
        )
    )


    overall_rank.to_csv(
        OUT_DIR
        / "experiment_method_mean_ranks.csv",
        index=False,
    )


    print()
    print("=" * 110)
    print(
        "MEAN MEDIAN-RANK SUMMARY"
    )
    print("=" * 110)

    print(
        overall_rank.to_string(
            index=False
        )
    )


    manifest = {
        "methods":
            [
                x[1]
                for x in METHODS
            ],

        "experiments":
            [
                x[0]
                for x in EXPERIMENTS
            ],

        "metrics":
            METRICS,

        "expression_strata":
            STRATA,

        "genes_per_stratum":
            100,

        "ranking_basis":
            (
                "Median metric within each "
                "experiment × expression stratum; "
                "direction-aware rank. "
                "Mean rank is descriptive only."
            ),

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
    print("=" * 110)
    print("DONE")
    print("=" * 110)


if __name__ == "__main__":
    main()
