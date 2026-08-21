#!/usr/bin/env python3

from pathlib import Path
import json
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE = (
    ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

FINAL = (
    BASE
    / "final_visualizations"
)

OUT = (
    FINAL
    / "final_summary"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


GENE_RANK_FILE = (
    FINAL
    / "panel_level"
    / "across_methods_by_expression"
    / "cross_method_metric_ranks.csv"
)

TAIL_FILE = (
    FINAL
    / "panel_level"
    / "rmse_mae_logscale"
    / "rmse_mae_median_tail_comparison.csv"
)

POP_RANK_FILE = (
    FINAL
    / "cell_population"
    / "nmi_ari"
    / "nmi_ari_oof_ranks.csv"
)


METHODS = [
    ("vista", "VISTA"),
    ("gimvi", "gimVI"),
    ("tangram", "Tangram"),
    ("envi", "ENVI"),
    ("spage", "SpaGE"),
    ("transimpspa", "TransImpSpa"),
]

METHOD_LABEL = dict(METHODS)


GENE_METRICS = [
    "SCC",
    "SSIM",
    "RMSE",
    "MAE",
    "JSD",
    "Moran_error",
]

POP_METRICS = [
    "NMI",
    "ARI",
]

ALL_METRICS = (
    GENE_METRICS
    + POP_METRICS
)


# =============================================================================
# GENE-LEVEL RANKS
# =============================================================================

def summarize_gene_ranks():

    df = pd.read_csv(
        GENE_RANK_FILE
    )

    required = {
        "method",
        "method_label",
        "metric",
        "median_rank",
        "experiment",
        "expression_stratum",
    }

    missing = required - set(
        df.columns
    )

    if missing:
        raise KeyError(
            f"Missing gene-rank columns: {sorted(missing)}"
        )


    df = df[
        df["metric"].isin(
            GENE_METRICS
        )
    ].copy()


    # Expected:
    # 3 experiments x 3 expression strata
    # = 9 ranks per method per gene-level metric

    counts = (
        df.groupby(
            [
                "method",
                "metric",
            ]
        )
        .size()
    )

    bad = counts[
        counts != 9
    ]

    if len(bad):
        raise ValueError(
            "Expected 9 contexts per method/metric:\n"
            + bad.to_string()
        )


    summary = (
        df.groupby(
            [
                "method",
                "method_label",
                "metric",
            ],
            as_index=False,
        )
        .agg(
            mean_rank=(
                "median_rank",
                "mean",
            ),
            median_rank_across_contexts=(
                "median_rank",
                "median",
            ),
            best_rank=(
                "median_rank",
                "min",
            ),
            worst_rank=(
                "median_rank",
                "max",
            ),
        )
    )


    return summary


# =============================================================================
# NMI / ARI OOF RANKS
# =============================================================================

def summarize_population_ranks():

    df = pd.read_csv(
        POP_RANK_FILE
    )

    required = {
        "method",
        "method_label",
        "metric",
        "rank",
        "experiment",
        "value",
    }

    missing = required - set(
        df.columns
    )

    if missing:
        raise KeyError(
            f"Missing population-rank columns: {sorted(missing)}"
        )


    df = df[
        df["metric"].isin(
            POP_METRICS
        )
    ].copy()


    counts = (
        df.groupby(
            [
                "method",
                "metric",
            ]
        )
        .size()
    )

    bad = counts[
        counts != 3
    ]

    if len(bad):
        raise ValueError(
            "Expected 3 OOF experiment scores "
            "per method/population metric:\n"
            + bad.to_string()
        )


    summary = (
        df.groupby(
            [
                "method",
                "method_label",
                "metric",
            ],
            as_index=False,
        )
        .agg(
            mean_rank=(
                "rank",
                "mean",
            ),
            median_rank_across_contexts=(
                "rank",
                "median",
            ),
            mean_oof_score=(
                "value",
                "mean",
            ),
            min_oof_score=(
                "value",
                "min",
            ),
            max_oof_score=(
                "value",
                "max",
            ),
        )
    )


    return summary


# =============================================================================
# RMSE / MAE TAIL RISK
# =============================================================================

def summarize_tail_behavior():

    df = pd.read_csv(
        TAIL_FILE
    )


    for col in [
        "median",
        "q95",
        "q99",
        "maximum",
    ]:
        df[col] = pd.to_numeric(
            df[col],
            errors="raise",
        )


    # Protect against division by zero
    denom = df["median"].replace(
        0,
        np.nan,
    )


    df[
        "q99_to_median"
    ] = (
        df["q99"]
        / denom
    )


    df[
        "max_to_median"
    ] = (
        df["maximum"]
        / denom
    )


    summary = (
        df.groupby(
            [
                "method_label",
                "metric",
            ],
            as_index=False,
        )
        .agg(
            median_q99_to_median=(
                "q99_to_median",
                "median",
            ),
            worst_q99_to_median=(
                "q99_to_median",
                "max",
            ),
            worst_max_to_median=(
                "max_to_median",
                "max",
            ),
        )
    )


    return summary


# =============================================================================
# UMAP DISPLACEMENT
# =============================================================================

def summarize_umap():

    rows = []


    for method, method_label in METHODS:

        path = (
            FINAL
            / "cell_population"
            / "umap"
            / method
            / f"{method}_umap_displacement_summary.csv"
        )


        if not path.is_file():

            raise FileNotFoundError(
                path
            )


        df = pd.read_csv(
            path
        )


        required = {
            "experiment",
            "median_umap_displacement",
            "q95_umap_displacement",
        }

        missing = required - set(
            df.columns
        )

        if missing:

            raise KeyError(
                f"{path}: missing {sorted(missing)}"
            )


        if len(df) != 3:

            raise ValueError(
                f"{method}: expected 3 UMAP experiment rows, "
                f"found {len(df)}"
            )


        rows.append(
            {
                "method":
                    method,

                "method_label":
                    method_label,

                "mean_median_umap_displacement":
                    float(
                        df[
                            "median_umap_displacement"
                        ].mean()
                    ),

                "best_median_umap_displacement":
                    float(
                        df[
                            "median_umap_displacement"
                        ].min()
                    ),

                "worst_median_umap_displacement":
                    float(
                        df[
                            "median_umap_displacement"
                        ].max()
                    ),

                "mean_q95_umap_displacement":
                    float(
                        df[
                            "q95_umap_displacement"
                        ].mean()
                    ),
            }
        )


    return pd.DataFrame(
        rows
    )


# =============================================================================
# BUILD FINAL 8-METRIC RANK MATRIX
# =============================================================================

def build_rank_matrix(
    gene_summary,
    pop_summary,
):

    combined = pd.concat(
        [
            gene_summary[
                [
                    "method",
                    "method_label",
                    "metric",
                    "mean_rank",
                ]
            ],
            pop_summary[
                [
                    "method",
                    "method_label",
                    "metric",
                    "mean_rank",
                ]
            ],
        ],
        ignore_index=True,
    )


    matrix = combined.pivot(
        index=[
            "method",
            "method_label",
        ],
        columns="metric",
        values="mean_rank",
    )


    matrix = matrix.reindex(
        columns=ALL_METRICS
    )


    if matrix.isna().any().any():

        raise ValueError(
            "Missing metric ranks in final matrix:\n"
            + matrix.to_string()
        )


    # Equal weight for each of the eight metric families.
    matrix[
        "descriptive_mean_rank_8metrics"
    ] = matrix[
        ALL_METRICS
    ].mean(
        axis=1
    )


    matrix = (
        matrix
        .sort_values(
            "descriptive_mean_rank_8metrics"
        )
        .reset_index()
    )


    return matrix


# =============================================================================
# MASTER TABLE
# =============================================================================

def build_master_table(
    rank_matrix,
    tail_summary,
    umap_summary,
):

    master = rank_matrix.copy()


    # Tail summaries in wide form
    tail_wide = (
        tail_summary.pivot(
            index="method_label",
            columns="metric",
            values=[
                "median_q99_to_median",
                "worst_q99_to_median",
                "worst_max_to_median",
            ],
        )
    )


    tail_wide.columns = [
        f"{stat}_{metric}"
        for stat, metric
        in tail_wide.columns
    ]


    tail_wide = (
        tail_wide
        .reset_index()
    )


    master = master.merge(
        tail_wide,
        on="method_label",
        how="left",
        validate="one_to_one",
    )


    master = master.merge(
        umap_summary,
        on=[
            "method",
            "method_label",
        ],
        how="left",
        validate="one_to_one",
    )


    # Descriptive final order only.
    master[
        "descriptive_overall_rank"
    ] = np.arange(
        1,
        len(master) + 1,
    )


    return master


# =============================================================================
# HEATMAP
# =============================================================================

def make_heatmap(
    master,
):

    display_metrics = [
        "SCC",
        "SSIM",
        "RMSE",
        "MAE",
        "JSD",
        "Moran_error",
        "NMI",
        "ARI",
    ]


    display_labels = [
        "SCC",
        "SSIM",
        "RMSE",
        "MAE",
        "JSD",
        "Moran\nerror",
        "NMI",
        "ARI",
    ]


    values = (
        master[
            display_metrics
        ]
        .to_numpy(
            dtype=float
        )
    )


    methods = (
        master[
            "method_label"
        ]
        .tolist()
    )


    fig, ax = plt.subplots(
        figsize=(
            12,
            7,
        )
    )


    im = ax.imshow(
        values,
        aspect="auto",
        vmin=1,
        vmax=6,
        cmap="viridis_r",
    )


    ax.set_xticks(
        np.arange(
            len(display_labels)
        )
    )

    ax.set_xticklabels(
        display_labels,
        fontsize=11,
    )


    ax.set_yticks(
        np.arange(
            len(methods)
        )
    )

    ax.set_yticklabels(
        methods,
        fontsize=11,
    )


    # Annotate each cell
    for i in range(
        values.shape[0]
    ):

        for j in range(
            values.shape[1]
        ):

            ax.text(
                j,
                i,
                f"{values[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=10,
            )


    ax.set_title(
        (
            "Final imputation benchmark: "
            "mean method rank across eight metric families\n"
            "Lower rank is better"
        ),
        fontsize=15,
        fontweight="bold",
        pad=18,
    )


    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.035,
        pad=0.03,
    )

    cbar.set_label(
        "Mean rank",
        rotation=90,
    )


    fig.tight_layout()


    png = (
        OUT
        / "final_8metric_rank_heatmap.png"
    )

    pdf = (
        OUT
        / "final_8metric_rank_heatmap.pdf"
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


    return png, pdf


# =============================================================================
# MAIN
# =============================================================================

def main():

    print(
        "=" * 110
    )

    print(
        "FINAL IMPUTATION BENCHMARK SUMMARY"
    )

    print(
        "=" * 110
    )


    gene_summary = (
        summarize_gene_ranks()
    )


    pop_summary = (
        summarize_population_ranks()
    )


    tail_summary = (
        summarize_tail_behavior()
    )


    umap_summary = (
        summarize_umap()
    )


    rank_matrix = (
        build_rank_matrix(
            gene_summary,
            pop_summary,
        )
    )


    master = (
        build_master_table(
            rank_matrix,
            tail_summary,
            umap_summary,
        )
    )


    # -------------------------------------------------------------------------
    # Save all tables
    # -------------------------------------------------------------------------

    gene_summary.to_csv(
        OUT
        / "gene_metric_rank_summary.csv",
        index=False,
    )


    pop_summary.to_csv(
        OUT
        / "population_metric_rank_summary.csv",
        index=False,
    )


    tail_summary.to_csv(
        OUT
        / "rmse_mae_tail_robustness_summary.csv",
        index=False,
    )


    umap_summary.to_csv(
        OUT
        / "umap_displacement_summary_all_methods.csv",
        index=False,
    )


    rank_matrix.to_csv(
        OUT
        / "final_8metric_rank_matrix.csv",
        index=False,
    )


    master.to_csv(
        OUT
        / "FINAL_MASTER_BENCHMARK_TABLE.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Heatmap
    # -------------------------------------------------------------------------

    png, pdf = make_heatmap(
        master
    )


    # -------------------------------------------------------------------------
    # Print compact summary
    # -------------------------------------------------------------------------

    compact_cols = [
        "descriptive_overall_rank",
        "method_label",
        "SCC",
        "SSIM",
        "RMSE",
        "MAE",
        "JSD",
        "Moran_error",
        "NMI",
        "ARI",
        "descriptive_mean_rank_8metrics",
        "mean_median_umap_displacement",
    ]


    print()
    print(
        "FINAL DESCRIPTIVE RANK SUMMARY"
    )

    print()


    print(
        master[
            compact_cols
        ].to_string(
            index=False
        )
    )


    print()
    print(
        "IMPORTANT:"
    )

    print(
        "The overall rank is descriptive, not an "
        "inferential statistical claim."
    )

    print(
        "All eight metric families receive equal weight."
    )

    print(
        "UMAP displacement and RMSE/MAE tail behavior "
        "are reported separately and are NOT included "
        "in the 8-metric mean rank."
    )


    manifest = {
        "gene_level_metrics": GENE_METRICS,

        "population_level_metrics": POP_METRICS,

        "gene_metric_rank_definition":
            (
                "Mean of direction-aware median ranks "
                "across 3 experiments x 3 expression "
                "strata = 9 contexts."
            ),

        "nmi_ari_rank_definition":
            (
                "Mean rank of combined OOF NMI/ARI "
                "across the 3 experiments."
            ),

        "overall_descriptive_rank":
            (
                "Equal-weight mean of the 8 metric-family "
                "mean ranks."
            ),

        "not_in_overall_rank": [
            "UMAP displacement",
            "RMSE/MAE tail ratios",
        ],

        "heatmap_png": str(png),

        "heatmap_pdf": str(pdf),
    }


    with (
        OUT
        / "final_summary_manifest.json"
    ).open(
        "w"
    ) as handle:

        json.dump(
            manifest,
            handle,
            indent=2,
        )


    (
        OUT
        / "summary_complete.flag"
    ).write_text(
        "PASS\n"
    )


    print()
    print(
        "=" * 110
    )

    print(
        "DONE"
    )

    print(
        "=" * 110
    )

    print(
        "Master table:",
        OUT
        / "FINAL_MASTER_BENCHMARK_TABLE.csv"
    )

    print(
        "Heatmap:",
        png
    )


if __name__ == "__main__":
    main()
