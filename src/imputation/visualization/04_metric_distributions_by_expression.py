#!/usr/bin/env python3

from pathlib import Path
import argparse
import json

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# =============================================================================
# CONFIG
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

OUT_ROOT = (
    BASE_DIR
    / "final_visualizations"
    / "panel_level"
    / "expression_stratified"
)

COMBINED_DIR = "combined_v2"


METHODS = {
    "vista": "VISTA",
    "gimvi": "gimVI",
    "tangram": "Tangram",
    "envi": "ENVI",
    "spage": "SpaGE",
    "transimpspa": "TransImpSpa",
}


EXPERIMENTS = [
    ("ex5", "Ex5"),
    ("ex5_1", "Ex5.1"),
    ("ex5_3", "Ex5.3"),
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
# DEFINE EXPRESSION STRATA ONCE
# =============================================================================

def build_expression_strata():

    """
    Define gene-expression bins ONCE from observed Xenium expression.

    Canonical source:
        Ex5 / VISTA / combined_v2 / gene_level_metrics_300genes.csv

    The aggregate file already contains observed_mean.

    Genes are ranked by observed mean expression and divided into
    three equal-sized groups of 100 genes.
    """

    path = (
        BASE_DIR
        / "ex5"
        / "vista"
        / COMBINED_DIR
        / "gene_level_metrics_300genes.csv"
    )

    if not path.is_file():
        raise FileNotFoundError(path)

    print("=" * 100)
    print("BUILDING FIXED EXPRESSION STRATA")
    print("=" * 100)
    print("Canonical source:")
    print(path)

    df = pd.read_csv(path)

    required = {
        "gene",
        "observed_mean",
    }

    missing = required - set(df.columns)

    if missing:
        raise KeyError(
            f"Missing required columns: {sorted(missing)}"
        )

    gene_df = (
        df[
            [
                "gene",
                "observed_mean",
            ]
        ]
        .copy()
    )

    gene_df["gene"] = (
        gene_df["gene"]
        .astype(str)
    )

    gene_df["observed_mean"] = (
        pd.to_numeric(
            gene_df["observed_mean"],
            errors="raise",
        )
    )

    if len(gene_df) != 300:
        raise ValueError(
            f"Expected 300 genes, found {len(gene_df)}"
        )

    if gene_df["gene"].duplicated().any():
        raise ValueError(
            "Duplicate genes found in canonical metric table."
        )

    if not np.isfinite(
        gene_df["observed_mean"]
        .to_numpy(dtype=float)
    ).all():
        raise ValueError(
            "Non-finite observed_mean values detected."
        )


    # -------------------------------------------------------------------------
    # Use log1p observed mean for interpretability.
    # Ranking is identical to raw observed_mean because log1p is monotonic.
    # -------------------------------------------------------------------------

    gene_df[
        "log1p_observed_mean"
    ] = np.log1p(
        gene_df[
            "observed_mean"
        ]
    )


    # -------------------------------------------------------------------------
    # Exactly 100 genes per stratum.
    #
    # Ranking first avoids qcut failures if several genes have identical means.
    # -------------------------------------------------------------------------

    ranks = (
        gene_df[
            "log1p_observed_mean"
        ]
        .rank(
            method="first",
            ascending=True,
        )
    )

    gene_df[
        "expression_stratum"
    ] = pd.qcut(
        ranks,
        q=3,
        labels=STRATA,
    )


    counts = (
        gene_df[
            "expression_stratum"
        ]
        .value_counts()
        .reindex(
            STRATA
        )
    )

    print()
    print("Genes per expression stratum:")
    print(counts)


    if not (
        counts.to_numpy()
        == 100
    ).all():

        raise ValueError(
            "Expected exactly 100 genes "
            "in each expression stratum."
        )


    # -------------------------------------------------------------------------
    # Save ranges
    # -------------------------------------------------------------------------

    ranges = (
        gene_df
        .groupby(
            "expression_stratum",
            observed=True,
        )
        .agg(
            n_genes=(
                "gene",
                "size",
            ),
            min_observed_mean=(
                "observed_mean",
                "min",
            ),
            median_observed_mean=(
                "observed_mean",
                "median",
            ),
            max_observed_mean=(
                "observed_mean",
                "max",
            ),
            min_log1p_mean=(
                "log1p_observed_mean",
                "min",
            ),
            median_log1p_mean=(
                "log1p_observed_mean",
                "median",
            ),
            max_log1p_mean=(
                "log1p_observed_mean",
                "max",
            ),
        )
        .reset_index()
    )


    print()
    print("Expression-stratum ranges:")
    print(
        ranges.to_string(
            index=False
        )
    )


    return (
        gene_df,
        ranges,
    )


# =============================================================================
# LOAD METHOD METRICS
# =============================================================================

def load_method_metrics(
    method,
    strata,
):

    frames = []


    for (
        experiment,
        experiment_label,
    ) in EXPERIMENTS:

        path = (
            BASE_DIR
            / experiment
            / method
            / COMBINED_DIR
            / "gene_level_metrics_300genes.csv"
        )

        if not path.is_file():
            raise FileNotFoundError(path)


        print()
        print(
            f"Loading {method} / "
            f"{experiment_label}"
        )

        print(path)


        source = pd.read_csv(
            path
        )


        if "gene" not in source.columns:
            raise KeyError(
                f"'gene' missing from {path}"
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
            }
        )


        for (
            display_name,
            source_name,
        ) in SOURCE_COLUMNS.items():

            if source_name not in source.columns:

                raise KeyError(
                    f"{source_name} missing from {path}\n"
                    f"Available columns: "
                    f"{list(source.columns)}"
                )

            out[
                display_name
            ] = pd.to_numeric(
                source[
                    source_name
                ],
                errors="raise",
            )


        # ---------------------------------------------------------------------
        # Merge FIXED observed-expression stratum
        # ---------------------------------------------------------------------

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

            missing_genes = (
                out.loc[
                    out[
                        "expression_stratum"
                    ].isna(),
                    "gene",
                ]
                .tolist()
            )

            raise ValueError(
                "Genes missing expression strata: "
                + ", ".join(
                    missing_genes
                )
            )


        frames.append(
            out
        )


    combined = pd.concat(
        frames,
        ignore_index=True,
    )


    return combined


# =============================================================================
# VALIDATE
# =============================================================================

def validate(
    df,
):

    print()
    print("=" * 100)
    print("VALIDATION")
    print("=" * 100)


    for experiment, _ in EXPERIMENTS:

        subset = df[
            df[
                "experiment"
            ]
            == experiment
        ]

        if len(subset) != 300:
            raise ValueError(
                f"{experiment}: expected 300 genes, "
                f"found {len(subset)}"
            )


        counts = (
            subset[
                "expression_stratum"
            ]
            .value_counts()
            .reindex(
                STRATA
            )
        )

        print()
        print(experiment)
        print(counts)

        if not (
            counts.to_numpy()
            == 100
        ).all():

            raise ValueError(
                f"{experiment}: expected "
                "100 genes/stratum."
            )


    for metric in METRICS:

        values = (
            df[
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
                f"{metric}: non-finite values."
            )


# =============================================================================
# SUMMARY
# =============================================================================

def create_summary(
    df,
):

    rows = []


    for metric in METRICS:

        for stratum in STRATA:

            for (
                experiment,
                experiment_label,
            ) in EXPERIMENTS:

                values = (
                    df.loc[
                        (
                            df[
                                "expression_stratum"
                            ].astype(str)
                            == stratum
                        )
                        &
                        (
                            df[
                                "experiment"
                            ]
                            == experiment
                        ),
                        metric,
                    ]
                    .to_numpy(
                        dtype=float
                    )
                )


                rows.append(
                    {
                        "metric":
                            metric,

                        "expression_stratum":
                            stratum,

                        "experiment":
                            experiment,

                        "experiment_label":
                            experiment_label,

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
                    }
                )


    return pd.DataFrame(
        rows
    )


# =============================================================================
# FIGURE
# =============================================================================

def create_plot(
    df,
    method,
    method_label,
    out_dir,
):

    fig, axes = plt.subplots(
        3,
        6,
        figsize=(
            25,
            12,
        ),
        squeeze=False,
    )


    positions = np.array(
        [
            1,
            2,
            3,
        ],
        dtype=float,
    )


    experiment_labels = [
        "Ex5",
        "Ex5.1",
        "Ex5.3",
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
                experiment,
                experiment_label,
            ) in EXPERIMENTS:

                values = (
                    df.loc[
                        (
                            df[
                                "experiment"
                            ]
                            == experiment
                        )
                        &
                        (
                            df[
                                "expression_stratum"
                            ].astype(str)
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
                widths=0.75,
                showmeans=False,
                showmedians=False,
                showextrema=False,
            )


            for body in violin[
                "bodies"
            ]:

                body.set_alpha(
                    0.35
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
            # Individual genes
            # -----------------------------------------------------------------

            for (
                position,
                values,
            ) in zip(
                positions,
                distributions,
            ):

                jitter = (
                    rng.uniform(
                        -0.10,
                        0.10,
                        size=len(values),
                    )
                )

                ax.scatter(
                    np.full(
                        len(values),
                        position,
                    )
                    + jitter,
                    values,
                    s=7,
                    alpha=0.25,
                    linewidths=0,
                )


            # -----------------------------------------------------------------
            # Formatting
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
                    (
                        f"{stratum} expression\n"
                        "Raw metric value"
                    ),
                    fontsize=10,
                    fontweight="bold",
                )


            ax.set_xticks(
                positions
            )

            ax.set_xticklabels(
                experiment_labels,
                fontsize=8,
            )


            ax.grid(
                axis="y",
                alpha=0.20,
            )


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

                max_value = max(
                    np.max(x)
                    for x in distributions
                )

                min_value = min(
                    np.min(x)
                    for x in distributions
                )

                if (
                    min_value >= 0
                    and
                    max_value <= 1
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
            f"{method_label}: Imputation performance "
            "stratified by observed mean expression\n"
            "300 Xenium genes divided into fixed "
            "low / medium / high expression tertiles"
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
        out_dir
        / f"{method}_metrics_by_expression.png"
    )

    pdf = (
        out_dir
        / f"{method}_metrics_by_expression.pdf"
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

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--method",
        required=True,
        choices=list(
            METHODS.keys()
        ),
    )

    args = parser.parse_args()


    method = args.method

    method_label = (
        METHODS[
            method
        ]
    )


    out_dir = (
        OUT_ROOT
        / method
    )

    out_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    print("=" * 100)
    print(
        "EXPRESSION-STRATIFIED "
        "GENE-LEVEL PERFORMANCE"
    )
    print("=" * 100)

    print(
        "Method:",
        method_label,
    )


    # =========================================================================
    # Fixed expression strata
    # =========================================================================

    strata, ranges = (
        build_expression_strata()
    )


    strata.to_csv(
        out_dir
        / "gene_expression_strata.csv",
        index=False,
    )


    ranges.to_csv(
        out_dir
        / "expression_stratum_ranges.csv",
        index=False,
    )


    # =========================================================================
    # Method metrics
    # =========================================================================

    metrics = load_method_metrics(
        method,
        strata,
    )


    validate(
        metrics
    )


    metrics.to_csv(
        out_dir
        / "gene_metrics_with_expression_strata.csv",
        index=False,
    )


    # =========================================================================
    # Summary
    # =========================================================================

    summary = create_summary(
        metrics
    )


    summary.to_csv(
        out_dir
        / "metric_summary_by_expression.csv",
        index=False,
    )


    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        summary[
            [
                "metric",
                "expression_stratum",
                "experiment_label",
                "mean",
                "median",
            ]
        ]
        .to_string(
            index=False
        )
    )


    # =========================================================================
    # Plot
    # =========================================================================

    png, pdf = create_plot(
        metrics,
        method,
        method_label,
        out_dir,
    )


    manifest = {
        "method":
            method,

        "method_label":
            method_label,

        "expression_definition":
            (
                "Observed Xenium mean expression; "
                "fixed tertiles from canonical "
                "Ex5/VISTA ground-truth metric table"
            ),

        "strata":
            STRATA,

        "genes_per_stratum":
            100,

        "metrics":
            METRICS,

        "experiments":
            [
                x[0]
                for x in EXPERIMENTS
            ],

        "png":
            str(png),

        "pdf":
            str(pdf),
    }


    with (
        out_dir
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
        out_dir
        / "plot_complete.flag"
    ).write_text(
        "PASS\n"
    )


    print()
    print("=" * 100)
    print("DONE")
    print("=" * 100)

    print(
        "Figure:",
        png,
    )


if __name__ == "__main__":
    main()
