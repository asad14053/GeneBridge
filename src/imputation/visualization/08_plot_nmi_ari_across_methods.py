#!/usr/bin/env python3

from pathlib import Path
import json
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BASE = (
    ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

OUT_DIR = (
    BASE
    / "final_visualizations"
    / "cell_population"
    / "nmi_ari"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


EXPERIMENTS = [
    ("ex5", "Experiment 5"),
    ("ex5_1", "Experiment 5.1"),
    ("ex5_3", "Experiment 5.3"),
]

METHODS = [
    ("vista", "VISTA"),
    ("gimvi", "gimVI"),
    ("tangram", "Tangram"),
    ("envi", "ENVI"),
    ("spage", "SpaGE"),
    ("transimpspa", "TransImpSpa"),
]

METRICS = [
    ("NMI", "nmi"),
    ("ARI", "ari"),
]


def load_data():

    oof_rows = []
    fold_rows = []

    for exp, exp_label in EXPERIMENTS:

        for method, method_label in METHODS:

            run_dir = (
                BASE
                / exp
                / method
                / "combined_v2"
            )

            summary_file = (
                run_dir
                / "model_experiment_summary.csv"
            )

            fold_file = (
                run_dir
                / "fold_level_metrics_3folds.csv"
            )

            if not summary_file.is_file():
                raise FileNotFoundError(summary_file)

            if not fold_file.is_file():
                raise FileNotFoundError(fold_file)


            summary = pd.read_csv(summary_file)

            if len(summary) != 1:
                raise ValueError(
                    f"Expected 1 row: {summary_file}"
                )

            folds = pd.read_csv(fold_file)

            if len(folds) != 3:
                raise ValueError(
                    f"Expected 3 folds: {fold_file}"
                )


            for metric_label, metric_col in METRICS:

                if metric_col not in summary.columns:
                    raise KeyError(
                        f"{metric_col} missing: {summary_file}"
                    )

                if metric_col not in folds.columns:
                    raise KeyError(
                        f"{metric_col} missing: {fold_file}"
                    )


                oof_value = float(
                    summary.loc[0, metric_col]
                )

                oof_rows.append(
                    {
                        "experiment": exp,
                        "experiment_label": exp_label,
                        "method": method,
                        "method_label": method_label,
                        "metric": metric_label,
                        "value": oof_value,
                    }
                )


                for fold_idx, value in enumerate(
                    folds[metric_col].astype(float),
                    start=1,
                ):

                    fold_rows.append(
                        {
                            "experiment": exp,
                            "experiment_label": exp_label,
                            "method": method,
                            "method_label": method_label,
                            "metric": metric_label,
                            "fold": fold_idx,
                            "value": float(value),
                        }
                    )


    return (
        pd.DataFrame(oof_rows),
        pd.DataFrame(fold_rows),
    )


def rank_oof(oof):

    ranked = []

    for (
        experiment,
        experiment_label,
        metric,
    ), block in oof.groupby(
        [
            "experiment",
            "experiment_label",
            "metric",
        ],
        sort=False,
    ):

        block = (
            block
            .sort_values(
                "value",
                ascending=False,
            )
            .copy()
        )

        block["rank"] = np.arange(
            1,
            len(block) + 1,
        )

        ranked.append(block)

    return pd.concat(
        ranked,
        ignore_index=True,
    )


def make_figure(
    oof,
    folds,
):

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(17, 7),
        squeeze=False,
    )

    axes = axes[0]

    method_positions = np.arange(
        len(METHODS),
        dtype=float,
    )

    offsets = {
        "ex5": -0.22,
        "ex5_1": 0.0,
        "ex5_3": 0.22,
    }

    markers = {
        "ex5": "o",
        "ex5_1": "s",
        "ex5_3": "^",
    }


    for ax, (metric_label, _) in zip(
        axes,
        METRICS,
    ):

        for exp, exp_label in EXPERIMENTS:

            for method_idx, (
                method,
                method_label,
            ) in enumerate(METHODS):

                x = (
                    method_positions[method_idx]
                    + offsets[exp]
                )

                fold_values = (
                    folds.loc[
                        (
                            folds["experiment"]
                            == exp
                        )
                        &
                        (
                            folds["method"]
                            == method
                        )
                        &
                        (
                            folds["metric"]
                            == metric_label
                        ),
                        "value",
                    ]
                    .to_numpy(dtype=float)
                )

                oof_value = float(
                    oof.loc[
                        (
                            oof["experiment"]
                            == exp
                        )
                        &
                        (
                            oof["method"]
                            == method
                        )
                        &
                        (
                            oof["metric"]
                            == metric_label
                        ),
                        "value",
                    ]
                    .iloc[0]
                )


                # fold-level values
                jitter = np.linspace(
                    -0.035,
                    0.035,
                    len(fold_values),
                )

                ax.scatter(
                    np.full(
                        len(fold_values),
                        x,
                    )
                    + jitter,
                    fold_values,
                    s=24,
                    alpha=0.28,
                    linewidths=0,
                )


                # OOF primary score
                label = (
                    exp_label
                    if method_idx == 0
                    else None
                )

                ax.scatter(
                    x,
                    oof_value,
                    s=90,
                    marker=markers[exp],
                    edgecolors="black",
                    linewidths=0.7,
                    label=label,
                    zorder=4,
                )


        ax.set_title(
            f"{metric_label}\nHigher is better",
            fontsize=14,
            fontweight="bold",
        )

        ax.set_xticks(
            method_positions
        )

        ax.set_xticklabels(
            [
                label
                for _, label in METHODS
            ],
            rotation=30,
            ha="right",
        )

        ax.set_ylim(
            0,
            0.8,
        )

        ax.set_ylabel(
            metric_label
        )

        ax.grid(
            axis="y",
            alpha=0.20,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


    axes[0].legend(
        title="OOF score",
        frameon=False,
        loc="lower left",
    )


    fig.suptitle(
        (
            "Cell-population preservation across "
            "six imputation methods\n"
            "Large symbols = combined out-of-fold score; "
            "small points = three gene-CV folds"
        ),
        fontsize=16,
        fontweight="bold",
        y=1.02,
    )


    fig.tight_layout()


    png = OUT_DIR / "nmi_ari_across_methods.png"
    pdf = OUT_DIR / "nmi_ari_across_methods.pdf"

    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    plt.close(fig)


def main():

    print("=" * 90)
    print("NMI / ARI CROSS-METHOD COMPARISON")
    print("=" * 90)

    oof, folds = load_data()

    ranked = rank_oof(oof)


    oof.to_csv(
        OUT_DIR / "nmi_ari_oof_scores.csv",
        index=False,
    )

    folds.to_csv(
        OUT_DIR / "nmi_ari_fold_scores.csv",
        index=False,
    )

    ranked.to_csv(
        OUT_DIR / "nmi_ari_oof_ranks.csv",
        index=False,
    )


    print()
    print("OOF RANKING")
    print()

    print(
        ranked[
            [
                "experiment_label",
                "metric",
                "rank",
                "method_label",
                "value",
            ]
        ].to_string(index=False)
    )


    make_figure(
        oof,
        folds,
    )


    (
        OUT_DIR
        / "plot_complete.flag"
    ).write_text("PASS\n")


    manifest = {
        "primary_score":
            "combined out-of-fold NMI/ARI",

        "secondary_display":
            "three fold-level values",

        "interpretation":
            (
                "OOF values are primary. Fold points "
                "show variability and are not treated "
                "as gene-level observations."
            ),
    }

    with (
        OUT_DIR
        / "nmi_ari_manifest.json"
    ).open("w") as handle:

        json.dump(
            manifest,
            handle,
            indent=2,
        )


    print()
    print("DONE")
    print(OUT_DIR)


if __name__ == "__main__":
    main()
