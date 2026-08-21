#!/usr/bin/env python3

from pathlib import Path
import argparse
import json

import anndata as ad
import numpy as np
from scipy import sparse

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


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
    / "gene_level"
    / "across_experiments"
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
    ("ex5", "Experiment 5"),
    ("ex5_1", "Experiment 5.1"),
    ("ex5_3", "Experiment 5.3"),
]


GENE_PANEL = [
    ("MBP", "WM"),
    ("RELN", "Layer 1"),
    ("CUX2", "Layer 2"),
    ("ADCYAP1", "Layer 3"),
    ("VAMP1", "Layer 4"),
    ("RORB", "Layer 4"),
    ("PCP4", "Layer 5"),
    ("NPTX1", "Layer 6"),
    ("CCK", "Interneuron"),
    ("PVALB", "Interneuron"),
    ("ENC1", "Cortical"),
]


GENES = [
    x[0]
    for x in GENE_PANEL
]


# Same scale as finalized across-method figures.
VMIN = 0.0
VMAX = 7.0


def dense_vector(matrix, idx):

    x = matrix[:, idx]

    if hasattr(x, "to_memory"):
        x = x.to_memory()

    if sparse.issparse(x):
        x = x.toarray()

    return (
        np.asarray(x)
        .reshape(-1)
        .astype(np.float64)
    )


def close_adata(adata):

    try:
        adata.file.close()
    except Exception:
        pass


def extract_expression(
    adata,
    prediction=False,
):

    missing = [
        gene
        for gene in GENES
        if gene not in adata.var_names
    ]

    if missing:
        raise ValueError(
            "Missing genes: "
            + ", ".join(missing)
        )

    if prediction:

        if "count_scale" not in adata.layers:
            raise KeyError(
                "Prediction H5AD is missing "
                "layers['count_scale']"
            )

        matrix = adata.layers[
            "count_scale"
        ]

    else:
        matrix = adata.X


    output = {}

    for gene in GENES:

        idx = int(
            adata.var_names.get_loc(
                gene
            )
        )

        values = dense_vector(
            matrix,
            idx,
        )

        if not np.isfinite(values).all():
            raise ValueError(
                f"Non-finite values for {gene}"
            )

        if np.min(values) < 0:
            raise ValueError(
                f"Negative values for {gene}"
            )

        output[gene] = values


    return output


def load_observed():

    # Ground truth is identical across all 18 audited runs.
    path = (
        BASE_DIR
        / "ex5"
        / "vista"
        / COMBINED_DIR
        / "oof_ground_truth_300genes.h5ad"
    )

    if not path.is_file():
        raise FileNotFoundError(path)

    adata = ad.read_h5ad(
        path,
        backed="r",
    )

    if "spatial" not in adata.obsm:
        raise KeyError(
            "'spatial' missing from observed H5AD"
        )

    obs_names = (
        adata.obs_names
        .astype(str)
        .to_numpy()
        .copy()
    )

    var_names = (
        adata.var_names
        .astype(str)
        .to_numpy()
        .copy()
    )

    coords = (
        np.asarray(
            adata.obsm["spatial"]
        )[:, :2]
        .copy()
    )

    expression = extract_expression(
        adata,
        prediction=False,
    )

    close_adata(
        adata
    )

    return {
        "path": path,
        "obs_names": obs_names,
        "var_names": var_names,
        "coords": coords,
        "expression": expression,
    }


def load_prediction(
    method,
    experiment,
    observed,
):

    path = (
        BASE_DIR
        / experiment
        / method
        / COMBINED_DIR
        / "oof_predictions_300genes.h5ad"
    )

    if not path.is_file():
        raise FileNotFoundError(path)

    print(
        "Loading:",
        path,
    )

    adata = ad.read_h5ad(
        path,
        backed="r",
    )

    obs_names = (
        adata.obs_names
        .astype(str)
        .to_numpy()
    )

    var_names = (
        adata.var_names
        .astype(str)
        .to_numpy()
    )

    if not np.array_equal(
        obs_names,
        observed["obs_names"],
    ):
        raise ValueError(
            f"{method}/{experiment}: "
            "cell order mismatch"
        )

    if not np.array_equal(
        var_names,
        observed["var_names"],
    ):
        raise ValueError(
            f"{method}/{experiment}: "
            "gene order mismatch"
        )

    if "spatial" not in adata.obsm:
        raise KeyError(
            f"{method}/{experiment}: "
            "'spatial' missing"
        )

    coords = (
        np.asarray(
            adata.obsm["spatial"]
        )[:, :2]
    )

    if not np.array_equal(
        coords,
        observed["coords"],
    ):
        raise ValueError(
            f"{method}/{experiment}: "
            "coordinate mismatch"
        )

    expression = extract_expression(
        adata,
        prediction=True,
    )

    close_adata(
        adata
    )

    return {
        "path": path,
        "expression": expression,
    }


def make_plot(
    method,
    method_label,
    observed,
    predictions,
    out_dir,
):

    n_rows = len(
        GENE_PANEL
    )

    column_labels = [
        "Observed Xenium",
        "Experiment 5",
        "Experiment 5.1",
        "Experiment 5.3",
    ]

    expression_sets = [
        observed["expression"],
        predictions["ex5"]["expression"],
        predictions["ex5_1"]["expression"],
        predictions["ex5_3"]["expression"],
    ]


    fig, axes = plt.subplots(
        n_rows,
        4,
        figsize=(
            11.5,
            25.5,
        ),
        squeeze=False,
    )


    x = observed["coords"][:, 0]
    y = observed["coords"][:, 1]

    norm = Normalize(
        vmin=VMIN,
        vmax=VMAX,
        clip=True,
    )


    for row_idx, (
        gene,
        region,
    ) in enumerate(
        GENE_PANEL
    ):

        for col_idx in range(4):

            ax = axes[
                row_idx,
                col_idx,
            ]

            values = np.log1p(
                np.maximum(
                    expression_sets[
                        col_idx
                    ][gene],
                    0,
                )
            )


            ax.scatter(
                x,
                y,
                c=values,
                cmap="viridis",
                norm=norm,
                s=0.30,
                linewidths=0,
                rasterized=True,
            )

            ax.set_aspect(
                "equal",
                adjustable="box",
            )

            ax.set_xticks([])
            ax.set_yticks([])

            for spine in ax.spines.values():
                spine.set_visible(False)


            if row_idx == 0:

                ax.set_title(
                    column_labels[
                        col_idx
                    ],
                    fontsize=12,
                    fontweight="bold",
                    pad=8,
                )


            if col_idx == 0:

                ax.set_ylabel(
                    f"{gene}\n{region}",
                    fontsize=11,
                    fontweight="bold",
                    rotation=0,
                    ha="right",
                    va="center",
                    labelpad=18,
                )


    # =========================================================================
    # Header
    # =========================================================================

    fig.text(
        0.5,
        0.985,
        (
            f"{method_label}: Spatial expression "
            "of the fixed 11-gene panel"
        ),
        ha="center",
        va="top",
        fontsize=16,
        fontweight="bold",
    )

    fig.text(
        0.5,
        0.962,
        (
            "Observed Xenium vs Experiments "
            "5, 5.1, and 5.3"
        ),
        ha="center",
        va="top",
        fontsize=14,
        fontweight="bold",
    )


    # =========================================================================
    # Exact same numeric legend style
    # =========================================================================

    shared_sm = ScalarMappable(
        norm=norm,
        cmap="viridis",
    )

    shared_sm.set_array([])


    cax = fig.add_axes(
        [
            0.385,
            0.920,
            0.230,
            0.010,
        ]
    )


    cbar = fig.colorbar(
        shared_sm,
        cax=cax,
        orientation="horizontal",
    )

    cbar.set_ticks(
        [
            0,
            2,
            4,
            6,
        ]
    )

    cbar.ax.tick_params(
        labelsize=9,
        length=3,
        pad=2,
    )

    cbar.set_label(
        "log1p(count)",
        fontsize=9,
        labelpad=3,
    )


    fig.subplots_adjust(
        left=0.13,
        right=0.985,
        top=0.865,
        bottom=0.02,
        wspace=0.035,
        hspace=0.08,
    )


    png = (
        out_dir
        / f"{method}_spatial_11genes_across_experiments.png"
    )

    pdf = (
        out_dir
        / f"{method}_spatial_11genes_across_experiments.pdf"
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

    plt.close(fig)

    return png, pdf


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
    method_label = METHODS[
        method
    ]


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
        "11-GENE SPATIAL PANEL "
        "ACROSS EXPERIMENTS"
    )
    print("=" * 100)

    print(
        "Method:",
        method_label,
    )


    observed = load_observed()

    predictions = {}

    for experiment, label in EXPERIMENTS:

        predictions[
            experiment
        ] = load_prediction(
            method,
            experiment,
            observed,
        )


    png, pdf = make_plot(
        method,
        method_label,
        observed,
        predictions,
        out_dir,
    )


    manifest = {
        "method": method,
        "method_label": method_label,
        "experiments": [
            "ex5",
            "ex5_1",
            "ex5_3",
        ],
        "genes": GENES,
        "display": "log1p(count_scale)",
        "display_vmin": VMIN,
        "display_vmax": VMAX,
        "png": str(png),
        "pdf": str(pdf),
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
    print("DONE")
    print("PNG:", png)
    print("PDF:", pdf)


if __name__ == "__main__":
    main()
