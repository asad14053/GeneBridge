#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# Configuration
# =============================================================================

PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

OUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "celltype_comparison"
    / "02_celltype_proportions"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# -----------------------------------------------------------------------------
# Input datasets
# -----------------------------------------------------------------------------

HUUKI_FILE = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
    "huuki_snrna_reference_full_allgenes.h5ad"
)

BAITUK_CANDIDATES = [
    PROJECT_ROOT
    / "data/processed/snrnaseq/Baituk/"
    "baituk_snrna_reference_full_allgenes.h5ad",

    PROJECT_ROOT
    / "data/processed/snrnaseq/Batiuk/"
    "baituk_snrna_reference_full_allgenes.h5ad",
]

XENIUM_CANDIDATES = [
    PROJECT_ROOT
    / "data/processed/xenium/"
    "xenium_N24_layer_celltype_annotated.h5ad",

    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/"
    "xenium_Br8667_annotated.h5ad",

    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/"
    "spatial_data_xenium_Br8667_vista_qc.h5ad",
]


# Baituk med labels exported previously from annotations_final.RDS
BAITUK_MED_LABEL_FILE = (
    PROJECT_ROOT
    / "outputs/celltype_comparison/"
    "01_annotation_inventory/"
    "baituk_published_rds/"
    "cell_labels__med.csv"
)


# =============================================================================
# Canonical annotation columns
# =============================================================================

HUUKI_COLUMN = "cellType_broad_hc"
XENIUM_COLUMN = "cell_type_annotation"


# =============================================================================
# Harmonization dictionaries
# =============================================================================

HUUKI_MAP = {

    "Excit":
        "Excitatory neuron",

    "Inhib":
        "Inhibitory neuron",

    "Astro":
        "Non-neuronal",

    "Oligo":
        "Non-neuronal",

    "OPC":
        "Non-neuronal",

    "Micro":
        "Non-neuronal",

    "EndoMural":
        "Non-neuronal",

    "Ambiguous":
        "Other/Ambiguous",
}


BAITUK_MAP = {

    # -------------------------------------------------------------------------
    # Excitatory neurons
    # -------------------------------------------------------------------------

    "L2_CUX2_LAMP5":
        "Excitatory neuron",

    "L2_3_CUX2_FREM3":
        "Excitatory neuron",

    "L3_CUX2_PRSS12":
        "Excitatory neuron",

    "L4_RORB_SCHLAP1":
        "Excitatory neuron",

    "L4_5_FEZF2_LRRK1":
        "Excitatory neuron",

    "L5_6_FEZF2_TLE4":
        "Excitatory neuron",

    "L5_6_THEMIS":
        "Excitatory neuron",

    "L5_FEZF2_ADRA1A":
        "Excitatory neuron",


    # -------------------------------------------------------------------------
    # Inhibitory neurons
    # -------------------------------------------------------------------------

    "PVALB":
        "Inhibitory neuron",

    "SST":
        "Inhibitory neuron",

    "VIP":
        "Inhibitory neuron",

    "ID2_LAMP5":
        "Inhibitory neuron",

    "ID2_NCKAP5":
        "Inhibitory neuron",

    "ID2_PAX6":
        "Inhibitory neuron",


    # -------------------------------------------------------------------------
    # Broad non-neuronal / unresolved
    # -------------------------------------------------------------------------

    "Glia":
        "Non-neuronal",

    "Other":
        "Other/Ambiguous",
}


XENIUM_MAP = {

    # -------------------------------------------------------------------------
    # Excitatory neurons
    # -------------------------------------------------------------------------

    "L2/3 Ex":
        "Excitatory neuron",

    "L4/5 Ex":
        "Excitatory neuron",

    "L5 Ex":
        "Excitatory neuron",

    "L6 Ex":
        "Excitatory neuron",


    # -------------------------------------------------------------------------
    # Inhibitory neurons
    # -------------------------------------------------------------------------

    "CGE":
        "Inhibitory neuron",

    "MGE":
        "Inhibitory neuron",


    # -------------------------------------------------------------------------
    # Non-neuronal
    # -------------------------------------------------------------------------

    "Ast":
        "Non-neuronal",

    "Oligo":
        "Non-neuronal",

    "Mic":
        "Non-neuronal",

    "Endo":
        "Non-neuronal",


    # -------------------------------------------------------------------------
    # Ambiguous
    # -------------------------------------------------------------------------

    "Ambig/In/Endo":
        "Other/Ambiguous",

    "Ambig/Oligo":
        "Other/Ambiguous",
}


CATEGORY_ORDER = [
    "Excitatory neuron",
    "Inhibitory neuron",
    "Non-neuronal",
    "Other/Ambiguous",
]


# =============================================================================
# Helpers
# =============================================================================

def resolve_file(candidates, dataset):

    for path in candidates:
        if path.is_file():
            print(
                f"{dataset} file: {path}"
            )
            return path

    raise FileNotFoundError(
        f"No valid {dataset} file found.\n"
        + "\n".join(
            str(x)
            for x in candidates
        )
    )


def check_mapping(
    dataset,
    labels,
    mapping,
):

    observed = set(
        labels.dropna().astype(str).unique()
    )

    expected = set(
        mapping.keys()
    )

    unmapped = sorted(
        observed - expected
    )

    unused = sorted(
        expected - observed
    )

    print()
    print(
        f"{dataset} unique original labels: "
        f"{len(observed)}"
    )

    print(
        "Observed labels:"
    )

    for x in sorted(observed):
        print(
            "  ",
            x,
        )

    if unused:

        print()
        print(
            f"{dataset}: mapping entries not "
            "present in this dataset:"
        )

        for x in unused:
            print(
                "  ",
                x,
            )

    if unmapped:

        print()
        print(
            f"ERROR: {dataset} has unmapped labels:"
        )

        for x in unmapped:
            print(
                "  ",
                x,
            )

        raise ValueError(
            f"{dataset} contains unmapped "
            "cell-type labels."
        )


def summarize(
    dataset,
    original_labels,
    harmonized_labels,
):

    df = pd.DataFrame(
        {
            "dataset":
                dataset,

            "original_annotation":
                original_labels.astype(str),

            "harmonized_cell_type":
                harmonized_labels.astype(str),
        }
    )


    # -------------------------------------------------------------------------
    # Native/original annotation proportions
    # -------------------------------------------------------------------------

    original = (
        df.groupby(
            [
                "dataset",
                "original_annotation",
            ],
            observed=True,
        )
        .size()
        .reset_index(
            name="n_cells"
        )
    )

    original["proportion"] = (
        original["n_cells"]
        / original["n_cells"].sum()
    )

    original["percent"] = (
        original["proportion"]
        * 100
    )


    # -------------------------------------------------------------------------
    # Harmonized proportions
    # -------------------------------------------------------------------------

    harmonized = (
        df.groupby(
            [
                "dataset",
                "harmonized_cell_type",
            ],
            observed=True,
        )
        .size()
        .reset_index(
            name="n_cells"
        )
    )

    harmonized["proportion"] = (
        harmonized["n_cells"]
        / harmonized["n_cells"].sum()
    )

    harmonized["percent"] = (
        harmonized["proportion"]
        * 100
    )

    return (
        df,
        original,
        harmonized,
    )


# =============================================================================
# Huuki
# =============================================================================

def process_huuki():

    print()
    print("=" * 100)
    print("HUUKI")
    print("=" * 100)

    if not HUUKI_FILE.is_file():

        raise FileNotFoundError(
            HUUKI_FILE
        )

    adata = ad.read_h5ad(
        HUUKI_FILE,
        backed="r",
    )

    print(
        "Shape:",
        adata.shape,
    )

    if HUUKI_COLUMN not in adata.obs.columns:

        raise KeyError(
            f"{HUUKI_COLUMN} not found.\n"
            f"Available columns:\n"
            f"{list(adata.obs.columns)}"
        )

    labels = (
        adata.obs[
            HUUKI_COLUMN
        ]
        .copy()
    )

    try:
        adata.file.close()
    except Exception:
        pass


    check_mapping(
        "Huuki",
        labels,
        HUUKI_MAP,
    )


    harmonized = (
        labels.astype(str)
        .map(
            HUUKI_MAP
        )
    )


    return summarize(
        "Huuki",
        labels,
        harmonized,
    )


# =============================================================================
# Baituk
# =============================================================================

def choose_baituk_cell_key(
    obs,
    annotation_ids,
):

    annotation_ids = set(
        annotation_ids.astype(str)
    )

    candidates = {}


    # -------------------------------------------------------------------------
    # Candidate 1: AnnData obs_names
    # -------------------------------------------------------------------------

    candidates["obs_names"] = pd.Series(
        obs.index.astype(str),
        index=obs.index,
    )


    # -------------------------------------------------------------------------
    # Candidate 2: original_cell_id
    # -------------------------------------------------------------------------

    if "original_cell_id" in obs.columns:

        candidates[
            "original_cell_id"
        ] = (
            obs[
                "original_cell_id"
            ]
            .astype(str)
        )


    # -------------------------------------------------------------------------
    # Candidate 3: sample + original_cell_id
    # -------------------------------------------------------------------------

    if (
        "sample" in obs.columns
        and
        "original_cell_id" in obs.columns
    ):

        candidates[
            "sample_original_cell_id"
        ] = (
            obs["sample"].astype(str)
            + "_"
            + obs[
                "original_cell_id"
            ].astype(str)
        )


    results = []

    for name, values in candidates.items():

        n_match = int(
            values.isin(
                annotation_ids
            ).sum()
        )

        fraction = (
            n_match
            / len(values)
        )

        results.append(
            (
                name,
                n_match,
                fraction,
            )
        )


    results = sorted(
        results,
        key=lambda x: x[2],
        reverse=True,
    )


    print()
    print(
        "Baituk cell-ID matching:"
    )

    for name, n_match, fraction in results:

        print(
            f"  {name:<28} "
            f"{n_match:>8,}/{len(obs):,} "
            f"({100*fraction:.2f}%)"
        )


    best_name = results[0][0]
    best_fraction = results[0][2]


    if best_fraction < 0.99:

        raise ValueError(
            "Could not confidently align Baituk "
            "H5AD cells with annotations_final.RDS. "
            f"Best matching fraction = "
            f"{100*best_fraction:.2f}%."
        )


    return (
        best_name,
        candidates[
            best_name
        ],
    )


def process_baituk():

    print()
    print("=" * 100)
    print("BAITUK")
    print("=" * 100)

    baituk_file = resolve_file(
        BAITUK_CANDIDATES,
        "Baituk",
    )

    if not BAITUK_MED_LABEL_FILE.is_file():

        raise FileNotFoundError(
            "Baituk med cell-label mapping "
            "does not exist:\n"
            f"{BAITUK_MED_LABEL_FILE}\n\n"
            "Run the Baituk annotation "
            "inspection/export first."
        )


    adata = ad.read_h5ad(
        baituk_file,
        backed="r",
    )

    print(
        "Shape:",
        adata.shape,
    )

    obs = (
        adata.obs
        .copy()
    )

    try:
        adata.file.close()
    except Exception:
        pass


    annotations = pd.read_csv(
        BAITUK_MED_LABEL_FILE,
        dtype=str,
    )


    required = {
        "cell_id",
        "label",
    }

    if not required.issubset(
        annotations.columns
    ):

        raise ValueError(
            "Baituk mapping file must contain "
            "columns: cell_id, label.\n"
            f"Found: "
            f"{list(annotations.columns)}"
        )


    if annotations[
        "cell_id"
    ].duplicated().any():

        raise ValueError(
            "Duplicate Baituk cell IDs found "
            "in med annotation mapping."
        )


    _, cell_ids = (
        choose_baituk_cell_key(
            obs,
            annotations["cell_id"],
        )
    )


    lookup = (
        annotations
        .set_index(
            "cell_id"
        )[
            "label"
        ]
    )


    labels = (
        cell_ids
        .astype(str)
        .map(
            lookup
        )
    )


    n_missing = int(
        labels.isna().sum()
    )

    print(
        "Missing Baituk annotations "
        "after alignment:",
        n_missing,
    )


    if n_missing > 0:

        raise ValueError(
            f"{n_missing} Baituk cells "
            "could not be annotated."
        )


    check_mapping(
        "Baituk",
        labels,
        BAITUK_MAP,
    )


    harmonized = (
        labels.astype(str)
        .map(
            BAITUK_MAP
        )
    )


    return summarize(
        "Baituk",
        labels,
        harmonized,
    )


# =============================================================================
# Xenium
# =============================================================================

def process_xenium():

    print()
    print("=" * 100)
    print("XENIUM")
    print("=" * 100)

    xenium_file = resolve_file(
        XENIUM_CANDIDATES,
        "Xenium",
    )


    adata = ad.read_h5ad(
        xenium_file,
        backed="r",
    )

    print(
        "Shape:",
        adata.shape,
    )


    if XENIUM_COLUMN not in adata.obs.columns:

        raise KeyError(
            f"{XENIUM_COLUMN} not found.\n"
            f"Available columns:\n"
            f"{list(adata.obs.columns)}"
        )


    labels = (
        adata.obs[
            XENIUM_COLUMN
        ]
        .copy()
    )


    try:
        adata.file.close()
    except Exception:
        pass


    check_mapping(
        "Xenium",
        labels,
        XENIUM_MAP,
    )


    harmonized = (
        labels.astype(str)
        .map(
            XENIUM_MAP
        )
    )


    return summarize(
        "Xenium",
        labels,
        harmonized,
    )


# =============================================================================
# Plot
# =============================================================================

def make_plot(
    harmonized_table,
):

    plot_df = (
        harmonized_table
        .pivot(
            index="dataset",
            columns="harmonized_cell_type",
            values="percent",
        )
        .fillna(0)
    )


    dataset_order = [
        "Huuki",
        "Baituk",
        "Xenium",
    ]

    plot_df = (
        plot_df
        .reindex(
            dataset_order
        )
    )


    for category in CATEGORY_ORDER:

        if category not in plot_df.columns:
            plot_df[
                category
            ] = 0


    plot_df = (
        plot_df[
            CATEGORY_ORDER
        ]
    )


    # -------------------------------------------------------------------------
    # Save wide comparison table
    # -------------------------------------------------------------------------

    plot_df.to_csv(
        OUT_DIR
        / "harmonized_celltype_percent_wide.csv"
    )


    # -------------------------------------------------------------------------
    # 100% stacked bar plot
    # -------------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )


    bottom = np.zeros(
        len(plot_df)
    )


    x = np.arange(
        len(plot_df)
    )


    for category in CATEGORY_ORDER:

        values = (
            plot_df[
                category
            ]
            .to_numpy()
        )

        ax.bar(
            x,
            values,
            bottom=bottom,
            label=category,
            width=0.65,
        )

        bottom += values


    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        plot_df.index,
        fontsize=12,
    )

    ax.set_ylabel(
        "Cells (%)",
        fontsize=12,
    )

    ax.set_xlabel(
        "Dataset",
        fontsize=12,
    )

    ax.set_ylim(
        0,
        100,
    )

    ax.set_title(
        "Cell-type composition across Huuki, Baituk, and Xenium\n"
        "Harmonized broad annotation",
        fontsize=13,
    )

    ax.legend(
        title="Harmonized cell type",
        bbox_to_anchor=(
            1.02,
            1,
        ),
        loc="upper left",
        frameon=False,
    )

    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)

    fig.tight_layout()


    fig.savefig(
        OUT_DIR
        / "harmonized_celltype_proportions.png",
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        OUT_DIR
        / "harmonized_celltype_proportions.pdf",
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# =============================================================================
# Main
# =============================================================================

def main():

    print("=" * 100)
    print(
        "HUUKI / BAITUK / XENIUM "
        "CELL-TYPE PROPORTION COMPARISON"
    )
    print("=" * 100)


    huuki = process_huuki()
    baituk = process_baituk()
    xenium = process_xenium()


    all_results = [
        huuki,
        baituk,
        xenium,
    ]


    # -------------------------------------------------------------------------
    # Combined cell-level table is intentionally NOT saved.
    # It would be unnecessarily large.
    # -------------------------------------------------------------------------

    original_tables = [
        x[1]
        for x in all_results
    ]

    harmonized_tables = [
        x[2]
        for x in all_results
    ]


    original = pd.concat(
        original_tables,
        ignore_index=True,
    )

    harmonized = pd.concat(
        harmonized_tables,
        ignore_index=True,
    )


    # -------------------------------------------------------------------------
    # Save original annotation proportions
    # -------------------------------------------------------------------------

    original.to_csv(
        OUT_DIR
        / "original_annotation_proportions.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Save harmonized proportions
    # -------------------------------------------------------------------------

    harmonized[
        "harmonized_cell_type"
    ] = pd.Categorical(
        harmonized[
            "harmonized_cell_type"
        ],
        categories=CATEGORY_ORDER,
        ordered=True,
    )


    harmonized = (
        harmonized
        .sort_values(
            [
                "dataset",
                "harmonized_cell_type",
            ]
        )
        .reset_index(
            drop=True
        )
    )


    harmonized.to_csv(
        OUT_DIR
        / "harmonized_celltype_proportions.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Total-cell summary
    # -------------------------------------------------------------------------

    totals = (
        harmonized
        .groupby(
            "dataset",
            observed=True,
        )[
            "n_cells"
        ]
        .sum()
        .reset_index(
            name="total_cells"
        )
    )


    totals.to_csv(
        OUT_DIR
        / "dataset_total_cells.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Plot
    # -------------------------------------------------------------------------

    make_plot(
        harmonized
    )


    # -------------------------------------------------------------------------
    # Print advisor-facing table
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        "HARMONIZED CELL-TYPE PROPORTIONS"
    )
    print("=" * 100)

    display = (
        harmonized[
            [
                "dataset",
                "harmonized_cell_type",
                "n_cells",
                "percent",
            ]
        ]
        .copy()
    )

    display[
        "percent"
    ] = (
        display[
            "percent"
        ]
        .round(2)
    )

    print(
        display.to_string(
            index=False
        )
    )


    print()
    print("=" * 100)
    print("OUTPUTS")
    print("=" * 100)

    for filename in [
        "dataset_total_cells.csv",
        "original_annotation_proportions.csv",
        "harmonized_celltype_proportions.csv",
        "harmonized_celltype_percent_wide.csv",
        "harmonized_celltype_proportions.png",
        "harmonized_celltype_proportions.pdf",
    ]:

        print(
            OUT_DIR
            / filename
        )


if __name__ == "__main__":
    main()
