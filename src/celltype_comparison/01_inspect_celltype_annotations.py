#!/usr/bin/env python3

from pathlib import Path
import json

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

OUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "celltype_comparison"
    / "01_annotation_inventory"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


DATASETS = {

    "Huuki": [
        PROJECT_ROOT
        / "data/processed/snrnaseq/sce_DLPFC_annotated/"
        "huuki_snrna_reference_full_allgenes.h5ad",

        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667/"
        "seq_data_huuki_snrna_9brains_excluding_Br8667.h5ad",

        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667/"
        "seq_data_huuki_snrna_Br8667_vista.h5ad",
    ],

    "Xenium": [
        PROJECT_ROOT
        / "data/processed/xenium/"
        "xenium_N24_layer_celltype_annotated.h5ad",

        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667/"
        "xenium_Br8667_annotated.h5ad",

        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667/"
        "spatial_data_xenium_Br8667_vista_qc.h5ad",
    ],
}


ANNOTATION_KEYWORDS = [
    "celltype",
    "cell_type",
    "cell type",
    "cellclass",
    "cell_class",
    "class",
    "subclass",
    "annotation",
    "annot",
    "label",
    "cluster",
    "type",
    "broad",
    "major",
    "subtype",
]


LAYER_KEYWORDS = [
    "layer",
]


def find_existing_file(paths):

    for path in paths:
        if path.is_file():
            return path

    return None


def looks_like_celltype(column):

    name = str(column).lower()

    return any(
        keyword in name
        for keyword in ANNOTATION_KEYWORDS
    )


def looks_like_layer(column):

    name = str(column).lower()

    return any(
        keyword in name
        for keyword in LAYER_KEYWORDS
    )


def count_labels(
    dataset,
    column,
    series,
):

    values = (
        series
        .astype("string")
        .fillna("<NA>")
    )

    counts = (
        values
        .value_counts(
            dropna=False
        )
        .rename_axis("label")
        .reset_index(
            name="n_cells"
        )
    )

    counts["percent"] = (
        counts["n_cells"]
        / counts["n_cells"].sum()
        * 100
    )

    counts.insert(
        0,
        "annotation_column",
        str(column),
    )

    counts.insert(
        0,
        "dataset",
        dataset,
    )

    return counts


def inspect_dataset(
    dataset,
    path,
):

    print()
    print("=" * 110)
    print(dataset.upper())
    print("=" * 110)

    print("File:")
    print(path)

    adata = ad.read_h5ad(
        path,
        backed="r",
    )

    n_cells = int(
        adata.n_obs
    )

    n_genes = int(
        adata.n_vars
    )

    obs = (
        adata.obs
        .copy()
    )

    try:
        adata.file.close()
    except Exception:
        pass

    print()
    print(
        "Shape:",
        f"{n_cells:,} cells x {n_genes:,} genes",
    )

    dataset_dir = (
        OUT_DIR
        / dataset.lower()
    )

    dataset_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    # -------------------------------------------------------------------------
    # All obs columns
    # -------------------------------------------------------------------------

    print()
    print("OBS COLUMNS")
    print("-" * 110)

    inventory_rows = []

    for index, column in enumerate(
        obs.columns,
        start=1,
    ):

        series = obs[column]

        n_unique = int(
            series.nunique(
                dropna=True
            )
        )

        n_missing = int(
            series.isna().sum()
        )

        is_celltype = (
            looks_like_celltype(
                column
            )
        )

        is_layer = (
            looks_like_layer(
                column
            )
        )

        print(
            f"{index:3d}. "
            f"{column:<40} "
            f"unique={n_unique:<8} "
            f"missing={n_missing}"
        )

        inventory_rows.append(
            {
                "dataset":
                    dataset,

                "column":
                    str(column),

                "dtype":
                    str(series.dtype),

                "n_unique_nonmissing":
                    n_unique,

                "n_missing":
                    n_missing,

                "candidate_celltype":
                    is_celltype,

                "candidate_layer":
                    is_layer,
            }
        )


    inventory = pd.DataFrame(
        inventory_rows
    )

    inventory.to_csv(
        dataset_dir
        / "obs_column_inventory.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Candidate annotation columns
    # -------------------------------------------------------------------------

    candidate_columns = [
        column
        for column in obs.columns
        if looks_like_celltype(
            column
        )
    ]

    layer_columns = [
        column
        for column in obs.columns
        if looks_like_layer(
            column
        )
    ]


    print()
    print("CANDIDATE CELL-TYPE / ANNOTATION COLUMNS")
    print("-" * 110)

    if len(candidate_columns) == 0:
        print("NONE FOUND BY COLUMN NAME")
    else:
        for column in candidate_columns:
            print(column)


    print()
    print("CANDIDATE LAYER COLUMNS")
    print("-" * 110)

    if len(layer_columns) == 0:
        print("NONE FOUND BY COLUMN NAME")
    else:
        for column in layer_columns:
            print(column)


    all_counts = []

    for column in candidate_columns:

        counts = count_labels(
            dataset,
            column,
            obs[column],
        )

        all_counts.append(
            counts
        )

        safe_column = (
            str(column)
            .replace("/", "_")
            .replace(" ", "_")
        )

        counts.to_csv(
            dataset_dir
            / f"counts__{safe_column}.csv",
            index=False,
        )

        print()
        print(
            f"ANNOTATION COLUMN: {column}"
        )
        print("-" * 110)

        print(
            counts[
                [
                    "label",
                    "n_cells",
                    "percent",
                ]
            ]
            .head(100)
            .to_string(
                index=False
            )
        )


    if all_counts:

        pd.concat(
            all_counts,
            ignore_index=True,
        ).to_csv(
            dataset_dir
            / "all_candidate_annotation_counts.csv",
            index=False,
        )


    # -------------------------------------------------------------------------
    # Layers
    # -------------------------------------------------------------------------

    all_layer_counts = []

    for column in layer_columns:

        counts = count_labels(
            dataset,
            column,
            obs[column],
        )

        all_layer_counts.append(
            counts
        )


    if all_layer_counts:

        pd.concat(
            all_layer_counts,
            ignore_index=True,
        ).to_csv(
            dataset_dir
            / "all_candidate_layer_counts.csv",
            index=False,
        )


    return {
        "dataset":
            dataset,

        "file":
            str(path),

        "n_cells":
            n_cells,

        "n_genes":
            n_genes,

        "candidate_annotation_columns":
            ";".join(
                map(
                    str,
                    candidate_columns,
                )
            ),

        "candidate_layer_columns":
            ";".join(
                map(
                    str,
                    layer_columns,
                )
            ),
    }


def main():

    print("=" * 110)
    print("HUUKI + XENIUM CELL-TYPE ANNOTATION INVENTORY")
    print("=" * 110)

    summaries = []

    for dataset, candidates in DATASETS.items():

        path = find_existing_file(
            candidates
        )

        if path is None:

            print()
            print(
                f"ERROR: no candidate {dataset} H5AD found."
            )

            print(
                "Checked:"
            )

            for candidate in candidates:
                print(
                    "  ",
                    candidate,
                )

            continue

        result = inspect_dataset(
            dataset,
            path,
        )

        summaries.append(
            result
        )


    summary_df = pd.DataFrame(
        summaries
    )

    summary_path = (
        OUT_DIR
        / "huuki_xenium_inventory.csv"
    )

    summary_df.to_csv(
        summary_path,
        index=False,
    )

    print()
    print("=" * 110)
    print("DONE")
    print("=" * 110)

    print(
        summary_df.to_string(
            index=False
        )
    )

    print()
    print(
        "Summary:",
        summary_path,
    )


if __name__ == "__main__":
    main()
