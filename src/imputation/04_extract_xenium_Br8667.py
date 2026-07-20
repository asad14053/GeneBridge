#!/usr/bin/env python

"""
04_extract_xenium_Br8667.py

Extract Br8667 from the annotated Xenium N24 AnnData object.

Input:
    data/processed/xenium/
        xenium_N24_layer_celltype_annotated.h5ad

Outputs:
    data/processed/imputation_beta/Br8667/
        xenium_Br8667_annotated.h5ad
        xenium_Br8667_extraction_summary.csv
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

INPUT_H5AD = (
    PROJECT_ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
)

OUTPUT_H5AD = OUTPUT_DIR / "xenium_Br8667_annotated.h5ad"

SUMMARY_CSV = (
    OUTPUT_DIR
    / "xenium_Br8667_extraction_summary.csv"
)

TARGET_BRAIN = "Br8667"
BRAIN_COLUMN = "BrNum"


def prepare_output_directory(path: Path) -> None:
    """
    Create output directory and detect symbolic-link loops.
    """

    try:
        resolved = path.resolve(strict=False)
    except RuntimeError as exc:
        raise RuntimeError(
            f"Symbolic-link loop detected for output directory:\n{path}"
        ) from exc

    path.mkdir(parents=True, exist_ok=True)

    print("\nLogical output directory:")
    print(path)

    print("\nPhysical output directory:")
    print(resolved)


def print_created_file(label: str, path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Expected output was not created: {path}"
        )

    size_mb = path.stat().st_size / (1024**2)

    print(f"\nCreated {label}:")
    print(path)
    print(f"Size: {size_mb:.2f} MB")


def ensure_spatial_coordinates(adata: ad.AnnData) -> None:
    """
    Preserve existing coordinates or construct obsm['spatial']
    from centroid columns.
    """

    if "spatial" in adata.obsm:
        coordinates = np.asarray(adata.obsm["spatial"])

        if coordinates.ndim != 2 or coordinates.shape[1] < 2:
            raise ValueError(
                "Existing obsm['spatial'] is not an n_cells x 2 matrix."
            )

        adata.obsm["spatial"] = coordinates[:, :2].astype(
            np.float32,
            copy=False,
        )

        print("\nExisting obsm['spatial'] preserved.")
        print("Coordinate shape:", adata.obsm["spatial"].shape)
        return

    candidate_pairs = [
        ("x_centroid", "y_centroid"),
        ("cell_centroid_x", "cell_centroid_y"),
        ("x_location", "y_location"),
        ("x", "y"),
    ]

    for x_column, y_column in candidate_pairs:
        if (
            x_column in adata.obs.columns
            and y_column in adata.obs.columns
        ):
            adata.obsm["spatial"] = (
                adata.obs[[x_column, y_column]]
                .to_numpy(dtype=np.float32)
            )

            print(
                "\nCreated obsm['spatial'] using "
                f"obs['{x_column}'] and obs['{y_column}']."
            )
            print("Coordinate shape:", adata.obsm["spatial"].shape)
            return

    raise ValueError(
        "Could not find Xenium spatial coordinates. "
        "No obsm['spatial'] or centroid columns were found."
    )


def main() -> None:
    print("=" * 100)
    print("Task-01: Extract Xenium Br8667 from Xenium N24")
    print("=" * 100)

    if not INPUT_H5AD.exists():
        raise FileNotFoundError(
            f"Missing Xenium N24 input:\n{INPUT_H5AD}"
        )

    prepare_output_directory(OUTPUT_DIR)

    print("\nInput file:")
    print(INPUT_H5AD)

    print("\nLoading Xenium N24 in backed mode...")

    # Backed mode prevents the complete 1.2-million-cell
    # object from being loaded into RAM.
    xenium_n24 = ad.read_h5ad(
        INPUT_H5AD,
        backed="r",
    )

    print("\nXenium N24:")
    print(xenium_n24)

    if BRAIN_COLUMN not in xenium_n24.obs.columns:
        raise KeyError(
            f"Missing obs['{BRAIN_COLUMN}'].\n"
            f"Available obs columns:\n"
            f"{xenium_n24.obs.columns.tolist()}"
        )

    brain_ids = (
        xenium_n24.obs[BRAIN_COLUMN]
        .astype(str)
        .str.strip()
    )

    keep = brain_ids.eq(TARGET_BRAIN).to_numpy()

    selected_cells = int(keep.sum())

    print("\nBrain column:")
    print(BRAIN_COLUMN)

    print("\nTarget donor:")
    print(TARGET_BRAIN)

    print("\nSelected Xenium cells:")
    print(f"{selected_cells:,}")

    if selected_cells == 0:
        raise ValueError(
            f"No Xenium cells were found for {TARGET_BRAIN}."
        )

    print("\nLoading only Br8667 cells into memory...")

    xenium = xenium_n24[keep, :].to_memory()

    try:
        xenium_n24.file.close()
    except Exception:
        pass

    print("\nExtracted Xenium Br8667:")
    print(xenium)

    print("\nAvailable layers:")
    print(list(xenium.layers.keys()))

    if "counts" not in xenium.layers:
        raise KeyError(
            "Xenium N24 does not contain layers['counts']. "
            "Do not continue to VISTA until raw counts are located."
        )

    ensure_spatial_coordinates(xenium)

    # VISTA-compatible helper metadata.
    xenium.obs["batch"] = TARGET_BRAIN
    xenium.obs["names"] = xenium.obs_names.astype(str)

    if "cell_type_annotation" in xenium.obs.columns:
        xenium.obs["scClassify"] = (
            xenium.obs["cell_type_annotation"]
            .astype("string")
            .fillna("Unknown")
            .astype(str)
        )

        celltype_source = "cell_type_annotation"
    else:
        xenium.obs["scClassify"] = "Unknown"
        celltype_source = "Unknown"

    xenium.obs_names_make_unique()
    xenium.var_names_make_unique()

    # Avoid AnnData writing conflicts.
    xenium.obs.index.name = None
    xenium.var.index.name = None

    summary = pd.DataFrame(
        [
            {
                "metric": "input_file",
                "value": str(INPUT_H5AD),
            },
            {
                "metric": "output_file",
                "value": str(OUTPUT_H5AD),
            },
            {
                "metric": "brain_id",
                "value": TARGET_BRAIN,
            },
            {
                "metric": "brain_column",
                "value": BRAIN_COLUMN,
            },
            {
                "metric": "cells",
                "value": xenium.n_obs,
            },
            {
                "metric": "genes",
                "value": xenium.n_vars,
            },
            {
                "metric": "layers",
                "value": ",".join(xenium.layers.keys()),
            },
            {
                "metric": "has_counts_layer",
                "value": "counts" in xenium.layers,
            },
            {
                "metric": "has_log1p_norm_layer",
                "value": "log1p_norm" in xenium.layers,
            },
            {
                "metric": "has_spatial",
                "value": "spatial" in xenium.obsm,
            },
            {
                "metric": "spatial_shape",
                "value": str(xenium.obsm["spatial"].shape),
            },
            {
                "metric": "celltype_source",
                "value": celltype_source,
            },
            {
                "metric": "has_layer_annotation",
                "value": "layer_annotation" in xenium.obs.columns,
            },
        ]
    )

    summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    if OUTPUT_H5AD.exists():
        print("\nRemoving old output:")
        print(OUTPUT_H5AD)
        OUTPUT_H5AD.unlink()

    print("\nWriting Xenium Br8667 h5ad:")
    print(OUTPUT_H5AD)

    xenium.write_h5ad(
        OUTPUT_H5AD,
        compression="gzip",
    )

    print_created_file(
        "Xenium Br8667 h5ad",
        OUTPUT_H5AD,
    )

    print_created_file(
        "extraction summary",
        SUMMARY_CSV,
    )

    print("\nFinal Xenium Br8667 object:")
    print(xenium)

    print("\nFinal obs helper columns:")
    print(
        [
            column
            for column in [
                "BrNum",
                "batch",
                "names",
                "scClassify",
                "cell_type_annotation",
                "layer_annotation",
            ]
            if column in xenium.obs.columns
        ]
    )

    print("\nDONE: Task-01 completed.")


if __name__ == "__main__":
    main()
