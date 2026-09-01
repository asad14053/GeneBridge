#!/usr/bin/env python3

from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

META = ROOT / "data/metadata/xenium_DE_metadata_23.csv"

ANNOTATED = (
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)

EX1_DIR = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/targets"
)

EX2_DIR = (
    ROOT
    / "data/processed/imputation_full/ex2_scz/targets"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/"
      "paper_concordance/original_layer_pseudobulk"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


EXPECTED_LAYERS = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WMtz",
    "WM",
]

LAYER_TO_SPD = {
    "L1/M": "spd07",
    "L2/3": "spd06",
    "L3/4": "spd02",
    "L5": "spd05",
    "L6": "spd03",
    "WMtz": "spd01",
    "WM": "spd04",
}

MIN_NCELLS = 10


def section(title):
    print(
        "\n"
        + "=" * 100
        + f"\n{title}\n"
        + "=" * 100,
        flush=True,
    )


def sum_matrix(matrix):

    if sparse.issparse(matrix):

        x = np.asarray(
            matrix.sum(axis=0)
        ).ravel()

    else:

        x = np.asarray(
            matrix
        ).sum(
            axis=0,
            dtype=np.float64,
        )

    return np.asarray(
        x,
        dtype=np.float64,
    )


# =============================================================================
# Metadata
# =============================================================================

section("LOAD CANONICAL 23-DONOR METADATA")

meta = pd.read_csv(META)

required = {
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "run_date",
}

missing = required - set(meta.columns)

if missing:
    raise RuntimeError(
        f"Metadata missing columns: {sorted(missing)}"
    )

meta["BrNum"] = (
    meta["BrNum"]
    .astype(str)
    .str.strip()
)

if len(meta) != 23:
    raise RuntimeError(
        f"Expected 23 donors; found {len(meta)}"
    )

if meta["BrNum"].nunique() != 23:
    raise RuntimeError(
        "BrNum is not unique in canonical metadata."
    )

if "Br6432" in set(meta["BrNum"]):
    raise RuntimeError(
        "Br6432 must be excluded."
    )


print(
    meta[
        [
            "BrNum",
            "Dx",
            "Age",
            "Sex",
            "slide_id",
            "run_date",
        ]
    ].to_string(index=False)
)

print("\nDiagnosis counts:")
print(meta["Dx"].value_counts())


# =============================================================================
# Cell-level annotation source
# =============================================================================

section("LOAD CELL-LEVEL SPATIAL ANNOTATIONS")

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
)

required_obs = {
    "BrNum",
    "layer_annotation",
    "predictions_smooth",
}

missing_obs = required_obs - set(ann.obs.columns)

if missing_obs:
    raise RuntimeError(
        f"Annotated object missing: {sorted(missing_obs)}"
    )


annotation = ann.obs[
    [
        "BrNum",
        "layer_annotation",
        "predictions_smooth",
    ]
].copy()


annotation["BrNum"] = (
    annotation["BrNum"]
    .astype(str)
    .str.strip()
)

annotation["layer_annotation"] = (
    annotation["layer_annotation"]
    .astype(str)
    .str.strip()
)

annotation["predictions_smooth"] = (
    annotation["predictions_smooth"]
    .astype(str)
    .str.strip()
)


if not annotation.index.is_unique:
    raise RuntimeError(
        "Annotated N24 cell IDs are not globally unique."
    )


print("Annotated cells:", len(annotation))
print(
    "Annotated donors:",
    annotation["BrNum"].nunique(),
)

ann.file.close()


# =============================================================================
# Locate target files
# =============================================================================

section("LOCATE ORIGINAL XENIUM TARGET FILES")

records = []


for row in meta.itertuples(index=False):

    donor = str(row.BrNum)
    dx = str(row.Dx)

    if dx == "NTC":

        path = (
            EX1_DIR
            / f"spatial_data_xenium_{donor}_ex1_ntc.h5ad"
        )

    elif dx == "SCZ":

        path = (
            EX2_DIR
            / f"spatial_data_xenium_{donor}_ex2_scz.h5ad"
        )

    else:

        raise RuntimeError(
            f"{donor}: unexpected Dx={dx}"
        )


    if not path.is_file():
        raise FileNotFoundError(path)


    records.append(
        {
            "BrNum": donor,
            "Dx": dx,
            "Age": float(row.Age),
            "Sex": str(row.Sex),
            "slide_id": str(row.slide_id),
            "run_date": str(row.run_date),
            "file": path,
        }
    )


print(
    f"Found {len(records)} / 23 target files."
)


# =============================================================================
# Build donor × layer pseudobulk
# =============================================================================

section("BUILD ORIGINAL XENIUM DONOR × LAYER PSEUDOBULK")

gene_order = None

pb_counts = []
pb_rows = []

donor_validation_rows = []
excluded_cell_rows = []


for i, record in enumerate(
    records,
    start=1,
):

    donor = record["BrNum"]

    print(
        f"\n[{i:02d}/23] "
        f"{donor} | {record['Dx']}",
        flush=True,
    )


    # Only 300 genes, so donor-level objects are manageable in memory.
    a = ad.read_h5ad(
        record["file"]
    )


    # -------------------------------------------------------------------------
    # Gene validation
    # -------------------------------------------------------------------------

    genes = (
        a.var_names
        .astype(str)
        .tolist()
    )


    if gene_order is None:

        gene_order = genes

    elif genes != gene_order:

        raise RuntimeError(
            f"{donor}: gene order differs."
        )


    if len(genes) != 300:

        raise RuntimeError(
            f"{donor}: expected 300 genes, "
            f"found {len(genes)}"
        )


    # -------------------------------------------------------------------------
    # Use same count source as existing donor-level script
    # -------------------------------------------------------------------------

    if "counts" in a.layers:

        matrix = a.layers["counts"]
        matrix_source = "layers[counts]"

    else:

        matrix = a.X
        matrix_source = "X"


    print(
        "  expression source:",
        matrix_source,
    )


    # -------------------------------------------------------------------------
    # Match cells to annotated N24 object
    # -------------------------------------------------------------------------

    target_ids = pd.Index(
        a.obs_names.astype(str)
    )


    ann_donor = annotation.loc[
        annotation["BrNum"] == donor
    ]


    ann_ids = pd.Index(
        ann_donor.index.astype(str)
    )


    target_only = target_ids.difference(
        ann_ids
    )

    annotated_only = ann_ids.difference(
        target_ids
    )


    if len(annotated_only) > 0:

        raise RuntimeError(
            f"{donor}: {len(annotated_only)} "
            "annotated cells missing from target."
        )


    if len(target_only) > 0:

        print(
            f"  excluding {len(target_only)} "
            "cells without layer annotation"
        )

        for cell_id in target_only:

            excluded_cell_rows.append(
                {
                    "BrNum": donor,
                    "cell_id": str(cell_id),
                    "reason":
                        "missing_layer_annotation",
                }
            )


    matched_mask = target_ids.isin(
        ann_ids
    )

    matched_positions = np.flatnonzero(
        matched_mask
    )

    matched_ids = target_ids[
        matched_mask
    ]


    # Cell annotation in exact expression-matrix order.
    cell_ann = annotation.loc[
        matched_ids,
        [
            "BrNum",
            "layer_annotation",
            "predictions_smooth",
        ],
    ].copy()


    if not (
        cell_ann["BrNum"] == donor
    ).all():

        raise RuntimeError(
            f"{donor}: BrNum mismatch after cell join."
        )


    n_target = a.n_obs
    n_matched = len(matched_ids)
    n_excluded = len(target_only)


    print(
        f"  target cells: {n_target}"
    )

    print(
        f"  layer-annotated cells: {n_matched}"
    )

    print(
        f"  excluded cells: {n_excluded}"
    )


    # -------------------------------------------------------------------------
    # Validate full and matched donor totals
    # -------------------------------------------------------------------------

    full_donor_sum = sum_matrix(
        matrix
    )


    if sparse.issparse(matrix):

        matched_matrix = matrix[
            matched_positions,
            :
        ]

    else:

        matched_matrix = np.asarray(
            matrix[
                matched_positions,
                :
            ]
        )


    matched_donor_sum = sum_matrix(
        matched_matrix
    )


    excluded_library = float(
        full_donor_sum.sum()
        - matched_donor_sum.sum()
    )


    # -------------------------------------------------------------------------
    # Seven spatial-layer pseudobulks
    # -------------------------------------------------------------------------

    layer_sum_check = np.zeros(
        len(gene_order),
        dtype=np.float64,
    )


    for layer in EXPECTED_LAYERS:

        spd = LAYER_TO_SPD[layer]

        layer_mask = (
            cell_ann["layer_annotation"]
            .eq(layer)
            .to_numpy()
        )

        # Check equivalent SpD label.
        if not (
            cell_ann.loc[
                layer_mask,
                "predictions_smooth",
            ]
            == spd
        ).all():

            raise RuntimeError(
                f"{donor} {layer}: "
                "SpD/layer mapping mismatch."
            )


        layer_positions_in_matched = np.flatnonzero(
            layer_mask
        )


        n_cells = len(
            layer_positions_in_matched
        )


        if n_cells < MIN_NCELLS:

            raise RuntimeError(
                f"{donor} {layer}: "
                f"{n_cells} cells < min_ncells={MIN_NCELLS}"
            )


        if sparse.issparse(
            matched_matrix
        ):

            layer_matrix = matched_matrix[
                layer_positions_in_matched,
                :
            ]

        else:

            layer_matrix = matched_matrix[
                layer_positions_in_matched,
                :
            ]


        summed = sum_matrix(
            layer_matrix
        )


        if (
            not np.isfinite(summed).all()
            or np.any(summed < 0)
        ):

            raise RuntimeError(
                f"{donor} {layer}: "
                "invalid pseudobulk counts."
            )


        pb_id = (
            f"{donor}__{spd}__"
            f"{layer.replace('/', '-')}"
        )


        pb_counts.append(
            summed
        )


        pb_rows.append(
            {
                "pseudobulk_id": pb_id,
                "BrNum": donor,
                "Dx": record["Dx"],
                "Age": record["Age"],
                "Sex": record["Sex"],
                "slide_id": record["slide_id"],
                "run_date": record["run_date"],
                "predictions_smooth": spd,
                "layer_annotation": layer,
                "n_cells": n_cells,
                "library_size":
                    float(summed.sum()),
                "matrix_source":
                    matrix_source,
            }
        )


        layer_sum_check += summed


        print(
            f"    {spd:5s} | "
            f"{layer:4s} | "
            f"cells={n_cells:6d} | "
            f"library={summed.sum():.0f}"
        )


    # -------------------------------------------------------------------------
    # Critical validation:
    # Sum of 7 layers must equal matched donor pseudobulk.
    # -------------------------------------------------------------------------

    if not np.allclose(
        layer_sum_check,
        matched_donor_sum,
        rtol=0,
        atol=1e-8,
    ):

        diff = np.max(
            np.abs(
                layer_sum_check
                - matched_donor_sum
            )
        )

        raise RuntimeError(
            f"{donor}: sum of 7 layers does not equal "
            f"matched donor counts. max_diff={diff}"
        )


    donor_validation_rows.append(
        {
            "BrNum": donor,
            "Dx": record["Dx"],
            "target_cells": n_target,
            "matched_layer_cells": n_matched,
            "excluded_unannotated_cells":
                n_excluded,
            "full_target_library":
                float(full_donor_sum.sum()),
            "layer_matched_library":
                float(matched_donor_sum.sum()),
            "excluded_library":
                excluded_library,
            "layer_sum_matches_donor":
                True,
        }
    )


    del (
        a,
        matrix,
        matched_matrix,
    )


# =============================================================================
# Assemble outputs
# =============================================================================

section("ASSEMBLE PSEUDOBULK MATRICES")


counts = np.vstack(
    pb_counts
)


pb_meta = pd.DataFrame(
    pb_rows
)


donor_validation = pd.DataFrame(
    donor_validation_rows
)


excluded_cells = pd.DataFrame(
    excluded_cell_rows,
    columns=[
        "BrNum",
        "cell_id",
        "reason",
    ],
)


print(
    "Pseudobulk count matrix:",
    counts.shape,
)

print(
    "Pseudobulk metadata:",
    pb_meta.shape,
)


if counts.shape != (161, 300):

    raise RuntimeError(
        f"Expected 161 × 300 matrix; "
        f"found {counts.shape}"
    )


if len(pb_meta) != 161:

    raise RuntimeError(
        f"Expected 161 metadata rows; "
        f"found {len(pb_meta)}"
    )


if pb_meta["BrNum"].nunique() != 23:

    raise RuntimeError(
        "Expected 23 donors."
    )


if pb_meta[
    "layer_annotation"
].nunique() != 7:

    raise RuntimeError(
        "Expected 7 spatial layers."
    )


if (
    pb_meta["n_cells"]
    < MIN_NCELLS
).any():

    raise RuntimeError(
        "A donor × layer pseudobulk "
        "contains <10 cells."
    )


# =============================================================================
# Save
# =============================================================================

section("SAVE OUTPUTS")


counts_df = pd.DataFrame(
    counts,
    index=pb_meta["pseudobulk_id"],
    columns=gene_order,
)


counts_file = (
    OUT
    / "original_xenium_300gene_donor_layer_pseudobulk_counts.csv.gz"
)

metadata_file = (
    OUT
    / "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)

validation_file = (
    OUT
    / "original_xenium_300gene_donor_layer_validation.csv"
)

excluded_file = (
    OUT
    / "original_xenium_cells_excluded_missing_layer.csv"
)

genes_file = (
    OUT
    / "original_xenium_300gene_order.txt"
)


counts_df.to_csv(
    counts_file
)

pb_meta.to_csv(
    metadata_file,
    index=False,
)

donor_validation.to_csv(
    validation_file,
    index=False,
)

excluded_cells.to_csv(
    excluded_file,
    index=False,
)

genes_file.write_text(
    "\n".join(gene_order)
    + "\n"
)


# =============================================================================
# Final summary
# =============================================================================

section("FINAL SUMMARY")


print(
    "Donors:",
    pb_meta["BrNum"].nunique(),
)

print(
    "Layers:",
    pb_meta["layer_annotation"].nunique(),
)

print(
    "Pseudobulk samples:",
    len(pb_meta),
)

print(
    "Genes:",
    len(gene_order),
)

print(
    "Minimum cells in any donor × layer:",
    pb_meta["n_cells"].min(),
)

print(
    "Maximum cells in any donor × layer:",
    pb_meta["n_cells"].max(),
)

print(
    "Excluded target cells:",
    len(excluded_cells),
)


print("\nPseudobulks per layer:")

print(
    pb_meta[
        "layer_annotation"
    ]
    .value_counts()
    .reindex(EXPECTED_LAYERS)
    .to_string()
)


print("\nPseudobulks per diagnosis:")

print(
    pb_meta[
        [
            "BrNum",
            "Dx",
        ]
    ]
    .drop_duplicates()
    ["Dx"]
    .value_counts()
    .to_string()
)


print(
    "\nAll layer sums equal "
    "their matched donor-level counts:",
    donor_validation[
        "layer_sum_matches_donor"
    ].all(),
)


print("\nExcluded cells:")

if excluded_cells.empty:

    print("NONE")

else:

    print(
        excluded_cells.to_string(
            index=False
        )
    )


print("\nOutputs:")
print(counts_file)
print(metadata_file)
print(validation_file)
print(excluded_file)
print(genes_file)


print(
    "\nFINAL STATUS: PASS — "
    "original Xenium 23-donor × 7-layer "
    "300-gene pseudobulk successfully created."
)
