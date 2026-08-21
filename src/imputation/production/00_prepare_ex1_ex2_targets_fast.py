#!/usr/bin/env python3

"""
Prepare production Xenium targets for Ex1 / Ex2.

Ex1:
    all-10 Huuki snRNA reference -> NTC Xenium targets

Ex2:
    all-10 Huuki snRNA reference -> SCZ Xenium targets

Design
------
Expression / QC source:
    data/processed/xenium/xenium_N24_imputation_ready.h5ad

Annotation source:
    data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad

Join key:
    donor_cell_id

Filtering rule:
    Keep only cells that:
        1. are present in the annotated Xenium object
        2. have non-missing layer_annotation
        3. have non-missing cell_type_annotation

Br6432 is intentionally excluded.

The expression matrix always comes from the imputation-ready object.
The annotation object contributes metadata only.
"""

from pathlib import Path
import argparse

import anndata as ad
import numpy as np
import pandas as pd


# =============================================================================
# PATHS
# =============================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

EXPRESSION_SOURCE = (
    ROOT
    / "data/processed/xenium/xenium_N24_imputation_ready.h5ad"
)

ANNOTATION_SOURCE = (
    ROOT
    / "data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad"
)

METADATA = (
    ROOT
    / "data/metadata/patient_xenium_visium_24_common_with_dx.csv"
)

EX1_DIR = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/targets"
)

EX2_DIR = (
    ROOT
    / "data/processed/imputation_full/ex2_scz/targets"
)

SUMMARY_OUT = (
    ROOT
    / "data/metadata/imputation_ex1_ex2_target_summary.csv"
)


# =============================================================================
# ANNOTATION COLUMNS
# =============================================================================

ANNOTATION_COLUMNS = [
    "predictions_smooth",

    "layer_annotation",
    "layer_annotation_source",
    "layer_annotation_stage",

    "Banksy",

    "cell_type_annotation_banksy",
    "cell_type_annotation",
    "cell_type_annotation_source",
    "cell_type_annotation_stage",
    "cell_type_annotation_confidence",
]

REQUIRED_FINAL_ANNOTATIONS = [
    "layer_annotation",
    "cell_type_annotation",
]


# =============================================================================
# HELPERS
# =============================================================================

def section(title):

    print(
        "\n" + "=" * 110
    )

    print(title)

    print(
        "=" * 110
    )


def clean_string(x):

    return (
        x.astype(str)
        .str.strip()
    )


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing donor target H5AD files.",
    )

    return parser.parse_args()


# =============================================================================
# MAIN
# =============================================================================

def main():

    args = parse_args()

    section(
        "PREPARE EX1 / EX2 TARGETS "
        "WITH CELL-TYPE + LAYER FILTERING"
    )

    print(
        "Expression source :",
        EXPRESSION_SOURCE
    )

    print(
        "Annotation source :",
        ANNOTATION_SOURCE
    )

    print(
        "Metadata          :",
        METADATA
    )

    print(
        "Overwrite         :",
        args.overwrite
    )


    # =========================================================================
    # INPUT CHECK
    # =========================================================================

    for path in [
        EXPRESSION_SOURCE,
        ANNOTATION_SOURCE,
        METADATA,
    ]:

        if not path.exists():

            raise FileNotFoundError(
                path
            )


    EX1_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    EX2_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


    # =========================================================================
    # AUTHORITATIVE DONOR METADATA
    # =========================================================================

    section(
        "LOAD DONOR METADATA"
    )

    meta = pd.read_csv(
        METADATA
    )


    if "patient_id" not in meta.columns:

        raise RuntimeError(
            "patient_id missing from metadata."
        )


    if "Dx" not in meta.columns:

        raise RuntimeError(
            "Dx missing from metadata."
        )


    meta["patient_id"] = clean_string(
        meta["patient_id"]
    )

    meta["Dx"] = (
        clean_string(
            meta["Dx"]
        )
        .str.upper()
    )


    # Corrupted donor.
    meta = meta.loc[
        meta["patient_id"] != "Br6432"
    ].copy()


    meta = meta.loc[
        meta["Dx"].isin(
            [
                "NTC",
                "SCZ",
            ]
        )
    ].copy()


    if len(meta) != 23:

        raise RuntimeError(
            f"Expected 23 usable donors; found {len(meta)}"
        )


    if (
        meta["Dx"]
        .eq("NTC")
        .sum()
        != 11
    ):

        raise RuntimeError(
            "Expected 11 NTC donors."
        )


    if (
        meta["Dx"]
        .eq("SCZ")
        .sum()
        != 12
    ):

        raise RuntimeError(
            "Expected 12 SCZ donors."
        )


    print(
        meta[
            [
                "patient_id",
                "Dx",
            ]
        ]
        .sort_values(
            [
                "Dx",
                "patient_id",
            ]
        )
        .to_string(
            index=False
        )
    )


    # =========================================================================
    # LOAD ANNOTATION TABLE
    # =========================================================================

    section(
        "LOAD CELL / LAYER ANNOTATION"
    )

    ann = ad.read_h5ad(
        ANNOTATION_SOURCE,
        backed="r"
    )


    required_columns = [
        "donor_cell_id",
        "BrNum",
        *ANNOTATION_COLUMNS,
    ]


    missing_columns = [
        c
        for c in required_columns
        if c not in ann.obs.columns
    ]


    if missing_columns:

        raise RuntimeError(
            "Missing annotation columns:\n"
            + "\n".join(
                missing_columns
            )
        )


    annotation = (
        ann.obs[
            required_columns
        ]
        .copy()
    )


    annotation[
        "donor_cell_id"
    ] = clean_string(
        annotation[
            "donor_cell_id"
        ]
    )


    annotation[
        "BrNum"
    ] = clean_string(
        annotation[
            "BrNum"
        ]
    )


    if annotation[
        "donor_cell_id"
    ].duplicated().any():

        raise RuntimeError(
            "donor_cell_id is not globally unique "
            "in annotation source."
        )


    annotation = annotation.set_index(
        "donor_cell_id",
        drop=True
    )


    print(
        "Annotated cells :",
        f"{len(annotation):,}"
    )

    print(
        "Annotated donors:",
        annotation[
            "BrNum"
        ].nunique()
    )


    ann.file.close()


    # =========================================================================
    # LOAD EXPRESSION SOURCE
    # =========================================================================

    section(
        "LOAD IMPUTATION-READY EXPRESSION SOURCE"
    )

    source = ad.read_h5ad(
        EXPRESSION_SOURCE,
        backed="r"
    )


    print(
        "Shape:",
        source.shape
    )

    print(
        "Layers:",
        list(
            source.layers.keys()
        )
    )

    print(
        "obsm:",
        list(
            source.obsm.keys()
        )
    )


    if source.n_vars != 300:

        raise RuntimeError(
            f"Expected 300 Xenium genes; "
            f"found {source.n_vars}"
        )


    for col in [
        "BrNum",
        "donor_cell_id",
    ]:

        if col not in source.obs.columns:

            raise RuntimeError(
                f"{col} missing from expression source."
            )


    # =========================================================================
    # DONOR-BY-DONOR TARGET PREPARATION
    # =========================================================================

    section(
        "PREPARE 23 DONOR TARGETS"
    )

    summary_rows = []


    for i, row in enumerate(
        meta.itertuples(
            index=False
        ),
        start=1,
    ):

        donor = str(
            row.patient_id
        )

        dx = str(
            row.Dx
        )


        if dx == "NTC":

            experiment = "ex1_ntc"
            out_dir = EX1_DIR

        else:

            experiment = "ex2_scz"
            out_dir = EX2_DIR


        out_path = (
            out_dir
            / (
                f"spatial_data_xenium_"
                f"{donor}_"
                f"{experiment}.h5ad"
            )
        )


        print(
            f"\n[{i:02d}/23] "
            f"{donor} | {dx} | {experiment}",
            flush=True
        )


        # ---------------------------------------------------------------------
        # Extract donor from expression/QC object
        # ---------------------------------------------------------------------

        donor_mask = (
            source.obs[
                "BrNum"
            ]
            .astype(str)
            .eq(
                donor
            )
            .to_numpy()
        )


        original_n_cells = int(
            donor_mask.sum()
        )


        if original_n_cells == 0:

            raise RuntimeError(
                f"{donor}: no cells found."
            )


        target = (
            source[
                donor_mask,
                :
            ]
            .to_memory()
        )


        print(
            "  source cells       :",
            f"{target.n_obs:,}"
        )


        # ---------------------------------------------------------------------
        # Build annotation alignment
        # ---------------------------------------------------------------------

        target_ids = clean_string(
            target.obs[
                "donor_cell_id"
            ]
        )


        if target_ids.duplicated().any():

            raise RuntimeError(
                f"{donor}: duplicate donor_cell_id "
                "in expression target."
            )


        matched = annotation.reindex(
            target_ids.to_numpy()
        )


        # Cell exists in annotation object.
        has_annotation_record = (
            matched[
                "BrNum"
            ]
            .notna()
        )


        # Verify matched records belong to the correct donor.
        wrong_donor = (
            has_annotation_record
            & (
                matched[
                    "BrNum"
                ]
                .astype(str)
                != donor
            )
        )


        if wrong_donor.any():

            raise RuntimeError(
                f"{donor}: {wrong_donor.sum()} "
                "annotation records map to wrong donor."
            )


        # ---------------------------------------------------------------------
        # Require final layer + cell-type annotation
        # ---------------------------------------------------------------------

        has_layer = (
            matched[
                "layer_annotation"
            ]
            .notna()
        )


        has_celltype = (
            matched[
                "cell_type_annotation"
            ]
            .notna()
        )


        keep = (
            has_annotation_record
            & has_layer
            & has_celltype
        ).to_numpy()


        n_keep = int(
            keep.sum()
        )

        n_drop = int(
            (~keep).sum()
        )


        n_missing_record = int(
            (~has_annotation_record)
            .sum()
        )


        n_missing_layer = int(
            (
                has_annotation_record
                & ~has_layer
            )
            .sum()
        )


        n_missing_celltype = int(
            (
                has_annotation_record
                & ~has_celltype
            )
            .sum()
        )


        print(
            "  matched annotation :",
            f"{int(has_annotation_record.sum()):,}"
            f"/{original_n_cells:,}"
        )

        print(
            "  missing record     :",
            n_missing_record
        )

        print(
            "  missing layer      :",
            n_missing_layer
        )

        print(
            "  missing cell type  :",
            n_missing_celltype
        )

        print(
            "  cells removed      :",
            n_drop
        )

        print(
            "  final cells        :",
            f"{n_keep:,}"
        )


        if n_keep == 0:

            raise RuntimeError(
                f"{donor}: no annotated cells remain."
            )


        # ---------------------------------------------------------------------
        # REMOVE unannotated cells
        # ---------------------------------------------------------------------

        target = (
            target[
                keep,
                :
            ]
            .copy()
        )


        matched = (
            matched.iloc[
                np.flatnonzero(
                    keep
                )
            ]
            .copy()
        )


        # ---------------------------------------------------------------------
        # Add annotation metadata
        # ---------------------------------------------------------------------

        for col in ANNOTATION_COLUMNS:

            target.obs[
                col
            ] = (
                matched[
                    col
                ]
                .to_numpy()
            )


        # ---------------------------------------------------------------------
        # Final strict annotation QC
        # ---------------------------------------------------------------------

        for col in REQUIRED_FINAL_ANNOTATIONS:

            if target.obs[
                col
            ].isna().any():

                raise RuntimeError(
                    f"{donor}: missing values remain in {col}"
                )


        if target.obs[
            "donor_cell_id"
        ].duplicated().any():

            raise RuntimeError(
                f"{donor}: duplicate donor_cell_id after filtering."
            )


        # ---------------------------------------------------------------------
        # Authoritative diagnosis
        # ---------------------------------------------------------------------

        if "Dx" in target.obs.columns:

            target.obs[
                "Dx_original"
            ] = (
                target.obs[
                    "Dx"
                ]
                .astype(str)
                .to_numpy()
            )


        target.obs[
            "Dx"
        ] = dx


        # ---------------------------------------------------------------------
        # Expression / spatial QC
        # ---------------------------------------------------------------------

        if target.n_vars != 300:

            raise RuntimeError(
                f"{donor}: expected 300 genes; "
                f"found {target.n_vars}"
            )


        if "spatial" not in target.obsm:

            raise RuntimeError(
                f"{donor}: obsm['spatial'] missing."
            )


        if (
            target.obsm[
                "spatial"
            ].shape
            != (
                target.n_obs,
                2
            )
        ):

            raise RuntimeError(
                f"{donor}: unexpected spatial shape "
                f"{target.obsm['spatial'].shape}"
            )


        # ---------------------------------------------------------------------
        # Output
        # ---------------------------------------------------------------------

        if (
            out_path.exists()
            and not args.overwrite
        ):

            status = (
                "EXISTS_NOT_OVERWRITTEN"
            )

            print(
                "  EXISTS — not overwriting"
            )


        else:

            target.write_h5ad(
                out_path,
                compression="gzip"
            )

            status = "WRITTEN"

            print(
                "  WROTE:",
                out_path
            )


        summary_rows.append({
            "donor":
                donor,

            "Dx":
                dx,

            "experiment":
                experiment,

            "source_cells":
                original_n_cells,

            "cells_removed":
                n_drop,

            "final_cells":
                target.n_obs,

            "missing_annotation_record":
                n_missing_record,

            "missing_layer_annotation":
                n_missing_layer,

            "missing_cell_type_annotation":
                n_missing_celltype,

            "n_genes":
                target.n_vars,

            "target_file":
                str(
                    out_path
                ),

            "status":
                status,
        })


        del target


    source.file.close()


    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================

    summary = pd.DataFrame(
        summary_rows
    )


    summary.to_csv(
        SUMMARY_OUT,
        index=False
    )


    section(
        "FINAL SUMMARY"
    )


    print(
        summary[
            [
                "donor",
                "Dx",
                "source_cells",
                "cells_removed",
                "final_cells",
                "status",
            ]
        ]
        .to_string(
            index=False
        )
    )


    print(
        "\nNTC donors:",
        int(
            summary[
                "Dx"
            ]
            .eq("NTC")
            .sum()
        )
    )


    print(
        "SCZ donors:",
        int(
            summary[
                "Dx"
            ]
            .eq("SCZ")
            .sum()
        )
    )


    print(
        "Total donors:",
        len(
            summary
        )
    )


    print(
        "\nTotal source cells :",
        f"{summary['source_cells'].sum():,}"
    )


    print(
        "Total cells removed:",
        int(
            summary[
                "cells_removed"
            ].sum()
        )
    )


    print(
        "Final total cells  :",
        f"{summary['final_cells'].sum():,}"
    )


    print(
        "\nMissing annotation records:",
        int(
            summary[
                "missing_annotation_record"
            ].sum()
        )
    )


    print(
        "Missing layer among matched cells:",
        int(
            summary[
                "missing_layer_annotation"
            ].sum()
        )
    )


    print(
        "Missing cell type among matched cells:",
        int(
            summary[
                "missing_cell_type_annotation"
            ].sum()
        )
    )


    print(
        "\nSummary:"
    )

    print(
        SUMMARY_OUT
    )


    section(
        "DONE"
    )


if __name__ == "__main__":
    main()
