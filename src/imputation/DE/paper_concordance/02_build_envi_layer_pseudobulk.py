#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = (
    ROOT
    / "data/metadata/envi_production_23donors.tsv"
)

DONOR_META = (
    ROOT
    / "data/metadata/xenium_DE_metadata_23.csv"
)

ANNOTATED = (
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)

ORIGINAL_PB_META = (
    ROOT
    / "outputs/imputation_full/DE/"
      "paper_concordance/original_layer_pseudobulk/"
      "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)

ORIGINAL_GENES = (
    ROOT
    / "outputs/imputation_full/DE/"
      "paper_concordance/original_layer_pseudobulk/"
      "original_xenium_300gene_order.txt"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/"
      "paper_concordance/envi_layer_pseudobulk"
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

# Keep chunk size consistent with existing full-transcriptome builder.
CHUNK = 512


def section(title):
    print(
        "\n"
        + "=" * 100
        + f"\n{title}\n"
        + "=" * 100,
        flush=True,
    )


def sum_block(x):

    if sparse.issparse(x):

        return np.asarray(
            x.sum(axis=0)
        ).ravel().astype(
            np.float64,
            copy=False,
        )

    return np.asarray(
        x,
        dtype=np.float64,
    ).sum(
        axis=0,
        dtype=np.float64,
    )


# =============================================================================
# Metadata
# =============================================================================

section("LOAD MANIFEST AND CANONICAL DONOR METADATA")

manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
)

manifest["donor"] = (
    manifest["donor"]
    .astype(str)
    .str.strip()
)

manifest["Dx"] = (
    manifest["Dx"]
    .astype(str)
    .str.strip()
    .str.upper()
)


donor_meta = pd.read_csv(
    DONOR_META
)

donor_meta["BrNum"] = (
    donor_meta["BrNum"]
    .astype(str)
    .str.strip()
)


if len(manifest) != 23:
    raise RuntimeError(
        f"Expected 23 manifest donors; found {len(manifest)}"
    )

if len(donor_meta) != 23:
    raise RuntimeError(
        f"Expected 23 canonical donors; found {len(donor_meta)}"
    )


if set(manifest["donor"]) != set(donor_meta["BrNum"]):

    raise RuntimeError(
        "Manifest donors and canonical metadata donors differ."
    )


if (manifest["Dx"] == "NTC").sum() != 11:
    raise RuntimeError("Expected 11 NTC donors.")

if (manifest["Dx"] == "SCZ").sum() != 12:
    raise RuntimeError("Expected 12 SCZ donors.")


meta = manifest.merge(
    donor_meta[
        [
            "BrNum",
            "Dx",
            "Age",
            "Sex",
            "slide_id",
            "run_date",
        ]
    ],
    left_on="donor",
    right_on="BrNum",
    suffixes=("_manifest", "_canonical"),
    validate="one_to_one",
)


if not (
    meta["Dx_manifest"]
    == meta["Dx_canonical"]
).all():

    raise RuntimeError(
        "Diagnosis disagreement between manifest and canonical metadata."
    )


print("Donors:", len(meta))
print(
    "NTC:",
    (meta["Dx_manifest"] == "NTC").sum(),
)
print(
    "SCZ:",
    (meta["Dx_manifest"] == "SCZ").sum(),
)


# =============================================================================
# Original Step-2 sample universe
# =============================================================================

section("LOAD ORIGINAL XENIUM DONOR × LAYER SAMPLE UNIVERSE")

original_pb_meta = pd.read_csv(
    ORIGINAL_PB_META
)


if len(original_pb_meta) != 161:
    raise RuntimeError(
        f"Expected 161 original pseudobulks; "
        f"found {len(original_pb_meta)}"
    )


canonical_pb_ids = (
    original_pb_meta[
        "pseudobulk_id"
    ]
    .astype(str)
    .tolist()
)


original_gene_order = [
    x.strip()
    for x in ORIGINAL_GENES.read_text().splitlines()
    if x.strip()
]


if len(original_gene_order) != 300:

    raise RuntimeError(
        f"Expected 300 original genes; "
        f"found {len(original_gene_order)}"
    )


print(
    "Canonical pseudobulk samples:",
    len(canonical_pb_ids),
)

print(
    "Original measured genes:",
    len(original_gene_order),
)


# =============================================================================
# Load cell-level annotations
# =============================================================================

section("LOAD CELL-LEVEL SPATIAL ANNOTATIONS")

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
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
        "Annotated cell IDs are not globally unique."
    )


print(
    "Annotated cells:",
    len(annotation),
)

ann.file.close()


# =============================================================================
# Build ENVI donor × layer pseudobulk
# =============================================================================

section("BUILD ENVI DONOR × LAYER PSEUDOBULK")


canonical_full_genes = None
canonical_full_source = None

measured_rows = []
full_rows = []

pb_meta_rows = []
donor_qc_rows = []
excluded_rows = []


for i, row in enumerate(
    meta.itertuples(index=False),
    start=1,
):

    donor = str(row.donor)
    experiment = str(row.experiment)

    print(
        f"\n[{i:02d}/23] "
        f"{donor} | {row.Dx_manifest}",
        flush=True,
    )


    target_path = Path(
        row.target
    )


    after_path = (
        ROOT
        / f"data/processed/imputation_full/{experiment}/envi/{donor}"
        / f"spatial_data_xenium_{donor}_ENVI_full_transcriptome.h5ad"
    )


    if not target_path.exists():
        raise FileNotFoundError(target_path)

    if not after_path.exists():
        raise FileNotFoundError(after_path)


    # -------------------------------------------------------------------------
    # Open target and ENVI
    # -------------------------------------------------------------------------

    target = ad.read_h5ad(
        target_path,
        backed="r",
    )

    after = ad.read_h5ad(
        after_path,
        backed="r",
    )


    # Production ENVI must preserve cell universe and order.
    target_ids = pd.Index(
        target.obs_names.astype(str)
    )

    after_ids = pd.Index(
        after.obs_names.astype(str)
    )


    if len(target_ids) != len(after_ids):

        raise RuntimeError(
            f"{donor}: target/ENVI cell number differs: "
            f"{len(target_ids)} vs {len(after_ids)}"
        )


    if not np.array_equal(
        target_ids.to_numpy(),
        after_ids.to_numpy(),
    ):

        raise RuntimeError(
            f"{donor}: target and ENVI cell names/order differ."
        )


    # -------------------------------------------------------------------------
    # Cell annotation matching
    # -------------------------------------------------------------------------

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
            f"{donor}: {len(annotated_only)} annotated "
            "cells are absent from ENVI target."
        )


    for cell_id in target_only:

        excluded_rows.append(
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

    matched_ids = target_ids[
        matched_mask
    ]


    cell_ann = annotation.loc[
        matched_ids,
        [
            "BrNum",
            "layer_annotation",
            "predictions_smooth",
        ],
    ].copy()


    # -------------------------------------------------------------------------
    # Full transcriptome gene validation
    # -------------------------------------------------------------------------

    genes = (
        after.var_names
        .astype(str)
        .tolist()
    )


    if after.n_vars != 34987:

        raise RuntimeError(
            f"{donor}: expected 34,987 genes; "
            f"found {after.n_vars}"
        )


    if "expression_source" not in after.var.columns:

        raise RuntimeError(
            f"{donor}: expression_source missing from ENVI var."
        )


    source = (
        after.var[
            "expression_source"
        ]
        .astype(str)
        .tolist()
    )


    n_measured = sum(
        x == "measured_xenium"
        for x in source
    )

    n_imputed = sum(
        x == "envi_imputed"
        for x in source
    )


    if n_measured != 300:

        raise RuntimeError(
            f"{donor}: expected 300 measured genes; "
            f"found {n_measured}"
        )


    if n_imputed != 34687:

        raise RuntimeError(
            f"{donor}: expected 34,687 imputed genes; "
            f"found {n_imputed}"
        )


    # -------------------------------------------------------------------------
    # Establish canonical full-transcriptome gene order
    # -------------------------------------------------------------------------

    if canonical_full_genes is None:

        canonical_full_genes = genes
        canonical_full_source = source

        full_reorder_idx = np.arange(
            len(genes)
        )

    else:

        if set(genes) != set(
            canonical_full_genes
        ):

            raise RuntimeError(
                f"{donor}: full-transcriptome gene set differs."
            )


        donor_gene_index = pd.Index(
            genes
        )


        full_reorder_idx = (
            donor_gene_index
            .get_indexer(
                canonical_full_genes
            )
        )


        if np.any(
            full_reorder_idx < 0
        ):

            raise RuntimeError(
                f"{donor}: full gene alignment failed."
            )


        source_aligned = [
            source[j]
            for j in full_reorder_idx
        ]


        if (
            source_aligned
            != canonical_full_source
        ):

            raise RuntimeError(
                f"{donor}: expression_source mismatch "
                "after full gene alignment."
            )


        if genes != canonical_full_genes:

            print(
                "  NOTE: full gene order differs; "
                "canonical alignment will be applied."
            )


    # -------------------------------------------------------------------------
    # Measured-300 gene indices
    # -------------------------------------------------------------------------

    measured_native_idx = np.flatnonzero(
        np.asarray(source)
        == "measured_xenium"
    )


    measured_native_genes = [
        genes[j]
        for j in measured_native_idx
    ]


    if set(measured_native_genes) != set(
        original_gene_order
    ):

        missing = sorted(
            set(original_gene_order)
            - set(measured_native_genes)
        )

        extra = sorted(
            set(measured_native_genes)
            - set(original_gene_order)
        )

        raise RuntimeError(
            f"{donor}: measured-300 gene set mismatch.\n"
            f"Missing: {missing[:20]}\n"
            f"Extra: {extra[:20]}"
        )


    native_gene_index = pd.Index(
        genes
    )


    measured_idx_original_order = (
        native_gene_index
        .get_indexer(
            original_gene_order
        )
    )


    if np.any(
        measured_idx_original_order < 0
    ):

        raise RuntimeError(
            f"{donor}: measured gene alignment failed."
        )


    # -------------------------------------------------------------------------
    # Matrix sources
    # -------------------------------------------------------------------------

    # SAME source as existing 01_after_envi_measured300_qc.py
    measured_matrix = after.X

    # SAME source as existing 02_build_envi_full_pseudobulk.py
    if "count_scale" not in after.layers:

        raise RuntimeError(
            f"{donor}: count_scale layer missing."
        )

    full_matrix = after.layers[
        "count_scale"
    ]


    # -------------------------------------------------------------------------
    # Allocate seven layer totals
    # -------------------------------------------------------------------------

    measured_layer_totals = {
        layer: np.zeros(
            300,
            dtype=np.float64,
        )
        for layer in EXPECTED_LAYERS
    }


    full_layer_totals_native = {
        layer: np.zeros(
            34987,
            dtype=np.float64,
        )
        for layer in EXPECTED_LAYERS
    }


    layer_cell_counts = {
        layer: 0
        for layer in EXPECTED_LAYERS
    }


    matched_measured_total = np.zeros(
        300,
        dtype=np.float64,
    )


    matched_full_total_native = np.zeros(
        34987,
        dtype=np.float64,
    )


    # -------------------------------------------------------------------------
    # Chunk through cells
    # -------------------------------------------------------------------------

    n_cells = after.n_obs


    for start in range(
        0,
        n_cells,
        CHUNK,
    ):

        stop = min(
            start + CHUNK,
            n_cells,
        )


        chunk_ids = target_ids[
            start:stop
        ]


        chunk_is_annotated = chunk_ids.isin(
            ann_ids
        )


        if not chunk_is_annotated.any():

            continue


        # Cell positions inside this chunk.
        chunk_keep_positions = np.flatnonzero(
            chunk_is_annotated
        )


        kept_ids = chunk_ids[
            chunk_is_annotated
        ]


        kept_ann = annotation.loc[
            kept_ids,
            [
                "layer_annotation",
                "predictions_smooth",
            ],
        ]


        # -------------------------------------------------------------
        # Read measured 300 from after.X
        # -------------------------------------------------------------

        measured_chunk_all = measured_matrix[
            start:stop,
            measured_idx_original_order,
        ]


        if sparse.issparse(
            measured_chunk_all
        ):

            measured_chunk = measured_chunk_all[
                chunk_keep_positions,
                :
            ]

        else:

            measured_chunk = np.asarray(
                measured_chunk_all,
                dtype=np.float64,
            )[
                chunk_keep_positions,
                :
            ]


        # -------------------------------------------------------------
        # Read all 34,987 count-scale genes
        # -------------------------------------------------------------

        full_chunk_all = full_matrix[
            start:stop,
            :
        ]


        if sparse.issparse(
            full_chunk_all
        ):

            full_chunk = full_chunk_all[
                chunk_keep_positions,
                :
            ]

        else:

            full_chunk = np.asarray(
                full_chunk_all,
                dtype=np.float64,
            )[
                chunk_keep_positions,
                :
            ]


        matched_measured_total += sum_block(
            measured_chunk
        )

        matched_full_total_native += sum_block(
            full_chunk
        )


        # -------------------------------------------------------------
        # Aggregate chunk into seven layers
        # -------------------------------------------------------------

        for layer in EXPECTED_LAYERS:

            expected_spd = LAYER_TO_SPD[
                layer
            ]


            layer_mask = (
                kept_ann[
                    "layer_annotation"
                ]
                .eq(layer)
                .to_numpy()
            )


            if not layer_mask.any():

                continue


            if not (
                kept_ann.loc[
                    layer_mask,
                    "predictions_smooth",
                ]
                == expected_spd
            ).all():

                raise RuntimeError(
                    f"{donor} {layer}: "
                    "SpD/layer mapping mismatch."
                )


            layer_cell_counts[
                layer
            ] += int(
                layer_mask.sum()
            )


            if sparse.issparse(
                measured_chunk
            ):

                m_sub = measured_chunk[
                    layer_mask,
                    :
                ]

            else:

                m_sub = measured_chunk[
                    layer_mask,
                    :
                ]


            if sparse.issparse(
                full_chunk
            ):

                f_sub = full_chunk[
                    layer_mask,
                    :
                ]

            else:

                f_sub = full_chunk[
                    layer_mask,
                    :
                ]


            measured_layer_totals[
                layer
            ] += sum_block(
                m_sub
            )


            full_layer_totals_native[
                layer
            ] += sum_block(
                f_sub
            )


        if (
            start == 0
            or stop == n_cells
            or start % (CHUNK * 20) == 0
        ):

            print(
                f"  cells {stop:,}/{n_cells:,}",
                flush=True,
            )


    # -------------------------------------------------------------------------
    # Validation: seven layers sum to all matched cells
    # -------------------------------------------------------------------------

    measured_layer_sum = np.sum(
        np.stack(
            [
                measured_layer_totals[x]
                for x in EXPECTED_LAYERS
            ],
            axis=0,
        ),
        axis=0,
    )


    full_layer_sum_native = np.sum(
        np.stack(
            [
                full_layer_totals_native[x]
                for x in EXPECTED_LAYERS
            ],
            axis=0,
        ),
        axis=0,
    )


    if not np.allclose(
        measured_layer_sum,
        matched_measured_total,
        rtol=1e-10,
        atol=1e-6,
    ):

        raise RuntimeError(
            f"{donor}: measured-300 layer sum mismatch."
        )


    if not np.allclose(
        full_layer_sum_native,
        matched_full_total_native,
        rtol=1e-10,
        atol=1e-6,
    ):

        raise RuntimeError(
            f"{donor}: full-transcriptome layer sum mismatch."
        )


    # -------------------------------------------------------------------------
    # Store seven pseudobulk samples
    # -------------------------------------------------------------------------

    for layer in EXPECTED_LAYERS:

        spd = LAYER_TO_SPD[
            layer
        ]


        n_layer_cells = (
            layer_cell_counts[
                layer
            ]
        )


        if n_layer_cells < MIN_NCELLS:

            raise RuntimeError(
                f"{donor} {layer}: "
                f"{n_layer_cells} cells <10."
            )


        pb_id = (
            f"{donor}__{spd}__"
            f"{layer.replace('/', '-')}"
        )


        measured_values = (
            measured_layer_totals[
                layer
            ]
        )


        # Convert native donor order → canonical full gene order.
        full_values = (
            full_layer_totals_native[
                layer
            ][
                full_reorder_idx
            ]
        )


        if not np.isfinite(
            measured_values
        ).all():

            raise RuntimeError(
                f"{pb_id}: non-finite measured values."
            )


        if not np.isfinite(
            full_values
        ).all():

            raise RuntimeError(
                f"{pb_id}: non-finite full values."
            )


        if np.any(
            measured_values < 0
        ):

            raise RuntimeError(
                f"{pb_id}: negative measured values."
            )


        if np.any(
            full_values < 0
        ):

            raise RuntimeError(
                f"{pb_id}: negative full values."
            )


        measured_rows.append(
            pd.Series(
                measured_values,
                index=original_gene_order,
                name=pb_id,
            )
        )


        # IMPORTANT:
        # Use canonical_full_genes after reordering.
        # Do NOT use donor-native 'genes' here.
        full_rows.append(
            pd.Series(
                full_values,
                index=canonical_full_genes,
                name=pb_id,
            )
        )


        pb_meta_rows.append(
            {
                "pseudobulk_id": pb_id,
                "BrNum": donor,
                "Dx": row.Dx_manifest,
                "Age": float(row.Age),
                "Sex": str(row.Sex),
                "slide_id":
                    str(row.slide_id),
                "run_date":
                    str(row.run_date),
                "predictions_smooth":
                    spd,
                "layer_annotation":
                    layer,
                "n_cells":
                    n_layer_cells,
                "measured300_library":
                    float(
                        measured_values.sum()
                    ),
                "full34987_library":
                    float(
                        full_values.sum()
                    ),
            }
        )


    donor_qc_rows.append(
        {
            "BrNum": donor,
            "Dx": row.Dx_manifest,
            "target_cells":
                len(target_ids),
            "annotated_layer_cells":
                len(matched_ids),
            "excluded_unannotated_cells":
                len(target_only),
            "measured_genes":
                n_measured,
            "imputed_genes":
                n_imputed,
            "full_genes":
                after.n_vars,
            "measured_layer_sum_matches":
                True,
            "full_layer_sum_matches":
                True,
        }
    )


    print(
        "  layer cells:",
        layer_cell_counts,
    )

    print(
        "  excluded unannotated cells:",
        len(target_only),
    )


    target.file.close()
    after.file.close()


# =============================================================================
# Assemble matrices
# =============================================================================

section("ASSEMBLE MATRICES")


measured_pb = pd.DataFrame(
    measured_rows
)


full_pb = pd.DataFrame(
    full_rows
)


pb_meta = pd.DataFrame(
    pb_meta_rows
)


donor_qc = pd.DataFrame(
    donor_qc_rows
)


excluded = pd.DataFrame(
    excluded_rows,
    columns=[
        "BrNum",
        "cell_id",
        "reason",
    ],
)


print(
    "Measured-300 matrix:",
    measured_pb.shape,
)

print(
    "Full-transcriptome matrix:",
    full_pb.shape,
)

print(
    "Metadata:",
    pb_meta.shape,
)


# =============================================================================
# Force exact Step-2 pseudobulk sample order
# =============================================================================

if set(measured_pb.index) != set(
    canonical_pb_ids
):

    raise RuntimeError(
        "ENVI measured-300 pseudobulk IDs "
        "do not match original Xenium Step-2 IDs."
    )


if set(full_pb.index) != set(
    canonical_pb_ids
):

    raise RuntimeError(
        "ENVI full-transcriptome pseudobulk IDs "
        "do not match original Xenium Step-2 IDs."
    )


measured_pb = measured_pb.loc[
    canonical_pb_ids
]

full_pb = full_pb.loc[
    canonical_pb_ids
]


pb_meta = (
    pb_meta
    .set_index("pseudobulk_id")
    .loc[canonical_pb_ids]
    .reset_index()
)


if measured_pb.shape != (
    161,
    300,
):

    raise RuntimeError(
        f"Unexpected measured matrix shape: "
        f"{measured_pb.shape}"
    )


if full_pb.shape != (
    161,
    34987,
):

    raise RuntimeError(
        f"Unexpected full matrix shape: "
        f"{full_pb.shape}"
    )


if len(pb_meta) != 161:

    raise RuntimeError(
        "Expected 161 ENVI pseudobulk metadata rows."
    )


# =============================================================================
# Compare sample structure with original Xenium
# =============================================================================

section("VALIDATE AGAINST ORIGINAL XENIUM STEP 2")


compare_cols = [
    "pseudobulk_id",
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "run_date",
    "predictions_smooth",
    "layer_annotation",
    "n_cells",
]


original_compare = (
    original_pb_meta[
        compare_cols
    ]
    .set_index(
        "pseudobulk_id"
    )
    .loc[
        canonical_pb_ids
    ]
)


envi_compare = (
    pb_meta[
        compare_cols
    ]
    .set_index(
        "pseudobulk_id"
    )
    .loc[
        canonical_pb_ids
    ]
)


for col in [
    "BrNum",
    "Dx",
    "Sex",
    "slide_id",
    "run_date",
    "predictions_smooth",
    "layer_annotation",
    "n_cells",
]:

    if not (
        original_compare[col]
        .astype(str)
        .to_numpy()
        ==
        envi_compare[col]
        .astype(str)
        .to_numpy()
    ).all():

        raise RuntimeError(
            f"Original-vs-ENVI mismatch in {col}"
        )


if not np.allclose(
    original_compare[
        "Age"
    ].astype(float),
    envi_compare[
        "Age"
    ].astype(float),
):

    raise RuntimeError(
        "Original-vs-ENVI mismatch in Age."
    )


print(
    "Original and ENVI pseudobulk IDs:",
    "IDENTICAL",
)

print(
    "Original and ENVI donor/layer cell counts:",
    "IDENTICAL",
)

print(
    "Original and ENVI covariates:",
    "IDENTICAL",
)


# =============================================================================
# Gene information
# =============================================================================

gene_info = pd.DataFrame(
    {
        "gene":
            canonical_full_genes,

        "expression_source":
            canonical_full_source,
    }
)


if (
    gene_info[
        "expression_source"
    ]
    .eq("measured_xenium")
    .sum()
    != 300
):

    raise RuntimeError(
        "Canonical full matrix does not contain "
        "exactly 300 measured genes."
    )


if (
    gene_info[
        "expression_source"
    ]
    .eq("envi_imputed")
    .sum()
    != 34687
):

    raise RuntimeError(
        "Canonical full matrix does not contain "
        "exactly 34,687 imputed genes."
    )


# =============================================================================
# Save
# =============================================================================

section("SAVE OUTPUTS")


measured_file = (
    OUT
    / "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
)

full_file = (
    OUT
    / "ENVI_full34987_donor_layer_pseudobulk_countscale.csv.gz"
)

metadata_file = (
    OUT
    / "ENVI_donor_layer_pseudobulk_metadata.csv"
)

gene_info_file = (
    OUT
    / "ENVI_full34987_gene_info.csv"
)

donor_qc_file = (
    OUT
    / "ENVI_donor_layer_validation.csv"
)

excluded_file = (
    OUT
    / "ENVI_cells_excluded_missing_layer.csv"
)


measured_pb.to_csv(
    measured_file,
    compression="gzip",
)


full_pb.to_csv(
    full_file,
    compression="gzip",
)


pb_meta.to_csv(
    metadata_file,
    index=False,
)


gene_info.to_csv(
    gene_info_file,
    index=False,
)


donor_qc.to_csv(
    donor_qc_file,
    index=False,
)


excluded.to_csv(
    excluded_file,
    index=False,
)


# =============================================================================
# Final
# =============================================================================

section("FINAL SUMMARY")


print(
    "Donors:",
    pb_meta["BrNum"].nunique(),
)

print(
    "Layers:",
    pb_meta[
        "layer_annotation"
    ].nunique(),
)

print(
    "Pseudobulk samples:",
    len(pb_meta),
)

print(
    "Measured genes:",
    measured_pb.shape[1],
)

print(
    "Full transcriptome genes:",
    full_pb.shape[1],
)

print(
    "ENVI-imputed genes:",
    (
        gene_info[
            "expression_source"
        ]
        == "envi_imputed"
    ).sum(),
)

print(
    "Minimum donor-layer cells:",
    pb_meta[
        "n_cells"
    ].min(),
)

print(
    "Maximum donor-layer cells:",
    pb_meta[
        "n_cells"
    ].max(),
)

print(
    "Excluded cells:",
    len(excluded),
)


print("\nExcluded cells:")

if excluded.empty:

    print("NONE")

else:

    print(
        excluded.to_string(
            index=False
        )
    )


print(
    "\nMeasured layer sums validated for all donors:",
    donor_qc[
        "measured_layer_sum_matches"
    ].all(),
)

print(
    "Full layer sums validated for all donors:",
    donor_qc[
        "full_layer_sum_matches"
    ].all(),
)


print("\nOutput files:")
print(measured_file)
print(full_file)
print(metadata_file)
print(gene_info_file)
print(donor_qc_file)
print(excluded_file)


print(
    "\nFINAL STATUS: PASS — "
    "ENVI measured-300 and full-transcriptome "
    "donor × spatial-layer pseudobulks created "
    "on the exact original Xenium cell/layer universe."
)
