#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = ROOT / "data/metadata/envi_production_23donors.tsv"

IN_DIR = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "envi_layer_pseudobulk/per_donor"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "envi_layer_pseudobulk"
)

ORIGINAL_META = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "original_layer_pseudobulk/"
      "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)

ORIGINAL_GENES = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "original_layer_pseudobulk/"
      "original_xenium_300gene_order.txt"
)


print("=" * 100)
print("AGGREGATE 23 PARALLEL ENVI DONOR × LAYER RESULTS")
print("=" * 100)


manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
)

manifest["donor"] = manifest["donor"].astype(str)

if len(manifest) != 23:
    raise RuntimeError("Expected 23 manifest donors")


orig_meta = pd.read_csv(
    ORIGINAL_META
)

if len(orig_meta) != 161:
    raise RuntimeError("Expected 161 original pseudobulk rows")


canonical_ids = (
    orig_meta["pseudobulk_id"]
    .astype(str)
    .tolist()
)


original_genes = [
    x.strip()
    for x in ORIGINAL_GENES.read_text().splitlines()
    if x.strip()
]


measured_parts = []
full_parts = []
meta_parts = []
qc_parts = []
excluded_parts = []

canonical_full_genes = None
canonical_sources = None


for i, donor in enumerate(
    manifest["donor"],
    start=1,
):

    print(
        f"[{i:02d}/23] {donor}",
        flush=True,
    )

    measured_file = (
        IN_DIR
        / f"{donor}_ENVI_measured300_layer_pb.csv.gz"
    )

    full_file = (
        IN_DIR
        / f"{donor}_ENVI_full34987_layer_pb.csv.gz"
    )

    meta_file = (
        IN_DIR
        / f"{donor}_layer_metadata.csv"
    )

    gene_file = (
        IN_DIR
        / f"{donor}_gene_info.csv.gz"
    )

    qc_file = (
        IN_DIR
        / f"{donor}_qc.csv"
    )

    excluded_file = (
        IN_DIR
        / f"{donor}_excluded_cells.csv"
    )


    for f in [
        measured_file,
        full_file,
        meta_file,
        gene_file,
        qc_file,
        excluded_file,
    ]:
        if not f.exists():
            raise FileNotFoundError(f)


    measured = pd.read_csv(
        measured_file,
        index_col=0,
    )

    full = pd.read_csv(
        full_file,
        index_col=0,
    )

    md = pd.read_csv(
        meta_file
    )

    gi = pd.read_csv(
        gene_file
    )

    qc = pd.read_csv(
        qc_file
    )


    if measured.shape != (7, 300):
        raise RuntimeError(
            f"{donor}: measured shape {measured.shape}"
        )

    if full.shape != (7, 34987):
        raise RuntimeError(
            f"{donor}: full shape {full.shape}"
        )

    if measured.columns.tolist() != original_genes:
        raise RuntimeError(
            f"{donor}: measured gene order mismatch"
        )


    donor_genes = gi["gene"].astype(str).tolist()
    donor_source = gi["expression_source"].astype(str).tolist()


    if canonical_full_genes is None:

        canonical_full_genes = donor_genes
        canonical_sources = donor_source

    else:

        if set(donor_genes) != set(
            canonical_full_genes
        ):
            raise RuntimeError(
                f"{donor}: full gene set mismatch"
            )


        gi_map = (
            gi
            .set_index("gene")
            .loc[
                canonical_full_genes,
                "expression_source",
            ]
            .astype(str)
            .tolist()
        )

        if gi_map != canonical_sources:
            raise RuntimeError(
                f"{donor}: expression_source mismatch"
            )


        if donor_genes != canonical_full_genes:

            print(
                "    reordering full genes to canonical order"
            )

            full = full.loc[
                :,
                canonical_full_genes,
            ]


    measured_parts.append(
        measured
    )

    full_parts.append(
        full
    )

    meta_parts.append(
        md
    )

    qc_parts.append(
        qc
    )


    ex = pd.read_csv(
        excluded_file
    )

    if len(ex) > 0:
        excluded_parts.append(
            ex
        )


# =============================================================================
# Combine
# =============================================================================

measured = pd.concat(
    measured_parts,
    axis=0,
)

full = pd.concat(
    full_parts,
    axis=0,
)

meta = pd.concat(
    meta_parts,
    axis=0,
    ignore_index=True,
)

qc = pd.concat(
    qc_parts,
    axis=0,
    ignore_index=True,
)


if excluded_parts:
    excluded = pd.concat(
        excluded_parts,
        axis=0,
        ignore_index=True,
    )
else:
    excluded = pd.DataFrame(
        columns=[
            "BrNum",
            "cell_id",
            "reason",
        ]
    )


# =============================================================================
# Exact Step-2 sample universe
# =============================================================================

if set(measured.index) != set(canonical_ids):
    raise RuntimeError(
        "Measured pseudobulk IDs differ from original Xenium"
    )

if set(full.index) != set(canonical_ids):
    raise RuntimeError(
        "Full pseudobulk IDs differ from original Xenium"
    )


measured = measured.loc[
    canonical_ids
]

full = full.loc[
    canonical_ids
]

meta = (
    meta
    .set_index("pseudobulk_id")
    .loc[canonical_ids]
    .reset_index()
)


if measured.shape != (161, 300):
    raise RuntimeError(
        f"Measured final shape: {measured.shape}"
    )

if full.shape != (161, 34987):
    raise RuntimeError(
        f"Full final shape: {full.shape}"
    )

if len(meta) != 161:
    raise RuntimeError(
        f"Metadata rows: {len(meta)}"
    )


# =============================================================================
# Compare against Step 2
# =============================================================================

compare_cols = [
    "pseudobulk_id",
    "BrNum",
    "Dx",
    "Sex",
    "slide_id",
    "run_date",
    "predictions_smooth",
    "layer_annotation",
    "n_cells",
]


a = (
    orig_meta[compare_cols]
    .set_index("pseudobulk_id")
    .loc[canonical_ids]
)

b = (
    meta[compare_cols]
    .set_index("pseudobulk_id")
    .loc[canonical_ids]
)


for col in compare_cols[1:]:

    if not np.array_equal(
        a[col].astype(str).to_numpy(),
        b[col].astype(str).to_numpy(),
    ):
        raise RuntimeError(
            f"Original-vs-ENVI mismatch: {col}"
        )


if not np.allclose(
    orig_meta.set_index("pseudobulk_id")
    .loc[canonical_ids, "Age"]
    .astype(float),
    meta.set_index("pseudobulk_id")
    .loc[canonical_ids, "Age"]
    .astype(float),
):
    raise RuntimeError(
        "Original-vs-ENVI Age mismatch"
    )


# =============================================================================
# Gene information
# =============================================================================

gene_info = pd.DataFrame(
    {
        "gene": canonical_full_genes,
        "expression_source": canonical_sources,
    }
)


if (
    gene_info["expression_source"]
    .eq("measured_xenium")
    .sum()
    != 300
):
    raise RuntimeError(
        "Expected 300 measured genes"
    )

if (
    gene_info["expression_source"]
    .eq("envi_imputed")
    .sum()
    != 34687
):
    raise RuntimeError(
        "Expected 34687 imputed genes"
    )


# =============================================================================
# Save
# =============================================================================

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


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

qc_file = (
    OUT
    / "ENVI_donor_layer_validation.csv"
)

excluded_file = (
    OUT
    / "ENVI_cells_excluded_missing_layer.csv"
)


measured.to_csv(
    measured_file,
    compression="gzip",
)

full.to_csv(
    full_file,
    compression="gzip",
)

meta.to_csv(
    metadata_file,
    index=False,
)

gene_info.to_csv(
    gene_info_file,
    index=False,
)

qc.to_csv(
    qc_file,
    index=False,
)

excluded.to_csv(
    excluded_file,
    index=False,
)


# =============================================================================
# Final
# =============================================================================

print("\n" + "=" * 100)
print("FINAL SUMMARY")
print("=" * 100)

print("Donors:", meta["BrNum"].nunique())
print("Layers:", meta["layer_annotation"].nunique())
print("Pseudobulk samples:", len(meta))
print("Measured matrix:", measured.shape)
print("Full matrix:", full.shape)

print(
    "ENVI-imputed genes:",
    gene_info["expression_source"]
    .eq("envi_imputed")
    .sum(),
)

print(
    "Minimum donor-layer cells:",
    meta["n_cells"].min(),
)

print(
    "Maximum donor-layer cells:",
    meta["n_cells"].max(),
)

print(
    "Excluded cells:",
    len(excluded),
)

print(
    "All measured donor validations:",
    qc["measured_layer_sum_matches"].all(),
)

print(
    "All full donor validations:",
    qc["full_layer_sum_matches"].all(),
)

print(
    "\nOriginal and ENVI pseudobulk IDs: IDENTICAL"
)

print(
    "Original and ENVI donor/layer cell counts: IDENTICAL"
)

print(
    "Original and ENVI covariates: IDENTICAL"
)

print(
    "\nFINAL STATUS: PASS — "
    "23 parallel donor outputs successfully aggregated."
)
