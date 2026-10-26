#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

SCZ_FILE = ROOT / "data/processed/snrnaseq/Baituk/baituk_snrna_reference_SCZ_full_allgenes.h5ad"
NTC_FILE = ROOT / "data/processed/snrnaseq/Baituk/baituk_snrna_reference_NTC_full_allgenes.h5ad"

ANNOT_FILE = (
    ROOT /
    "outputs/snrnaseq_branch/04_baituk_celltype_annotation/"
    "baituk_selected_cell_annotations.csv"
)

OUT = ROOT / "outputs/snrnaseq_branch/05_baituk_scz_ntc_de"
OUT.mkdir(parents=True, exist_ok=True)

DONOR_COUNTS = OUT / "baituk_donor_pseudobulk_counts.csv.gz"
DONOR_META   = OUT / "baituk_donor_pseudobulk_metadata.csv"

CELL_COUNTS = OUT / "baituk_medclean_donor_celltype_pseudobulk_counts.csv.gz"
CELL_META   = OUT / "baituk_medclean_donor_celltype_pseudobulk_metadata.csv"


def section(x):
    print("\n" + "=" * 110)
    print(x)
    print("=" * 110)


def aggregate_sparse(X, labels):
    """
    Aggregate cells × genes sparse matrix into group × genes sums.
    """

    labels = pd.Series(labels, dtype="string")

    codes, groups = pd.factorize(
        labels,
        sort=True
    )

    if np.any(codes < 0):
        raise RuntimeError("Missing aggregation labels found.")

    indicator = sparse.csr_matrix(
        (
            np.ones(len(codes), dtype=np.int64),
            (
                codes,
                np.arange(len(codes))
            )
        ),
        shape=(len(groups), len(codes))
    )

    out = indicator @ X

    if sparse.issparse(out):
        out = out.toarray()

    out = np.asarray(out)

    n_cells = np.bincount(
        codes,
        minlength=len(groups)
    )

    return list(groups), out, n_cells


# ==================================================================================================
# LOAD MED.CLEAN ANNOTATION
# ==================================================================================================

section("LOAD BAITUK MED.CLEAN CELL ANNOTATIONS")

annot = pd.read_csv(
    ANNOT_FILE,
    low_memory=False
)

print("Annotation shape:", annot.shape)
print("Columns:", list(annot.columns))

if "cell_type" not in annot.columns:
    raise RuntimeError(
        "Expected column 'cell_type' in selected annotation."
    )

print("\nCell types:")
print(
    annot["cell_type"]
    .value_counts()
    .to_string()
)

if annot["cell_type"].nunique() != 14:
    raise RuntimeError(
        f"Expected 14 med.clean cell types, "
        f"found {annot['cell_type'].nunique()}."
    )

bad = {"Glia", "Other"}

if set(annot["cell_type"].dropna()) & bad:
    raise RuntimeError(
        "Glia/Other unexpectedly present in med.clean annotation."
    )


# ==================================================================================================
# DETECT CELL-ID COLUMN
# ==================================================================================================

section("DETECT ANNOTATION CELL-ID COLUMN")

# Use SCZ IDs as a reference to identify which annotation column is the barcode/cell ID.
tmp = sc.read_h5ad(
    SCZ_FILE,
    backed="r"
)

target_ids = set(
    tmp.obs["original_cell_id"]
    .astype(str)
)

del tmp

best_col = None
best_overlap = -1

for col in annot.columns:

    if col == "cell_type":
        continue

    vals = annot[col].astype(str)

    overlap = vals.isin(target_ids).sum()

    print(
        f"{col:40s} overlap with SCZ original_cell_id = {overlap}"
    )

    if overlap > best_overlap:
        best_overlap = overlap
        best_col = col


if best_overlap <= 0:
    raise RuntimeError(
        "Could not identify annotation cell-ID column."
    )

print("\nSelected annotation ID column:", best_col)
print("SCZ overlap:", best_overlap)


# Validate duplicate IDs
dup = annot.duplicated(
    subset=best_col,
    keep=False
)

if dup.any():

    x = (
        annot.loc[dup, [best_col, "cell_type"]]
        .groupby(best_col)["cell_type"]
        .nunique()
    )

    conflicting = x[x > 1]

    if len(conflicting):
        raise RuntimeError(
            f"{len(conflicting)} annotation IDs have conflicting cell types."
        )

    annot = annot.drop_duplicates(
        subset=best_col
    )


annotation_map = pd.Series(
    annot["cell_type"].astype(str).values,
    index=annot[best_col].astype(str).values
)


# ==================================================================================================
# PROCESS EACH DIAGNOSIS FILE
# ==================================================================================================

donor_count_frames = []
donor_meta_frames = []

cell_count_frames = []
cell_meta_frames = []

gene_names_reference = None


for diagnosis, infile in [
    ("SCZ", SCZ_FILE),
    ("NTC", NTC_FILE),
]:

    section(f"PROCESS {diagnosis}")

    a = sc.read_h5ad(
        infile
    )

    print("Input:", infile)
    print("Shape:", a.shape)

    if not sparse.issparse(a.X):
        raise RuntimeError(
            f"{diagnosis}: expected sparse raw-count .X"
        )

    # ----------------------------------------------------------------------
    # Gene names
    # ----------------------------------------------------------------------

    if "gene_name" in a.var.columns:
        genes = (
            a.var["gene_name"]
            .astype(str)
            .tolist()
        )
    else:
        genes = (
            a.var_names
            .astype(str)
            .tolist()
        )

    if len(set(genes)) != len(genes):
        raise RuntimeError(
            f"{diagnosis}: gene names are not unique."
        )

    if gene_names_reference is None:
        gene_names_reference = genes
    else:
        if genes != gene_names_reference:
            raise RuntimeError(
                "SCZ and NTC gene universes/order are not identical."
            )

    print("Genes:", len(genes))

    # ----------------------------------------------------------------------
    # Validate diagnosis
    # ----------------------------------------------------------------------

    observed_dx = set(
        a.obs["diagnosis"]
        .astype(str)
        .unique()
    )

    if observed_dx != {diagnosis}:
        raise RuntimeError(
            f"{diagnosis}: diagnosis mismatch: {observed_dx}"
        )

    obs = a.obs.copy()

    obs["original_cell_id"] = (
        obs["original_cell_id"]
        .astype(str)
    )

    obs["donor_id"] = (
        obs["donor_id"]
        .astype(str)
    )

    obs["sample"] = (
        obs["sample"]
        .astype(str)
    )

    obs["diagnosis"] = diagnosis


    # ==============================================================================================
    # DONOR-LEVEL PSEUDOBULK
    # ==============================================================================================

    section(f"{diagnosis}: DONOR-LEVEL PSEUDOBULK")

    donor_pb_id = (
        "donor::"
        + obs["donor_id"]
    )

    donor_groups, donor_counts, donor_ncells = aggregate_sparse(
        a.X,
        donor_pb_id
    )

    donor_meta = []

    for pb_id, ncells in zip(
        donor_groups,
        donor_ncells
    ):

        donor = pb_id.replace(
            "donor::",
            "",
            1
        )

        z = obs.loc[
            obs["donor_id"] == donor
        ]

        donor_meta.append({
            "pseudobulk_id": pb_id,
            "donor_id": donor,
            "diagnosis": diagnosis,
            "n_cells": int(ncells),
            "samples": ";".join(
                sorted(
                    z["sample"]
                    .unique()
                )
            )
        })


    donor_meta = pd.DataFrame(
        donor_meta
    )

    donor_df = pd.DataFrame(
        donor_counts.T,
        index=genes,
        columns=donor_groups
    )

    donor_df.index.name = "gene"

    donor_count_frames.append(
        donor_df
    )

    donor_meta_frames.append(
        donor_meta
    )


    print(donor_meta.to_string(index=False))


    # ==============================================================================================
    # ATTACH MED.CLEAN CELL TYPE
    # ==============================================================================================

    section(f"{diagnosis}: ATTACH MED.CLEAN CELL TYPES")

    obs["cell_type"] = (
        obs["original_cell_id"]
        .map(annotation_map)
    )

    matched = obs["cell_type"].notna().sum()
    frac = matched / len(obs)

    print("Cells total   :", len(obs))
    print("Cells matched :", matched)
    print("Match fraction:", f"{frac:.4%}")

    if frac < 0.99:
        raise RuntimeError(
            f"{diagnosis}: only {frac:.2%} cells matched annotation. "
            "Expected >=99%."
        )

    keep = obs["cell_type"].notna()

    obs_cell = obs.loc[
        keep
    ].copy()

    X_cell = a.X[
        np.asarray(keep)
    ]


    print("\nCell-type counts:")
    print(
        obs_cell["cell_type"]
        .value_counts()
        .to_string()
    )


    # ==============================================================================================
    # DONOR × CELL-TYPE PSEUDOBULK
    # ==============================================================================================

    section(
        f"{diagnosis}: DONOR × MED.CLEAN CELL-TYPE PSEUDOBULK"
    )

    cell_pb_id = (
        "celltype::"
        + obs_cell["donor_id"].astype(str)
        + "::"
        + obs_cell["cell_type"].astype(str)
    )

    cell_groups, cell_counts, cell_ncells = aggregate_sparse(
        X_cell,
        cell_pb_id
    )

    cell_meta = []

    for pb_id, ncells in zip(
        cell_groups,
        cell_ncells
    ):

        parts = pb_id.split(
            "::",
            2
        )

        donor = parts[1]
        cell_type = parts[2]

        z = obs_cell.loc[
            (obs_cell["donor_id"] == donor)
            &
            (obs_cell["cell_type"] == cell_type)
        ]

        cell_meta.append({
            "pseudobulk_id": pb_id,
            "donor_id": donor,
            "diagnosis": diagnosis,
            "cell_type": cell_type,
            "n_cells": int(ncells),
            "samples": ";".join(
                sorted(
                    z["sample"]
                    .unique()
                )
            )
        })


    cell_meta = pd.DataFrame(
        cell_meta
    )

    cell_df = pd.DataFrame(
        cell_counts.T,
        index=genes,
        columns=cell_groups
    )

    cell_df.index.name = "gene"

    cell_count_frames.append(
        cell_df
    )

    cell_meta_frames.append(
        cell_meta
    )


    print(
        cell_meta[
            [
                "donor_id",
                "diagnosis",
                "cell_type",
                "n_cells"
            ]
        ]
        .sort_values(
            ["cell_type", "donor_id"]
        )
        .to_string(index=False)
    )

    del a


# ==================================================================================================
# COMBINE SCZ + NTC
# ==================================================================================================

section("COMBINE DONOR-LEVEL PSEUDOBULK")

donor_counts = pd.concat(
    donor_count_frames,
    axis=1
)

donor_meta = pd.concat(
    donor_meta_frames,
    ignore_index=True
)

print(
    donor_meta["diagnosis"]
    .value_counts()
)

print(
    "Biological donors:",
    donor_meta["donor_id"].nunique()
)

if donor_meta["donor_id"].nunique() != 23:
    raise RuntimeError(
        "Expected 23 biological donors."
    )

if (
    donor_meta.loc[
        donor_meta["diagnosis"] == "SCZ",
        "donor_id"
    ].nunique()
    != 9
):
    raise RuntimeError("Expected 9 SCZ donors.")

if (
    donor_meta.loc[
        donor_meta["diagnosis"] == "NTC",
        "donor_id"
    ].nunique()
    != 14
):
    raise RuntimeError("Expected 14 NTC donors.")


# Save genes × pseudobulk
donor_counts.to_csv(
    DONOR_COUNTS,
    compression="gzip"
)

donor_meta.to_csv(
    DONOR_META,
    index=False
)


# ==================================================================================================
# COMBINE CELL-TYPE PSEUDOBULK
# ==================================================================================================

section("COMBINE DONOR × CELL-TYPE PSEUDOBULK")

cell_counts = pd.concat(
    cell_count_frames,
    axis=1
)

cell_meta = pd.concat(
    cell_meta_frames,
    ignore_index=True
)

print(
    "Pseudobulks:",
    len(cell_meta)
)

print(
    "Cell types:",
    cell_meta["cell_type"].nunique()
)

print(
    "Donors:",
    cell_meta["donor_id"].nunique()
)

print("\nMinimum cells per donor×celltype:")
print(
    cell_meta["n_cells"]
    .describe()
)

print("\nCoverage by cell type:")
coverage = (
    cell_meta
    .groupby(
        ["cell_type", "diagnosis"]
    )["donor_id"]
    .nunique()
    .unstack(
        fill_value=0
    )
)

print(coverage.to_string())


cell_counts.to_csv(
    CELL_COUNTS,
    compression="gzip"
)

cell_meta.to_csv(
    CELL_META,
    index=False
)


section("SAVED")

print("Donor counts:")
print(DONOR_COUNTS)

print("\nDonor metadata:")
print(DONOR_META)

print("\nCell-type counts:")
print(CELL_COUNTS)

print("\nCell-type metadata:")
print(CELL_META)

print(
    "\nFINAL STATUS: PASS — BAITUK SCZ/NTC DONOR + MED.CLEAN CELL-TYPE PSEUDOBULK COMPLETE"
)
