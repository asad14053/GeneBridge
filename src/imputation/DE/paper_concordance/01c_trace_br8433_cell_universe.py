#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

DONOR = "Br8433"

PAPER_LABELS = (
    ROOT
    / "data/reference/paper_xenium/"
      "label_transfer_N24_k50_smoothed_labels.csv"
)

PAPER_OUTLIERS = (
    ROOT
    / "data/reference/paper_xenium/"
      "outlier_ids.csv"
)

MANIFEST = (
    ROOT
    / "data/metadata/"
      "envi_production_23donors.tsv"
)

ANNOTATED_N24 = (
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)


def section(x):
    print("\n" + "=" * 100)
    print(x)
    print("=" * 100)


# =============================================================================
# PAPER CELL LABELS
# =============================================================================

section("LOAD PAPER CELL-LEVEL SPATIAL LABELS")

labels = pd.read_csv(PAPER_LABELS)

print("Columns:", labels.columns.tolist())
print("Shape:", labels.shape)

cell_col = labels.columns[0]

if "predictions_smooth" in labels.columns:
    pred_col = "predictions_smooth"
else:
    pred_col = labels.columns[1]

labels[cell_col] = labels[cell_col].astype(str)

paper_br = labels[
    labels[cell_col].str.startswith(f"{DONOR}_")
].copy()

paper_br = paper_br.rename(
    columns={
        cell_col: "cell_id",
        pred_col: "paper_predictions_smooth",
    }
)

print("\nPaper Br8433 cells:", len(paper_br))

print("\nPaper Br8433 cells by domain:")
print(
    paper_br["paper_predictions_smooth"]
    .value_counts()
    .sort_index()
    .to_string()
)


# =============================================================================
# FIND TARGET H5AD
# =============================================================================

section("LOAD GENEBRIDGE TARGET")

manifest = pd.read_csv(
    MANIFEST,
    sep="\t"
)

row = manifest[
    manifest["donor"].astype(str) == DONOR
]

if len(row) != 1:
    raise RuntimeError(
        f"Expected exactly one {DONOR} row in manifest; found {len(row)}."
    )

target_str = str(
    row.iloc[0]["target"]
)

TARGET = Path(target_str)

if not TARGET.is_absolute():
    TARGET = ROOT / TARGET

print("Target:", TARGET)

if not TARGET.exists():
    raise FileNotFoundError(TARGET)

target = ad.read_h5ad(
    TARGET,
    backed="r"
)

target_ids = set(
    target.obs_names.astype(str)
)

print("Target cells:", len(target_ids))


# =============================================================================
# ANNOTATED N24
# =============================================================================

section("LOAD GENEBRIDGE ANNOTATED N24")

annot = ad.read_h5ad(
    ANNOTATED_N24,
    backed="r"
)

if "BrNum" not in annot.obs.columns:
    raise RuntimeError(
        "BrNum missing from annotated N24 object."
    )

mask = (
    annot.obs["BrNum"]
    .astype(str)
    == DONOR
)

annot_ids = set(
    annot.obs_names[
        mask
    ].astype(str)
)

print("Annotated N24 Br8433 cells:", len(annot_ids))


# =============================================================================
# CELL-ID SET COMPARISON
# =============================================================================

section("COMPARE CELL-ID UNIVERSES")

paper_ids = set(
    paper_br["cell_id"]
)

paper_not_target = sorted(
    paper_ids - target_ids
)

target_not_paper = sorted(
    target_ids - paper_ids
)

paper_not_annot = sorted(
    paper_ids - annot_ids
)

annot_not_paper = sorted(
    annot_ids - paper_ids
)


print("Paper cells              :", len(paper_ids))
print("Target cells             :", len(target_ids))
print("Annotated N24 cells      :", len(annot_ids))

print()
print("Paper but NOT target     :", len(paper_not_target))
print("Target but NOT paper     :", len(target_not_paper))
print("Paper but NOT annotation :", len(paper_not_annot))
print("Annotation but NOT paper :", len(annot_not_paper))


# =============================================================================
# DETAILS OF PAPER CELLS ABSENT FROM TARGET
# =============================================================================

section("PAPER CELLS MISSING FROM TARGET")

if paper_not_target:

    x = paper_br[
        paper_br["cell_id"].isin(
            paper_not_target
        )
    ].copy()

    print(
        x[
            [
                "cell_id",
                "paper_predictions_smooth",
            ]
        ].to_string(
            index=False
        )
    )

    print("\nMissing target cells by domain:")

    print(
        x["paper_predictions_smooth"]
        .value_counts()
        .sort_index()
        .to_string()
    )

else:
    print("NONE")


# =============================================================================
# DETAILS OF PAPER CELLS ABSENT FROM ANNOTATION
# =============================================================================

section("PAPER CELLS MISSING FROM GENEBRIDGE ANNOTATION")

if paper_not_annot:

    x2 = paper_br[
        paper_br["cell_id"].isin(
            paper_not_annot
        )
    ].copy()

    print(
        x2[
            [
                "cell_id",
                "paper_predictions_smooth",
            ]
        ].to_string(
            index=False
        )
    )

    print("\nMissing annotation cells by domain:")

    print(
        x2["paper_predictions_smooth"]
        .value_counts()
        .sort_index()
        .to_string()
    )

else:
    print("NONE")


# =============================================================================
# CHECK PAPER QC OUTLIER LIST
# =============================================================================

section("CHECK AGAINST PAPER OUTLIER LIST")

outliers = pd.read_csv(
    PAPER_OUTLIERS
)

# Flatten all values because the published file may use
# an arbitrary column name.
outlier_ids = set(
    outliers.astype(str)
    .stack()
    .values
)


candidate_missing = sorted(
    paper_ids - annot_ids
)

print(
    "Candidate missing cells:",
    len(candidate_missing)
)

for cid in candidate_missing:

    print(
        cid,
        "paper_domain=",
        paper_br.loc[
            paper_br["cell_id"] == cid,
            "paper_predictions_smooth"
        ].iloc[0],
        "in_paper_outlier_ids=",
        cid in outlier_ids,
        "in_target=",
        cid in target_ids,
        "in_annotation=",
        cid in annot_ids,
    )


# =============================================================================
# DECISION
# =============================================================================

section("INTERPRETATION")

if (
    len(paper_not_annot) == 5
    and
    len(paper_not_target) == 0
):

    print(
        "RESULT: The five paper cells exist in the target H5AD "
        "but were lost from the GeneBridge annotation object."
    )

    print(
        "ACTION: repair the annotation join; do NOT change expression data."
    )

elif (
    len(paper_not_target) == 5
    and
    len(paper_not_annot) == 5
):

    print(
        "RESULT: The five paper cells are absent from both "
        "the target expression H5AD and annotation object."
    )

    print(
        "ACTION: trace the GeneBridge target extraction/QC step."
    )

elif (
    len(paper_not_target) == 0
    and
    len(paper_not_annot) == 0
):

    print(
        "RESULT: Cell IDs are present in both objects. "
        "The mismatch must come from grouping/filtering logic."
    )

else:

    print(
        "RESULT: Mixed cell-universe discrepancy. "
        "Inspect the cell lists printed above."
    )


print("\nFINAL STATUS: TRACE COMPLETE")
