#!/usr/bin/env python3

from pathlib import Path
import pickle

import anndata as ad
import pandas as pd


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

OUT_DIR = ROOT / "outputs/celltype_comparison/raw_proportions"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# INPUTS
# =============================================================================

HUUKI = (
    ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
      "huuki_snrna_reference_full_allgenes.h5ad"
)

BAITUK_RDS = (
    ROOT
    / "data/raw/Baituk/annotations_final.RDS"
)

XENIUM_CANDIDATES = [
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad",

    ROOT
    / "data/processed/imputation_beta/Br8667/"
      "xenium_Br8667_annotated.h5ad",

    ROOT
    / "data/processed/imputation_beta/Br8667/"
      "spatial_data_xenium_Br8667_vista_qc.h5ad",
]


HUUKI_COLUMN = "cellType_broad_hc"
XENIUM_COLUMN = "cell_type_annotation"


# =============================================================================
# HELPER
# =============================================================================

def summarize(dataset, annotation_column, labels):

    labels = pd.Series(labels).dropna().astype(str)

    counts = (
        labels
        .value_counts()
        .rename_axis("cell_type")
        .reset_index(name="n_cells")
    )

    total = counts["n_cells"].sum()

    counts["percentage"] = (
        counts["n_cells"]
        / total
        * 100
    )

    counts.insert(
        0,
        "annotation_column",
        annotation_column,
    )

    counts.insert(
        0,
        "dataset",
        dataset,
    )

    counts["total_cells"] = total

    return counts


# =============================================================================
# HUUKI
# =============================================================================

print("=" * 80)
print("HUUKI")
print("=" * 80)

if not HUUKI.exists():
    raise FileNotFoundError(HUUKI)

adata = ad.read_h5ad(
    HUUKI,
    backed="r",
)

if HUUKI_COLUMN not in adata.obs.columns:
    raise KeyError(
        f"{HUUKI_COLUMN} not found.\n"
        f"Available columns:\n{list(adata.obs.columns)}"
    )

huuki = summarize(
    dataset="Huuki",
    annotation_column=HUUKI_COLUMN,
    labels=adata.obs[HUUKI_COLUMN],
)

print(huuki.to_string(index=False))

try:
    adata.file.close()
except Exception:
    pass


# =============================================================================
# XENIUM
# =============================================================================

print()
print("=" * 80)
print("XENIUM")
print("=" * 80)

XENIUM = next(
    (p for p in XENIUM_CANDIDATES if p.exists()),
    None,
)

if XENIUM is None:
    raise FileNotFoundError(
        "No Xenium annotated H5AD found."
    )

print("Using:", XENIUM)

adata = ad.read_h5ad(
    XENIUM,
    backed="r",
)

if XENIUM_COLUMN not in adata.obs.columns:
    raise KeyError(
        f"{XENIUM_COLUMN} not found.\n"
        f"Available columns:\n{list(adata.obs.columns)}"
    )

xenium = summarize(
    dataset="Xenium",
    annotation_column=XENIUM_COLUMN,
    labels=adata.obs[XENIUM_COLUMN],
)

print(xenium.to_string(index=False))

try:
    adata.file.close()
except Exception:
    pass


# =============================================================================
# SAVE TEMPORARY HUUKI + XENIUM
# Baituk will be appended from R because annotations_final.RDS is an R object.
# =============================================================================

hx = pd.concat(
    [huuki, xenium],
    ignore_index=True,
)

hx.to_csv(
    OUT_DIR / "huuki_xenium_raw_celltype_proportions.csv",
    index=False,
)

print()
print("Saved Huuki/Xenium temporary table.")
