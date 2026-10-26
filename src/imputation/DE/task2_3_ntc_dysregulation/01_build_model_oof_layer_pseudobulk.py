#!/usr/bin/env python3

from pathlib import Path
import argparse
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

META = ROOT / "data/metadata/patient_xenium_visium_24_common_with_dx.csv"

OUT_ROOT = (
    ROOT
    / "outputs/imputation_full/DE/task2_3_ntc_dysregulation/pseudobulk"
)

MODEL_CONFIG = {
    "VISTA": {
        "root": "VISTA_N24",
        "method": "vista",
    },
    "gimVI": {
        "root": "GIMVI_N24",
        "method": "gimvi",
    },
    "ENVI": {
        "root": "ENVI_N24",
        "method": "envi",
    },
    "TransImp": {
        "root": "TRANSIMP_N24",
        "method": "transimpspa",
    },
    "Tangram": {
        "root": "TANGRAM_N24",
        "method": "tangram",
    },
    "SpaGE": {
        "root": "SPAGE_N24",
        "method": "spage",
    },
}

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
    print("\n" + "=" * 110)
    print(title)
    print("=" * 110, flush=True)


def sum_matrix(x):
    if sparse.issparse(x):
        return np.asarray(x.sum(axis=0)).ravel().astype(np.float64)

    return np.asarray(x).sum(
        axis=0,
        dtype=np.float64,
    )


parser = argparse.ArgumentParser()

parser.add_argument(
    "--model",
    required=True,
    choices=list(MODEL_CONFIG),
)

args = parser.parse_args()

MODEL = args.model
CFG = MODEL_CONFIG[MODEL]

OUT = OUT_ROOT / MODEL
OUT.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Metadata
# =============================================================================

section(f"LOAD METADATA — {MODEL}")

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
        f"Missing metadata columns: {sorted(missing)}"
    )

meta["BrNum"] = meta["BrNum"].astype(str).str.strip()
meta["Dx"] = meta["Dx"].astype(str).str.strip()

if len(meta) != 24:
    raise RuntimeError(
        f"Expected 24 donors; found {len(meta)}"
    )

if meta["BrNum"].nunique() != 24:
    raise RuntimeError("Expected 24 unique donors.")

dx_counts = meta["Dx"].value_counts().to_dict()

if dx_counts.get("NTC", 0) != 12:
    raise RuntimeError(
        f"Expected 12 NTC; found {dx_counts.get('NTC', 0)}"
    )

if dx_counts.get("SCZ", 0) != 12:
    raise RuntimeError(
        f"Expected 12 SCZ; found {dx_counts.get('SCZ', 0)}"
    )

print("Donors:", len(meta))
print("Dx:", dx_counts)


# =============================================================================
# Build pseudobulk
# =============================================================================

section(
    f"BUILD {MODEL} NTC-TRAINED OOF DONOR × SPATIAL-DOMAIN PSEUDOBULK"
)

gene_order = None

pb_arrays = []
pb_rows = []
validation_rows = []

for i, row in enumerate(
    meta.itertuples(index=False),
    start=1,
):

    donor = str(row.BrNum)

    path = (
        ROOT
        / "outputs/imputation_full"
        / CFG["root"]
        / donor
        / "ex5_1"
        / CFG["method"]
        / "oof_300"
        / "oof_predicted_300genes.h5ad"
    )

    if not path.is_file():
        raise FileNotFoundError(
            f"{MODEL} {donor}: {path}"
        )

    print(
        f"\n[{i:02d}/24] "
        f"{MODEL:8s} | {donor} | {row.Dx}",
        flush=True,
    )

    a = ad.read_h5ad(path)


    # -------------------------------------------------------------------------
    # OOF structure
    # -------------------------------------------------------------------------

    if a.n_vars != 300:
        raise RuntimeError(
            f"{MODEL} {donor}: "
            f"expected 300 genes, found {a.n_vars}"
        )

    if "count_scale" not in a.layers:
        raise RuntimeError(
            f"{MODEL} {donor}: count_scale missing."
        )

    required_obs = {
        "layer_annotation",
        "predictions_smooth",
    }

    missing_obs = required_obs - set(a.obs.columns)

    if missing_obs:
        raise RuntimeError(
            f"{MODEL} {donor}: "
            f"missing obs columns {sorted(missing_obs)}"
        )


    # -------------------------------------------------------------------------
    # Fold sanity check
    # -------------------------------------------------------------------------

    if "oof_fold" not in a.var.columns:
        raise RuntimeError(
            f"{MODEL} {donor}: oof_fold missing."
        )

    fold_counts = (
        a.var["oof_fold"]
        .value_counts()
        .sort_index()
        .to_dict()
    )

    if fold_counts != {
        1: 100,
        2: 100,
        3: 100,
    }:
        raise RuntimeError(
            f"{MODEL} {donor}: "
            f"invalid OOF folds {fold_counts}"
        )


    # -------------------------------------------------------------------------
    # Genes
    # -------------------------------------------------------------------------

    genes = a.var_names.astype(str).tolist()

    if gene_order is None:
        gene_order = genes

    elif genes != gene_order:
        raise RuntimeError(
            f"{MODEL} {donor}: gene order differs."
        )


    # -------------------------------------------------------------------------
    # Validate donor / Dx
    # -------------------------------------------------------------------------

    if "BrNum" in a.obs.columns:

        x = set(
            a.obs["BrNum"]
            .astype(str)
            .str.strip()
        )

        if x != {donor}:
            raise RuntimeError(
                f"{MODEL} {donor}: BrNum mismatch {x}"
            )


    if "Dx" in a.obs.columns:

        x = set(
            a.obs["Dx"]
            .astype(str)
            .str.strip()
        )

        # OOF files can contain placeholder Dx="unknown".
        # Canonical N24 patient metadata is the source of truth.
        known_dx = {
            v for v in x
            if v in {"NTC", "SCZ"}
        }

        if known_dx and known_dx != {str(row.Dx)}:
            raise RuntimeError(
                f"{MODEL} {donor}: "
                f"OOF Dx={known_dx}, metadata Dx={row.Dx}"
            )

        if not known_dx:
            print(
                f"  OOF Dx={x}; using canonical metadata Dx={row.Dx}"
            )


    # -------------------------------------------------------------------------
    # Prediction matrix ONLY
    # -------------------------------------------------------------------------

    X = a.layers["count_scale"]

    layers = (
        a.obs["layer_annotation"]
        .astype(str)
        .str.strip()
    )

    spds = (
        a.obs["predictions_smooth"]
        .astype(str)
        .str.strip()
    )


    valid = layers.isin(
        EXPECTED_LAYERS
    ).to_numpy()

    valid_pos = np.flatnonzero(valid)

    excluded = a.n_obs - len(valid_pos)

    if sparse.issparse(X):
        donor_valid_total = sum_matrix(
            X[valid_pos, :]
        )
    else:
        donor_valid_total = sum_matrix(
            np.asarray(X)[valid_pos, :]
        )


    layer_sum_check = np.zeros(
        300,
        dtype=np.float64,
    )


    # -------------------------------------------------------------------------
    # Seven spatial domains
    # -------------------------------------------------------------------------

    for layer in EXPECTED_LAYERS:

        expected_spd = LAYER_TO_SPD[layer]

        mask = layers.eq(layer).to_numpy()

        pos = np.flatnonzero(mask)

        n_cells = len(pos)

        if n_cells < MIN_NCELLS:
            raise RuntimeError(
                f"{MODEL} {donor} {layer}: "
                f"{n_cells} cells < {MIN_NCELLS}"
            )


        found_spd = set(
            spds.iloc[pos]
        )

        if found_spd != {expected_spd}:
            raise RuntimeError(
                f"{MODEL} {donor} {layer}: "
                f"expected {expected_spd}; "
                f"found {found_spd}"
            )


        if sparse.issparse(X):
            layer_matrix = X[pos, :]
        else:
            layer_matrix = np.asarray(X)[pos, :]


        summed = sum_matrix(
            layer_matrix
        )


        if not np.isfinite(summed).all():
            raise RuntimeError(
                f"{MODEL} {donor} {layer}: "
                "non-finite pseudobulk values."
            )


        if np.any(summed < 0):
            raise RuntimeError(
                f"{MODEL} {donor} {layer}: "
                "negative pseudobulk values."
            )


        pb_id = (
            f"{donor}__{expected_spd}__"
            f"{layer.replace('/', '-')}"
        )


        pb_arrays.append(
            summed
        )


        pb_rows.append(
            {
                "pseudobulk_id": pb_id,
                "BrNum": donor,
                "Dx": str(row.Dx),
                "Age": float(row.Age),
                "Sex": str(row.Sex),
                "slide_id": str(row.slide_id),
                "run_date": str(row.run_date),
                "predictions_smooth": expected_spd,
                "domain_annotations": layer,
                "n_cells": n_cells,
                "model": MODEL,
                "experiment": "ex5_1",
                "training": "NTC_only",
                "expression_source":
                    "layers[count_scale]",
                "library_size":
                    float(summed.sum()),
            }
        )


        layer_sum_check += summed


        print(
            f"    {expected_spd:5s} | "
            f"{layer:4s} | "
            f"cells={n_cells:6d} | "
            f"library={summed.sum():.3f}"
        )


    # -------------------------------------------------------------------------
    # Critical layer-sum validation
    # -------------------------------------------------------------------------

    matches = np.allclose(
        layer_sum_check,
        donor_valid_total,
        rtol=1e-8,
        atol=1e-5,
    )

    if not matches:

        max_diff = np.max(
            np.abs(
                layer_sum_check
                - donor_valid_total
            )
        )

        raise RuntimeError(
            f"{MODEL} {donor}: "
            f"7-layer sum mismatch; max_diff={max_diff}"
        )


    validation_rows.append(
        {
            "model": MODEL,
            "BrNum": donor,
            "Dx": str(row.Dx),
            "oof_cells": a.n_obs,
            "layer_valid_cells": len(valid_pos),
            "excluded_layer_cells": excluded,
            "n_genes": a.n_vars,
            "fold1_genes": fold_counts[1],
            "fold2_genes": fold_counts[2],
            "fold3_genes": fold_counts[3],
            "layer_sum_matches_donor": True,
        }
    )


    del a, X


# =============================================================================
# Assemble
# =============================================================================

section(f"ASSEMBLE — {MODEL}")

pb = np.vstack(
    pb_arrays
)

pb_meta = pd.DataFrame(
    pb_rows
)

validation = pd.DataFrame(
    validation_rows
)


print("Matrix  :", pb.shape)
print("Metadata:", pb_meta.shape)


if pb.shape != (168, 300):
    raise RuntimeError(
        f"{MODEL}: expected 168 × 300; "
        f"found {pb.shape}"
    )


if pb_meta["BrNum"].nunique() != 24:
    raise RuntimeError(
        f"{MODEL}: expected 24 donors."
    )


if pb_meta[
    "domain_annotations"
].nunique() != 7:
    raise RuntimeError(
        f"{MODEL}: expected 7 spatial domains."
    )


per_donor = (
    pb_meta
    .groupby("BrNum")
    ["domain_annotations"]
    .nunique()
)

if not (per_donor == 7).all():
    raise RuntimeError(
        f"{MODEL}: not every donor has 7 domains."
    )


dx = (
    pb_meta[
        ["BrNum", "Dx"]
    ]
    .drop_duplicates()
    ["Dx"]
    .value_counts()
    .to_dict()
)

if dx != {
    "NTC": 12,
    "SCZ": 12,
}:
    raise RuntimeError(
        f"{MODEL}: unexpected Dx counts {dx}"
    )


# =============================================================================
# Save
# =============================================================================

section(f"SAVE — {MODEL}")

pb_df = pd.DataFrame(
    pb,
    index=pb_meta["pseudobulk_id"],
    columns=gene_order,
)


PB_FILE = (
    OUT
    / f"{MODEL}_NTCtrained_OOF_donor_domain_pseudobulk_300genes.csv.gz"
)

META_FILE = (
    OUT
    / f"{MODEL}_NTCtrained_OOF_donor_domain_pseudobulk_metadata.csv"
)

VALID_FILE = (
    OUT
    / f"{MODEL}_NTCtrained_OOF_donor_domain_validation.csv"
)

GENE_FILE = (
    OUT
    / f"{MODEL}_gene_order_300.txt"
)


pb_df.to_csv(
    PB_FILE
)

pb_meta.to_csv(
    META_FILE,
    index=False,
)

validation.to_csv(
    VALID_FILE,
    index=False,
)

GENE_FILE.write_text(
    "\n".join(gene_order) + "\n"
)


# =============================================================================
# Summary
# =============================================================================

section(f"FINAL SUMMARY — {MODEL}")

print("Model              :", MODEL)
print("Donors             :", pb_meta["BrNum"].nunique())
print("NTC donors         :", dx["NTC"])
print("SCZ donors         :", dx["SCZ"])
print("Spatial domains    :", pb_meta["domain_annotations"].nunique())
print("Pseudobulk samples :", len(pb_meta))
print("Genes              :", len(gene_order))
print("Minimum cells / PB :", pb_meta["n_cells"].min())
print("Maximum cells / PB :", pb_meta["n_cells"].max())

print(
    "All donor layer sums match:",
    validation["layer_sum_matches_donor"].all(),
)

print("\nOutput:")
print(PB_FILE)
print(META_FILE)
print(VALID_FILE)

print(
    f"\nFINAL STATUS: PASS — "
    f"{MODEL} NTC-trained OOF "
    "24-donor × 7-domain pseudobulk created."
)
