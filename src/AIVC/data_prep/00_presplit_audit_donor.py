#!/usr/bin/env python3

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

ANNOTATED = (
    ROOT
    / "data/processed/xenium/"
      "xenium_N24_layer_celltype_annotated.h5ad"
)

SNRNA_REFERENCE = (
    ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
      "huuki_snrna_reference_full_allgenes.h5ad"
)

OUT = (
    ROOT
    / "outputs/AIVC_beta/presplit_audit/per_donor"
)

OUT.mkdir(parents=True, exist_ok=True)


EXPECTED_CELL_TYPES = [
    "Oligo",
    "Endo",
    "Ast",
    "L2/3 Ex",
    "Mic",
    "In: SST, PVALB",
    "L6 Ex",
    "Ambig/In/Endo",
    "L4/5 Ex",
    "In: VIP, LAMP5",
    "Ambig/Oligo",
    "L5 Ex",
]

EXPECTED_LAYERS = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WMtz",
    "WM",
]


def clean_series(x):
    return x.astype(str).str.strip()


def unique_nonempty(x):
    s = clean_series(x)

    s = s.replace(
        {
            "nan": np.nan,
            "None": np.nan,
            "": np.nan,
        }
    )

    return sorted(
        s.dropna().unique().tolist()
    )


def collapse(x):
    v = unique_nonempty(x)

    if len(v) == 0:
        return np.nan

    return "|".join(v)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--envi",
    required=True,
)

parser.add_argument(
    "--donor",
    required=True,
)

args = parser.parse_args()

path = Path(args.envi)
donor = args.donor


print("=" * 100, flush=True)
print("AIVC PRE-SPLIT AUDIT — ONE DONOR", flush=True)
print("=" * 100, flush=True)

print("Donor:", donor, flush=True)
print("ENVI :", path, flush=True)


# ======================================================================================
# ENVI
# ======================================================================================

print("\nLoading ENVI...", flush=True)

e = ad.read_h5ad(
    path,
    backed="r",
)

print(
    "ENVI shape:",
    e.shape,
    flush=True,
)

e_ids = pd.Index(
    e.obs_names.astype(str)
)

duplicate_envi_ids = int(
    e_ids.duplicated().sum()
)


# ======================================================================================
# DONOR / DX
# ======================================================================================

br_values = (
    unique_nonempty(e.obs["BrNum"])
    if "BrNum" in e.obs.columns
    else []
)

dx_values = (
    unique_nonempty(e.obs["Dx"])
    if "Dx" in e.obs.columns
    else []
)

donor_actual = (
    br_values[0]
    if len(br_values) == 1
    else "|".join(br_values)
)

dx = (
    dx_values[0]
    if len(dx_values) == 1
    else "|".join(dx_values)
)

print(
    "File donor:",
    donor_actual,
    " Dx:",
    dx,
    flush=True,
)


# ======================================================================================
# ANNOTATED XENIUM
# ======================================================================================

print("\nLoading annotated Xenium metadata...", flush=True)

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
)

ann_br = clean_series(
    ann.obs["BrNum"]
)

mask = (
    ann_br.to_numpy()
    == donor
)

ann_d = ann.obs.loc[
    mask
].copy()

ann_d.index = ann_d.index.astype(str)

a_ids = pd.Index(
    ann_d.index
)

duplicate_annotation_ids = int(
    a_ids.duplicated().sum()
)

overlap = e_ids.intersection(
    a_ids
)

join_n = int(
    len(overlap)
)

join_pct = (
    100.0
    * join_n
    / len(e_ids)
    if len(e_ids)
    else np.nan
)

envi_missing_annotation = int(
    len(
        e_ids.difference(
            a_ids
        )
    )
)

annotation_missing_envi = int(
    len(
        a_ids.difference(
            e_ids
        )
    )
)

print(
    f"Annotation join: {join_n}/{len(e_ids)} "
    f"({join_pct:.4f}%)",
    flush=True,
)


# ======================================================================================
# COORDINATES
# ======================================================================================

has_spatial = (
    "spatial"
    in e.obsm.keys()
)

spatial_finite = False
spatial_shape = ""

if has_spatial:

    coords = np.asarray(
        e.obsm["spatial"]
    )

    spatial_shape = str(
        tuple(coords.shape)
    )

    spatial_finite = bool(
        np.isfinite(coords).all()
    )

coordinate_max_abs_diff = np.nan

if (
    join_n == len(e_ids)
    and "x_centroid" in e.obs.columns
    and "y_centroid" in e.obs.columns
    and "x_centroid" in ann_d.columns
    and "y_centroid" in ann_d.columns
):

    aligned = ann_d.loc[
        e_ids
    ]

    ex = pd.to_numeric(
        e.obs["x_centroid"],
        errors="coerce",
    ).to_numpy()

    ey = pd.to_numeric(
        e.obs["y_centroid"],
        errors="coerce",
    ).to_numpy()

    ax = pd.to_numeric(
        aligned["x_centroid"],
        errors="coerce",
    ).to_numpy()

    ay = pd.to_numeric(
        aligned["y_centroid"],
        errors="coerce",
    ).to_numpy()

    coordinate_max_abs_diff = float(
        np.nanmax(
            np.abs(
                np.concatenate(
                    [
                        ex - ax,
                        ey - ay,
                    ]
                )
            )
        )
    )


# ======================================================================================
# CELL TYPES
# ======================================================================================

if "cell_type_annotation" not in ann_d.columns:
    raise RuntimeError(
        "cell_type_annotation missing"
    )

ct = clean_series(
    ann_d[
        "cell_type_annotation"
    ]
)

ct_counts = (
    ct.value_counts()
)

celltype_row = {
    "BrNum": donor,
    "Dx": dx,
}

for label in EXPECTED_CELL_TYPES:

    celltype_row[label] = int(
        ct_counts.get(
            label,
            0,
        )
    )

missing_celltypes = sorted(
    set(EXPECTED_CELL_TYPES)
    - set(ct.unique())
)

unexpected_celltypes = sorted(
    set(ct.unique())
    - set(EXPECTED_CELL_TYPES)
)


# ======================================================================================
# LAYERS
# ======================================================================================

layer_row = {
    "BrNum": donor,
    "Dx": dx,
}

missing_layers = []
unexpected_layers = []

if "layer_annotation" in ann_d.columns:

    la = clean_series(
        ann_d[
            "layer_annotation"
        ]
    )

    la_counts = (
        la.value_counts()
    )

    for layer in EXPECTED_LAYERS:

        layer_row[layer] = int(
            la_counts.get(
                layer,
                0,
            )
        )

    missing_layers = sorted(
        set(EXPECTED_LAYERS)
        - set(la.unique())
    )

    unexpected_layers = sorted(
        set(la.unique())
        - set(EXPECTED_LAYERS)
    )


# ======================================================================================
# GENE CONSISTENCY INFORMATION
# ======================================================================================

genes = pd.Index(
    e.var_names.astype(str)
)

gene_hash = str(
    hash(
        tuple(
            genes.tolist()
        )
    )
)


# ======================================================================================
# GENE PROVENANCE
# ======================================================================================

measured_genes = np.nan
imputed_genes = np.nan

if "measured_in_xenium" in e.var.columns:

    measured = (
        e.var[
            "measured_in_xenium"
        ]
        .astype(bool)
        .to_numpy()
    )

    measured_genes = int(
        measured.sum()
    )

    imputed_genes = int(
        len(measured)
        - measured.sum()
    )


# ======================================================================================
# EXPRESSION QC — SAMPLE ONLY
# ======================================================================================

expression_rows = []

n_cells = e.n_obs
n_genes = e.n_vars

row_idx = np.unique(
    np.linspace(
        0,
        n_cells - 1,
        min(
            128,
            n_cells,
        ),
        dtype=int,
    )
)

col_idx = np.unique(
    np.linspace(
        0,
        n_genes - 1,
        min(
            1024,
            n_genes,
        ),
        dtype=int,
    )
)

for layer in [
    "count_scale",
    "log1p",
    "native",
]:

    if layer not in e.layers.keys():
        continue

    try:

        x = e.layers[layer][
            row_idx,
            :
        ]

        x = np.asarray(x)[
            :,
            col_idx
        ]

        expression_rows.append(
            {
                "BrNum": donor,
                "Dx": dx,
                "layer": layer,
                "sample_cells":
                    int(x.shape[0]),
                "sample_genes":
                    int(x.shape[1]),
                "all_finite":
                    bool(
                        np.isfinite(
                            x
                        ).all()
                    ),
                "nan_n":
                    int(
                        np.isnan(
                            x
                        ).sum()
                    ),
                "inf_n":
                    int(
                        np.isinf(
                            x
                        ).sum()
                    ),
                "min":
                    float(
                        np.nanmin(
                            x
                        )
                    ),
                "max":
                    float(
                        np.nanmax(
                            x
                        )
                    ),
                "mean":
                    float(
                        np.nanmean(
                            x
                        )
                    ),
                "negative_n":
                    int(
                        (x < 0).sum()
                    ),
            }
        )

    except Exception as exc:

        expression_rows.append(
            {
                "BrNum": donor,
                "Dx": dx,
                "layer": layer,
                "all_finite": False,
                "error": repr(exc),
            }
        )


# ======================================================================================
# DONOR METADATA
# ======================================================================================

metadata = {}

for col in [
    "Age",
    "Sex",
    "slide_id",
    "run_date",
    "sample_id",
    "CaptureArea",
    "PNN",
    "tear",
]:

    if col in ann_d.columns:

        metadata[col] = collapse(
            ann_d[col]
        )

    elif col in e.obs.columns:

        metadata[col] = collapse(
            e.obs[col]
        )

    else:

        metadata[col] = np.nan


# ======================================================================================
# snRNA REFERENCE OVERLAP
# ======================================================================================

print(
    "\nChecking snRNA donor overlap...",
    flush=True,
)

ref = ad.read_h5ad(
    SNRNA_REFERENCE,
    backed="r",
)

ref_donors = set(
    clean_series(
        ref.obs[
            "BrNum"
        ]
    )
)

reference_overlap = (
    donor
    in ref_donors
)


# ======================================================================================
# MAIN DONOR ROW
# ======================================================================================

ambiguous_n = int(
    ct.isin(
        [
            "Ambig/In/Endo",
            "Ambig/Oligo",
        ]
    ).sum()
)

donor_row = {
    "BrNum": donor,
    "Dx": dx,
    "file": str(path),

    "n_cells": int(
        e.n_obs
    ),

    "n_genes": int(
        e.n_vars
    ),

    "donor_in_file":
        donor_actual,

    "donor_matches_path":
        donor_actual == donor,

    "n_dx_values":
        len(dx_values),

    "duplicate_envi_ids":
        duplicate_envi_ids,

    "annotated_cells":
        int(
            len(ann_d)
        ),

    "duplicate_annotation_ids":
        duplicate_annotation_ids,

    "join_n":
        join_n,

    "join_pct":
        join_pct,

    "envi_missing_annotation":
        envi_missing_annotation,

    "annotation_missing_envi":
        annotation_missing_envi,

    "has_spatial":
        has_spatial,

    "spatial_shape":
        spatial_shape,

    "spatial_finite":
        spatial_finite,

    "coordinate_max_abs_diff":
        coordinate_max_abs_diff,

    "n_cell_types":
        int(
            ct.nunique()
        ),

    "missing_cell_types":
        ";".join(
            missing_celltypes
        ),

    "unexpected_cell_types":
        ";".join(
            unexpected_celltypes
        ),

    "ambiguous_cells":
        ambiguous_n,

    "ambiguous_pct":
        100.0
        * ambiguous_n
        / len(ct),

    "missing_layers":
        ";".join(
            missing_layers
        ),

    "unexpected_layers":
        ";".join(
            unexpected_layers
        ),

    "gene_names_unique":
        bool(
            genes.nunique()
            == len(genes)
        ),

    "gene_hash":
        gene_hash,

    "measured_genes":
        measured_genes,

    "imputed_genes":
        imputed_genes,

    "reference_donor_overlap":
        reference_overlap,
}

donor_row.update(
    metadata
)


# ======================================================================================
# SAVE
# ======================================================================================

pd.DataFrame(
    [donor_row]
).to_csv(
    OUT
    / f"{donor}_donor_audit.csv",
    index=False,
)

pd.DataFrame(
    [celltype_row]
).to_csv(
    OUT
    / f"{donor}_celltype_counts.csv",
    index=False,
)

pd.DataFrame(
    [layer_row]
).to_csv(
    OUT
    / f"{donor}_layer_counts.csv",
    index=False,
)

pd.DataFrame(
    expression_rows
).to_csv(
    OUT
    / f"{donor}_expression_qc.csv",
    index=False,
)


with open(
    OUT
    / f"{donor}_status.json",
    "w",
) as f:

    json.dump(
        {
            "BrNum": donor,
            "Dx": dx,
            "completed": True,
            "join_pct":
                join_pct,
            "reference_overlap":
                reference_overlap,
        },
        f,
        indent=2,
    )


e.file.close()
ann.file.close()
ref.file.close()

print("\nDONE:", donor, flush=True)

