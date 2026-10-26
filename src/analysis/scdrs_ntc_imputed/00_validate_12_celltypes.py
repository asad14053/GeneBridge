#!/usr/bin/env python3

from pathlib import Path
import anndata as ad
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

ENVI_ROOT = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/envi"
)

OUTDIR = (
    ROOT
    / "results/analysis/scdrs_ntc_imputed"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)

files = sorted(
    ENVI_ROOT.glob(
        "Br*/spatial_data_xenium_*_ENVI_full_transcriptome.h5ad"
    )
)

if len(files) != 12:
    raise RuntimeError(
        f"Expected exactly 12 NTC ENVI files, found {len(files)}"
    )

print("=" * 100)
print("12-NTC ENVI CELL-TYPE VALIDATION")
print("=" * 100)
print("Files:", len(files))


# ============================================================
# Collect candidate columns
# ============================================================

rows = []
label_sets = {}

for f in files:

    br = f.parent.name

    a = ad.read_h5ad(
        f,
        backed="r"
    )

    print()
    print("=" * 100)
    print(br)
    print("=" * 100)
    print("shape:", a.shape)

    for col in a.obs.columns:

        try:
            s = (
                a.obs[col]
                .dropna()
                .astype(str)
            )

            n = s.nunique()

        except Exception:
            continue

        if n == 12:

            labels = tuple(
                sorted(
                    s.unique().tolist()
                )
            )

            print()
            print(
                f"CANDIDATE COLUMN: {col}"
            )

            print(
                a.obs[col]
                .astype(str)
                .value_counts()
                .to_string()
            )

            rows.append({
                "BrNum": br,
                "column": col,
                "n_cells": a.n_obs,
                "n_categories": n,
                "labels": " | ".join(labels),
            })

            label_sets[
                (br, col)
            ] = labels

    a.file.close()


cand = pd.DataFrame(rows)

OUT_CAND = (
    OUTDIR
    / "00_celltype_12category_candidates.csv"
)

cand.to_csv(
    OUT_CAND,
    index=False
)


# ============================================================
# Which columns work in ALL 12 donors?
# ============================================================

if cand.empty:
    raise RuntimeError(
        "No obs column with exactly 12 categories was found."
    )

counts = (
    cand.groupby("column")["BrNum"]
    .nunique()
    .sort_values(
        ascending=False
    )
)

print()
print("=" * 100)
print("NUMBER OF DONORS WITH EXACTLY 12 CATEGORIES")
print("=" * 100)

print(
    counts.to_string()
)


common = counts[
    counts == 12
].index.tolist()


print()
print("=" * 100)
print("COLUMNS WITH 12 CATEGORIES IN ALL 12 DONORS")
print("=" * 100)

if not common:
    print("NONE")

else:

    for col in common:

        print()
        print("COLUMN:", col)

        sub = cand[
            cand["column"] == col
        ]

        unique_label_sets = (
            sub["labels"]
            .drop_duplicates()
            .tolist()
        )

        print(
            "different label sets across donors:",
            len(unique_label_sets)
        )

        for x in unique_label_sets:
            print(" ", x)


# ============================================================
# Strict candidates:
# exactly 12 categories AND identical labels in all donors
# ============================================================

strict = []

for col in common:

    sub = cand[
        cand["column"] == col
    ]

    if (
        sub["labels"].nunique()
        ==
        1
    ):

        strict.append(col)


print()
print("=" * 100)
print("STRICT CONSISTENT 12-CELL-TYPE COLUMNS")
print("=" * 100)

if strict:

    for col in strict:
        print(col)

else:
    print("NONE")


# ============================================================
# Save donor manifest
# ============================================================

manifest = pd.DataFrame({
    "BrNum": [
        f.parent.name
        for f in files
    ],
    "Dx": "NTC",
    "envi_file": [
        str(f)
        for f in files
    ],
})

OUT_MANIFEST = (
    OUTDIR
    / "00_12NTC_manifest.csv"
)

manifest.to_csv(
    OUT_MANIFEST,
    index=False
)


# ============================================================
# Summary
# ============================================================

summary = pd.DataFrame([
    {
        "n_donors": len(files),
        "expected_genes_per_ENVI": 34987,
        "measured_xenium_genes": 300,
        "target_imputed_only_genes": 34687,
        "common_12category_columns":
            ";".join(common),
        "strict_consistent_columns":
            ";".join(strict),
    }
])

OUT_SUMMARY = (
    OUTDIR
    / "00_validation_summary.csv"
)

summary.to_csv(
    OUT_SUMMARY,
    index=False
)


print()
print("=" * 100)
print("OUTPUT")
print("=" * 100)

print("Candidates :", OUT_CAND)
print("Manifest   :", OUT_MANIFEST)
print("Summary    :", OUT_SUMMARY)
