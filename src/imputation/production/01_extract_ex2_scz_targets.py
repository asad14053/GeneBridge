#!/usr/bin/env python3

from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

INPUT_H5AD = (
    ROOT
    / "data/processed/xenium/xenium_N24_imputation_ready.h5ad"
)

MANIFEST = (
    ROOT
    / "data/metadata/imputation_ex2_scz_samples.csv"
)

OUT_DIR = (
    ROOT
    / "data/processed/imputation_full/ex2_scz/targets"
)

SUMMARY_CSV = (
    ROOT
    / "data/metadata/imputation_ex2_scz_target_summary.csv"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


def matrix_stats(x):
    if sparse.issparse(x):
        values = x.data
        total = x.shape[0] * x.shape[1]
        nnz = x.nnz
        zero_fraction = 1.0 - nnz / total

        if values.size == 0:
            return {
                "min": 0.0,
                "max": 0.0,
                "mean": 0.0,
                "zero_fraction": 1.0,
            }

        matrix_sum = float(x.sum())

        return {
            "min": min(0.0, float(values.min())),
            "max": float(values.max()),
            "mean": matrix_sum / total,
            "zero_fraction": zero_fraction,
        }

    arr = np.asarray(x)

    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "zero_fraction": float(np.mean(arr == 0)),
    }


def is_integer_like(x, atol=1e-6):
    if sparse.issparse(x):
        values = x.data
    else:
        values = np.asarray(x).ravel()

    if len(values) == 0:
        return True

    return bool(
        np.all(
            np.abs(values - np.rint(values)) <= atol
        )
    )


print("=" * 100)
print("EXPERIMENT 2 — EXTRACT 12 SCZ XENIUM TARGETS")
print("=" * 100)

print("Input:")
print(INPUT_H5AD)

print("\nManifest:")
print(MANIFEST)

manifest = pd.read_csv(MANIFEST)

expected = (
    manifest.loc[
        manifest["Dx"].astype(str).str.upper().eq("SCZ"),
        "patient_id",
    ]
    .astype(str)
    .tolist()
)

if len(expected) != 12:
    raise RuntimeError(
        f"Expected 12 SCZ donors in manifest; found {len(expected)}"
    )

print("\nSCZ donors:")
print(expected)

print("\nReading N24 Xenium object...")
adata = ad.read_h5ad(INPUT_H5AD)

print("N24 shape:", adata.shape)
print("Layers:", list(adata.layers.keys()))
print("obsm:", list(adata.obsm.keys()))

# ------------------------------------------------------------------
# Required checks
# ------------------------------------------------------------------

for required in ["BrNum", "Dx"]:
    if required not in adata.obs.columns:
        raise RuntimeError(
            f"Required obs column missing: {required}"
        )

if "counts" not in adata.layers:
    raise RuntimeError(
        "Required layer['counts'] is missing."
    )

if "spatial" not in adata.obsm:
    raise RuntimeError(
        "Required obsm['spatial'] is missing."
    )

if adata.n_vars != 300:
    raise RuntimeError(
        f"Expected 300 Xenium genes; found {adata.n_vars}"
    )

# ------------------------------------------------------------------
# Verify diagnosis labels against manifest
# ------------------------------------------------------------------

for donor in expected:

    mask = (
        adata.obs["BrNum"]
        .astype(str)
        .eq(donor)
    )

    if mask.sum() == 0:
        raise RuntimeError(
            f"{donor}: no cells found in N24 object."
        )

    dx = (
        adata.obs.loc[mask, "Dx"]
        .astype(str)
        .str.upper()
        .unique()
        .tolist()
    )

    if dx != ["SCZ"]:
        raise RuntimeError(
            f"{donor}: unexpected Dx values: {dx}"
        )

# ------------------------------------------------------------------
# Extract each donor
# ------------------------------------------------------------------

summary_rows = []

for i, donor in enumerate(expected, start=1):

    print("\n" + "-" * 100)
    print(f"[{i:02d}/12] {donor}")
    print("-" * 100)

    mask = (
        adata.obs["BrNum"]
        .astype(str)
        .eq(donor)
        .to_numpy()
    )

    target = adata[mask].copy()

    # --------------------------------------------------------------
    # CRITICAL:
    # Production ENVI target uses raw/count-scale Xenium expression.
    # --------------------------------------------------------------

    counts = target.layers["counts"].copy()

    target.X = counts.copy()
    target.layers["counts"] = counts.copy()

    # Keep normalized layer only as auxiliary information.
    # ENVI production input will explicitly use count scale.
    if "log1p_norm" in target.layers:
        pass

    # Store provenance
    target.uns["production_imputation"] = {
        "experiment": "ex2_scz",
        "diagnosis": "SCZ",
        "reference": "all_10_Huuki_snRNA",
        "target_donor": donor,
        "gene_holdout": False,
        "n_xenium_input_genes": int(target.n_vars),
        "X_scale": "raw_counts",
        "counts_layer": "counts",
        "spatial_key": "spatial",
    }

    stats = matrix_stats(target.X)
    integer_like = is_integer_like(target.X)

    coords = np.asarray(target.obsm["spatial"])

    finite_spatial = bool(np.isfinite(coords).all())

    print("Shape             :", target.shape)
    print("X source          : layers['counts']")
    print("X scale           : raw/count scale")
    print("X min             :", stats["min"])
    print("X max             :", stats["max"])
    print("X mean            :", stats["mean"])
    print("X zero fraction   :", stats["zero_fraction"])
    print("Integer-like      :", integer_like)
    print("Spatial shape     :", coords.shape)
    print("Spatial finite    :", finite_spatial)

    if stats["min"] < 0:
        raise RuntimeError(
            f"{donor}: negative count values detected."
        )

    if not integer_like:
        raise RuntimeError(
            f"{donor}: count layer is not integer-like."
        )

    if coords.shape != (target.n_obs, 2):
        raise RuntimeError(
            f"{donor}: unexpected spatial shape {coords.shape}"
        )

    if not finite_spatial:
        raise RuntimeError(
            f"{donor}: non-finite spatial coordinates."
        )

    out_file = (
        OUT_DIR
        / f"spatial_data_xenium_{donor}_ex2_scz.h5ad"
    )

    target.write_h5ad(
        out_file,
        compression="gzip",
    )

    print("Saved:")
    print(out_file)

    summary_rows.append(
        {
            "experiment": "ex2_scz",
            "patient_id": donor,
            "Dx": "SCZ",
            "n_cells": target.n_obs,
            "n_xenium_genes": target.n_vars,
            "x_scale": "raw_counts",
            "counts_integer_like": integer_like,
            "spatial_available": True,
            "spatial_finite": finite_spatial,
            "x_min": stats["min"],
            "x_max": stats["max"],
            "x_mean": stats["mean"],
            "x_zero_fraction": stats["zero_fraction"],
            "target_h5ad": str(out_file),
        }
    )

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    SUMMARY_CSV,
    index=False,
)

print("\n" + "=" * 100)
print("EXPERIMENT 2 TARGET SUMMARY")
print("=" * 100)

print(
    summary[
        [
            "patient_id",
            "Dx",
            "n_cells",
            "n_xenium_genes",
            "x_scale",
            "counts_integer_like",
            "spatial_available",
        ]
    ].to_string(index=False)
)

print("\nTotal SCZ donors:", len(summary))
print("Total SCZ cells :", summary["n_cells"].sum())

if len(summary) != 12:
    raise RuntimeError(
        f"Expected 12 output targets; produced {len(summary)}"
    )

print("\nSummary CSV:")
print(SUMMARY_CSV)

print("\nTarget directory:")
print(OUT_DIR)

print("\n" + "=" * 100)
print("SUCCESS: ALL 12 SCZ TARGETS EXTRACTED")
print("=" * 100)
