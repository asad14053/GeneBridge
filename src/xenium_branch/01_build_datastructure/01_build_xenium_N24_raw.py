#!/usr/bin/env python3

"""
01_build_xenium_N24_raw_h5ad.py

Goal
----
Load all 24 raw Xenium samples and create:

    data/processed/xenium/xenium_N24_raw.h5ad

This is the AnnData/Python equivalent of the PI's:

    processed-data/01_build_spe/raw_spe_N24.RDS

Input assumed
-------------
Only raw Xenium folders:

    data/raw/xenium/Br2039/
    data/raw/xenium/Br2719/
    ...

Each folder should contain:

    cell_feature_matrix.h5
    cells.csv.gz OR cells.parquet.gz

Optional but recorded if present:

    cell_boundaries.csv.gz
    nucleus_boundaries.csv.gz

Outputs
-------
Main h5ad:

    data/processed/xenium/xenium_N24_raw.h5ad

Reusable CSVs:

    outputs/xenium_branch/01_build_datastructure/tables/sample_paths.csv
    outputs/xenium_branch/01_build_datastructure/tables/xenium_N24_sample_summary.csv
    outputs/xenium_branch/01_build_datastructure/tables/xenium_N24_cell_metadata.csv.gz
    outputs/xenium_branch/01_build_datastructure/tables/xenium_N24_feature_metadata.csv
    outputs/xenium_branch/01_build_datastructure/tables/xenium_N24_build_status.csv

Plots:

    outputs/xenium_branch/01_build_datastructure/figures/xenium_N24_cells_per_sample.png
    outputs/xenium_branch/01_build_datastructure/figures/xenium_N24_features_per_sample.png
    outputs/xenium_branch/01_build_datastructure/figures/xenium_N24_raw_spatial_locations.png
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import shutil
import tempfile
import traceback
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[3]

RAW_XENIUM_DIR = PROJECT_ROOT / "data" / "raw" / "xenium"

OUT_DIR = PROJECT_ROOT / "outputs" / "xenium_branch" / "01_build_datastructure"
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "xenium"
OUT_H5AD = PROCESSED_DIR / "xenium_N24_raw.h5ad"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def section(title: str) -> None:
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def extract_brnum(path: Path) -> str | None:
    """
    Extract Br donor ID from folder/file name.

    Works for:
        Br2039
        GSM9223468_Br2039
        GSM9223468_Br2039-cell_feature_matrix.h5
    """
    text = str(path)
    match = re.search(r"(Br\d+)", text)
    if match:
        return match.group(1)
    return None


def find_first(folder: Path, patterns: list[str]) -> str:
    """
    Find first file matching any pattern.
    """
    for pattern in patterns:
        hits = sorted(folder.glob(pattern))
        if hits:
            return str(hits[0])
    return ""


def scan_xenium_folders(raw_root: Path) -> pd.DataFrame:
    """
    Scan data/raw/xenium/Br*/ folders and create sample_paths.csv.
    """
    if not raw_root.exists():
        raise FileNotFoundError(f"Raw Xenium folder does not exist: {raw_root}")

    folders = sorted([p for p in raw_root.iterdir() if p.is_dir()])

    records = []

    for folder in folders:
        brnum = extract_brnum(folder)

        if brnum is None:
            print(f"Skipping folder without BrNum: {folder}")
            continue

        matrix_path = find_first(
            folder,
            [
                "cell_feature_matrix.h5",
                "*cell_feature_matrix.h5",
            ],
        )

        cells_path = find_first(
            folder,
            [
                "cells.csv.gz",
                "*cells.csv.gz",
                "cells.parquet.gz",
                "*cells.parquet.gz",
                "cells.parquet",
                "*cells.parquet",
                "cells.csv",
                "*cells.csv",
            ],
        )

        cell_boundaries_path = find_first(
            folder,
            [
                "cell_boundaries.csv.gz",
                "*cell_boundaries.csv.gz",
                "cell_boundaries.csv",
                "*cell_boundaries.csv",
            ],
        )

        nucleus_boundaries_path = find_first(
            folder,
            [
                "nucleus_boundaries.csv.gz",
                "*nucleus_boundaries.csv.gz",
                "nucleus_boundaries.csv",
                "*nucleus_boundaries.csv",
            ],
        )

        records.append(
            {
                "BrNum": brnum,
                "path": str(folder),
                "matrix_path": matrix_path,
                "cells_path": cells_path,
                "cell_boundaries_path": cell_boundaries_path,
                "nucleus_boundaries_path": nucleus_boundaries_path,
                "has_matrix": bool(matrix_path),
                "has_cells": bool(cells_path),
                "has_cell_boundaries": bool(cell_boundaries_path),
                "has_nucleus_boundaries": bool(nucleus_boundaries_path),
            }
        )

    sample_paths = pd.DataFrame(records)

    if sample_paths.empty:
        raise RuntimeError(f"No Xenium Br folders found inside: {raw_root}")

    sample_paths = sample_paths.sort_values("BrNum").reset_index(drop=True)

    out_csv = TABLE_DIR / "sample_paths.csv"
    sample_paths.to_csv(out_csv, index=False)

    print(f"Saved sample paths: {out_csv}")
    print(sample_paths[["BrNum", "has_matrix", "has_cells"]].to_string(index=False))

    return sample_paths


def read_cells_metadata(path: str | Path) -> pd.DataFrame:
    """
    Read Xenium cells metadata.

    Supports:
        cells.csv.gz
        cells.csv
        cells.parquet
        cells.parquet.gz

    Some GEO downloads can be gzip-wrapped parquet.
    """
    path = Path(path)
    name = path.name.lower()

    if name.endswith(".csv.gz") or name.endswith(".csv"):
        return pd.read_csv(path)

    if name.endswith(".parquet") or name.endswith(".parquet.gz"):
        try:
            return pd.read_parquet(path)
        except Exception:
            pass

        with open(path, "rb") as f:
            first_bytes = f.read(2)

        if first_bytes == b"\x1f\x8b":
            temp_path = None
            try:
                with gzip.open(path, "rb") as gz:
                    with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
                        shutil.copyfileobj(gz, tmp)
                        temp_path = Path(tmp.name)

                return pd.read_parquet(temp_path)

            finally:
                if temp_path is not None and temp_path.exists():
                    temp_path.unlink()

    raise ValueError(f"Unsupported or unreadable cells metadata file: {path}")


def clean_var_metadata(adata: ad.AnnData) -> ad.AnnData:
    """
    Make AnnData var similar to SpatialExperiment rowData.

    PI logic:
        rownames(spe) <- rowData(spe)$Symbol

    AnnData logic:
        adata.var_names = gene symbols
        adata.var contains ID, Symbol, Type
    """
    adata.var["Symbol"] = adata.var_names.astype(str)

    if "gene_ids" in adata.var.columns:
        adata.var["ID"] = adata.var["gene_ids"].astype(str)
    elif "ID" not in adata.var.columns:
        adata.var["ID"] = adata.var_names.astype(str)

    if "feature_types" in adata.var.columns:
        adata.var["Type"] = adata.var["feature_types"].astype(str)
    elif "Type" not in adata.var.columns:
        adata.var["Type"] = "Gene Expression"

    adata.var_names = adata.var["Symbol"].astype(str)
    adata.var_names_make_unique()

    adata.var["is_gene_expression"] = adata.var["Type"].astype(str).eq("Gene Expression")

    return adata


def attach_cells_metadata(
    adata: ad.AnnData,
    cells: pd.DataFrame,
    brnum: str,
) -> ad.AnnData:
    """
    Attach cells.csv.gz / cells.parquet.gz metadata to AnnData obs.

    PI logic:
        colData(sce) <- cbind(colData(sce), cell_info)
        colnames(spe) <- paste(brnum, rownames(cell_info), sep = "_")

    Python logic:
        obs gets cell metadata
        obs_names become BrNum_1, BrNum_2, ...
        original Xenium cell_id is preserved
    """
    required = ["cell_id", "x_centroid", "y_centroid"]

    missing = [c for c in required if c not in cells.columns]
    if missing:
        raise ValueError(f"cells metadata is missing required columns: {missing}")

    cells = cells.copy()
    cells["cell_id"] = cells["cell_id"].astype(str)

    original_h5_cell_ids = pd.Index(adata.obs_names.astype(str))
    cells_cell_ids = pd.Index(cells["cell_id"].astype(str))

    common_ids = original_h5_cell_ids.intersection(cells_cell_ids)

    print(f"Cells in H5 matrix:      {adata.n_obs}")
    print(f"Cells in cells metadata: {cells.shape[0]}")
    print(f"Matched cell IDs:        {len(common_ids)}")

    if len(common_ids) > 0:
        # Best case: align by real Xenium cell_id.
        ordered_ids = [cid for cid in cells["cell_id"].tolist() if cid in set(common_ids)]

        cells = cells.set_index("cell_id", drop=False).loc[ordered_ids].copy()
        adata = adata[ordered_ids, :].copy()

    else:
        # Fallback: PI R code cbinds by order. Use order only if lengths match.
        if adata.n_obs != cells.shape[0]:
            raise ValueError(
                "Could not match H5 cell IDs to cells metadata, and row counts differ. "
                f"H5 cells={adata.n_obs}, cells metadata={cells.shape[0]}"
            )

        print("WARNING: No direct cell_id match. Attaching cells metadata by row order.")
        cells.index = adata.obs_names

    # Create PI-style donor-specific cell IDs.
    cells["r_cell_index"] = np.arange(1, cells.shape[0] + 1)
    cells["donor_cell_id"] = brnum + "_" + cells["r_cell_index"].astype(str)

    # Attach all cell metadata.
    adata.obs = adata.obs.join(cells, how="left")

    adata.obs["original_cell_id"] = adata.obs["cell_id"].astype(str)
    adata.obs["BrNum"] = brnum
    adata.obs["sample_id"] = brnum
    adata.obs["technology"] = "Xenium"

    # These columns existed in PI object from Excel.
    # Raw Xenium alone cannot provide them, so we keep placeholders.
    adata.obs["Dx"] = "unknown"
    adata.obs["CaptureArea"] = "unknown"
    adata.obs["PNN"] = "unknown"
    adata.obs["tear"] = "unknown"
    adata.obs["Age"] = np.nan
    adata.obs["Sex"] = "unknown"

    # AnnData equivalent of SpatialExperiment spatial coordinates.
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    # Set final unique cell names.
    adata.obs_names = adata.obs["donor_cell_id"].astype(str)
    adata.obs_names_make_unique()

    return adata


def load_one_xenium_sample(row: pd.Series) -> tuple[ad.AnnData, dict]:
    """
    Load one raw Xenium sample.
    """
    brnum = str(row["BrNum"])
    matrix_path = Path(row["matrix_path"])
    cells_path = Path(row["cells_path"])

    section(f"Loading {brnum}")

    print(f"Matrix: {matrix_path}")
    print(f"Cells:  {cells_path}")

    if not matrix_path.exists():
        raise FileNotFoundError(f"Missing cell_feature_matrix.h5 for {brnum}: {matrix_path}")

    if not cells_path.exists():
        raise FileNotFoundError(f"Missing cells metadata for {brnum}: {cells_path}")

    # Important: gex_only=False preserves non-gene features.
    # Later QC may need negative probes, negative codewords, unassigned counts, etc.
    adata = sc.read_10x_h5(matrix_path, gex_only=False)

    adata = clean_var_metadata(adata)

    cells = read_cells_metadata(cells_path)

    adata = attach_cells_metadata(adata, cells, brnum)

    # Keep raw counts.
    adata.layers["counts"] = adata.X.copy()

    # Store raw path metadata.
    adata.obs["raw_matrix_path"] = str(matrix_path)
    adata.obs["raw_cells_path"] = str(cells_path)

    summary = {
        "BrNum": brnum,
        "n_cells": int(adata.n_obs),
        "n_features": int(adata.n_vars),
        "n_gene_expression_features": int(adata.var["is_gene_expression"].sum()),
        "n_non_gene_features": int((~adata.var["is_gene_expression"]).sum()),
        "has_spatial": "spatial" in adata.obsm,
        "has_transcript_counts": "transcript_counts" in adata.obs.columns,
        "has_cell_area": "cell_area" in adata.obs.columns,
        "has_nucleus_area": "nucleus_area" in adata.obs.columns,
    }

    if "transcript_counts" in adata.obs.columns:
        tc = pd.to_numeric(adata.obs["transcript_counts"], errors="coerce")
        summary["median_transcript_counts"] = float(tc.median())
        summary["mean_transcript_counts"] = float(tc.mean())
    else:
        summary["median_transcript_counts"] = np.nan
        summary["mean_transcript_counts"] = np.nan

    print(adata)
    return adata, summary


def save_plots(adata: ad.AnnData, sample_summary: pd.DataFrame) -> None:
    """
    Build-stage diagnostic plots.
    These are not QC filtering plots yet.
    """

    section("Saving diagnostic plots")

    # Plot 1: cells per sample
    plt.figure(figsize=(12, 5))
    plt.bar(sample_summary["BrNum"], sample_summary["n_cells"])
    plt.xticks(rotation=90)
    plt.ylabel("Number of cells")
    plt.title("Xenium raw cells per sample")
    plt.tight_layout()
    out = FIGURE_DIR / "xenium_N24_cells_per_sample.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")

    # Plot 2: features per sample
    plt.figure(figsize=(12, 5))
    plt.bar(sample_summary["BrNum"], sample_summary["n_features"])
    plt.xticks(rotation=90)
    plt.ylabel("Number of features")
    plt.title("Xenium raw features per sample")
    plt.tight_layout()
    out = FIGURE_DIR / "xenium_N24_features_per_sample.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")

    # Plot 3: spatial locations panel
    brnums = sorted(adata.obs["BrNum"].unique())
    n = len(brnums)
    ncols = 4
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.array(axes).reshape(-1)

    for ax, brnum in zip(axes, brnums):
        obs = adata.obs.loc[adata.obs["BrNum"] == brnum, ["x_centroid", "y_centroid"]].copy()

        # Downsample for plotting speed.
        if obs.shape[0] > 25000:
            obs = obs.sample(25000, random_state=0)

        ax.scatter(obs["x_centroid"], obs["y_centroid"], s=0.1)
        ax.invert_yaxis()
        ax.set_title(brnum)
        ax.set_xticks([])
        ax.set_yticks([])

    for ax in axes[len(brnums):]:
        ax.axis("off")

    plt.suptitle("Raw Xenium cell spatial locations", y=1.0)
    plt.tight_layout()
    out = FIGURE_DIR / "xenium_N24_raw_spatial_locations.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--raw-xenium-dir",
        default=str(RAW_XENIUM_DIR),
        help="Folder containing raw Xenium Br* folders.",
    )

    parser.add_argument(
        "--sample",
        default=None,
        help="Optional single sample test, e.g. Br2039.",
    )

    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional maximum number of samples for testing.",
    )

    args = parser.parse_args()

    section("Build Xenium N24 raw AnnData")

    raw_xenium_dir = Path(args.raw_xenium_dir)
    print(f"Project root:    {PROJECT_ROOT}")
    print(f"Raw Xenium dir:  {raw_xenium_dir}")
    print(f"Output h5ad:     {OUT_H5AD}")

    sample_paths = scan_xenium_folders(raw_xenium_dir)

    # Keep only samples with required files.
    sample_paths = sample_paths[
        sample_paths["has_matrix"].astype(bool) &
        sample_paths["has_cells"].astype(bool)
    ].copy()

    if args.sample is not None:
        sample_paths = sample_paths[sample_paths["BrNum"].eq(args.sample)].copy()

    if args.max_samples is not None:
        sample_paths = sample_paths.head(args.max_samples).copy()

    if sample_paths.empty:
        raise RuntimeError("No valid samples to load.")

    print("\nSamples selected:")
    print(sample_paths[["BrNum", "matrix_path", "cells_path"]].to_string(index=False))

    adata_list = []
    summaries = []
    status_records = []

    for _, row in sample_paths.iterrows():
        brnum = row["BrNum"]

        try:
            adata_one, summary = load_one_xenium_sample(row)
            adata_list.append(adata_one)
            summaries.append(summary)

            status_records.append(
                {
                    "BrNum": brnum,
                    "status": "ok",
                    "error": "",
                    "n_cells": adata_one.n_obs,
                    "n_features": adata_one.n_vars,
                }
            )

        except Exception as e:
            print(f"\nFAILED: {brnum}")
            print(e)
            traceback.print_exc()

            status_records.append(
                {
                    "BrNum": brnum,
                    "status": "failed",
                    "error": str(e),
                    "n_cells": 0,
                    "n_features": 0,
                }
            )

    build_status = pd.DataFrame(status_records)
    build_status.to_csv(TABLE_DIR / "xenium_N24_build_status.csv", index=False)

    if len(adata_list) == 0:
        raise RuntimeError("No samples loaded successfully. Check xenium_N24_build_status.csv.")

    section("Combining samples")

    # AnnData rows are cells. Combining donors means stacking rows.
    adata_all = ad.concat(
        adata_list,
        axis=0,
        join="outer",
        merge="first",
        uns_merge="first",
        index_unique=None,
    )

    # Add useful global metadata.
    adata_all.uns["stage"] = "01_build_datastructure"
    adata_all.uns["object_name"] = "xenium_N24_raw"
    adata_all.uns["description"] = (
        "Raw Xenium AnnData built from cell_feature_matrix.h5 and cells metadata. "
        "Equivalent purpose to PI raw_spe_N24.RDS."
    )
    adata_all.uns["counts_layer"] = "counts"
    adata_all.uns["spatial_key"] = "spatial"

    print(adata_all)

    section("Saving outputs")

    # Main h5ad
    adata_all.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"Saved h5ad: {OUT_H5AD}")

    # CSV outputs
    sample_summary = pd.DataFrame(summaries).sort_values("BrNum").reset_index(drop=True)

    sample_summary.to_csv(
        TABLE_DIR / "xenium_N24_sample_summary.csv",
        index=False,
    )

    adata_all.obs.to_csv(
        TABLE_DIR / "xenium_N24_cell_metadata.csv.gz",
        compression="gzip",
    )

    adata_all.var.to_csv(
        TABLE_DIR / "xenium_N24_feature_metadata.csv",
    )

    save_plots(adata_all, sample_summary)

    section("Done")

    print(f"Final h5ad: {OUT_H5AD}")
    print(f"Cells:      {adata_all.n_obs}")
    print(f"Features:   {adata_all.n_vars}")
    print(f"Donors:     {adata_all.obs['BrNum'].nunique()}")

    if adata_all.obs["BrNum"].nunique() != 24 and args.sample is None and args.max_samples is None:
        print("\nWARNING: Final object does not contain 24 donors.")
        print("Check data/raw/xenium and sample_paths.csv.")


if __name__ == "__main__":
    main()