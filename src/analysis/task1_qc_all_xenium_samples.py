"""
Task 1 helper: Run Xenium QC and basic spatial plots for ALL Xenium Br folders.

This script is a generalized version of your one-sample Br2039 code.

It loops over:
    data/raw/xenium/Br*/

For each Br sample, it expects:
    *cell_feature_matrix.h5
    *cells.parquet.gz

For each sample, it saves:
    outputs/task1_xenium_all_samples/<BrID>/tables/
    outputs/task1_xenium_all_samples/<BrID>/figures/

It also saves:
    outputs/task1_xenium_all_samples/tables/all_xenium_qc_summary.csv
    outputs/task1_xenium_all_samples/tables/all_xenium_processing_status.csv

Run from project root:
    python src/analysis/task1_qc_all_xenium_samples.py

Optional:
    python src/analysis/task1_qc_all_xenium_samples.py --sample Br2039
    python src/analysis/task1_qc_all_xenium_samples.py --max-samples 3
"""

from __future__ import annotations

from pathlib import Path
import argparse
import gzip
import shutil
import tempfile
import os
import traceback

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

XENIUM_ROOT = PROJECT_ROOT / "data" / "raw" / "xenium"

OUT_ROOT = PROJECT_ROOT / "outputs" / "task1_xenium_all_samples"
ALL_TABLE_DIR = OUT_ROOT / "tables"
ALL_FIGURE_DIR = OUT_ROOT / "figures"

ALL_TABLE_DIR.mkdir(parents=True, exist_ok=True)
ALL_FIGURE_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Utility functions
# =============================================================================

def print_section(title: str):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def find_one_file(folder: Path, pattern: str):
    """
    Find one file in a sample folder.
    Uses recursive search because files may be nested.
    """
    files = sorted(folder.rglob(pattern))

    if len(files) == 0:
        return None

    if len(files) > 1:
        print(f"WARNING: Multiple files found for {pattern} in {folder}. Using:")
        print(files[0])

    return files[0]


def read_parquet_maybe_gzip(path: Path) -> pd.DataFrame:
    """
    Read cells.parquet.gz robustly.

    Handles:
        1. normal parquet
        2. gzip-wrapped parquet

    This avoids:
        pyarrow.lib.ArrowInvalid: Parquet magic bytes not found in footer
    """
    path = Path(path)

    print(f"Reading cells metadata: {path}")

    # First try normal pandas reader.
    try:
        return pd.read_parquet(path)
    except Exception as normal_error:
        print("Normal pd.read_parquet() failed.")
        print(f"Reason: {normal_error}")

    # Check first bytes.
    with open(path, "rb") as f:
        first_bytes = f.read(4)

    # Normal parquet magic bytes.
    if first_bytes == b"PAR1":
        return pd.read_parquet(path)

    # Gzip magic bytes.
    if first_bytes[:2] == b"\x1f\x8b":
        print("Detected gzip-compressed parquet. Decompressing temporarily...")

        temp_path = None

        try:
            with gzip.open(path, "rb") as gz:
                with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
                    shutil.copyfileobj(gz, tmp)
                    temp_path = Path(tmp.name)

            df = pd.read_parquet(temp_path)
            return df

        finally:
            if temp_path is not None and temp_path.exists():
                try:
                    os.remove(temp_path)
                    print("Temporary parquet removed.")
                except Exception as cleanup_error:
                    print(f"WARNING: could not remove temporary file: {temp_path}")
                    print(cleanup_error)

    raise ValueError(
        f"Could not read file as parquet or gzip-compressed parquet: {path}\n"
        f"First bytes: {first_bytes}\n"
        "The file may be corrupted or incompletely downloaded."
    )


def safe_quantile(series: pd.Series, q: float):
    clean = series.replace([np.inf, -np.inf], np.nan).dropna()

    if clean.empty:
        return np.nan

    return float(clean.quantile(q))


def save_hist(cells: pd.DataFrame, column: str, out_path: Path, threshold=None):
    plt.figure(figsize=(7, 5))
    plt.hist(cells[column].replace([np.inf, -np.inf], np.nan).dropna(), bins=60)

    if threshold is not None and not pd.isna(threshold):
        plt.axvline(threshold, linestyle="--")

    plt.xlabel(column)
    plt.ylabel("Number of cells")
    plt.title(column)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def get_sample_dirs(xenium_root: Path, requested_samples=None, max_samples=None):
    """
    Return all Br* sample directories under data/raw/xenium.
    """
    if not xenium_root.exists():
        raise FileNotFoundError(f"Xenium root does not exist: {xenium_root}")

    sample_dirs = sorted([p for p in xenium_root.iterdir() if p.is_dir() and p.name.startswith("Br")])

    if requested_samples:
        requested = set(requested_samples)
        sample_dirs = [p for p in sample_dirs if p.name in requested]

    if max_samples is not None:
        sample_dirs = sample_dirs[:max_samples]

    return sample_dirs


# =============================================================================
# Core sample processing
# =============================================================================

def process_one_xenium_sample(sample_dir: Path):
    """
    Process one Xenium Br folder.
    """
    sample_id = sample_dir.name

    print_section(f"Processing Xenium sample: {sample_id}")

    sample_out = OUT_ROOT / sample_id
    table_dir = sample_out / "tables"
    figure_dir = sample_out / "figures"

    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    xenium_h5 = find_one_file(sample_dir, "*cell_feature_matrix.h5")
    cell_parquet_gz = find_one_file(sample_dir, "*cells.parquet.gz")

    if xenium_h5 is None:
        raise FileNotFoundError(f"No *cell_feature_matrix.h5 found in {sample_dir}")

    if cell_parquet_gz is None:
        raise FileNotFoundError(f"No *cells.parquet.gz found in {sample_dir}")

    print("Project root:", PROJECT_ROOT)
    print("Sample ID:", sample_id)
    print("Sample dir:", sample_dir)
    print("Xenium H5:", xenium_h5)
    print("Cells parquet gz:", cell_parquet_gz)

    # -------------------------------------------------------------------------
    # Read expression matrix
    # -------------------------------------------------------------------------

    xenium = sc.read_10x_h5(xenium_h5)
    xenium.var_names_make_unique()

    print("\n===== Xenium expression matrix =====")
    print(xenium)
    print("Cells:", xenium.n_obs)
    print("Features/genes:", xenium.n_vars)

    total_cells_h5_before_match = xenium.n_obs
    total_features = xenium.n_vars

    # -------------------------------------------------------------------------
    # Read cells metadata
    # -------------------------------------------------------------------------

    cells = read_parquet_maybe_gzip(cell_parquet_gz)

    print("\n===== Xenium cell metadata =====")
    print("Cells parquet shape:", cells.shape)
    print("Columns:", cells.columns.tolist())

    # -------------------------------------------------------------------------
    # Required columns
    # -------------------------------------------------------------------------

    required_columns = [
        "cell_id",
        "x_centroid",
        "y_centroid",
        "transcript_counts",
        "control_probe_counts",
        "control_codeword_counts",
        "unassigned_codeword_counts",
        "total_counts",
        "cell_area",
        "nucleus_area",
    ]

    missing_columns = [col for col in required_columns if col not in cells.columns]

    if len(missing_columns) > 0:
        raise ValueError(f"Missing required columns in cells.parquet: {missing_columns}")

    # -------------------------------------------------------------------------
    # Match cell IDs
    # -------------------------------------------------------------------------

    cells["cell_id"] = cells["cell_id"].astype(str)
    xenium.obs_names = xenium.obs_names.astype(str)

    common_cell_ids = sorted(set(cells["cell_id"]) & set(xenium.obs_names))

    print("\n===== Cell ID matching =====")
    print("Cells in H5:", xenium.n_obs)
    print("Cells in parquet:", cells.shape[0])
    print("Matched cells:", len(common_cell_ids))

    if len(common_cell_ids) == 0:
        print("\nFirst 5 H5 cell IDs:")
        print(list(xenium.obs_names[:5]))
        print("\nFirst 5 parquet cell IDs:")
        print(cells["cell_id"].head().tolist())

        raise ValueError("Cell IDs do not match between H5 and cells.parquet.")

    # Keep only matched cells
    xenium = xenium[common_cell_ids, :].copy()

    cells = (
        cells
        .set_index("cell_id")
        .loc[common_cell_ids]
        .copy()
    )

    cells["sample_id"] = sample_id

    # -------------------------------------------------------------------------
    # Expression-based QC
    # -------------------------------------------------------------------------

    xenium.var["mt"] = xenium.var_names.str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        xenium,
        qc_vars=["mt"],
        inplace=True,
        percent_top=None,
    )

    cells["n_genes_by_counts"] = xenium.obs["n_genes_by_counts"].values
    cells["matrix_total_counts"] = xenium.obs["total_counts"].values
    cells["pct_counts_mt"] = xenium.obs["pct_counts_mt"].values

    # -------------------------------------------------------------------------
    # Metadata-based QC percentages
    # -------------------------------------------------------------------------

    safe_total = cells["total_counts"].replace(0, np.nan)

    cells["pct_control_probe"] = 100 * cells["control_probe_counts"] / safe_total
    cells["pct_control_codeword"] = 100 * cells["control_codeword_counts"] / safe_total
    cells["pct_unassigned"] = 100 * cells["unassigned_codeword_counts"] / safe_total

    # -------------------------------------------------------------------------
    # Outlier thresholds, per sample
    # -------------------------------------------------------------------------

    cells["empty_cell_out"] = cells["matrix_total_counts"] == 0

    probe_thresh = safe_quantile(cells["pct_control_probe"], 0.99)
    codeword_thresh = safe_quantile(cells["pct_control_codeword"], 0.99)
    unassigned_thresh = safe_quantile(cells["pct_unassigned"], 0.99)

    mito_thresh = safe_quantile(cells["pct_counts_mt"], 0.99)
    low_detected_thresh = safe_quantile(cells["n_genes_by_counts"], 0.01)

    low_total_counts_thresh = safe_quantile(cells["matrix_total_counts"], 0.01)
    high_total_counts_thresh = safe_quantile(cells["matrix_total_counts"], 0.99)

    cells["neg_probe_out"] = cells["pct_control_probe"] >= probe_thresh
    cells["neg_codeword_out"] = cells["pct_control_codeword"] >= codeword_thresh
    cells["unassigned_out"] = cells["pct_unassigned"] >= unassigned_thresh

    cells["mito_out"] = cells["pct_counts_mt"] >= mito_thresh
    cells["detected_out"] = cells["n_genes_by_counts"] <= low_detected_thresh

    cells["total_counts_out"] = (
        (cells["matrix_total_counts"] <= low_total_counts_thresh)
        |
        (cells["matrix_total_counts"] >= high_total_counts_thresh)
    )

    # Lieber-like final union did NOT include mito_out.
    cells["is_outlier_lieber_like"] = (
        cells["empty_cell_out"]
        |
        cells["neg_probe_out"]
        |
        cells["neg_codeword_out"]
        |
        cells["unassigned_out"]
        |
        cells["detected_out"]
        |
        cells["total_counts_out"]
    )

    cells["is_outlier_with_mito"] = (
        cells["is_outlier_lieber_like"]
        |
        cells["mito_out"]
    )

    # -------------------------------------------------------------------------
    # Save QC tables
    # -------------------------------------------------------------------------

    qc_metrics_file = table_dir / f"{sample_id}_xenium_cell_qc_metrics.csv"
    outlier_file = table_dir / f"{sample_id}_xenium_outlier_cells.csv"
    summary_file = table_dir / f"{sample_id}_xenium_qc_summary.csv"

    cells.to_csv(qc_metrics_file)

    outlier_cells = (
        cells[cells["is_outlier_lieber_like"]]
        .reset_index()[["cell_id"]]
    )

    outlier_cells.to_csv(outlier_file, index=False)

    summary_dict = {
        "sample_id": sample_id,
        "sample_dir": str(sample_dir),
        "xenium_h5": str(xenium_h5),
        "cells_parquet": str(cell_parquet_gz),
        "total_cells_h5": int(total_cells_h5_before_match),
        "total_cells_parquet": int(cells.shape[0]),
        "matched_cells": int(len(common_cell_ids)),
        "total_features": int(total_features),
        "empty_cell_out": int(cells["empty_cell_out"].sum()),
        "neg_probe_threshold_99pct": probe_thresh,
        "neg_probe_out": int(cells["neg_probe_out"].sum()),
        "neg_codeword_threshold_99pct": codeword_thresh,
        "neg_codeword_out": int(cells["neg_codeword_out"].sum()),
        "unassigned_threshold_99pct": unassigned_thresh,
        "unassigned_out": int(cells["unassigned_out"].sum()),
        "mito_threshold_99pct": mito_thresh,
        "mito_out": int(cells["mito_out"].sum()),
        "low_detected_threshold_1pct": low_detected_thresh,
        "detected_out": int(cells["detected_out"].sum()),
        "low_total_counts_threshold_1pct": low_total_counts_thresh,
        "high_total_counts_threshold_99pct": high_total_counts_thresh,
        "total_counts_out": int(cells["total_counts_out"].sum()),
        "final_outliers_lieber_like": int(cells["is_outlier_lieber_like"].sum()),
        "final_outliers_with_mito": int(cells["is_outlier_with_mito"].sum()),
        "kept_cells_lieber_like": int((~cells["is_outlier_lieber_like"]).sum()),
        "pct_outlier_lieber_like": float(100 * cells["is_outlier_lieber_like"].mean()),
        "pct_outlier_with_mito": float(100 * cells["is_outlier_with_mito"].mean()),
        "median_total_counts": float(cells["matrix_total_counts"].median()),
        "median_detected_genes": float(cells["n_genes_by_counts"].median()),
        "median_pct_control_probe": float(cells["pct_control_probe"].median(skipna=True)),
        "median_pct_control_codeword": float(cells["pct_control_codeword"].median(skipna=True)),
        "median_pct_unassigned": float(cells["pct_unassigned"].median(skipna=True)),
        "median_pct_mito": float(cells["pct_counts_mt"].median(skipna=True)),
    }

    summary = pd.DataFrame(
        [{"metric": k, "value": v} for k, v in summary_dict.items()]
    )

    summary.to_csv(summary_file, index=False)

    # -------------------------------------------------------------------------
    # Basic plots
    # -------------------------------------------------------------------------

    save_hist(
        cells,
        "pct_control_probe",
        figure_dir / f"{sample_id}_xenium_pct_control_probe_hist.png",
        probe_thresh,
    )

    save_hist(
        cells,
        "pct_control_codeword",
        figure_dir / f"{sample_id}_xenium_pct_control_codeword_hist.png",
        codeword_thresh,
    )

    save_hist(
        cells,
        "pct_unassigned",
        figure_dir / f"{sample_id}_xenium_pct_unassigned_hist.png",
        unassigned_thresh,
    )

    save_hist(
        cells,
        "pct_counts_mt",
        figure_dir / f"{sample_id}_xenium_pct_mito_hist.png",
        mito_thresh,
    )

    save_hist(
        cells,
        "n_genes_by_counts",
        figure_dir / f"{sample_id}_xenium_detected_genes_hist.png",
        low_detected_thresh,
    )

    save_hist(
        cells,
        "matrix_total_counts",
        figure_dir / f"{sample_id}_xenium_total_counts_hist.png",
        high_total_counts_thresh,
    )

    # Spatial plot: all cells
    plt.figure(figsize=(7, 7))
    plt.scatter(
        cells["x_centroid"],
        cells["y_centroid"],
        s=0.2,
    )
    plt.gca().invert_yaxis()
    plt.xlabel("x_centroid")
    plt.ylabel("y_centroid")
    plt.title(f"{sample_id} Xenium cell locations")
    plt.tight_layout()
    cell_locations_file = figure_dir / f"{sample_id}_xenium_cell_locations.png"
    plt.savefig(cell_locations_file, dpi=300)
    plt.close()

    # Spatial plot: outliers
    plt.figure(figsize=(7, 7))

    normal = cells[~cells["is_outlier_lieber_like"]]
    outliers = cells[cells["is_outlier_lieber_like"]]

    plt.scatter(
        normal["x_centroid"],
        normal["y_centroid"],
        s=0.2,
        alpha=0.4,
        label="Kept",
    )

    plt.scatter(
        outliers["x_centroid"],
        outliers["y_centroid"],
        s=1,
        alpha=0.8,
        label="Outlier",
    )

    plt.gca().invert_yaxis()
    plt.xlabel("x_centroid")
    plt.ylabel("y_centroid")
    plt.title(f"{sample_id} Xenium QC outliers")
    plt.legend(frameon=False, markerscale=5)
    plt.tight_layout()
    outlier_plot_file = figure_dir / f"{sample_id}_xenium_spatial_outliers.png"
    plt.savefig(outlier_plot_file, dpi=300)
    plt.close()

    print("\n===== Xenium QC summary =====")
    for k, v in summary_dict.items():
        print(f"{k}: {v}")

    print("\nSaved tables:")
    print(qc_metrics_file)
    print(outlier_file)
    print(summary_file)

    print("\nSaved figures in:")
    print(figure_dir)

    return summary_dict


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample",
        nargs="*",
        default=None,
        help="Optional one or more sample IDs to process, e.g. --sample Br2039 Br2719",
    )

    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional maximum number of samples to process for testing.",
    )

    args = parser.parse_args()

    print_section("Task 1: Run Xenium QC for all Br folders")

    print("Project root:", PROJECT_ROOT)
    print("Xenium root:", XENIUM_ROOT)
    print("Output root:", OUT_ROOT)

    sample_dirs = get_sample_dirs(
        XENIUM_ROOT,
        requested_samples=args.sample,
        max_samples=args.max_samples,
    )

    print(f"Found {len(sample_dirs)} Xenium sample folders:")
    for p in sample_dirs:
        print("  ", p.name)

    all_summaries = []
    status_records = []

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name

        try:
            summary = process_one_xenium_sample(sample_dir)
            all_summaries.append(summary)

            status_records.append({
                "sample_id": sample_id,
                "status": "ok",
                "error": "",
            })

        except Exception as e:
            print_section(f"ERROR while processing {sample_id}")
            print(e)
            traceback.print_exc()

            status_records.append({
                "sample_id": sample_id,
                "status": "failed",
                "error": str(e),
            })

    status_df = pd.DataFrame(status_records)
    status_out = ALL_TABLE_DIR / "all_xenium_processing_status.csv"
    status_df.to_csv(status_out, index=False)

    if all_summaries:
        all_summary_df = pd.DataFrame(all_summaries)
    else:
        all_summary_df = pd.DataFrame()

    all_summary_out = ALL_TABLE_DIR / "all_xenium_qc_summary.csv"
    all_summary_df.to_csv(all_summary_out, index=False)

    print_section("Finished all Xenium samples")

    print("Saved all-sample summary:")
    print(all_summary_out)

    print("Saved processing status:")
    print(status_out)

    if not all_summary_df.empty:
        print("\nPreview all-sample summary:")
        preview_cols = [
            "sample_id",
            "matched_cells",
            "total_features",
            "final_outliers_lieber_like",
            "kept_cells_lieber_like",
            "pct_outlier_lieber_like",
            "median_total_counts",
            "median_detected_genes",
        ]
        preview_cols = [c for c in preview_cols if c in all_summary_df.columns]
        print(all_summary_df[preview_cols])


if __name__ == "__main__":
    main()
