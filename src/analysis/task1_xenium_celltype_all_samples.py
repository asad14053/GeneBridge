"""
Task 1: Xenium cell-type analysis for ALL Xenium Br folders.

This is the all-sample version of your Br2039 marker-based cell-type script.
It loops over data/raw/xenium/Br*/ and generates one cell-type map plus
annotation/count tables for each sample.

Expected per sample:
    *cell_feature_matrix.h5
    *cells.parquet.gz

Run from project root:
    python src/analysis/task1_xenium_celltype_all_samples.py

Test one sample:
    python src/analysis/task1_xenium_celltype_all_samples.py --sample Br2039

Test first 3 samples:
    python src/analysis/task1_xenium_celltype_all_samples.py --max-samples 3

Make dots bigger:
    python src/analysis/task1_xenium_celltype_all_samples.py --point-size 6
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
from anndata import AnnData


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
XENIUM_ROOT = PROJECT_ROOT / "data" / "raw" / "xenium"

OUT_ROOT = PROJECT_ROOT / "outputs" / "task1_xenium_celltype_all_samples"
ALL_TABLE_DIR = OUT_ROOT / "tables"
ALL_TABLE_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# PI-provided categorical color map
# =============================================================================

cat_color = [
    "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
    "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
    "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
    "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785"
]


# =============================================================================
# Broad cell-type markers
# =============================================================================

CELLTYPE_MARKERS = {
    "Excitatory_neuron": ["SLC17A7", "SATB2", "RORB", "PCP4"],
    "Inhibitory_neuron": ["GAD1", "GAD2", "PVALB", "SST", "VIP"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG"],
    "OPC": ["PDGFRA", "CSPG4"],
    "Microglia": ["CX3CR1", "P2RY12", "C3"],
    "Endothelial": ["CLDN5", "FLT1", "VWF"],
    "Mural": ["PDGFRB", "RGS5"],
}

CELLTYPE_ORDER = [
    "Excitatory_neuron",
    "Inhibitory_neuron",
    "Astrocyte",
    "Oligodendrocyte",
    "OPC",
    "Microglia",
    "Endothelial",
    "Mural",
    "Unknown",
]

CELLTYPE_COLORS = {
    celltype: cat_color[i]
    for i, celltype in enumerate(CELLTYPE_ORDER)
}


# =============================================================================
# Utility functions
# =============================================================================

def print_section(title: str):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def find_one_file(folder: Path, pattern: str):
    """Find one matching file recursively."""
    files = sorted(folder.rglob(pattern))

    if len(files) == 0:
        return None

    if len(files) > 1:
        print(f"WARNING: Multiple files found for pattern '{pattern}' in {folder}. Using:")
        print(files[0])

    return files[0]


def read_parquet_maybe_gzip(path: Path) -> pd.DataFrame:
    """
    Read cells.parquet.gz robustly.

    Handles normal Parquet and gzip-wrapped Parquet.
    This fixes: Parquet magic bytes not found in footer.
    """
    path = Path(path)
    print(f"Reading parquet-like file: {path}")

    try:
        return pd.read_parquet(path)
    except Exception as normal_error:
        print("Normal pd.read_parquet() failed.")
        print(f"Reason: {normal_error}")

    with open(path, "rb") as f:
        first_bytes = f.read(4)

    if first_bytes == b"PAR1":
        return pd.read_parquet(path)

    if first_bytes[:2] == b"\x1f\x8b":
        print("Detected gzip-compressed parquet. Decompressing temporarily...")
        temp_path = None
        try:
            with gzip.open(path, "rb") as gz:
                with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
                    shutil.copyfileobj(gz, tmp)
                    temp_path = Path(tmp.name)
            return pd.read_parquet(temp_path)
        finally:
            if temp_path is not None and temp_path.exists():
                try:
                    os.remove(temp_path)
                    print("Temporary parquet removed.")
                except Exception as cleanup_error:
                    print(f"WARNING: Could not remove temporary file: {temp_path}")
                    print(cleanup_error)

    raise ValueError(
        f"Could not read file as Parquet or gzip-compressed Parquet: {path}\n"
        f"First bytes were: {first_bytes}\n"
        "Possible causes: incomplete/corrupted download, HTML file, or wrong extension."
    )


def sparse_mean_by_rows(X):
    """Calculate mean expression per row for dense or sparse matrices."""
    if sp.issparse(X):
        return np.asarray(X.mean(axis=1)).ravel()
    return np.mean(X, axis=1)


def get_sample_dirs(xenium_root: Path, requested_samples=None, max_samples=None):
    """Return Br* directories under data/raw/xenium."""
    if not xenium_root.exists():
        raise FileNotFoundError(f"Xenium root does not exist: {xenium_root}")

    sample_dirs = sorted(
        [p for p in xenium_root.iterdir() if p.is_dir() and p.name.startswith("Br")]
    )

    if requested_samples:
        requested = set(requested_samples)
        sample_dirs = [p for p in sample_dirs if p.name in requested]

    if max_samples is not None:
        sample_dirs = sample_dirs[:max_samples]

    return sample_dirs


def normalize_log1p(adata: AnnData) -> AnnData:
    """Normalize total counts and log-transform."""
    adata = adata.copy()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


# =============================================================================
# Loading and annotation
# =============================================================================

def load_xenium_sample(sample_dir: Path) -> AnnData:
    """
    Load one Xenium sample and add x/y centroid coordinates to adata.obs.
    """
    sample_id = sample_dir.name
    print_section(f"Loading Xenium sample: {sample_id}")

    h5_path = find_one_file(sample_dir, "*cell_feature_matrix.h5")
    cells_path = find_one_file(sample_dir, "*cells.parquet.gz")

    if h5_path is None:
        raise FileNotFoundError(f"No *cell_feature_matrix.h5 found in: {sample_dir}")
    if cells_path is None:
        raise FileNotFoundError(f"No *cells.parquet.gz found in: {sample_dir}")

    print(f"Xenium matrix file: {h5_path}")
    print(f"Xenium cells metadata file: {cells_path}")

    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    cells = read_parquet_maybe_gzip(cells_path)

    for col in ["cell_id", "x_centroid", "y_centroid"]:
        if col not in cells.columns:
            raise ValueError(f"cells.parquet.gz must contain column: {col}")

    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    adata.obs_names = adata.obs_names.astype(str)

    common_cells = adata.obs_names.intersection(cells.index)

    print(f"Xenium matrix cells: {adata.n_obs}")
    print(f"Xenium metadata cells: {cells.shape[0]}")
    print(f"Matched cells: {len(common_cells)}")
    print(f"Xenium genes/features: {adata.n_vars}")

    if len(common_cells) == 0:
        print("\nFirst 5 matrix cell IDs:")
        print(list(adata.obs_names[:5]))
        print("\nFirst 5 metadata cell IDs:")
        print(list(cells.index[:5]))
        raise ValueError("No matching cell IDs between H5 matrix and cells.parquet.gz.")

    adata = adata[common_cells].copy()
    adata.obs = adata.obs.join(cells.loc[common_cells], how="left")
    adata.obs["sample_id"] = sample_id
    adata.obs["technology"] = "Xenium"

    return adata


def marker_score_annotation(adata: AnnData, sample_id: str, table_dir: Path) -> AnnData:
    """
    First-pass broad cell-type annotation using mean expression of marker genes.
    """
    print_section(f"Marker-based Xenium cell-type annotation: {sample_id}")

    adata = normalize_log1p(adata)
    marker_summary = []

    for celltype, marker_genes in CELLTYPE_MARKERS.items():
        available_genes = [g for g in marker_genes if g in adata.var_names]
        score_col = f"score_{celltype}"

        if len(available_genes) == 0:
            adata.obs[score_col] = 0.0
            print(f"{celltype}: no markers found")
        else:
            X_sub = adata[:, available_genes].X
            adata.obs[score_col] = sparse_mean_by_rows(X_sub)
            print(f"{celltype}: using {len(available_genes)} markers: {available_genes}")

        marker_summary.append({
            "sample_id": sample_id,
            "celltype": celltype,
            "requested_markers": ",".join(marker_genes),
            "available_markers": ",".join(available_genes),
            "n_available_markers": len(available_genes),
        })

    score_cols = [f"score_{ct}" for ct in CELLTYPE_MARKERS.keys()]
    score_df = adata.obs[score_cols].copy()

    adata.obs["predicted_celltype"] = (
        score_df.idxmax(axis=1).str.replace("score_", "", regex=False)
    )
    adata.obs["max_celltype_score"] = score_df.max(axis=1)

    adata.obs.loc[
        adata.obs["max_celltype_score"] <= 0,
        "predicted_celltype"
    ] = "Unknown"

    final_order = [
        ct for ct in CELLTYPE_ORDER
        if ct in list(adata.obs["predicted_celltype"].unique())
    ]

    adata.obs["predicted_celltype"] = pd.Categorical(
        adata.obs["predicted_celltype"],
        categories=final_order,
        ordered=True,
    )

    marker_summary_df = pd.DataFrame(marker_summary)
    marker_out = table_dir / f"{sample_id}_xenium_marker_availability.csv"
    marker_summary_df.to_csv(marker_out, index=False)
    print(f"Saved marker availability: {marker_out}")

    print(f"\n{sample_id} predicted cell-type counts:")
    print(adata.obs["predicted_celltype"].value_counts(dropna=False))

    return adata


# =============================================================================
# Output functions
# =============================================================================

def save_celltype_tables(adata: AnnData, sample_id: str, table_dir: Path):
    """Save cell-level annotation and per-celltype count table."""
    keep_cols = [
        "sample_id",
        "technology",
        "x_centroid",
        "y_centroid",
        "predicted_celltype",
        "max_celltype_score",
    ]
    keep_cols = [c for c in keep_cols if c in adata.obs.columns]

    anno = adata.obs[keep_cols].copy()
    anno.index.name = "cell_id"

    anno_out = table_dir / f"{sample_id}_xenium_celltype_annotations.csv"
    anno.to_csv(anno_out)

    counts = (
        adata.obs["predicted_celltype"]
        .value_counts(dropna=False)
        .rename_axis("predicted_celltype")
        .reset_index(name="n_cells")
    )
    counts["sample_id"] = sample_id
    counts["pct_cells"] = 100 * counts["n_cells"] / counts["n_cells"].sum()
    counts = counts[["sample_id", "predicted_celltype", "n_cells", "pct_cells"]]

    counts_out = table_dir / f"{sample_id}_xenium_celltype_counts.csv"
    counts.to_csv(counts_out, index=False)

    print(f"Saved: {anno_out}")
    print(f"Saved: {counts_out}")

    return counts


def plot_xenium_celltype_map(
    adata: AnnData,
    sample_id: str,
    figure_dir: Path,
    point_size: float = 5,
):
    """Plot Xenium cell-level spatial map colored by predicted cell type."""
    print_section(f"Plotting Xenium cell-type map: {sample_id}")

    plot_df = adata.obs.copy()

    plt.figure(figsize=(8, 8))

    for celltype in CELLTYPE_ORDER:
        sub = plot_df[plot_df["predicted_celltype"].astype(str) == celltype]
        if sub.shape[0] == 0:
            continue

        plt.scatter(
            sub["x_centroid"],
            sub["y_centroid"],
            s=point_size,
            alpha=0.75,
            color=CELLTYPE_COLORS.get(celltype, "#787878"),
            label=f"{celltype} ({sub.shape[0]})",
            linewidths=0,
        )

    plt.gca().invert_yaxis()
    plt.xlabel("X centroid")
    plt.ylabel("Y centroid")
    plt.title(f"{sample_id} Xenium cell-level map\nBroad marker-based cell type")

    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        frameon=False,
        markerscale=5,
        fontsize=8,
    )

    plt.tight_layout()

    out_path = figure_dir / f"{sample_id}_xenium_celltype_map.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def save_color_legend():
    """Save global color legend."""
    legend = pd.DataFrame({
        "celltype": CELLTYPE_ORDER,
        "color": [CELLTYPE_COLORS[ct] for ct in CELLTYPE_ORDER],
    })
    out = ALL_TABLE_DIR / "xenium_celltype_color_legend.csv"
    legend.to_csv(out, index=False)
    print(f"Saved color legend: {out}")


def process_one_sample(sample_dir: Path, point_size: float = 5):
    """Full cell-type workflow for one Xenium sample."""
    sample_id = sample_dir.name

    sample_out = OUT_ROOT / sample_id
    table_dir = sample_out / "tables"
    figure_dir = sample_out / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    adata = load_xenium_sample(sample_dir)
    adata = marker_score_annotation(adata, sample_id=sample_id, table_dir=table_dir)

    counts = save_celltype_tables(adata, sample_id=sample_id, table_dir=table_dir)
    plot_xenium_celltype_map(
        adata,
        sample_id=sample_id,
        figure_dir=figure_dir,
        point_size=point_size,
    )

    summary = {
        "sample_id": sample_id,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_predicted_celltypes": int(adata.obs["predicted_celltype"].nunique()),
        "n_unknown": int((adata.obs["predicted_celltype"].astype(str) == "Unknown").sum()),
        "pct_unknown": float(100 * (adata.obs["predicted_celltype"].astype(str) == "Unknown").mean()),
        "median_max_celltype_score": float(adata.obs["max_celltype_score"].median()),
    }

    summary_out = table_dir / f"{sample_id}_xenium_celltype_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_out, index=False)
    print(f"Saved: {summary_out}")

    return summary, counts


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample",
        nargs="*",
        default=None,
        help="Optional sample IDs, e.g. --sample Br2039 Br2719",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Optional maximum number of samples to process for testing.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=5,
        help="Xenium cell dot size. Try 3, 5, or 8.",
    )

    args = parser.parse_args()

    print_section("Task 1: Xenium cell-type analysis for all samples")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Xenium root: {XENIUM_ROOT}")
    print(f"Output root: {OUT_ROOT}")

    save_color_legend()

    sample_dirs = get_sample_dirs(
        XENIUM_ROOT,
        requested_samples=args.sample,
        max_samples=args.max_samples,
    )

    print(f"Found {len(sample_dirs)} Xenium sample folders:")
    for p in sample_dirs:
        print("  ", p.name)

    all_summaries = []
    all_counts = []
    status_records = []

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        try:
            summary, counts = process_one_sample(
                sample_dir=sample_dir,
                point_size=args.point_size,
            )
            all_summaries.append(summary)
            all_counts.append(counts)
            status_records.append({"sample_id": sample_id, "status": "ok", "error": ""})
        except Exception as e:
            print_section(f"ERROR while processing {sample_id}")
            print(e)
            traceback.print_exc()
            status_records.append({"sample_id": sample_id, "status": "failed", "error": str(e)})

    status_df = pd.DataFrame(status_records)
    status_out = ALL_TABLE_DIR / "all_xenium_celltype_processing_status.csv"
    status_df.to_csv(status_out, index=False)

    summary_df = pd.DataFrame(all_summaries) if all_summaries else pd.DataFrame()
    summary_out = ALL_TABLE_DIR / "all_xenium_celltype_summary.csv"
    summary_df.to_csv(summary_out, index=False)

    counts_df = pd.concat(all_counts, ignore_index=True) if all_counts else pd.DataFrame()
    counts_out = ALL_TABLE_DIR / "all_xenium_celltype_counts.csv"
    counts_df.to_csv(counts_out, index=False)

    print_section("Done")
    print("Saved all-sample summary:")
    print(summary_out)
    print("Saved all-sample counts:")
    print(counts_out)
    print("Saved processing status:")
    print(status_out)

    if not summary_df.empty:
        print("\nPreview:")
        print(summary_df)


if __name__ == "__main__":
    main()
