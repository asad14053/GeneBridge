#!/usr/bin/env python3

"""
02_banksy_celltype_annotation.py

Stage 03 Annotation Task 2:
BANKSY-based cell-type annotation for Xenium N24.

Purpose
-------
Use PI-generated BANKSY clustering output:

    banksy_clustering_lambda0.1_res0.7.csv

and attach BANKSY clusters + manual cell-type annotation to:

    xenium_N24_layer_annotated.h5ad

This is a first-pass / BANKSY-based cell-type annotation.

Important
---------
This is NOT the supervised snRNA-seq label-transfer result.
This uses BANKSY clusters and the PI/manual cluster-to-cell-type mapping.

Input
-----
1. data/processed/xenium/xenium_N24_layer_annotated.h5ad
2. data/metadata/celltype_annotations/banksy_clustering_lambda0.1_res0.7.csv

Output h5ad
-----------
data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad

CSV outputs
-----------
outputs/xenium_branch/03_annotation/celltype_banksy/tables/
    xenium_N24_banksy_celltype_annotations.csv
    xenium_N24_banksy_match_report.csv
    xenium_N24_banksy_cluster_to_celltype_mapping.csv
    xenium_N24_banksy_cluster_counts_by_BrNum.csv
    xenium_N24_banksy_cluster_proportions_by_BrNum.csv
    xenium_N24_celltype_counts_by_BrNum.csv
    xenium_N24_celltype_proportions_by_BrNum.csv
    xenium_N24_celltype_summary_by_BrNum.csv

Plot outputs
------------
outputs/xenium_branch/03_annotation/celltype_banksy/figures/
    xenium_N24_banksy_celltype_maps_panel.png
    xenium_N24_banksy_cluster_maps_panel.png
    xenium_N24_banksy_celltype_proportions_by_BrNum.png
    xenium_N24_banksy_cluster_proportions_by_BrNum.png
    xenium_N24_banksy_celltype_overall_distribution.png

Per-sample plot outputs
-----------------------
outputs/xenium_branch/03_annotation/celltype_banksy/figures/per_sample_celltype_maps/
    Br1113_banksy_celltype_map.png
    Br2039_banksy_celltype_map.png
    ...

outputs/xenium_branch/03_annotation/celltype_banksy/figures/per_sample_cluster_maps/
    Br1113_banksy_cluster_map.png
    Br2039_banksy_cluster_map.png
    ...

Notes on inhibitory mapping
---------------------------
The PI scripts contain slightly different interpretations for BANKSY clusters 9 and 12.

Default here uses the final tissue-plotting style:
    cluster 12 -> CGE
    cluster 9  -> MGE

You can switch with:
    --inhibitory-mode marker_heatmap
or:
    --inhibitory-mode explore_subtypes
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_INPUT_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "xenium"
    / "xenium_N24_layer_annotated.h5ad"
)

DEFAULT_BANKSY_CSV = (
    PROJECT_ROOT
    / "data"
    / "metadata"
    / "celltype_annotations"
    / "banksy_clustering_lambda0.1_res0.7.csv"
)

DEFAULT_OUTPUT_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "xenium"
    / "xenium_N24_layer_celltype_annotated.h5ad"
)

OUT_DIR = PROJECT_ROOT / "outputs" / "xenium_branch" / "03_annotation" / "celltype_banksy"
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"
PER_SAMPLE_CELLTYPE_DIR = FIGURE_DIR / "per_sample_celltype_maps"
PER_SAMPLE_CLUSTER_DIR = FIGURE_DIR / "per_sample_cluster_maps"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
PER_SAMPLE_CELLTYPE_DIR.mkdir(parents=True, exist_ok=True)
PER_SAMPLE_CLUSTER_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# BANKSY cluster-to-cell-type mappings
# =============================================================================

# This is the final tissue-plot style mapping from PI code.
# It maps BANKSY cluster IDs to broad/interpreted cell types.
BANKSY_TO_CELLTYPE_PLOT_CLUSTERS = {
    "1": "Oligo",
    "6": "Oligo",
    "8": "Oligo",
    "10": "Ambig/Oligo",
    "7": "Mic",
    "15": "Ambig/In/Endo",
    "18": "L5 Ex",
    "14": "L6 Ex",
    "11": "L4/5 Ex",
    "2": "L2/3 Ex",
    "17": "Ast",
    "16": "Ast",
    "4": "Ast",
    "5": "Endo",
    "13": "Endo",
    "3": "Endo",
    "12": "CGE",
    "9": "MGE",
}

# This is the marker-heatmap script version.
# Main difference:
#     12 -> MGE
#     9  -> CGE
BANKSY_TO_CELLTYPE_MARKER_HEATMAP = {
    **BANKSY_TO_CELLTYPE_PLOT_CLUSTERS,
    "12": "MGE",
    "9": "CGE",
}

# This is the earlier exploratory subtype naming.
BANKSY_TO_CELLTYPE_EXPLORE_SUBTYPES = {
    **BANKSY_TO_CELLTYPE_PLOT_CLUSTERS,
    "12": "In: VIP+, LAMP5",
    "9": "In: SST+, PVALB+",
}

CELLTYPE_ORDER = [
    "L2/3 Ex",
    "L4/5 Ex",
    "L5 Ex",
    "L6 Ex",
    "MGE",
    "CGE",
    "In: VIP+, LAMP5",
    "In: SST+, PVALB+",
    "Ast",
    "Oligo",
    "Ambig/Oligo",
    "Mic",
    "Endo",
    "Ambig/In/Endo",
    "Unknown",
]

# Use stable colors for plots.
CELLTYPE_COLORS = {
    "L2/3 Ex": "#1f77b4",
    "L4/5 Ex": "#aec7e8",
    "L5 Ex": "#2ca02c",
    "L6 Ex": "#d62728",
    "MGE": "#9467bd",
    "CGE": "#ff7f0e",
    "In: VIP+, LAMP5": "#ffbb78",
    "In: SST+, PVALB+": "#c5b0d5",
    "Ast": "#8c564b",
    "Oligo": "#17becf",
    "Ambig/Oligo": "#9edae5",
    "Mic": "#bcbd22",
    "Endo": "#e377c2",
    "Ambig/In/Endo": "#7f7f7f",
    "Unknown": "#999999",
}


# =============================================================================
# Helper functions
# =============================================================================

def section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def ensure_brnum(adata: ad.AnnData) -> ad.AnnData:
    """
    Ensure adata.obs['BrNum'] exists.
    If missing, derive it from cell IDs like Br1113_1.
    """
    if "BrNum" not in adata.obs.columns:
        print("BrNum not found. Deriving BrNum from cell IDs.")
        adata.obs["BrNum"] = pd.Index(adata.obs_names.astype(str)).str.split("_").str[0]
    else:
        adata.obs["BrNum"] = adata.obs["BrNum"].astype(str)

    return adata


def get_spatial_xy(adata: ad.AnnData) -> pd.DataFrame:
    """
    Return x/y coordinates from obs x_centroid/y_centroid or obsm['spatial'].
    """
    if "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        return pd.DataFrame(
            {
                "x": pd.to_numeric(adata.obs["x_centroid"], errors="coerce").values,
                "y": pd.to_numeric(adata.obs["y_centroid"], errors="coerce").values,
            },
            index=adata.obs_names,
        )

    if "spatial" in adata.obsm:
        spatial = np.asarray(adata.obsm["spatial"])
        if spatial.shape[1] < 2:
            raise ValueError("adata.obsm['spatial'] exists but has fewer than 2 columns.")

        return pd.DataFrame(
            {
                "x": spatial[:, 0],
                "y": spatial[:, 1],
            },
            index=adata.obs_names,
        )

    raise ValueError(
        "No spatial coordinates found. Need obs x_centroid/y_centroid or obsm['spatial']."
    )


def normalize_cluster_id(x) -> str:
    """
    Normalize cluster IDs:
        1.0 -> "1"
        "01" -> "1"
        "cluster_1" stays as "cluster_1" unless numeric-like.
    """
    if pd.isna(x):
        return "Unknown"

    x = str(x).strip().replace('"', "").replace("'", "")

    if x in ["", "nan", "None", "NA", "NaN"]:
        return "Unknown"

    try:
        f = float(x)
        if f.is_integer():
            return str(int(f))
    except Exception:
        pass

    return x


def sanitize_for_h5ad(adata: ad.AnnData) -> ad.AnnData:
    """
    Prevent h5ad writing errors from mixed object columns.
    """
    for col in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[col]):
            adata.obs[col] = (
                adata.obs[col]
                .astype(str)
                .replace({"nan": "Unknown", "None": "Unknown"})
            )

    for col in adata.var.columns:
        if pd.api.types.is_object_dtype(adata.var[col]):
            adata.var[col] = (
                adata.var[col]
                .astype(str)
                .replace({"nan": "Unknown", "None": "Unknown"})
            )

    return adata


def get_mapping(inhibitory_mode: str) -> dict[str, str]:
    """
    Choose which PI/manual mapping version to use.
    """
    if inhibitory_mode == "plot_clusters":
        return BANKSY_TO_CELLTYPE_PLOT_CLUSTERS

    if inhibitory_mode == "marker_heatmap":
        return BANKSY_TO_CELLTYPE_MARKER_HEATMAP

    if inhibitory_mode == "explore_subtypes":
        return BANKSY_TO_CELLTYPE_EXPLORE_SUBTYPES

    raise ValueError(f"Unknown inhibitory mode: {inhibitory_mode}")


# =============================================================================
# BANKSY CSV reading and matching
# =============================================================================

def read_banksy_csv_robust(
    banksy_csv: Path,
    adata_obs_names: pd.Index,
    allow_order_match: bool,
) -> pd.DataFrame:
    """
    Read PI BANKSY CSV robustly.

    PI R code usually saves:
        clusts <- cbind(colData(spe_joint)[, cnm], rownames(colData(spe_joint)))
        write.csv(clusts, "banksy_clustering_lambda0.1_res0.7.csv")

    Typical CSV may look like:
        X,V1,V2
        1,6,Br1113_1
        2,6,Br1113_2

    where:
        V1 = BANKSY cluster
        V2 = cell ID

    This function detects the cell ID column by overlap with adata.obs_names.
    Then it detects the cluster column.
    """
    if not banksy_csv.exists():
        raise FileNotFoundError(f"BANKSY CSV not found: {banksy_csv}")

    raw = pd.read_csv(banksy_csv)
    raw.columns = [str(c) for c in raw.columns]

    print("Raw BANKSY CSV shape:", raw.shape)
    print("Raw BANKSY CSV columns:", list(raw.columns))
    print(raw.head())

    adata_set = set(adata_obs_names.astype(str))

    # Step 1: detect cell ID column by maximum overlap with AnnData cell names.
    overlap_by_col = {}
    for col in raw.columns:
        vals = raw[col].astype(str)
        overlap_by_col[col] = int(vals.isin(adata_set).sum())

    cell_col = max(overlap_by_col, key=overlap_by_col.get)
    max_overlap = overlap_by_col[cell_col]

    print("Column overlap with h5ad cells:")
    for col, n in overlap_by_col.items():
        print(f"  {col}: {n}")

    if max_overlap == 0:
        if not allow_order_match:
            raise ValueError(
                "Could not identify a cell ID column in BANKSY CSV. "
                "No column overlaps with adata.obs_names. "
                "Use --allow-order-match only if the CSV row order exactly matches h5ad cell order."
            )

        if raw.shape[0] != len(adata_obs_names):
            raise ValueError(
                "Order-match requested, but BANKSY CSV row count does not match h5ad cells."
            )

        print(
            "WARNING: No cell ID column detected. "
            "Using row order to match BANKSY clusters to h5ad cells."
        )
        cell_ids = adata_obs_names.astype(str).to_numpy()
        cell_col = "__cell_id_from_h5ad_order__"
        raw[cell_col] = cell_ids
    else:
        print(f"Detected cell ID column: {cell_col}")

    # Step 2: detect cluster column.
    # Prefer columns named V1, Banksy, cluster, or the non-cell column with numeric cluster-like values.
    candidate_cols = [c for c in raw.columns if c != cell_col]

    priority_names = ["Banksy", "banksy", "cluster", "Cluster", "V1"]
    cluster_col = None

    for name in priority_names:
        if name in candidate_cols:
            cluster_col = name
            break

    if cluster_col is None:
        scores = {}
        for col in candidate_cols:
            vals = raw[col].map(normalize_cluster_id)
            unique_n = vals.nunique(dropna=True)
            non_unknown = (vals != "Unknown").sum()

            # Cluster column should have fewer unique values than cell column.
            # It should also not be just a row index.
            scores[col] = {
                "unique_n": unique_n,
                "non_unknown": non_unknown,
                "overlap": overlap_by_col.get(col, 0),
            }

        # Choose the column with small unique_n and low cell overlap.
        sorted_candidates = sorted(
            scores.keys(),
            key=lambda c: (
                scores[c]["overlap"],
                scores[c]["unique_n"],
                -scores[c]["non_unknown"],
            ),
        )

        cluster_col = sorted_candidates[0]

    print(f"Detected BANKSY cluster column: {cluster_col}")

    banksy = pd.DataFrame(
        {
            "cell_id": raw[cell_col].astype(str).values,
            "Banksy": raw[cluster_col].map(normalize_cluster_id).values,
        }
    )

    banksy = banksy[banksy["cell_id"].notna()].copy()
    banksy = banksy[banksy["cell_id"].astype(str).str.len() > 0].copy()

    if banksy["cell_id"].duplicated().any():
        dup_n = int(banksy["cell_id"].duplicated().sum())
        print(f"WARNING: duplicated cell IDs in BANKSY CSV: {dup_n}. Keeping first.")
        banksy = banksy.drop_duplicates("cell_id", keep="first")

    banksy = banksy.set_index("cell_id")

    return banksy


def add_banksy_celltype_annotations(
    adata: ad.AnnData,
    banksy: pd.DataFrame,
    cluster_to_celltype: dict[str, str],
    keep_unlabeled: bool,
    strict: bool,
    inhibitory_mode: str,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Match BANKSY clusters to AnnData cells and add cell-type annotations.
    """
    section("Matching BANKSY clusters to AnnData cells")

    adata.obs_names = adata.obs_names.astype(str)

    adata_cells = pd.Index(adata.obs_names.astype(str))
    banksy_cells = pd.Index(banksy.index.astype(str))

    banksy_set = set(banksy_cells)
    adata_set = set(adata_cells)

    common_cells = [cell for cell in adata_cells if cell in banksy_set]
    missing_in_banksy = [cell for cell in adata_cells if cell not in banksy_set]
    extra_in_banksy = [cell for cell in banksy_cells if cell not in adata_set]

    match_report = pd.DataFrame(
        [
            {"metric": "n_cells_in_h5ad", "value": len(adata_cells)},
            {"metric": "n_cells_in_banksy_csv", "value": len(banksy_cells)},
            {"metric": "n_common_cells", "value": len(common_cells)},
            {"metric": "n_h5ad_cells_missing_banksy", "value": len(missing_in_banksy)},
            {"metric": "n_extra_banksy_cells_not_in_h5ad", "value": len(extra_in_banksy)},
            {
                "metric": "n_BrNum_in_h5ad_before_matching",
                "value": adata.obs["BrNum"].astype(str).nunique(),
            },
            {"metric": "inhibitory_mapping_mode", "value": inhibitory_mode},
        ]
    )

    match_report.to_csv(TABLE_DIR / "xenium_N24_banksy_match_report.csv", index=False)

    pd.Series(missing_in_banksy, name="cell_id").to_csv(
        TABLE_DIR / "xenium_N24_cells_missing_banksy_clusters.csv",
        index=False,
    )

    pd.Series(extra_in_banksy, name="cell_id").to_csv(
        TABLE_DIR / "xenium_N24_extra_banksy_clusters_not_in_h5ad.csv",
        index=False,
    )

    print(match_report.to_string(index=False))

    if strict and (len(missing_in_banksy) > 0 or len(extra_in_banksy) > 0):
        raise ValueError(
            "Strict mode failed: h5ad cells and BANKSY CSV cells do not match exactly. "
            "Check xenium_N24_banksy_match_report.csv."
        )

    if len(common_cells) == 0:
        raise ValueError(
            "No matching cell IDs between h5ad and BANKSY CSV. "
            "Check whether the BANKSY CSV contains cell IDs like Br1113_1."
        )

    if keep_unlabeled:
        print("Keeping all h5ad cells. Unlabeled cells will be marked Unknown.")
        banksy_ordered = banksy.reindex(adata.obs_names)
    else:
        if len(missing_in_banksy) > 0:
            print(f"WARNING: dropping {len(missing_in_banksy)} h5ad cells without BANKSY clusters.")
        adata = adata[common_cells, :].copy()
        banksy_ordered = banksy.loc[adata.obs_names]

    adata.obs["Banksy"] = (
        banksy_ordered["Banksy"]
        .fillna("Unknown")
        .map(normalize_cluster_id)
        .values
    )

    adata.obs["cell_type_annotation_banksy"] = (
        adata.obs["Banksy"]
        .map(cluster_to_celltype)
        .fillna("Unknown")
        .astype(str)
    )

    # Current final cell type column for downstream imputation.
    # Later, if supervised snRNA-seq transfer becomes available, you can replace or compare this.
    adata.obs["cell_type_annotation"] = adata.obs["cell_type_annotation_banksy"].astype(str)

    adata.obs["cell_type_annotation_source"] = (
        "BANKSY_manual_mapping_lambda0.1_res0.7"
    )

    adata.obs["cell_type_annotation_stage"] = (
        "03_annotation_celltype_banksy"
    )

    adata.obs["cell_type_annotation_confidence"] = (
        "first_pass_Banksy_manual_mapping_not_snRNA_transfer"
    )

    # Category ordering
    banksy_categories = sorted(
        [x for x in adata.obs["Banksy"].astype(str).unique() if x != "Unknown"],
        key=lambda x: int(x) if str(x).isdigit() else 9999,
    )
    if "Unknown" in adata.obs["Banksy"].astype(str).unique():
        banksy_categories.append("Unknown")

    adata.obs["Banksy"] = pd.Categorical(
        adata.obs["Banksy"],
        categories=banksy_categories,
        ordered=True,
    )

    adata.obs["cell_type_annotation_banksy"] = pd.Categorical(
        adata.obs["cell_type_annotation_banksy"],
        categories=CELLTYPE_ORDER,
        ordered=True,
    )

    adata.obs["cell_type_annotation"] = pd.Categorical(
        adata.obs["cell_type_annotation"],
        categories=CELLTYPE_ORDER,
        ordered=True,
    )

    adata.uns["cell_type_annotation_source"] = "BANKSY_manual_mapping_lambda0.1_res0.7"
    adata.uns["cell_type_annotation_mapping"] = cluster_to_celltype
    adata.uns["cell_type_annotation_order"] = CELLTYPE_ORDER
    adata.uns["banksy_lambda"] = 0.1
    adata.uns["banksy_resolution"] = 0.7
    adata.uns["banksy_inhibitory_mapping_mode"] = inhibitory_mode

    return adata, match_report


# =============================================================================
# CSV outputs
# =============================================================================

def save_csv_outputs(
    adata: ad.AnnData,
    cluster_to_celltype: dict[str, str],
    inhibitory_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    section("Saving CSV outputs")

    xy = get_spatial_xy(adata)

    annotation_cols = [
        "BrNum",
        "Banksy",
        "cell_type_annotation_banksy",
        "cell_type_annotation",
        "cell_type_annotation_source",
        "cell_type_annotation_confidence",
    ]

    optional_cols = [
        "layer_annotation",
        "predictions_smooth",
        "Dx",
        "Sex",
        "Age",
        "PNN",
        "CaptureArea",
        "tear",
        "qc_total_counts",
        "qc_detected_genes",
        "total_counts",
        "detected",
    ]

    keep_cols = annotation_cols + [c for c in optional_cols if c in adata.obs.columns]

    annotation_table = adata.obs[keep_cols].copy()
    annotation_table.insert(0, "cell_id", adata.obs_names)
    annotation_table["x"] = xy.loc[adata.obs_names, "x"].values
    annotation_table["y"] = xy.loc[adata.obs_names, "y"].values

    annotation_table.to_csv(
        TABLE_DIR / "xenium_N24_banksy_celltype_annotations.csv",
        index=False,
    )

    mapping = pd.DataFrame(
        {
            "Banksy": list(cluster_to_celltype.keys()),
            "cell_type_annotation": list(cluster_to_celltype.values()),
            "inhibitory_mapping_mode": inhibitory_mode,
            "source": "BANKSY_manual_mapping_lambda0.1_res0.7",
        }
    )

    mapping.to_csv(
        TABLE_DIR / "xenium_N24_banksy_cluster_to_celltype_mapping.csv",
        index=False,
    )

    cluster_counts = (
        adata.obs
        .groupby(["BrNum", "Banksy"], observed=False)
        .size()
        .reset_index(name="n_cells")
    )

    cluster_counts.to_csv(
        TABLE_DIR / "xenium_N24_banksy_cluster_counts_by_BrNum.csv",
        index=False,
    )

    cluster_props = cluster_counts.copy()
    cluster_props["proportion"] = cluster_props.groupby("BrNum")["n_cells"].transform(
        lambda x: x / x.sum() if x.sum() > 0 else np.nan
    )

    cluster_props.to_csv(
        TABLE_DIR / "xenium_N24_banksy_cluster_proportions_by_BrNum.csv",
        index=False,
    )

    celltype_counts = (
        adata.obs
        .groupby(["BrNum", "cell_type_annotation"], observed=False)
        .size()
        .reset_index(name="n_cells")
    )

    celltype_counts.to_csv(
        TABLE_DIR / "xenium_N24_celltype_counts_by_BrNum.csv",
        index=False,
    )

    celltype_props = celltype_counts.copy()
    celltype_props["proportion"] = celltype_props.groupby("BrNum")["n_cells"].transform(
        lambda x: x / x.sum() if x.sum() > 0 else np.nan
    )

    celltype_props.to_csv(
        TABLE_DIR / "xenium_N24_celltype_proportions_by_BrNum.csv",
        index=False,
    )

    summary = (
        adata.obs
        .groupby("BrNum")
        .agg(
            n_cells=("cell_type_annotation", "size"),
            n_banksy_clusters=("Banksy", lambda x: x.astype(str).nunique()),
            n_cell_types=("cell_type_annotation", lambda x: x.astype(str).nunique()),
            most_common_cell_type=(
                "cell_type_annotation",
                lambda x: x.astype(str).value_counts().idxmax(),
            ),
            most_common_cell_type_n=(
                "cell_type_annotation",
                lambda x: int(x.astype(str).value_counts().max()),
            ),
        )
        .reset_index()
    )

    summary["most_common_cell_type_fraction"] = (
        summary["most_common_cell_type_n"] / summary["n_cells"]
    )

    summary.to_csv(
        TABLE_DIR / "xenium_N24_celltype_summary_by_BrNum.csv",
        index=False,
    )

    print("Saved:", TABLE_DIR / "xenium_N24_banksy_celltype_annotations.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_banksy_cluster_to_celltype_mapping.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_celltype_counts_by_BrNum.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_celltype_proportions_by_BrNum.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_celltype_summary_by_BrNum.csv")

    return cluster_counts, cluster_props, celltype_counts, celltype_props


# =============================================================================
# Plotting helpers
# =============================================================================

def maybe_downsample(df: pd.DataFrame, max_cells: int) -> pd.DataFrame:
    """
    max_cells <= 0 means plot all cells.
    """
    if max_cells is None or max_cells <= 0:
        return df

    if df.shape[0] > max_cells:
        return df.sample(max_cells, random_state=0)

    return df


def present_celltypes(adata: ad.AnnData) -> list[str]:
    found = adata.obs["cell_type_annotation"].astype(str).unique().tolist()
    ordered = [x for x in CELLTYPE_ORDER if x in found]

    extras = [x for x in found if x not in ordered]
    return ordered + sorted(extras)


def present_clusters(adata: ad.AnnData) -> list[str]:
    found = adata.obs["Banksy"].astype(str).unique().tolist()
    ordered = sorted(
        [x for x in found if x != "Unknown"],
        key=lambda x: int(x) if str(x).isdigit() else 9999,
    )

    if "Unknown" in found:
        ordered.append("Unknown")

    return ordered


def cluster_color_map(clusters: list[str]) -> dict[str, tuple]:
    cmap = plt.get_cmap("tab20")
    return {cluster: cmap(i % 20) for i, cluster in enumerate(clusters)}


def plot_spatial_panel(
    adata: ad.AnnData,
    value_col: str,
    output_file: Path,
    title: str,
    max_cells_per_sample: int,
    color_map: dict | None = None,
    order: list[str] | None = None,
) -> None:
    section(f"Creating panel plot: {title}")

    xy = get_spatial_xy(adata)

    plot_df = adata.obs[["BrNum", value_col]].copy()
    plot_df["x"] = xy.loc[adata.obs_names, "x"].values
    plot_df["y"] = xy.loc[adata.obs_names, "y"].values

    brnums = sorted(plot_df["BrNum"].astype(str).unique())

    if order is None:
        order = sorted(plot_df[value_col].astype(str).unique())

    if color_map is None:
        color_map = cluster_color_map(order)

    n = len(brnums)
    ncols = 4
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.2 * nrows),
        squeeze=False,
    )
    axes = axes.ravel()

    for ax, br in zip(axes, brnums):
        sub = plot_df[plot_df["BrNum"].astype(str).eq(br)].copy()
        sub = maybe_downsample(sub, max_cells=max_cells_per_sample)

        for val in order:
            one = sub[sub[value_col].astype(str).eq(str(val))]
            if one.empty:
                continue

            ax.scatter(
                one["x"],
                one["y"],
                s=0.25,
                alpha=0.8,
                color=color_map.get(val, "#999999"),
                linewidths=0,
                label=val,
            )

        ax.set_title(br)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")

    for ax in axes[len(brnums):]:
        ax.axis("off")

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=6,
            color=color_map.get(val, "#999999"),
            label=val,
        )
        for val in order
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(len(order), 8),
        frameon=False,
    )

    fig.suptitle(title, y=0.995)
    plt.tight_layout(rect=[0, 0.06, 1, 0.98])
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def plot_spatial_per_sample(
    adata: ad.AnnData,
    value_col: str,
    output_dir: Path,
    file_suffix: str,
    title_prefix: str,
    max_cells_per_sample: int,
    color_map: dict | None = None,
    order: list[str] | None = None,
) -> None:
    section(f"Creating separate per-BrNum maps: {title_prefix}")

    output_dir.mkdir(parents=True, exist_ok=True)

    xy = get_spatial_xy(adata)

    plot_df = adata.obs[["BrNum", value_col]].copy()
    plot_df["x"] = xy.loc[adata.obs_names, "x"].values
    plot_df["y"] = xy.loc[adata.obs_names, "y"].values

    brnums = sorted(plot_df["BrNum"].astype(str).unique())

    if order is None:
        order = sorted(plot_df[value_col].astype(str).unique())

    if color_map is None:
        color_map = cluster_color_map(order)

    for br in brnums:
        sub_all = plot_df[plot_df["BrNum"].astype(str).eq(br)].copy()
        n_total = sub_all.shape[0]

        sub = maybe_downsample(sub_all, max_cells=max_cells_per_sample)
        n_plotted = sub.shape[0]

        fig, ax = plt.subplots(figsize=(7.5, 7.5))

        for val in order:
            one = sub[sub[value_col].astype(str).eq(str(val))]
            if one.empty:
                continue

            ax.scatter(
                one["x"],
                one["y"],
                s=0.35,
                alpha=0.85,
                color=color_map.get(val, "#999999"),
                linewidths=0,
                label=val,
            )

        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")

        ax.set_title(
            f"{br} {title_prefix}\n"
            f"Plotted cells: {n_plotted:,} / Total cells: {n_total:,}"
        )

        ax.legend(
            title=value_col,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            markerscale=7,
            fontsize=8,
        )

        plt.tight_layout()

        out_file = output_dir / f"{br}_{file_suffix}.png"
        plt.savefig(out_file, dpi=300)
        plt.close()

        print("Saved:", out_file)


def plot_stacked_proportions(
    proportions: pd.DataFrame,
    category_col: str,
    value_col: str,
    output_file: Path,
    title: str,
    color_map: dict | None = None,
    order: list[str] | None = None,
) -> None:
    section(f"Creating proportion plot: {title}")

    prop = proportions.copy()

    if order is not None:
        prop[category_col] = pd.Categorical(
            prop[category_col].astype(str),
            categories=order,
            ordered=True,
        )

    pivot = prop.pivot_table(
        index="BrNum",
        columns=category_col,
        values=value_col,
        fill_value=0,
        observed=False,
    )

    if order is not None:
        ordered_cols = [x for x in order if x in pivot.columns]
        pivot = pivot[ordered_cols]

    colors = None
    if color_map is not None:
        colors = [color_map.get(c, "#999999") for c in pivot.columns]

    ax = pivot.plot(
        kind="bar",
        stacked=True,
        figsize=(14, 6),
        width=0.85,
        color=colors,
    )

    ax.set_ylabel("Proportion of cells")
    ax.set_xlabel("BrNum")
    ax.set_title(title)
    ax.legend(
        title=category_col,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def plot_overall_celltype_distribution(
    adata: ad.AnnData,
    output_file: Path,
) -> None:
    section("Creating overall cell-type distribution plot")

    vals = adata.obs["cell_type_annotation"].astype(str)
    counts = vals.value_counts()

    ordered = [x for x in CELLTYPE_ORDER if x in counts.index]
    counts = counts.loc[ordered]

    plt.figure(figsize=(10, 5))
    plt.bar(
        counts.index,
        counts.values,
        color=[CELLTYPE_COLORS.get(x, "#999999") for x in counts.index],
    )
    plt.ylabel("Number of cells")
    plt.xlabel("BANKSY-based cell type")
    plt.title("Overall Xenium BANKSY-based cell-type distribution")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def save_plots(
    adata: ad.AnnData,
    cluster_props: pd.DataFrame,
    celltype_props: pd.DataFrame,
    max_cells_per_sample: int,
) -> None:
    celltype_order = present_celltypes(adata)
    cluster_order = present_clusters(adata)
    cluster_colors = cluster_color_map(cluster_order)

    plot_spatial_panel(
        adata=adata,
        value_col="cell_type_annotation",
        output_file=FIGURE_DIR / "xenium_N24_banksy_celltype_maps_panel.png",
        title="Xenium BANKSY-based cell-type annotation by donor",
        max_cells_per_sample=max_cells_per_sample,
        color_map=CELLTYPE_COLORS,
        order=celltype_order,
    )

    plot_spatial_per_sample(
        adata=adata,
        value_col="cell_type_annotation",
        output_dir=PER_SAMPLE_CELLTYPE_DIR,
        file_suffix="banksy_celltype_map",
        title_prefix="BANKSY-based cell-type annotation",
        max_cells_per_sample=max_cells_per_sample,
        color_map=CELLTYPE_COLORS,
        order=celltype_order,
    )

    plot_spatial_panel(
        adata=adata,
        value_col="Banksy",
        output_file=FIGURE_DIR / "xenium_N24_banksy_cluster_maps_panel.png",
        title="Xenium BANKSY clusters by donor",
        max_cells_per_sample=max_cells_per_sample,
        color_map=cluster_colors,
        order=cluster_order,
    )

    plot_spatial_per_sample(
        adata=adata,
        value_col="Banksy",
        output_dir=PER_SAMPLE_CLUSTER_DIR,
        file_suffix="banksy_cluster_map",
        title_prefix="BANKSY cluster map",
        max_cells_per_sample=max_cells_per_sample,
        color_map=cluster_colors,
        order=cluster_order,
    )

    plot_stacked_proportions(
        proportions=celltype_props,
        category_col="cell_type_annotation",
        value_col="proportion",
        output_file=FIGURE_DIR / "xenium_N24_banksy_celltype_proportions_by_BrNum.png",
        title="BANKSY-based cell-type composition by donor",
        color_map=CELLTYPE_COLORS,
        order=celltype_order,
    )

    plot_stacked_proportions(
        proportions=cluster_props,
        category_col="Banksy",
        value_col="proportion",
        output_file=FIGURE_DIR / "xenium_N24_banksy_cluster_proportions_by_BrNum.png",
        title="BANKSY cluster composition by donor",
        color_map=cluster_colors,
        order=cluster_order,
    )

    plot_overall_celltype_distribution(
        adata=adata,
        output_file=FIGURE_DIR / "xenium_N24_banksy_celltype_overall_distribution.png",
    )


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-h5ad",
        default=str(DEFAULT_INPUT_H5AD),
        help="Input Xenium layer-annotated h5ad.",
    )

    parser.add_argument(
        "--banksy-csv",
        default=str(DEFAULT_BANKSY_CSV),
        help="BANKSY clustering CSV, usually banksy_clustering_lambda0.1_res0.7.csv.",
    )

    parser.add_argument(
        "--output-h5ad",
        default=str(DEFAULT_OUTPUT_H5AD),
        help="Output Xenium layer + BANKSY cell-type annotated h5ad.",
    )

    parser.add_argument(
        "--inhibitory-mode",
        default="plot_clusters",
        choices=["plot_clusters", "marker_heatmap", "explore_subtypes"],
        help=(
            "How to map BANKSY clusters 9 and 12. "
            "plot_clusters: 12=CGE, 9=MGE. "
            "marker_heatmap: 12=MGE, 9=CGE. "
            "explore_subtypes: 12=In: VIP/LAMP5, 9=In: SST/PVALB."
        ),
    )

    parser.add_argument(
        "--max-cells-per-sample",
        type=int,
        default=0,
        help=(
            "Maximum cells plotted per donor. "
            "Use 0 to plot all cells. Use 8000 for faster test plots."
        ),
    )

    parser.add_argument(
        "--keep-unlabeled",
        action="store_true",
        help="Keep h5ad cells without BANKSY clusters and mark them Unknown. Default drops unlabeled cells.",
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Require exact match between h5ad cell IDs and BANKSY CSV cell IDs.",
    )

    parser.add_argument(
        "--allow-order-match",
        action="store_true",
        help=(
            "Use only if BANKSY CSV has no cell_id column but rows are in exactly the same "
            "order as h5ad cells."
        ),
    )

    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plot generation.",
    )

    parser.add_argument(
        "--compression",
        default="lzf",
        choices=["lzf", "gzip", "none"],
        help="h5ad compression.",
    )

    args = parser.parse_args()

    input_h5ad = Path(args.input_h5ad)
    banksy_csv = Path(args.banksy_csv)
    output_h5ad = Path(args.output_h5ad)

    section("Stage 03: BANKSY-based cell-type annotation")

    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_h5ad}")

    if not banksy_csv.exists():
        raise FileNotFoundError(f"BANKSY CSV not found: {banksy_csv}")

    print("Input h5ad:", input_h5ad)
    print("BANKSY CSV:", banksy_csv)
    print("Output h5ad:", output_h5ad)
    print("Inhibitory mapping mode:", args.inhibitory_mode)

    section("Loading h5ad")

    adata = ad.read_h5ad(input_h5ad)
    adata.obs_names = adata.obs_names.astype(str)
    adata = ensure_brnum(adata)

    print(adata)

    section("Loading BANKSY CSV")

    banksy = read_banksy_csv_robust(
        banksy_csv=banksy_csv,
        adata_obs_names=pd.Index(adata.obs_names.astype(str)),
        allow_order_match=args.allow_order_match,
    )

    print("Processed BANKSY table shape:", banksy.shape)
    print(banksy.head())

    cluster_to_celltype = get_mapping(args.inhibitory_mode)

    adata, match_report = add_banksy_celltype_annotations(
        adata=adata,
        banksy=banksy,
        cluster_to_celltype=cluster_to_celltype,
        keep_unlabeled=args.keep_unlabeled,
        strict=args.strict,
        inhibitory_mode=args.inhibitory_mode,
    )

    cluster_counts, cluster_props, celltype_counts, celltype_props = save_csv_outputs(
        adata=adata,
        cluster_to_celltype=cluster_to_celltype,
        inhibitory_mode=args.inhibitory_mode,
    )

    if not args.skip_plots:
        save_plots(
            adata=adata,
            cluster_props=cluster_props,
            celltype_props=celltype_props,
            max_cells_per_sample=args.max_cells_per_sample,
        )
    else:
        print("Skipping plots.")

    section("Saving h5ad")

    adata = sanitize_for_h5ad(adata)

    if args.compression == "none":
        adata.write_h5ad(output_h5ad)
    else:
        adata.write_h5ad(output_h5ad, compression=args.compression)

    print("Saved h5ad:", output_h5ad)

    section("Final summary")

    print("Cells in output:", adata.n_obs)
    print("Genes/features in output:", adata.n_vars)
    print("Donors:", adata.obs["BrNum"].nunique())

    print("\nBANKSY cluster counts:")
    print(adata.obs["Banksy"].astype(str).value_counts().sort_index().to_string())

    print("\nBANKSY-based cell-type counts:")
    print(adata.obs["cell_type_annotation"].astype(str).value_counts().to_string())

    print("\nMain outputs:")
    print("H5AD:")
    print(" ", output_h5ad)
    print("CSV folder:")
    print(" ", TABLE_DIR)
    print("Figure folder:")
    print(" ", FIGURE_DIR)
    print("Per-sample cell-type PNG folder:")
    print(" ", PER_SAMPLE_CELLTYPE_DIR)
    print("Per-sample cluster PNG folder:")
    print(" ", PER_SAMPLE_CLUSTER_DIR)


if __name__ == "__main__":
    main()