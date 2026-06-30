#!/usr/bin/env python3

"""
01_layerwise_annotation.py

Stage 03 Annotation Task 1:
Layer-wise annotation for Xenium N24.

Purpose
-------
Take the PI smoothed spatial-domain label CSV:

    label_transfer_N24_k50_smoothed_labels.csv

and add those labels to:

    xenium_N24_imputation_ready.h5ad

Then save:

    xenium_N24_layer_annotated.h5ad

The key biological mapping is:

    spd07 -> L1/M
    spd06 -> L2/3
    spd02 -> L3/4
    spd05 -> L5
    spd03 -> L6
    spd01 -> WMtz
    spd04 -> WM

Main output h5ad
----------------
data/processed/xenium/xenium_N24_layer_annotated.h5ad

Main CSV outputs
----------------
outputs/xenium_branch/03_annotation/layerwise/tables/
    xenium_N24_layer_annotations.csv
    xenium_N24_layer_counts_by_BrNum.csv
    xenium_N24_layer_proportions_by_BrNum.csv
    xenium_N24_layer_mapping.csv
    xenium_N24_layer_match_report.csv
    xenium_N24_layer_summary.csv

Main plot outputs
-----------------
outputs/xenium_branch/03_annotation/layerwise/figures/
    xenium_N24_layer_maps_panel.png
    xenium_N24_layer_proportions_by_BrNum.png
    xenium_N24_layer_counts_by_BrNum.png
    xenium_N24_layer_overall_distribution.png

Per-sample plot outputs
-----------------------
outputs/xenium_branch/03_annotation/layerwise/figures/per_sample_layer_maps/
    Br1113_layer_map.png
    Br1139_layer_map.png
    ...
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
    / "xenium_N24_imputation_ready.h5ad"
)

DEFAULT_LABEL_CSV = (
    PROJECT_ROOT
    / "data"
    / "metadata"
    / "layer_annotations"
    / "label_transfer_N24_k50_smoothed_labels.csv"
)

DEFAULT_OUTPUT_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "xenium"
    / "xenium_N24_layer_annotated.h5ad"
)

OUT_DIR = PROJECT_ROOT / "outputs" / "xenium_branch" / "03_annotation" / "layerwise"
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"
PER_SAMPLE_FIGURE_DIR = FIGURE_DIR / "per_sample_layer_maps"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
PER_SAMPLE_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
DEFAULT_OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Layer mapping
# =============================================================================

SPD_TO_LAYER = {
    "spd07": "L1/M",
    "spd06": "L2/3",
    "spd02": "L3/4",
    "spd05": "L5",
    "spd03": "L6",
    "spd01": "WMtz",
    "spd04": "WM",
}

LAYER_ORDER = ["L1/M", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM"]
SPD_ORDER = ["spd07", "spd06", "spd02", "spd05", "spd03", "spd01", "spd04"]

# PI-style colors from the R layer plotting logic.
LAYER_COLORS = {
    "L1/M": "#FEAF16",
    "L2/3": "#3283FE",
    "L3/4": "#E4E1E3",
    "L5": "#16FF32",
    "L6": "#F6222E",
    "WMtz": "#5A5156",
    "WM": "#FE00FA",
    "Unknown": "#999999",
}


# =============================================================================
# Helpers
# =============================================================================

def section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def normalize_spd_label(x) -> str:
    """
    Normalize labels:
        spd6  -> spd06
        spd06 -> spd06
        SPD06 -> spd06
    """
    if pd.isna(x):
        return "Unknown"

    x = str(x).strip().replace('"', "").replace("'", "").lower()

    if x in ["", "nan", "none", "na"]:
        return "Unknown"

    if x.startswith("spd"):
        num = x.replace("spd", "")
        if num.isdigit():
            return f"spd{int(num):02d}"

    return x


def read_layer_label_csv(label_csv: Path) -> pd.DataFrame:
    """
    Read the PI smoothed-label CSV robustly.

    Example expected format:

        ,colData(spe_all)[, "predictions_smooth"]
        Br1113_1,spd06
        Br1113_2,spd07

    or:

        cell_id,predictions_smooth
        Br1113_1,spd06
    """
    if not label_csv.exists():
        raise FileNotFoundError(f"Label CSV not found: {label_csv}")

    labels = pd.read_csv(label_csv, sep=None, engine="python", index_col=0)

    if labels.shape[1] == 0:
        raise ValueError(f"No label column found in: {label_csv}")

    candidate_cols = [
        c for c in labels.columns
        if "predictions_smooth" in str(c) or "prediction" in str(c).lower()
    ]

    label_col = candidate_cols[0] if len(candidate_cols) > 0 else labels.columns[0]

    labels = labels[[label_col]].copy()
    labels.columns = ["predictions_smooth"]

    labels.index = labels.index.astype(str)
    labels["predictions_smooth"] = labels["predictions_smooth"].map(normalize_spd_label)

    if labels.index.duplicated().any():
        dup_n = int(labels.index.duplicated().sum())
        print(f"WARNING: duplicated cell IDs in label CSV: {dup_n}. Keeping first.")
        labels = labels[~labels.index.duplicated(keep="first")]

    return labels


def ensure_brnum(adata: ad.AnnData) -> ad.AnnData:
    """
    Ensure adata.obs['BrNum'] exists.
    """
    if "BrNum" not in adata.obs.columns:
        print("BrNum not found. Deriving BrNum from cell IDs.")
        adata.obs["BrNum"] = pd.Index(adata.obs_names.astype(str)).str.split("_").str[0]
    else:
        adata.obs["BrNum"] = adata.obs["BrNum"].astype(str)

    return adata


def get_spatial_xy(adata: ad.AnnData) -> pd.DataFrame:
    """
    Return x/y coordinates from:
        obs['x_centroid'], obs['y_centroid']
    or:
        obsm['spatial']
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


# =============================================================================
# Annotation
# =============================================================================

def add_layer_annotations(
    adata: ad.AnnData,
    labels: pd.DataFrame,
    keep_unlabeled: bool,
    strict: bool,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Match labels to adata.obs_names and add:
        predictions_smooth
        layer_annotation
        layer_annotation_source
    """
    section("Matching label CSV to AnnData cells")

    adata.obs_names = adata.obs_names.astype(str)

    adata_cells = pd.Index(adata.obs_names.astype(str))
    label_cells = pd.Index(labels.index.astype(str))
    label_set = set(label_cells)
    adata_set = set(adata_cells)

    common_cells = [cell for cell in adata_cells if cell in label_set]
    missing_in_labels = [cell for cell in adata_cells if cell not in label_set]
    extra_in_labels = [cell for cell in label_cells if cell not in adata_set]

    match_report = pd.DataFrame(
        [
            {"metric": "n_cells_in_h5ad", "value": len(adata_cells)},
            {"metric": "n_cells_in_label_csv", "value": len(label_cells)},
            {"metric": "n_common_cells", "value": len(common_cells)},
            {"metric": "n_h5ad_cells_missing_labels", "value": len(missing_in_labels)},
            {"metric": "n_extra_labels_not_in_h5ad", "value": len(extra_in_labels)},
            {
                "metric": "n_BrNum_in_h5ad_before_matching",
                "value": adata.obs["BrNum"].astype(str).nunique(),
            },
        ]
    )

    match_report.to_csv(TABLE_DIR / "xenium_N24_layer_match_report.csv", index=False)

    pd.Series(missing_in_labels, name="cell_id").to_csv(
        TABLE_DIR / "xenium_N24_cells_missing_layer_labels.csv",
        index=False,
    )

    pd.Series(extra_in_labels, name="cell_id").to_csv(
        TABLE_DIR / "xenium_N24_extra_layer_labels_not_in_h5ad.csv",
        index=False,
    )

    print(match_report.to_string(index=False))

    if strict and (len(missing_in_labels) > 0 or len(extra_in_labels) > 0):
        raise ValueError(
            "Strict mode failed: h5ad cells and label CSV cells do not match exactly. "
            "Check xenium_N24_layer_match_report.csv."
        )

    if len(common_cells) == 0:
        raise ValueError(
            "No matching cell IDs between h5ad and label CSV. "
            "Check whether your h5ad cell names look like Br1113_1."
        )

    if keep_unlabeled:
        print("Keeping all h5ad cells. Unlabeled cells will be marked Unknown.")
        labels_ordered = labels.reindex(adata.obs_names)
    else:
        if len(missing_in_labels) > 0:
            print(f"WARNING: dropping {len(missing_in_labels)} h5ad cells without labels.")
        adata = adata[common_cells, :].copy()
        labels_ordered = labels.loc[adata.obs_names]

    adata.obs["predictions_smooth"] = (
        labels_ordered["predictions_smooth"]
        .fillna("Unknown")
        .map(normalize_spd_label)
        .values
    )

    adata.obs["layer_annotation"] = (
        adata.obs["predictions_smooth"]
        .map(SPD_TO_LAYER)
        .fillna("Unknown")
        .astype(str)
    )

    adata.obs["predictions_smooth"] = pd.Categorical(
        adata.obs["predictions_smooth"],
        categories=SPD_ORDER + ["Unknown"],
        ordered=True,
    )

    adata.obs["layer_annotation"] = pd.Categorical(
        adata.obs["layer_annotation"],
        categories=LAYER_ORDER + ["Unknown"],
        ordered=True,
    )

    adata.obs["layer_annotation_source"] = "PI_label_transfer_N24_k50_smoothed_labels"
    adata.obs["layer_annotation_stage"] = "03_annotation_layerwise"

    adata.uns["layer_annotation_source"] = "PI_label_transfer_N24_k50_smoothed_labels"
    adata.uns["layer_annotation_mapping"] = SPD_TO_LAYER
    adata.uns["layer_annotation_order"] = LAYER_ORDER

    return adata, match_report


# =============================================================================
# CSV outputs
# =============================================================================

def save_csv_outputs(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    section("Saving CSV outputs")

    xy = get_spatial_xy(adata)

    annotation_cols = [
        "BrNum",
        "predictions_smooth",
        "layer_annotation",
        "layer_annotation_source",
    ]

    optional_cols = [
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
        TABLE_DIR / "xenium_N24_layer_annotations.csv",
        index=False,
    )

    mapping = pd.DataFrame(
        {
            "predictions_smooth": list(SPD_TO_LAYER.keys()),
            "layer_annotation": list(SPD_TO_LAYER.values()),
        }
    )

    mapping.to_csv(TABLE_DIR / "xenium_N24_layer_mapping.csv", index=False)

    counts = (
        adata.obs
        .groupby(["BrNum", "layer_annotation"], observed=False)
        .size()
        .reset_index(name="n_cells")
    )

    counts.to_csv(TABLE_DIR / "xenium_N24_layer_counts_by_BrNum.csv", index=False)

    proportions = counts.copy()
    proportions["proportion"] = proportions.groupby("BrNum")["n_cells"].transform(
        lambda x: x / x.sum() if x.sum() > 0 else np.nan
    )

    proportions.to_csv(
        TABLE_DIR / "xenium_N24_layer_proportions_by_BrNum.csv",
        index=False,
    )

    summary = (
        adata.obs
        .groupby("BrNum")
        .agg(
            n_cells=("layer_annotation", "size"),
            n_layers_detected=("layer_annotation", lambda x: x.astype(str).nunique()),
            most_common_layer=(
                "layer_annotation",
                lambda x: x.astype(str).value_counts().idxmax(),
            ),
            most_common_layer_n=(
                "layer_annotation",
                lambda x: int(x.astype(str).value_counts().max()),
            ),
        )
        .reset_index()
    )

    summary["most_common_layer_fraction"] = (
        summary["most_common_layer_n"] / summary["n_cells"]
    )

    summary.to_csv(TABLE_DIR / "xenium_N24_layer_summary.csv", index=False)

    print("Saved:", TABLE_DIR / "xenium_N24_layer_annotations.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_layer_counts_by_BrNum.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_layer_proportions_by_BrNum.csv")
    print("Saved:", TABLE_DIR / "xenium_N24_layer_summary.csv")

    return counts, proportions, summary


# =============================================================================
# Plotting
# =============================================================================

def layers_present(adata: ad.AnnData) -> list[str]:
    found = adata.obs["layer_annotation"].astype(str).unique().tolist()
    layers = [layer for layer in LAYER_ORDER if layer in found]

    if "Unknown" in found:
        layers.append("Unknown")

    return layers


def maybe_downsample(df: pd.DataFrame, max_cells: int) -> pd.DataFrame:
    """
    max_cells <= 0 means plot all cells.
    """
    if max_cells is None or max_cells <= 0:
        return df

    if df.shape[0] > max_cells:
        return df.sample(max_cells, random_state=0)

    return df


def plot_layer_maps_panel(
    adata: ad.AnnData,
    max_cells_per_sample: int,
    output_file: Path,
) -> None:
    """
    Combined panel with one subplot per BrNum.
    """
    section("Creating combined spatial layer map panel")

    xy = get_spatial_xy(adata)

    plot_df = adata.obs[["BrNum", "layer_annotation"]].copy()
    plot_df["x"] = xy.loc[adata.obs_names, "x"].values
    plot_df["y"] = xy.loc[adata.obs_names, "y"].values

    brnums = sorted(plot_df["BrNum"].astype(str).unique())
    layers = layers_present(adata)

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

        for layer in layers:
            one = sub[sub["layer_annotation"].astype(str).eq(layer)]
            if one.empty:
                continue

            ax.scatter(
                one["x"],
                one["y"],
                s=0.25,
                alpha=0.8,
                label=layer,
                color=LAYER_COLORS.get(layer, "#999999"),
                linewidths=0,
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
            color=LAYER_COLORS.get(layer, "#999999"),
            label=layer,
        )
        for layer in layers
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=min(len(layers), 8),
        frameon=False,
    )

    fig.suptitle("Xenium layer-wise annotation by donor", y=0.995)
    plt.tight_layout(rect=[0, 0.04, 1, 0.98])
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def plot_layer_maps_per_sample(
    adata: ad.AnnData,
    max_cells_per_sample: int,
    output_dir: Path,
) -> None:
    """
    Save one separate PNG per BrNum/sample.
    This is what your PI wants.

    Output example:
        Br1113_layer_map.png
        Br1139_layer_map.png
        Br2039_layer_map.png
    """
    section("Creating separate layer map PNGs for each BrNum")

    output_dir.mkdir(parents=True, exist_ok=True)

    xy = get_spatial_xy(adata)

    plot_df = adata.obs[["BrNum", "layer_annotation"]].copy()
    plot_df["x"] = xy.loc[adata.obs_names, "x"].values
    plot_df["y"] = xy.loc[adata.obs_names, "y"].values

    brnums = sorted(plot_df["BrNum"].astype(str).unique())
    layers = layers_present(adata)

    for br in brnums:
        sub_all = plot_df[plot_df["BrNum"].astype(str).eq(br)].copy()
        n_total = sub_all.shape[0]

        sub = maybe_downsample(sub_all, max_cells=max_cells_per_sample)
        n_plotted = sub.shape[0]

        fig, ax = plt.subplots(figsize=(7.5, 7.5))

        for layer in layers:
            one = sub[sub["layer_annotation"].astype(str).eq(layer)]
            if one.empty:
                continue

            ax.scatter(
                one["x"],
                one["y"],
                s=0.35,
                alpha=0.85,
                label=layer,
                color=LAYER_COLORS.get(layer, "#999999"),
                linewidths=0,
            )

        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal", adjustable="box")

        ax.set_title(
            f"{br} layer-wise annotation\n"
            f"Plotted cells: {n_plotted:,} / Total cells: {n_total:,}"
        )

        ax.legend(
            title="Layer",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            markerscale=7,
        )

        plt.tight_layout()

        out_file = output_dir / f"{br}_layer_map.png"
        plt.savefig(out_file, dpi=300)
        plt.close()

        print("Saved:", out_file)


def plot_layer_proportions(proportions: pd.DataFrame, output_file: Path) -> None:
    section("Creating layer proportion plot")

    prop = proportions.copy()
    prop["layer_annotation"] = pd.Categorical(
        prop["layer_annotation"],
        categories=LAYER_ORDER + ["Unknown"],
        ordered=True,
    )

    pivot = prop.pivot_table(
        index="BrNum",
        columns="layer_annotation",
        values="proportion",
        fill_value=0,
        observed=False,
    )

    ordered_cols = [c for c in LAYER_ORDER + ["Unknown"] if c in pivot.columns]
    pivot = pivot[ordered_cols]

    ax = pivot.plot(
        kind="bar",
        stacked=True,
        figsize=(14, 6),
        width=0.85,
        color=[LAYER_COLORS.get(c, "#999999") for c in pivot.columns],
    )

    ax.set_ylabel("Proportion of cells")
    ax.set_xlabel("BrNum")
    ax.set_title("Layer composition by donor")
    ax.legend(
        title="Layer",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def plot_layer_counts(counts: pd.DataFrame, output_file: Path) -> None:
    section("Creating layer count plot")

    cnt = counts.copy()
    cnt["layer_annotation"] = pd.Categorical(
        cnt["layer_annotation"],
        categories=LAYER_ORDER + ["Unknown"],
        ordered=True,
    )

    pivot = cnt.pivot_table(
        index="BrNum",
        columns="layer_annotation",
        values="n_cells",
        fill_value=0,
        observed=False,
    )

    ordered_cols = [c for c in LAYER_ORDER + ["Unknown"] if c in pivot.columns]
    pivot = pivot[ordered_cols]

    ax = pivot.plot(
        kind="bar",
        stacked=True,
        figsize=(14, 6),
        width=0.85,
        color=[LAYER_COLORS.get(c, "#999999") for c in pivot.columns],
    )

    ax.set_ylabel("Number of cells")
    ax.set_xlabel("BrNum")
    ax.set_title("Layer cell counts by donor")
    ax.legend(
        title="Layer",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )

    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def plot_overall_distribution(adata: ad.AnnData, output_file: Path) -> None:
    section("Creating overall layer distribution plot")

    vals = adata.obs["layer_annotation"].astype(str)
    counts = vals.value_counts()

    ordered = [x for x in LAYER_ORDER + ["Unknown"] if x in counts.index]
    counts = counts.loc[ordered]

    plt.figure(figsize=(8, 5))
    plt.bar(
        counts.index,
        counts.values,
        color=[LAYER_COLORS.get(x, "#999999") for x in counts.index],
    )
    plt.ylabel("Number of cells")
    plt.xlabel("Layer annotation")
    plt.title("Overall Xenium layer distribution")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print("Saved:", output_file)


def save_plots(
    adata: ad.AnnData,
    counts: pd.DataFrame,
    proportions: pd.DataFrame,
    max_cells_per_sample: int,
) -> None:
    """
    Save combined and per-sample plots.
    """
    plot_layer_maps_panel(
        adata,
        max_cells_per_sample=max_cells_per_sample,
        output_file=FIGURE_DIR / "xenium_N24_layer_maps_panel.png",
    )

    plot_layer_maps_per_sample(
        adata,
        max_cells_per_sample=max_cells_per_sample,
        output_dir=PER_SAMPLE_FIGURE_DIR,
    )

    plot_layer_proportions(
        proportions,
        output_file=FIGURE_DIR / "xenium_N24_layer_proportions_by_BrNum.png",
    )

    plot_layer_counts(
        counts,
        output_file=FIGURE_DIR / "xenium_N24_layer_counts_by_BrNum.png",
    )

    plot_overall_distribution(
        adata,
        output_file=FIGURE_DIR / "xenium_N24_layer_overall_distribution.png",
    )


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-h5ad",
        default=str(DEFAULT_INPUT_H5AD),
        help="Input Xenium imputation-ready h5ad.",
    )

    parser.add_argument(
        "--label-csv",
        default=str(DEFAULT_LABEL_CSV),
        help="PI label_transfer_N24_k50_smoothed_labels.csv.",
    )

    parser.add_argument(
        "--output-h5ad",
        default=str(DEFAULT_OUTPUT_H5AD),
        help="Output layer-annotated h5ad.",
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
        help="Keep h5ad cells without labels and mark them Unknown. Default drops unlabeled cells.",
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Require exact match between h5ad cell IDs and label CSV cell IDs.",
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
    label_csv = Path(args.label_csv)
    output_h5ad = Path(args.output_h5ad)

    section("Stage 03: Layer-wise annotation")

    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_h5ad}")

    if not label_csv.exists():
        raise FileNotFoundError(f"Label CSV not found: {label_csv}")

    print("Input h5ad:", input_h5ad)
    print("Label CSV:", label_csv)
    print("Output h5ad:", output_h5ad)

    section("Loading data")

    adata = ad.read_h5ad(input_h5ad)
    adata.obs_names = adata.obs_names.astype(str)

    print(adata)

    adata = ensure_brnum(adata)

    labels = read_layer_label_csv(label_csv)

    print("Label CSV shape:", labels.shape)
    print(labels.head())

    adata, match_report = add_layer_annotations(
        adata=adata,
        labels=labels,
        keep_unlabeled=args.keep_unlabeled,
        strict=args.strict,
    )

    counts, proportions, summary = save_csv_outputs(adata)

    if not args.skip_plots:
        save_plots(
            adata=adata,
            counts=counts,
            proportions=proportions,
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

    print("\nLayer counts:")
    print(adata.obs["layer_annotation"].astype(str).value_counts().to_string())

    print("\nMain outputs:")
    print("H5AD:")
    print(" ", output_h5ad)
    print("CSV folder:")
    print(" ", TABLE_DIR)
    print("Figure folder:")
    print(" ", FIGURE_DIR)
    print("Per-sample PNG folder:")
    print(" ", PER_SAMPLE_FIGURE_DIR)


if __name__ == "__main__":
    main()