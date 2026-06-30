#!/usr/bin/env python3

"""
01_xenium_qc_imputation_ready.py

Single-file Xenium QC pipeline for imputation preparation.

Input
-----
data/processed/xenium/xenium_N24_raw.h5ad

Main outputs
------------
data/processed/xenium/xenium_N24_qc_all_features.h5ad
data/processed/xenium/xenium_N24_imputation_ready.h5ad

The first output keeps all Xenium features after cell QC.
The second output keeps only Gene Expression features and adds a normalized layer
for imputation/modeling.

This script follows the biological logic of the PI Xenium QC stage:

1. Plot raw QC metrics on tissue before filtering.
2. Calculate per-cell QC metrics.
3. Detect outliers per donor/sample.
4. Remove empty cells and technical outliers.
5. Save outlier IDs and QC summaries.
6. Make kept/discarded supplemental plots.
7. Run donor-level pseudobulk PCA after QC.
8. Save an imputation-ready h5ad.

Important:
---------
Raw Xenium alone usually does not contain Dx, Sex, Age, PNN, CaptureArea, or tear.
If those columns are not present in the h5ad, this script keeps them as unknown
or skips metadata-based PCA coloring.
"""

from __future__ import annotations

import argparse
import math
import traceback
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_INPUT = PROJECT_ROOT / "data" / "processed" / "xenium" / "xenium_N24_raw.h5ad"

OUT_DIR = PROJECT_ROOT / "outputs" / "xenium_branch" / "02_qc"
TABLE_DIR = OUT_DIR / "tables"
FIGURE_DIR = OUT_DIR / "figures"

PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "xenium"

QC_ALL_FEATURES_H5AD = PROCESSED_DIR / "xenium_N24_qc_all_features.h5ad"
IMPUTATION_READY_H5AD = PROCESSED_DIR / "xenium_N24_imputation_ready.h5ad"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Utility functions
# =============================================================================

def section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def ensure_bool_array(x) -> np.ndarray:
    return np.asarray(x, dtype=bool)


def row_sum(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel()
    return np.asarray(X.sum(axis=1)).ravel()


def col_sum(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel()
    return np.asarray(X.sum(axis=0)).ravel()


def row_detected(X) -> np.ndarray:
    if sp.issparse(X):
        return np.asarray((X > 0).sum(axis=1)).ravel()
    return np.asarray((X > 0).sum(axis=1)).ravel()


def safe_percent(numer, denom) -> np.ndarray:
    numer = np.asarray(numer, dtype=float)
    denom = np.asarray(denom, dtype=float)

    out = np.full_like(numer, fill_value=np.nan, dtype=float)
    ok = denom > 0
    out[ok] = 100.0 * numer[ok] / denom[ok]
    return out


def get_counts(adata: ad.AnnData):
    if "counts" in adata.layers:
        return adata.layers["counts"]
    return adata.X


def robust_mad_outlier(x, side: str = "both", nmads: float = 3.0):
    """
    Approximate scuttle::isOutlier-style MAD filtering.

    side:
        "higher" = high outliers only
        "lower"  = low outliers only
        "both"   = low and high outliers
    """
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)

    flags = np.zeros(x.shape[0], dtype=bool)

    if finite.sum() == 0:
        return flags, np.nan, np.nan

    vals = x[finite]
    med = np.nanmedian(vals)
    mad = np.nanmedian(np.abs(vals - med)) * 1.4826

    if mad == 0 or not np.isfinite(mad):
        low = np.nanquantile(vals, 0.01)
        high = np.nanquantile(vals, 0.99)
    else:
        low = med - nmads * mad
        high = med + nmads * mad

    if side == "higher":
        flags[finite] = vals > high
    elif side == "lower":
        flags[finite] = vals < low
    elif side == "both":
        flags[finite] = (vals < low) | (vals > high)
    else:
        raise ValueError(f"Unknown side: {side}")

    return flags, low, high


def sanitize_for_h5ad(adata: ad.AnnData) -> ad.AnnData:
    """
    Avoid h5ad write errors from mixed object columns.
    """
    for df_name in ["obs", "var"]:
        df = getattr(adata, df_name)
        for col in df.columns:
            if pd.api.types.is_object_dtype(df[col]):
                df[col] = df[col].astype(str).replace({"nan": "NA", "None": "NA"})
        setattr(adata, df_name, df)
    return adata


def maybe_attach_sample_metadata(adata: ad.AnnData, metadata_csv: str | None) -> ad.AnnData:
    """
    Optional: attach donor metadata if the user later has a verified CSV.

    Expected columns:
        BrNum, Dx, Sex, Age, PNN, CaptureArea, tear, slide_id, run_date
    """
    if metadata_csv is None:
        return adata

    path = Path(metadata_csv)
    if not path.exists():
        print(f"WARNING: metadata CSV not found, skipping: {path}")
        return adata

    section("Attaching optional sample metadata")

    meta = pd.read_csv(path)
    if "BrNum" not in meta.columns:
        print("WARNING: metadata CSV does not contain BrNum. Skipping metadata merge.")
        return adata

    meta["BrNum"] = meta["BrNum"].astype(str)
    meta = meta.drop_duplicates("BrNum").set_index("BrNum")

    adata.obs["BrNum"] = adata.obs["BrNum"].astype(str)

    for col in meta.columns:
        mapping = meta[col].to_dict()
        adata.obs[col] = adata.obs["BrNum"].map(mapping).fillna(
            adata.obs[col] if col in adata.obs.columns else "unknown"
        )

    print("Attached columns:", list(meta.columns))
    return adata


# =============================================================================
# Feature masks
# =============================================================================

def infer_feature_masks(adata: ad.AnnData) -> dict[str, np.ndarray]:
    """
    Identify Xenium biological/technical feature groups.

    PI logic uses:
        ^NegControlProbe
        ^NegControlCodeword
        ^Unassigned
        MT-
        Type == Gene Expression
    """
    names = pd.Index(adata.var_names.astype(str))

    neg_probe = ensure_bool_array(names.str.startswith("NegControlProbe"))
    neg_codeword = ensure_bool_array(names.str.startswith("NegControlCodeword"))
    unassigned = ensure_bool_array(names.str.startswith("Unassigned"))
    mito = ensure_bool_array(names.str.upper().str.startswith("MT-"))

    if "Type" in adata.var.columns:
        gene_expression = ensure_bool_array(
            adata.var["Type"].astype(str).eq("Gene Expression")
        )
    elif "feature_types" in adata.var.columns:
        gene_expression = ensure_bool_array(
            adata.var["feature_types"].astype(str).eq("Gene Expression")
        )
    elif "is_gene_expression" in adata.var.columns:
        gene_expression = ensure_bool_array(adata.var["is_gene_expression"].astype(bool))
    else:
        gene_expression = ~(neg_probe | neg_codeword | unassigned)

    masks = {
        "gene_expression": ensure_bool_array(gene_expression),
        "neg_probe": ensure_bool_array(neg_probe),
        "neg_codeword": ensure_bool_array(neg_codeword),
        "unassigned": ensure_bool_array(unassigned),
        "mitochondrial": ensure_bool_array(mito),
    }

    feature_summary = pd.DataFrame({
        "feature_group": list(masks.keys()),
        "n_features": [int(v.sum()) for v in masks.values()],
    })

    feature_summary.to_csv(TABLE_DIR / "xenium_N24_qc_feature_summary.csv", index=False)

    print(feature_summary.to_string(index=False))

    return masks


# =============================================================================
# QC metrics
# =============================================================================

def add_qc_metrics(adata: ad.AnnData, masks: dict[str, np.ndarray]) -> ad.AnnData:
    """
    Add PI/scuttle-like QC fields.
    """
    section("Calculating QC metrics")

    X = get_counts(adata)

    total_counts_matrix = row_sum(X)
    detected_all = row_detected(X)

    gene_mask = masks["gene_expression"]
    neg_probe_mask = masks["neg_probe"]
    neg_codeword_mask = masks["neg_codeword"]
    unassigned_mask = masks["unassigned"]
    mito_mask = masks["mitochondrial"]

    gene_counts = row_sum(X[:, gene_mask]) if gene_mask.sum() else np.zeros(adata.n_obs)
    detected_genes = row_detected(X[:, gene_mask]) if gene_mask.sum() else np.zeros(adata.n_obs)

    neg_probe_counts = row_sum(X[:, neg_probe_mask]) if neg_probe_mask.sum() else np.zeros(adata.n_obs)
    neg_codeword_counts = row_sum(X[:, neg_codeword_mask]) if neg_codeword_mask.sum() else np.zeros(adata.n_obs)
    unassigned_counts = row_sum(X[:, unassigned_mask]) if unassigned_mask.sum() else np.zeros(adata.n_obs)
    mito_counts = row_sum(X[:, mito_mask]) if mito_mask.sum() else np.zeros(adata.n_obs)

    # Preserve existing Xenium metadata if present.
    if "total_counts" in adata.obs.columns:
        adata.obs["xenium_total_counts_metadata"] = adata.obs["total_counts"]

    # PI-like fields and explicit qc fields.
    adata.obs["qc_total_counts"] = total_counts_matrix
    adata.obs["qc_gene_counts"] = gene_counts
    adata.obs["qc_detected"] = detected_all
    adata.obs["qc_detected_genes"] = detected_genes

    # PI-style aliases.
    adata.obs["total_counts"] = total_counts_matrix
    adata.obs["detected"] = detected_all

    adata.obs["subsets_negProbe_sum"] = neg_probe_counts
    adata.obs["subsets_negCodeword_sum"] = neg_codeword_counts
    adata.obs["subsets_unassigned_sum"] = unassigned_counts
    adata.obs["subsets_mito_sum"] = mito_counts

    adata.obs["subsets_negProbe_percent"] = safe_percent(neg_probe_counts, total_counts_matrix)
    adata.obs["subsets_negCodeword_percent"] = safe_percent(neg_codeword_counts, total_counts_matrix)
    adata.obs["subsets_unassigned_percent"] = safe_percent(unassigned_counts, total_counts_matrix)
    adata.obs["subsets_mito_percent"] = safe_percent(mito_counts, total_counts_matrix)

    # More explicit duplicate names.
    adata.obs["qc_pct_neg_probe"] = adata.obs["subsets_negProbe_percent"]
    adata.obs["qc_pct_neg_codeword"] = adata.obs["subsets_negCodeword_percent"]
    adata.obs["qc_pct_unassigned"] = adata.obs["subsets_unassigned_percent"]
    adata.obs["qc_pct_mito"] = adata.obs["subsets_mito_percent"]

    adata.obs["qc_empty_cell"] = total_counts_matrix == 0

    return adata


# =============================================================================
# Plotting helpers
# =============================================================================

def plot_box_by_sample(adata: ad.AnnData, metric: str, output_name: str) -> None:
    if metric not in adata.obs.columns:
        print(f"Skipping {metric}: not found in obs")
        return

    brnums = sorted(adata.obs["BrNum"].astype(str).unique())
    data = []

    for br in brnums:
        vals = pd.to_numeric(
            adata.obs.loc[adata.obs["BrNum"].astype(str).eq(br), metric],
            errors="coerce",
        ).dropna()
        data.append(vals.to_numpy())

    plt.figure(figsize=(14, 5))
    plt.boxplot(data, labels=brnums, showfliers=False)
    plt.xticks(rotation=90)
    plt.ylabel(metric)
    plt.title(metric + " by BrNum")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / output_name, dpi=300)
    plt.close()


def plot_spatial_panel(
    adata: ad.AnnData,
    metric: str,
    output_name: str,
    max_cells_per_sample: int = 5000,
) -> None:
    if metric not in adata.obs.columns:
        print(f"Skipping spatial plot for {metric}: not found")
        return

    if "x_centroid" not in adata.obs.columns or "y_centroid" not in adata.obs.columns:
        print("Skipping spatial plot: x_centroid/y_centroid not found")
        return

    brnums = sorted(adata.obs["BrNum"].astype(str).unique())
    n = len(brnums)
    ncols = 4
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.asarray(axes).reshape(-1)

    for ax, br in zip(axes, brnums):
        sub = adata.obs.loc[
            adata.obs["BrNum"].astype(str).eq(br),
            ["x_centroid", "y_centroid", metric],
        ].copy()

        if sub.shape[0] > max_cells_per_sample:
            sub = sub.sample(max_cells_per_sample, random_state=0)

        values = pd.to_numeric(sub[metric], errors="coerce")
        sc = ax.scatter(
            sub["x_centroid"],
            sub["y_centroid"],
            c=values,
            s=0.2,
            alpha=0.8,
        )
        ax.invert_yaxis()
        ax.set_title(br)
        ax.set_xticks([])
        ax.set_yticks([])

    for ax in axes[len(brnums):]:
        ax.axis("off")

    fig.suptitle(metric, y=1.0)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / output_name, dpi=300)
    plt.close()


def plot_raw_qc_metrics(adata: ad.AnnData, max_cells_per_sample: int) -> None:
    """
    PI equivalent of 00_plot_metrics_on_tissue.R.
    """
    section("Plotting raw QC metrics before filtering")

    metrics = [
        "total_counts",
        "control_probe_counts",
        "unassigned_codeword_counts",
        "cell_area",
        "nucleus_area",
        "transcript_counts",
        "subsets_negProbe_percent",
        "subsets_negCodeword_percent",
        "subsets_unassigned_percent",
        "detected",
    ]

    for metric in metrics:
        if metric not in adata.obs.columns:
            continue

        safe_name = metric.replace("/", "_")
        plot_box_by_sample(
            adata,
            metric,
            f"raw_{safe_name}_by_BrNum.png",
        )
        plot_spatial_panel(
            adata,
            metric,
            f"raw_{safe_name}_on_tissue.png",
            max_cells_per_sample=max_cells_per_sample,
        )


def plot_one_sample_report(
    adata: ad.AnnData,
    brnum: str,
    thresholds: dict,
    max_cells: int,
) -> None:
    """
    Per-sample QC report: histograms + spatial final outliers.
    """
    sub = adata[adata.obs["BrNum"].astype(str).eq(brnum)].copy()
    obs = sub.obs.copy()

    if obs.shape[0] > max_cells:
        obs_plot = obs.sample(max_cells, random_state=0)
    else:
        obs_plot = obs

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    axes = axes.reshape(-1)

    hist_metrics = [
        ("subsets_negProbe_percent", thresholds.get("neg_probe_pct_q99")),
        ("subsets_negCodeword_percent", thresholds.get("neg_codeword_pct_q99")),
        ("subsets_unassigned_percent", thresholds.get("unassigned_pct_q99")),
        ("subsets_mito_percent", thresholds.get("mito_pct_high_threshold")),
        ("detected", thresholds.get("detected_low_threshold")),
        ("total_counts", thresholds.get("total_counts_high_threshold")),
    ]

    for ax, (metric, thr) in zip(axes[:6], hist_metrics):
        vals = pd.to_numeric(obs[metric], errors="coerce").dropna()
        ax.hist(vals, bins=50)
        if thr is not None and np.isfinite(thr):
            ax.axvline(thr, linestyle="--")
        ax.set_title(metric)

    if "x_centroid" in obs_plot.columns and "y_centroid" in obs_plot.columns:
        ax = axes[6]
        good = obs_plot[~obs_plot["qc_is_outlier"]]
        bad = obs_plot[obs_plot["qc_is_outlier"]]
        ax.scatter(good["x_centroid"], good["y_centroid"], s=0.2, alpha=0.4)
        ax.scatter(bad["x_centroid"], bad["y_centroid"], s=0.5, alpha=0.9)
        ax.invert_yaxis()
        ax.set_title("Kept vs discarded")
        ax.set_xticks([])
        ax.set_yticks([])

        ax = axes[7]
        ax.scatter(
            obs_plot["x_centroid"],
            obs_plot["y_centroid"],
            c=pd.to_numeric(obs_plot["subsets_unassigned_percent"], errors="coerce"),
            s=0.2,
        )
        ax.invert_yaxis()
        ax.set_title("unassigned %")
        ax.set_xticks([])
        ax.set_yticks([])

        ax = axes[8]
        ax.scatter(
            obs_plot["x_centroid"],
            obs_plot["y_centroid"],
            c=pd.to_numeric(obs_plot["total_counts"], errors="coerce"),
            s=0.2,
        )
        ax.invert_yaxis()
        ax.set_title("total counts")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"QC report: {brnum}", y=1.0)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / f"sample_qc_report_{brnum}.png", dpi=300)
    plt.close()


# =============================================================================
# Outlier calling
# =============================================================================

def call_outliers(
    adata: ad.AnnData,
    nmads: float,
    use_mito: bool,
    strict_pi_threshold: bool,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    section("Calling QC outliers per BrNum")

    for col in [
        "neg_probe_out",
        "neg_codeword_out",
        "unassigned_out",
        "mito_out",
        "detected_out",
        "total_counts_out",
    ]:
        adata.obs[col] = False

    threshold_records = []
    na_records = []

    brnums = sorted(adata.obs["BrNum"].astype(str).unique())

    for br in brnums:
        idx = adata.obs["BrNum"].astype(str).eq(br).to_numpy()
        sub = adata.obs.loc[idx].copy()

        neg_probe = pd.to_numeric(sub["subsets_negProbe_percent"], errors="coerce")
        neg_codeword = pd.to_numeric(sub["subsets_negCodeword_percent"], errors="coerce")
        unassigned = pd.to_numeric(sub["subsets_unassigned_percent"], errors="coerce")
        mito = pd.to_numeric(sub["subsets_mito_percent"], errors="coerce")
        detected = pd.to_numeric(sub["detected"], errors="coerce")
        total_counts = pd.to_numeric(sub["total_counts"], errors="coerce")

        neg_probe_q99 = float(np.nanquantile(neg_probe, 0.99))
        neg_codeword_q99 = float(np.nanquantile(neg_codeword, 0.99))
        unassigned_q99 = float(np.nanquantile(unassigned, 0.99))

        if strict_pi_threshold:
            neg_probe_out = neg_probe >= neg_probe_q99
            neg_codeword_out = neg_codeword >= neg_codeword_q99
            unassigned_out = unassigned >= unassigned_q99
        else:
            # Safer public-data mode:
            # keeps PI's q99 idea but avoids flagging all cells if q99 == 0.
            neg_probe_out = (neg_probe >= neg_probe_q99) & (neg_probe > 0)
            neg_codeword_out = (neg_codeword >= neg_codeword_q99) & (neg_codeword > 0)
            unassigned_out = (unassigned >= unassigned_q99) & (unassigned > 0)

        detected_out, detected_low, detected_high = robust_mad_outlier(
            detected.to_numpy(),
            side="lower",
            nmads=nmads,
        )

        total_counts_out, total_low, total_high = robust_mad_outlier(
            total_counts.to_numpy(),
            side="both",
            nmads=nmads,
        )

        mito_out, mito_low, mito_high = robust_mad_outlier(
            mito.to_numpy(),
            side="higher",
            nmads=nmads,
        )

        adata.obs.loc[idx, "neg_probe_out"] = np.asarray(neg_probe_out, dtype=bool)
        adata.obs.loc[idx, "neg_codeword_out"] = np.asarray(neg_codeword_out, dtype=bool)
        adata.obs.loc[idx, "unassigned_out"] = np.asarray(unassigned_out, dtype=bool)
        adata.obs.loc[idx, "mito_out"] = mito_out
        adata.obs.loc[idx, "detected_out"] = detected_out
        adata.obs.loc[idx, "total_counts_out"] = total_counts_out

        threshold_records.append({
            "BrNum": br,
            "n_cells_before_qc": int(idx.sum()),
            "neg_probe_pct_q99": neg_probe_q99,
            "neg_codeword_pct_q99": neg_codeword_q99,
            "unassigned_pct_q99": unassigned_q99,
            "detected_low_threshold": float(detected_low),
            "detected_high_threshold_unused": float(detected_high),
            "total_counts_low_threshold": float(total_low),
            "total_counts_high_threshold": float(total_high),
            "mito_pct_high_threshold": float(mito_high),
            "nmads": nmads,
            "strict_pi_threshold": strict_pi_threshold,
        })

        na_records.append({
            "BrNum": br,
            "n_na_neg_probe_percent": int(neg_probe.isna().sum()),
            "n_na_neg_codeword_percent": int(neg_codeword.isna().sum()),
            "n_na_unassigned_percent": int(unassigned.isna().sum()),
            "n_na_mito_percent": int(mito.isna().sum()),
            "n_empty_cells": int((total_counts == 0).sum()),
        })

    adata.obs["qc_is_outlier"] = (
        adata.obs["qc_empty_cell"].astype(bool)
        | adata.obs["neg_probe_out"].astype(bool)
        | adata.obs["neg_codeword_out"].astype(bool)
        | adata.obs["unassigned_out"].astype(bool)
        | adata.obs["detected_out"].astype(bool)
        | adata.obs["total_counts_out"].astype(bool)
    )

    if use_mito:
        adata.obs["qc_is_outlier"] = adata.obs["qc_is_outlier"] | adata.obs["mito_out"].astype(bool)

    adata.obs["qc_pass"] = ~adata.obs["qc_is_outlier"]

    thresholds = pd.DataFrame(threshold_records)
    na_diagnostics = pd.DataFrame(na_records)

    thresholds.to_csv(TABLE_DIR / "xenium_N24_qc_thresholds_by_sample.csv", index=False)
    na_diagnostics.to_csv(TABLE_DIR / "xenium_N24_qc_na_diagnostics_by_sample.csv", index=False)

    return adata, thresholds, na_diagnostics


def save_qc_summary_tables(adata: ad.AnnData) -> pd.DataFrame:
    section("Saving QC summary tables")

    summary_records = []

    for br, sub in adata.obs.groupby("BrNum"):
        n_before = sub.shape[0]
        n_out = int(sub["qc_is_outlier"].sum())
        n_after = int(sub["qc_pass"].sum())

        summary_records.append({
            "BrNum": br,
            "n_cells_before_qc": n_before,
            "n_outlier_cells": n_out,
            "n_cells_after_qc": n_after,
            "pct_outlier_cells": 100.0 * n_out / n_before if n_before else np.nan,
            "empty_cells": int(sub["qc_empty_cell"].sum()),
            "neg_probe_outliers": int(sub["neg_probe_out"].sum()),
            "neg_codeword_outliers": int(sub["neg_codeword_out"].sum()),
            "unassigned_outliers": int(sub["unassigned_out"].sum()),
            "mito_outliers_not_default_filter": int(sub["mito_out"].sum()),
            "detected_outliers": int(sub["detected_out"].sum()),
            "total_counts_outliers": int(sub["total_counts_out"].sum()),
            "median_total_counts": float(pd.to_numeric(sub["total_counts"], errors="coerce").median()),
            "median_detected": float(pd.to_numeric(sub["detected"], errors="coerce").median()),
            "median_neg_probe_percent": float(pd.to_numeric(sub["subsets_negProbe_percent"], errors="coerce").median()),
            "median_neg_codeword_percent": float(pd.to_numeric(sub["subsets_negCodeword_percent"], errors="coerce").median()),
            "median_unassigned_percent": float(pd.to_numeric(sub["subsets_unassigned_percent"], errors="coerce").median()),
            "median_mito_percent": float(pd.to_numeric(sub["subsets_mito_percent"], errors="coerce").median()),
        })

    summary = pd.DataFrame(summary_records).sort_values("BrNum").reset_index(drop=True)

    outlier_ids = pd.DataFrame({
        "cell_id": adata.obs.index[adata.obs["qc_is_outlier"].to_numpy()],
        "BrNum": adata.obs.loc[adata.obs["qc_is_outlier"], "BrNum"].values,
    })

    summary.to_csv(TABLE_DIR / "xenium_N24_qc_summary_by_sample.csv", index=False)
    outlier_ids.to_csv(TABLE_DIR / "xenium_N24_outlier_ids.csv", index=False)

    return summary


# =============================================================================
# Supplemental plots
# =============================================================================

def plot_qc_supplement(adata: ad.AnnData, summary: pd.DataFrame, max_cells_per_sample: int) -> None:
    section("Creating supplemental QC plots")

    # Before/after cells
    x = np.arange(summary.shape[0])
    width = 0.4

    plt.figure(figsize=(14, 5))
    plt.bar(x - width / 2, summary["n_cells_before_qc"], width, label="Before QC")
    plt.bar(x + width / 2, summary["n_cells_after_qc"], width, label="After QC")
    plt.xticks(x, summary["BrNum"], rotation=90)
    plt.ylabel("Number of cells")
    plt.title("Cells before and after QC")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "qc_cells_before_after.png", dpi=300)
    plt.close()

    # Outlier reason counts
    reason_cols = [
        "empty_cells",
        "neg_probe_outliers",
        "neg_codeword_outliers",
        "unassigned_outliers",
        "detected_outliers",
        "total_counts_outliers",
        "mito_outliers_not_default_filter",
    ]

    reason_counts = summary[reason_cols].sum().sort_values(ascending=False)

    plt.figure(figsize=(10, 5))
    plt.bar(reason_counts.index, reason_counts.values)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Number of cells")
    plt.title("QC outlier reason counts")
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / "qc_outlier_reason_counts.png", dpi=300)
    plt.close()

    # Spatial kept/discarded panel
    if "x_centroid" in adata.obs.columns and "y_centroid" in adata.obs.columns:
        brnums = sorted(adata.obs["BrNum"].astype(str).unique())
        n = len(brnums)
        ncols = 4
        nrows = math.ceil(n / ncols)

        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
        axes = np.asarray(axes).reshape(-1)

        for ax, br in zip(axes, brnums):
            obs = adata.obs.loc[
                adata.obs["BrNum"].astype(str).eq(br),
                ["x_centroid", "y_centroid", "qc_is_outlier"],
            ].copy()

            if obs.shape[0] > max_cells_per_sample:
                obs = obs.sample(max_cells_per_sample, random_state=0)

            kept = obs[~obs["qc_is_outlier"]]
            disc = obs[obs["qc_is_outlier"]]

            ax.scatter(kept["x_centroid"], kept["y_centroid"], s=0.2, alpha=0.3)
            ax.scatter(disc["x_centroid"], disc["y_centroid"], s=0.5, alpha=0.9)
            ax.invert_yaxis()
            ax.set_title(br)
            ax.set_xticks([])
            ax.set_yticks([])

        for ax in axes[len(brnums):]:
            ax.axis("off")

        fig.suptitle("QC kept vs discarded cells", y=1.0)
        plt.tight_layout()
        plt.savefig(FIGURE_DIR / "qc_kept_discarded_spatial_panel.png", dpi=300)
        plt.close()

    # Violin/box-style plot metrics
    metrics = [
        "detected",
        "total_counts",
        "subsets_negProbe_percent",
        "subsets_negCodeword_percent",
    ]

    for metric in metrics:
        plot_box_by_sample(adata, metric, f"qc_kept_discarded_{metric}_by_sample.png")


# =============================================================================
# Pseudobulk PCA
# =============================================================================

def pca_numpy(X: np.ndarray, n_components: int = 10):
    """
    PCA using numpy SVD. Rows are samples, columns are genes.
    """
    X = np.asarray(X, dtype=float)
    X_centered = X - X.mean(axis=0, keepdims=True)

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

    scores = U[:, :n_components] * S[:n_components]
    variance = (S ** 2) / max((X.shape[0] - 1), 1)
    explained = variance / variance.sum() if variance.sum() > 0 else np.zeros_like(variance)

    return scores, explained[:n_components]


def build_pseudobulk_pca(
    adata_qc: ad.AnnData,
    gene_mask: np.ndarray,
    target_sum: float,
) -> pd.DataFrame:
    section("Running donor-level pseudobulk PCA after QC")

    X = get_counts(adata_qc)
    Xg = X[:, gene_mask]

    brnums = sorted(adata_qc.obs["BrNum"].astype(str).unique())

    rows = []
    pb_counts = []

    for br in brnums:
        idx = adata_qc.obs["BrNum"].astype(str).eq(br).to_numpy()
        summed = col_sum(Xg[idx, :])
        pb_counts.append(summed)

        sub = adata_qc.obs.loc[idx]

        row = {
            "BrNum": br,
            "n_cells_after_qc": int(idx.sum()),
            "median_total_counts": float(pd.to_numeric(sub["total_counts"], errors="coerce").median()),
            "median_detected": float(pd.to_numeric(sub["detected"], errors="coerce").median()),
            "median_unassigned_percent": float(pd.to_numeric(sub["subsets_unassigned_percent"], errors="coerce").median()),
        }

        # Carry donor-level metadata if available.
        for col in ["Dx", "PNN", "Sex", "Age", "CaptureArea", "tear", "slide_id", "run_date"]:
            if col in sub.columns:
                vals = sub[col].dropna().astype(str).unique()
                row[col] = vals[0] if len(vals) else "unknown"

        rows.append(row)

    pb_counts = np.vstack(pb_counts)
    pb_meta = pd.DataFrame(rows)

    gene_names = adata_qc.var_names[gene_mask].astype(str).tolist()

    pb_counts_df = pd.DataFrame(pb_counts, index=brnums, columns=gene_names)
    pb_counts_df.to_csv(TABLE_DIR / "xenium_N24_pseudobulk_counts_by_BrNum.csv.gz", compression="gzip")

    lib = pb_counts.sum(axis=1)
    lib[lib == 0] = 1

    pb_log = np.log1p(pb_counts / lib[:, None] * target_sum)

    scores, explained = pca_numpy(pb_log, n_components=min(10, pb_log.shape[0], pb_log.shape[1]))

    for i in range(scores.shape[1]):
        pb_meta[f"PC{i+1}"] = scores[:, i]

    for i, val in enumerate(explained):
        pb_meta[f"PC{i+1}_variance_explained"] = val

    pb_meta.to_csv(TABLE_DIR / "xenium_N24_pseudobulk_pca_coordinates.csv", index=False)

    # PCA plot labeled by BrNum
    if "PC1" in pb_meta.columns and "PC2" in pb_meta.columns:
        plt.figure(figsize=(7, 6))
        plt.scatter(pb_meta["PC1"], pb_meta["PC2"], s=60)
        for _, row in pb_meta.iterrows():
            plt.text(row["PC1"], row["PC2"], row["BrNum"], fontsize=8)
        plt.xlabel(f"PC1 ({explained[0]*100:.1f}%)" if len(explained) > 0 else "PC1")
        plt.ylabel(f"PC2 ({explained[1]*100:.1f}%)" if len(explained) > 1 else "PC2")
        plt.title("Pseudobulk PCA by BrNum")
        plt.tight_layout()
        plt.savefig(FIGURE_DIR / "pseudobulk_pca_by_BrNum.png", dpi=300)
        plt.close()

    return pb_meta


# =============================================================================
# Imputation-ready object
# =============================================================================

def make_log1p_norm_layer(adata: ad.AnnData, target_sum: float) -> ad.AnnData:
    """
    Add normalized log1p layer for modeling.

    Keeps X as raw counts.
    Adds:
        layers["log1p_norm"]
    """
    X = get_counts(adata)

    lib = row_sum(X)
    lib_safe = lib.copy()
    lib_safe[lib_safe == 0] = 1.0

    scale = target_sum / lib_safe

    if sp.issparse(X):
        X_norm = X.multiply(scale[:, None]).tocsr()
        X_norm.data = np.log1p(X_norm.data)
    else:
        X_norm = np.log1p(X * scale[:, None])

    adata.layers["log1p_norm"] = X_norm
    adata.obs["imputation_size_factor"] = lib_safe / target_sum

    return adata


def save_qc_and_imputation_objects(
    adata: ad.AnnData,
    masks: dict[str, np.ndarray],
    target_sum: float,
    compression: str,
) -> tuple[ad.AnnData, ad.AnnData]:
    section("Saving QC and imputation-ready h5ad objects")

    adata_qc = adata[adata.obs["qc_pass"].to_numpy(), :].copy()

    adata_qc.uns["stage"] = "02_xenium_qc"
    adata_qc.uns["description"] = "QC-filtered Xenium object with all features retained."
    adata_qc.uns["outlier_logic"] = (
        "Empty cells OR high negative probe OR high negative codeword OR high unassigned "
        "OR low detected OR abnormal total counts. Mito calculated but optional."
    )

    adata_qc = sanitize_for_h5ad(adata_qc)

    if compression == "none":
        adata_qc.write_h5ad(QC_ALL_FEATURES_H5AD)
    else:
        adata_qc.write_h5ad(QC_ALL_FEATURES_H5AD, compression=compression)

    print(f"Saved: {QC_ALL_FEATURES_H5AD}")

    gene_mask = masks["gene_expression"]
    adata_imp = adata_qc[:, gene_mask].copy()

    adata_imp.uns["stage"] = "02_xenium_qc_imputation_ready"
    adata_imp.uns["description"] = (
        "QC-filtered, Gene Expression-only Xenium object prepared as imputation target."
    )
    adata_imp.uns["raw_counts_layer"] = "counts"
    adata_imp.uns["normalized_layer"] = "log1p_norm"
    adata_imp.uns["spatial_key"] = "spatial"
    adata_imp.uns["target_sum"] = target_sum

    adata_imp = make_log1p_norm_layer(adata_imp, target_sum=target_sum)
    adata_imp = sanitize_for_h5ad(adata_imp)

    if compression == "none":
        adata_imp.write_h5ad(IMPUTATION_READY_H5AD)
    else:
        adata_imp.write_h5ad(IMPUTATION_READY_H5AD, compression=compression)

    print(f"Saved: {IMPUTATION_READY_H5AD}")

    return adata_qc, adata_imp


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-h5ad",
        default=str(DEFAULT_INPUT),
        help="Input raw Xenium h5ad.",
    )

    parser.add_argument(
        "--sample-metadata-csv",
        default=None,
        help="Optional donor metadata CSV with BrNum column.",
    )

    parser.add_argument(
        "--nmads",
        type=float,
        default=3.0,
        help="MAD threshold for detected/total/mito outliers.",
    )

    parser.add_argument(
        "--use-mito",
        action="store_true",
        help="Include mitochondrial outliers in final filtering.",
    )

    parser.add_argument(
        "--strict-pi-threshold",
        action="store_true",
        help="Use exact PI >= q99 rule for control/unassigned outliers. Default is safer q99 and >0.",
    )

    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Target sum for log-normalized imputation layer.",
    )

    parser.add_argument(
        "--max-plot-cells-per-sample",
        type=int,
        default=5000,
        help="Maximum cells per sample for spatial plots.",
    )

    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip plots for faster execution.",
    )

    parser.add_argument(
        "--write-cell-qc-csv",
        action="store_true",
        help="Write per-cell QC CSV. Can be very large.",
    )

    parser.add_argument(
        "--compression",
        default="lzf",
        choices=["lzf", "gzip", "none"],
        help="h5ad compression. lzf is faster; gzip is smaller.",
    )

    args = parser.parse_args()

    section("Xenium QC pipeline for imputation preparation")

    input_h5ad = Path(args.input_h5ad)

    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_h5ad}")

    print(f"Input: {input_h5ad}")

    adata = ad.read_h5ad(input_h5ad)
    print(adata)

    if "BrNum" not in adata.obs.columns:
        raise ValueError("Input AnnData must contain adata.obs['BrNum'].")

    adata.obs["BrNum"] = adata.obs["BrNum"].astype(str)

    adata = maybe_attach_sample_metadata(adata, args.sample_metadata_csv)

    # Placeholder metadata columns for future compatibility.
    for col in ["Dx", "PNN", "Sex", "Age", "CaptureArea", "tear"]:
        if col not in adata.obs.columns:
            adata.obs[col] = "unknown"

    masks = infer_feature_masks(adata)

    adata = add_qc_metrics(adata, masks)

    if not args.skip_plots:
        plot_raw_qc_metrics(
            adata,
            max_cells_per_sample=args.max_plot_cells_per_sample,
        )

    adata, thresholds, na_diagnostics = call_outliers(
        adata=adata,
        nmads=args.nmads,
        use_mito=args.use_mito,
        strict_pi_threshold=args.strict_pi_threshold,
    )

    summary = save_qc_summary_tables(adata)

    if args.write_cell_qc_csv:
        cell_qc_cols = [
            "BrNum",
            "Dx",
            "total_counts",
            "detected",
            "subsets_negProbe_percent",
            "subsets_negCodeword_percent",
            "subsets_unassigned_percent",
            "subsets_mito_percent",
            "qc_empty_cell",
            "neg_probe_out",
            "neg_codeword_out",
            "unassigned_out",
            "mito_out",
            "detected_out",
            "total_counts_out",
            "qc_is_outlier",
            "qc_pass",
        ]
        existing = [c for c in cell_qc_cols if c in adata.obs.columns]
        adata.obs[existing].to_csv(
            TABLE_DIR / "xenium_N24_cell_qc_metrics.csv.gz",
            compression="gzip",
        )

    if not args.skip_plots:
        plot_qc_supplement(
            adata,
            summary,
            max_cells_per_sample=args.max_plot_cells_per_sample,
        )

        # Per-sample reports.
        thr_by_br = thresholds.set_index("BrNum").to_dict(orient="index")
        for br in sorted(adata.obs["BrNum"].astype(str).unique()):
            try:
                plot_one_sample_report(
                    adata=adata,
                    brnum=br,
                    thresholds=thr_by_br.get(br, {}),
                    max_cells=args.max_plot_cells_per_sample,
                )
            except Exception:
                print(f"WARNING: failed to make sample report for {br}")
                traceback.print_exc()

    adata_qc, adata_imp = save_qc_and_imputation_objects(
        adata=adata,
        masks=masks,
        target_sum=args.target_sum,
        compression=args.compression,
    )

    # Pseudobulk PCA should be based on QC-filtered cells and gene-expression features.
    build_pseudobulk_pca(
        adata_qc=adata_qc,
        gene_mask=masks["gene_expression"],
        target_sum=1e6,
    )

    section("Final summary")

    print(f"Cells before QC:         {adata.n_obs}")
    print(f"Cells removed:           {int(adata.obs['qc_is_outlier'].sum())}")
    print(f"Cells after QC:          {adata_qc.n_obs}")
    print(f"Features all QC object:  {adata_qc.n_vars}")
    print(f"Features imputation obj: {adata_imp.n_vars}")
    print(f"Donors:                  {adata_qc.obs['BrNum'].nunique()}")

    print("\nMain files:")
    print(f"  {QC_ALL_FEATURES_H5AD}")
    print(f"  {IMPUTATION_READY_H5AD}")
    print(f"  {TABLE_DIR / 'xenium_N24_outlier_ids.csv'}")
    print(f"  {TABLE_DIR / 'xenium_N24_qc_summary_by_sample.csv'}")
    print(f"  {TABLE_DIR / 'xenium_N24_pseudobulk_pca_coordinates.csv'}")
    print(f"  {FIGURE_DIR}")


if __name__ == "__main__":
    main()