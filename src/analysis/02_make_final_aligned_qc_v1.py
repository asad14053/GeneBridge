#!/usr/bin/env python3

"""
02_make_final_aligned_qc.py

QC of FINAL Xenium -> Visium alignment.

Uses:
    Xenium:
        x_aligned
        y_aligned
        layer
        inside_visium
        usable
        spot_dist_um
        visium_spot

    Visium:
        x
        y
        domain

NO additional rotation, reflection, translation or scaling is performed.
The delivered aligned coordinates are plotted directly.
"""

from pathlib import Path
import argparse
import traceback

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

DEFAULT_ROOT = (
    PROJECT_ROOT
    / "data"
    / "metadata"
    / "stalign_results_update"
)

DEFAULT_OUT = (
    PROJECT_ROOT
    / "results"
    / "analysis"
    / "stalign_final_aligned_qc"
)


# ============================================================
# Layer/domain ordering: pia -> WM
# ============================================================

LAYER_MAP = {
    "spd07": 1,
    "l1": 1,
    "l1/m": 1,
    "l1m": 1,

    "spd06": 2,
    "l2/3": 2,
    "l23": 2,

    "spd02": 3,
    "l3/4": 3,
    "l34": 3,

    "spd05": 4,
    "l5": 4,

    "spd03": 5,
    "l6": 5,

    "spd01": 6,
    "wmtz": 6,

    "spd04": 7,
    "wm": 7,
}


def normalize_label(x):
    if pd.isna(x):
        return None

    return (
        str(x)
        .strip()
        .lower()
        .replace("_", "")
        .replace("-", "")
        .replace(" ", "")
    )


def layer_to_ordinal(series):

    out = pd.Series(
        np.nan,
        index=series.index,
        dtype=float
    )

    numeric = pd.to_numeric(
        series,
        errors="coerce"
    )

    ok = numeric.between(1, 7)

    out.loc[ok] = numeric.loc[ok]

    for i, x in series.items():

        if pd.notna(out.loc[i]):
            continue

        key = normalize_label(x)

        if key in LAYER_MAP:
            out.loc[i] = LAYER_MAP[key]

    return out


def parse_bool(series):

    if series.dtype == bool:
        return series

    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({
            "true": True,
            "false": False,
            "1": True,
            "0": False,
            "yes": True,
            "no": False,
        })
        .fillna(False)
        .astype(bool)
    )


# ============================================================
# Process one donor
# ============================================================

def process_sample(
    br,
    aligned_file,
    visium_file,
    outdir
):

    print()
    print("=" * 80)
    print(br)
    print("=" * 80)

    print("Aligned:", aligned_file)
    print("Visium :", visium_file)

    xen = pd.read_csv(aligned_file)
    vis = pd.read_csv(visium_file)

    required_xen = [
        "x_aligned",
        "y_aligned",
        "layer"
    ]

    for c in required_xen:
        if c not in xen.columns:
            raise ValueError(
                f"{br}: missing Xenium column {c}"
            )

    required_vis = [
        "x",
        "y",
        "domain"
    ]

    for c in required_vis:
        if c not in vis.columns:
            raise ValueError(
                f"{br}: missing Visium column {c}"
            )

    # --------------------------------------------------------
    # Coordinates
    # --------------------------------------------------------

    xen_xy = np.c_[
        pd.to_numeric(
            xen["x_aligned"],
            errors="coerce"
        ),
        pd.to_numeric(
            xen["y_aligned"],
            errors="coerce"
        )
    ]

    vis_xy = np.c_[
        pd.to_numeric(
            vis["x"],
            errors="coerce"
        ),
        pd.to_numeric(
            vis["y"],
            errors="coerce"
        )
    ]

    xen_layer = (
        layer_to_ordinal(
            xen["layer"]
        ).to_numpy()
    )

    vis_layer = (
        layer_to_ordinal(
            vis["domain"]
        ).to_numpy()
    )

    xen_good = (
        np.isfinite(xen_xy).all(axis=1)
        & np.isfinite(xen_layer)
    )

    vis_good = (
        np.isfinite(vis_xy).all(axis=1)
        & np.isfinite(vis_layer)
    )

    # keep original indices for boolean masks
    xidx = np.where(xen_good)[0]

    xen_xy_good = xen_xy[xen_good]
    xen_layer_good = xen_layer[xen_good]

    vis_xy_good = vis_xy[vis_good]
    vis_layer_good = vis_layer[vis_good]

    # --------------------------------------------------------
    # Final-alignment masks
    # --------------------------------------------------------

    if "inside_visium" in xen.columns:
        inside_full = parse_bool(
            xen["inside_visium"]
        ).to_numpy()
    else:
        inside_full = np.ones(
            len(xen),
            dtype=bool
        )

    if "usable" in xen.columns:
        usable_full = parse_bool(
            xen["usable"]
        ).to_numpy()
    else:
        usable_full = inside_full.copy()

    inside = inside_full[xidx]
    usable = usable_full[xidx]

    # --------------------------------------------------------
    # Nearest Visium spot
    # --------------------------------------------------------

    tree = cKDTree(
        vis_xy_good
    )

    nearest_dist, nearest_idx = tree.query(
        xen_xy_good,
        k=1
    )

    nearest_vis_layer = (
        vis_layer_good[
            nearest_idx
        ]
    )

    # exact ordinal agreement
    layer_match = (
        xen_layer_good
        == nearest_vis_layer
    )

    # correlation across cells
    rho_all = spearmanr(
        xen_layer_good,
        nearest_vis_layer,
        nan_policy="omit"
    ).statistic

    if usable.sum() > 2:

        rho_usable = spearmanr(
            xen_layer_good[usable],
            nearest_vis_layer[usable],
            nan_policy="omit"
        ).statistic

        agreement_usable = (
            layer_match[usable].mean()
        )

    else:

        rho_usable = np.nan
        agreement_usable = np.nan

    # --------------------------------------------------------
    # Existing supplied distance metric
    # --------------------------------------------------------

    if "spot_dist_um" in xen.columns:

        supplied_dist = pd.to_numeric(
            xen.loc[
                xidx,
                "spot_dist_um"
            ],
            errors="coerce"
        ).to_numpy()

    else:

        supplied_dist = nearest_dist.copy()

    median_dist_all = np.nanmedian(
        supplied_dist
    )

    median_dist_usable = (
        np.nanmedian(
            supplied_dist[usable]
        )
        if usable.any()
        else np.nan
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    n_cells = len(xen)

    n_valid = len(
        xen_xy_good
    )

    n_inside = int(
        inside.sum()
    )

    n_usable = int(
        usable.sum()
    )

    inside_pct = (
        100 * n_inside / n_valid
    )

    usable_pct = (
        100 * n_usable / n_valid
    )

    assigned_pct = np.nan

    if "visium_spot" in xen.columns:

        assigned = (
            xen.loc[
                xidx,
                "visium_spot"
            ]
            .notna()
            .to_numpy()
        )

        # exclude blank strings
        vals = (
            xen.loc[
                xidx,
                "visium_spot"
            ]
            .astype(str)
            .str.strip()
            .to_numpy()
        )

        assigned &= (
            vals != ""
        )

        assigned &= (
            vals != "nan"
        )

        assigned_pct = (
            100 * assigned.mean()
        )

    print(
        f"Cells total             : {n_cells:,}"
    )

    print(
        f"Cells with layer+coords : {n_valid:,}"
    )

    print(
        f"Inside Visium           : "
        f"{n_inside:,} ({inside_pct:.2f}%)"
    )

    print(
        f"Usable                  : "
        f"{n_usable:,} ({usable_pct:.2f}%)"
    )

    print(
        f"Median spot distance    : "
        f"{median_dist_all:.2f} um"
    )

    print(
        f"Median usable distance  : "
        f"{median_dist_usable:.2f} um"
    )

    print(
        f"Layer rho, all          : "
        f"{rho_all:.3f}"
    )

    print(
        f"Layer rho, usable       : "
        f"{rho_usable:.3f}"
    )

    print(
        f"Exact layer agreement   : "
        f"{agreement_usable:.3f}"
    )

    # ========================================================
    # Figure
    # ========================================================

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 6)
    )

    fig.suptitle(
        f"{br} — Final Xenium → Visium alignment",
        fontsize=15,
        fontweight="bold"
    )

    # --------------------------------------------------------
    # PANEL 1: all aligned cells
    # --------------------------------------------------------

    ax = axes[0]

    ax.scatter(
        vis_xy_good[:, 0] / 1000,
        vis_xy_good[:, 1] / 1000,
        c=vis_layer_good,
        cmap="viridis",
        s=10,
        alpha=0.75,
        linewidths=0,
        rasterized=True
    )

    ax.scatter(
        xen_xy_good[:, 0] / 1000,
        xen_xy_good[:, 1] / 1000,
        c=xen_layer_good,
        cmap="autumn",
        s=1.2,
        alpha=0.30,
        linewidths=0,
        rasterized=True
    )

    ax.set_title(
        "All aligned Xenium cells\n"
        f"inside={inside_pct:.1f}%  "
        f"usable={usable_pct:.1f}%"
    )

    # --------------------------------------------------------
    # PANEL 2: usable cells only
    # --------------------------------------------------------

    ax = axes[1]

    ax.scatter(
        vis_xy_good[:, 0] / 1000,
        vis_xy_good[:, 1] / 1000,
        c=vis_layer_good,
        cmap="viridis",
        s=10,
        alpha=0.75,
        linewidths=0,
        rasterized=True
    )

    ax.scatter(
        xen_xy_good[usable, 0] / 1000,
        xen_xy_good[usable, 1] / 1000,
        c=xen_layer_good[usable],
        cmap="autumn",
        s=1.5,
        alpha=0.35,
        linewidths=0,
        rasterized=True
    )

    ax.set_title(
        "Usable Xenium cells\n"
        f"layer rho={rho_usable:.3f}  "
        f"median dist={median_dist_usable:.1f} µm"
    )

    # --------------------------------------------------------
    # PANEL 3: Visium reference
    # --------------------------------------------------------

    ax = axes[2]

    ax.scatter(
        vis_xy_good[:, 0] / 1000,
        vis_xy_good[:, 1] / 1000,
        c=vis_layer_good,
        cmap="viridis",
        s=12,
        linewidths=0,
        rasterized=True
    )

    ax.set_title(
        "Visium spatial-domain reference"
    )

    # --------------------------------------------------------
    # styling
    # --------------------------------------------------------

    combined = np.vstack([
        xen_xy_good,
        vis_xy_good
    ])

    xmin, ymin = np.nanmin(
        combined,
        axis=0
    )

    xmax, ymax = np.nanmax(
        combined,
        axis=0
    )

    padx = 0.03 * (
        xmax - xmin
    )

    pady = 0.03 * (
        ymax - ymin
    )

    for ax in axes:

        ax.set_xlabel("mm")
        ax.set_ylabel("mm")

        ax.set_xlim(
            (xmin - padx) / 1000,
            (xmax + padx) / 1000
        )

        ax.set_ylim(
            (ymin - pady) / 1000,
            (ymax + pady) / 1000
        )

        ax.invert_yaxis()
        ax.set_aspect("equal")

    footer = (
        f"Final delivered x_aligned/y_aligned coordinates; "
        f"no additional transformation applied. "
        f"Usable={usable_pct:.1f}%; "
        f"median usable spot distance={median_dist_usable:.1f} µm; "
        f"nearest-Visium-domain Spearman rho={rho_usable:.3f}; "
        f"exact ordinal domain agreement={agreement_usable:.1%}."
    )

    fig.text(
        0.5,
        0.015,
        footer,
        ha="center",
        fontsize=8
    )

    plt.tight_layout(
        rect=[
            0,
            0.06,
            1,
            0.94
        ]
    )

    out_png = (
        outdir
        / f"{br}_final_aligned_qc.png"
    )

    out_pdf = (
        outdir
        / f"{br}_final_aligned_qc.pdf"
    )

    fig.savefig(
        out_png,
        dpi=250,
        bbox_inches="tight"
    )

    fig.savefig(
        out_pdf,
        bbox_inches="tight"
    )

    plt.close(fig)

    return {
        "BrNum": br,
        "status": "PASS",

        "n_cells": n_cells,
        "n_valid": n_valid,

        "n_inside_visium": n_inside,
        "inside_visium_pct": inside_pct,

        "n_usable": n_usable,
        "usable_pct": usable_pct,

        "assigned_to_spot_pct": assigned_pct,

        "median_spot_dist_um": median_dist_all,
        "median_usable_spot_dist_um": median_dist_usable,

        "nearest_domain_rho_all": rho_all,
        "nearest_domain_rho_usable": rho_usable,

        "exact_layer_agreement_usable": agreement_usable,

        "figure_png": str(out_png),
        "figure_pdf": str(out_pdf)
    }


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT
    )

    parser.add_argument(
        "--br",
        type=str,
        default=None
    )

    args = parser.parse_args()

    root = args.root.resolve()
    outdir = args.out.resolve()

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    aligned_files = sorted(
        root.rglob(
            "Br*_aligned_coordinates.csv"
        )
    )

    samples = []

    for aligned_file in aligned_files:

        br = (
            aligned_file.name
            .replace(
                "_aligned_coordinates.csv",
                ""
            )
        )

        if (
            args.br is not None
            and br != args.br
        ):
            continue

        visium_file = (
            aligned_file.parent
            / f"{br}_visium_coordinates.csv"
        )

        if not visium_file.exists():

            print(
                f"WARNING {br}: missing "
                f"{visium_file.name}"
            )

            continue

        samples.append(
            (
                br,
                aligned_file,
                visium_file
            )
        )

    print(
        f"Samples found: {len(samples)}"
    )

    results = []

    for (
        br,
        aligned_file,
        visium_file
    ) in samples:

        try:

            r = process_sample(
                br,
                aligned_file,
                visium_file,
                outdir
            )

        except Exception as e:

            traceback.print_exc()

            r = {
                "BrNum": br,
                "status": "FAIL",
                "error": str(e)
            }

        results.append(r)

    summary = pd.DataFrame(
        results
    )

    summary_file = (
        outdir
        / "final_aligned_qc_summary.csv"
    )

    summary.to_csv(
        summary_file,
        index=False
    )

    print()
    print("=" * 80)
    print("FINAL ALIGNMENT SUMMARY")
    print("=" * 80)

    cols = [
        c
        for c in [
            "BrNum",
            "status",
            "inside_visium_pct",
            "usable_pct",
            "median_usable_spot_dist_um",
            "nearest_domain_rho_usable",
            "exact_layer_agreement_usable"
        ]
        if c in summary.columns
    ]

    print(
        summary[
            cols
        ].to_string(
            index=False
        )
    )

    print()
    print(
        "Summary:",
        summary_file
    )


if __name__ == "__main__":
    main()
