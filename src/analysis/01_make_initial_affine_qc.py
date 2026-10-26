#!/usr/bin/env python3

"""
01_make_initial_affine_qc.py

Create initial affine alignment QC plots for all Xenium -> Visium samples.

Expected files for each sample:
    BrXXXX_aligned_coordinates.csv
    visium_coordinate.csv

The script:
    1. Uses native Xenium x/y coordinates.
    2. Uses Xenium layer/domain annotations.
    3. Loads Visium coordinates + layer/domain annotations.
    4. Estimates the pia -> white-matter anatomical depth axis.
    5. Tests both:
           mirror=False
           mirror=True
    6. Scores each orientation using layer concordance.
    7. Generates a 3-panel QC figure:
           left   = mirror=False
           middle = mirror=True
           right  = Visium reference
    8. Writes a cohort-level summary CSV.

IMPORTANT:
    This is INITIAL AFFINE QC.
    It intentionally uses Xenium native x/y, NOT x_aligned/y_aligned.

Author: GeneBridge
"""

from pathlib import Path
import argparse
import traceback

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr
from scipy.spatial import cKDTree, ConvexHull
from matplotlib.path import Path as MplPath


# ============================================================
# DEFAULT PROJECT PATHS
# ============================================================

PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

DEFAULT_INPUT_ROOT = (
    PROJECT_ROOT / "data" / "metadata"
)

DEFAULT_OUTPUT_ROOT = (
    PROJECT_ROOT
    / "results"
    / "analysis"
    / "stalign_initial_affine_qc"
)


# ============================================================
# LAYER ORDER
#
# pia --> WM
#
# GeneBridge / paper mapping:
#
# spd07 -> L1/M
# spd06 -> L2/3
# spd02 -> L3/4
# spd05 -> L5
# spd03 -> L6
# spd01 -> WMtz
# spd04 -> WM
# ============================================================

LAYER_ORDINAL = {
    # paper spatial domains
    "spd07": 1,
    "spd06": 2,
    "spd02": 3,
    "spd05": 4,
    "spd03": 5,
    "spd01": 6,
    "spd04": 7,

    # textual layer names
    "l1": 1,
    "l1/m": 1,
    "l1m": 1,
    "layer1": 1,

    "l2/3": 2,
    "l23": 2,
    "layer2/3": 2,
    "layer23": 2,

    "l3/4": 3,
    "l34": 3,
    "layer3/4": 3,
    "layer34": 3,

    "l5": 4,
    "layer5": 4,

    "l6": 5,
    "layer6": 5,

    "wmtz": 6,
    "wmtransition": 6,
    "whitemattertransition": 6,

    "wm": 7,
    "whitematter": 7,
}


# ============================================================
# BASIC HELPERS
# ============================================================

def normalize_label(x):
    if pd.isna(x):
        return None

    s = str(x).strip().lower()

    s = s.replace("_", "")
    s = s.replace("-", "")
    s = s.replace(" ", "")

    return s


def layer_to_ordinal(series, name="layer"):
    """
    Convert layer/domain labels to ordered cortical-depth values.
    """

    out = pd.Series(
        np.nan,
        index=series.index,
        dtype=float
    )

    # allow already-numeric 1..7
    numeric = pd.to_numeric(
        series,
        errors="coerce"
    )

    mask = numeric.between(1, 7)

    out.loc[mask] = numeric.loc[mask]

    # map strings
    for idx, value in series.items():

        if pd.notna(out.loc[idx]):
            continue

        key = normalize_label(value)

        if key in LAYER_ORDINAL:
            out.loc[idx] = LAYER_ORDINAL[key]

    mapped_fraction = out.notna().mean()

    print(
        f"  {name}: mapped "
        f"{out.notna().sum():,}/{len(out):,} "
        f"({mapped_fraction:.1%})"
    )

    if mapped_fraction < 0.40:

        examples = (
            series.dropna()
            .astype(str)
            .unique()[:30]
        )

        raise ValueError(
            f"Could map only {mapped_fraction:.1%} "
            f"of {name}.\n"
            f"Observed labels: {examples}"
        )

    return out


# ============================================================
# VISIUM COLUMN DETECTION
# ============================================================

def detect_xy_columns(df):
    """
    Detect coordinate columns in visium_coordinate.csv.
    """

    lower = {
        str(c).lower(): c
        for c in df.columns
    }

    candidates = [
        ("x_um", "y_um"),
        ("x_micron", "y_micron"),
        ("x_microns", "y_microns"),
        ("x", "y"),
        ("spatial_x", "spatial_y"),
        ("imagecol", "imagerow"),
        ("pxl_col_in_fullres", "pxl_row_in_fullres"),
        ("pixel_x", "pixel_y"),
        ("col", "row"),
    ]

    for x_name, y_name in candidates:

        if (
            x_name.lower() in lower
            and y_name.lower() in lower
        ):
            return (
                lower[x_name.lower()],
                lower[y_name.lower()]
            )

    raise ValueError(
        "Could not determine Visium X/Y columns.\n"
        f"Available columns:\n{list(df.columns)}"
    )


def detect_layer_column(df):
    """
    Detect Visium cortical domain/layer column.
    """

    lower = {
        str(c).lower(): c
        for c in df.columns
    }

    candidates = [
        "predictions_smooth",
        "paper_predictions_smooth",
        "layer_annotation",
        "spatial_domain",
        "domain",
        "layer",
        "spd",
    ]

    for candidate in candidates:

        if candidate.lower() in lower:
            return lower[candidate.lower()]

    # more permissive fallback
    for c in df.columns:

        lc = str(c).lower()

        if (
            "prediction" in lc
            or "layer" in lc
            or "domain" in lc
        ):
            return c

    raise ValueError(
        "Could not determine Visium layer/domain column.\n"
        f"Available columns:\n{list(df.columns)}"
    )


# ============================================================
# VISIUM COORDINATE UNIT NORMALIZATION
# ============================================================

def normalize_visium_coordinates(
    xy,
    expected_pitch_um=100.0
):
    """
    Visium spots are approximately 100 um center-to-center.

    If coordinates are already microns:
        median nearest-neighbor ~= 100
        scale ~= 1

    If coordinates are full-resolution pixels:
        infer pixel->micron conversion from spot spacing.
    """

    xy = np.asarray(
        xy,
        dtype=float
    )

    good = np.isfinite(xy).all(axis=1)

    clean = xy[good]

    if len(clean) < 20:
        raise ValueError(
            "Too few valid Visium coordinates."
        )

    tree = cKDTree(clean)

    distances, _ = tree.query(
        clean,
        k=2
    )

    nearest = distances[:, 1]

    nearest = nearest[
        np.isfinite(nearest)
        & (nearest > 0)
    ]

    median_nn = np.median(nearest)

    scale = (
        expected_pitch_um
        / median_nn
    )

    return (
        xy * scale,
        median_nn,
        scale
    )


# ============================================================
# ANATOMICAL DEPTH AXIS
# ============================================================

def fit_depth_axis(
    xy,
    layer_ordinal
):
    """
    Fit

        layer_depth ~ x + y

    Regression gradient defines pia -> WM direction.
    """

    xy = np.asarray(
        xy,
        dtype=float
    )

    layer_ordinal = np.asarray(
        layer_ordinal,
        dtype=float
    )

    good = (
        np.isfinite(xy).all(axis=1)
        & np.isfinite(layer_ordinal)
    )

    xy = xy[good]
    layer_ordinal = layer_ordinal[good]

    if len(xy) < 30:

        raise ValueError(
            "Too few annotated coordinates "
            "to fit anatomical depth axis."
        )

    center = np.mean(
        xy,
        axis=0
    )

    centered = (
        xy - center
    )

    design = np.c_[
        centered[:, 0],
        centered[:, 1],
        np.ones(len(centered))
    ]

    beta, *_ = np.linalg.lstsq(
        design,
        layer_ordinal,
        rcond=None
    )

    depth = beta[:2]

    norm = np.linalg.norm(depth)

    if norm < 1e-12:

        raise ValueError(
            "Could not estimate depth direction."
        )

    depth = depth / norm

    tangent = np.array([
        -depth[1],
        depth[0]
    ])

    return (
        center,
        depth,
        tangent
    )


# ============================================================
# INITIAL AFFINE TRANSFORMATION
# ============================================================

def initial_affine(
    source_xy,
    source_center,
    source_depth,
    source_tangent,
    target_center,
    target_depth,
    target_tangent,
    mirror=False
):
    """
    Align anatomical depth axes.

    No scaling is introduced here.

    mirror=False:
        preserve tangent direction

    mirror=True:
        reverse tangent direction
    """

    p = (
        source_xy
        - source_center
    )

    depth_coord = (
        p @ source_depth
    )

    tangent_coord = (
        p @ source_tangent
    )

    if mirror:
        tangent_coord *= -1

    transformed = (
        target_center[None, :]
        + depth_coord[:, None]
        * target_depth[None, :]
        + tangent_coord[:, None]
        * target_tangent[None, :]
    )

    return transformed


def rotation_angle_deg(
    source_depth,
    target_depth
):
    """
    Signed angle needed to rotate source depth axis
    onto target depth axis.
    """

    cross = (
        source_depth[0] * target_depth[1]
        - source_depth[1] * target_depth[0]
    )

    dot = np.dot(
        source_depth,
        target_depth
    )

    return np.degrees(
        np.arctan2(
            cross,
            dot
        )
    )


# ============================================================
# LAYER-CONCORDANCE SCORE
# ============================================================

def make_spatial_bins(
    xy,
    layers,
    origin,
    bin_um=250
):
    """
    Mean cortical depth within fixed spatial bins.
    """

    xy = np.asarray(xy)
    layers = np.asarray(layers)

    good = (
        np.isfinite(xy).all(axis=1)
        & np.isfinite(layers)
    )

    xy = xy[good]
    layers = layers[good]

    ix = np.floor(
        (xy[:, 0] - origin[0])
        / bin_um
    ).astype(int)

    iy = np.floor(
        (xy[:, 1] - origin[1])
        / bin_um
    ).astype(int)

    d = pd.DataFrame({
        "ix": ix,
        "iy": iy,
        "layer": layers
    })

    return (
        d.groupby(
            ["ix", "iy"],
            as_index=False
        )["layer"]
        .mean()
    )


def layer_concordance(
    source_xy,
    source_layer,
    target_xy,
    target_layer,
    bin_um=250
):
    """
    Spearman correlation of binned layer-depth signals.
    """

    combined = np.vstack([
        source_xy,
        target_xy
    ])

    origin = (
        np.floor(
            np.nanmin(
                combined,
                axis=0
            )
            / bin_um
        )
        * bin_um
    )

    source_bins = make_spatial_bins(
        source_xy,
        source_layer,
        origin,
        bin_um
    ).rename(
        columns={
            "layer": "source_layer"
        }
    )

    target_bins = make_spatial_bins(
        target_xy,
        target_layer,
        origin,
        bin_um
    ).rename(
        columns={
            "layer": "target_layer"
        }
    )

    shared = source_bins.merge(
        target_bins,
        on=["ix", "iy"]
    )

    if len(shared) < 3:
        return np.nan, len(shared)

    rho = spearmanr(
        shared["source_layer"],
        shared["target_layer"],
        nan_policy="omit"
    ).statistic

    return (
        float(rho),
        len(shared)
    )


# ============================================================
# OUTLINE IoU
# ============================================================

def trimmed_points(
    xy,
    q=0.0025
):
    xy = np.asarray(xy)

    good = np.isfinite(
        xy
    ).all(axis=1)

    z = xy[good]

    if len(z) < 20:
        return z

    lo = np.quantile(
        z,
        q,
        axis=0
    )

    hi = np.quantile(
        z,
        1 - q,
        axis=0
    )

    keep = (
        (z[:, 0] >= lo[0])
        & (z[:, 0] <= hi[0])
        & (z[:, 1] >= lo[1])
        & (z[:, 1] <= hi[1])
    )

    return z[keep]


def convex_polygon(xy):

    xy = trimmed_points(xy)

    hull = ConvexHull(xy)

    return xy[
        hull.vertices
    ]


def outline_iou(
    source_xy,
    target_xy,
    grid_um=50
):
    """
    Approximate tissue-outline overlap with convex hulls.
    """

    try:

        p1 = convex_polygon(
            source_xy
        )

        p2 = convex_polygon(
            target_xy
        )

        combined = np.vstack([
            p1,
            p2
        ])

        minimum = np.floor(
            np.min(
                combined,
                axis=0
            )
            / grid_um
        ) * grid_um

        maximum = np.ceil(
            np.max(
                combined,
                axis=0
            )
            / grid_um
        ) * grid_um

        xs = np.arange(
            minimum[0] + grid_um / 2,
            maximum[0],
            grid_um
        )

        ys = np.arange(
            minimum[1] + grid_um / 2,
            maximum[1],
            grid_um
        )

        xx, yy = np.meshgrid(
            xs,
            ys
        )

        points = np.c_[
            xx.ravel(),
            yy.ravel()
        ]

        mask1 = MplPath(
            p1
        ).contains_points(
            points
        )

        mask2 = MplPath(
            p2
        ).contains_points(
            points
        )

        union = np.sum(
            mask1 | mask2
        )

        if union == 0:
            return np.nan

        return float(
            np.sum(
                mask1 & mask2
            )
            / union
        )

    except Exception:
        return np.nan


# ============================================================
# PLOTTING
# ============================================================

def draw_qc(
    br,
    visium_xy,
    visium_layer,
    xenium_false,
    xenium_true,
    xenium_layer,
    rho_false,
    rho_true,
    iou_false,
    iou_true,
    bins_false,
    bins_true,
    chosen_mirror,
    rotation_deg,
    output_png,
    output_pdf
):

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(17, 6)
    )

    fig.suptitle(
        (
            f"{br} — Initial affine derived "
            f"from anatomical depth axis"
        ),
        fontsize=14,
        fontweight="bold"
    )

    all_xy = np.vstack([
        visium_xy,
        xenium_false,
        xenium_true
    ])

    mins = np.nanmin(
        all_xy,
        axis=0
    )

    maxs = np.nanmax(
        all_xy,
        axis=0
    )

    pad = (
        maxs - mins
    ) * 0.04

    mins -= pad
    maxs += pad

    candidates = [
        (
            False,
            xenium_false,
            rho_false,
            iou_false,
            bins_false
        ),
        (
            True,
            xenium_true,
            rho_true,
            iou_true,
            bins_true
        )
    ]

    for ax, (
        mirror,
        xenium_xy,
        rho,
        iou,
        n_bins
    ) in zip(
        axes[:2],
        candidates
    ):

        # Visium spots
        ax.scatter(
            visium_xy[:, 0] / 1000,
            visium_xy[:, 1] / 1000,
            c=visium_layer,
            cmap="viridis",
            s=10,
            alpha=0.75,
            linewidths=0,
            rasterized=True
        )

        # Xenium cells
        ax.scatter(
            xenium_xy[:, 0] / 1000,
            xenium_xy[:, 1] / 1000,
            c=xenium_layer,
            cmap="autumn",
            s=1.5,
            alpha=0.35,
            linewidths=0,
            rasterized=True
        )

        chosen_text = (
            "  ← chosen"
            if mirror == chosen_mirror
            else ""
        )

        ax.set_title(
            (
                f"mirror={mirror}\n"
                f"rho={rho:.3f}   "
                f"IoU={iou:.3f}   "
                f"bins={n_bins}"
                f"{chosen_text}"
            ),
            fontsize=10
        )

        ax.set_xlim(
            mins[0] / 1000,
            maxs[0] / 1000
        )

        ax.set_ylim(
            mins[1] / 1000,
            maxs[1] / 1000
        )

        ax.set_xlabel("mm")
        ax.set_ylabel("mm")

        ax.invert_yaxis()
        ax.set_aspect("equal")

    # Visium reference only
    axes[2].scatter(
        visium_xy[:, 0] / 1000,
        visium_xy[:, 1] / 1000,
        c=visium_layer,
        cmap="viridis",
        s=12,
        linewidths=0,
        rasterized=True
    )

    axes[2].set_title(
        "Visium spatial-domain reference",
        fontsize=10
    )

    axes[2].set_xlabel("mm")
    axes[2].set_ylabel("mm")

    axes[2].invert_yaxis()
    axes[2].set_aspect("equal")

    fig.text(
        0.5,
        0.015,
        (
            f"Rotation = {rotation_deg:+.1f}°. "
            f"Anatomical depth direction is aligned "
            f"between Xenium and Visium. "
            f"Both reflection states are evaluated "
            f"using cortical-layer concordance "
            f"in 250-µm spatial bins."
        ),
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

    fig.savefig(
        output_png,
        dpi=250,
        bbox_inches="tight"
    )

    fig.savefig(
        output_pdf,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# FIND SAMPLE PAIRS
# ============================================================

def discover_samples(root):
    """
    Recursively find every Br*_aligned_coordinates.csv.

    The matching visium_coordinate.csv is expected
    in the same directory.
    """

    aligned_files = sorted(
        root.rglob(
            "Br*_aligned_coordinates.csv"
        )
    )

    samples = []

    for aligned_file in aligned_files:

        br = aligned_file.name.replace(
            "_aligned_coordinates.csv",
            ""
        )

        visium_file = (
            aligned_file.parent
            / f"{br}_visium_coordinates.csv"
        )

        if not visium_file.exists():

            print(
                f"WARNING: {br}: "
                f"{br}_visium_coordinates.csv not found "
                f"in {aligned_file.parent}"
            )

            continue

        samples.append(
            (
                br,
                aligned_file,
                visium_file
            )
        )

    return samples


# ============================================================
# PROCESS ONE DONOR
# ============================================================

def process_sample(
    br,
    aligned_file,
    visium_file,
    output_dir,
    visium_pitch_um
):

    print()
    print("=" * 80)
    print(br)
    print("=" * 80)

    print(
        f"Aligned CSV : {aligned_file}"
    )

    print(
        f"Visium CSV  : {visium_file}"
    )

    xenium = pd.read_csv(
        aligned_file
    )

    visium = pd.read_csv(
        visium_file
    )

    # -------------------------
    # Xenium
    # -------------------------

    required = {
        "x",
        "y",
        "layer"
    }

    missing = (
        required
        - set(xenium.columns)
    )

    if missing:

        raise ValueError(
            f"Missing Xenium columns: {missing}"
        )

    xenium_xy = np.c_[

        pd.to_numeric(
            xenium["x"],
            errors="coerce"
        ),

        pd.to_numeric(
            xenium["y"],
            errors="coerce"
        )
    ]

    xenium_layer = (
        layer_to_ordinal(
            xenium["layer"],
            name="Xenium layer"
        )
        .to_numpy()
    )

    # -------------------------
    # Visium
    # -------------------------

    x_col, y_col = (
        detect_xy_columns(
            visium
        )
    )

    layer_col = (
        detect_layer_column(
            visium
        )
    )

    print(
        f"  Visium X/Y : "
        f"{x_col}, {y_col}"
    )

    print(
        f"  Visium layer: {layer_col}"
    )

    visium_raw_xy = np.c_[

        pd.to_numeric(
            visium[x_col],
            errors="coerce"
        ),

        pd.to_numeric(
            visium[y_col],
            errors="coerce"
        )
    ]

    (
        visium_xy,
        raw_nn,
        visium_scale
    ) = normalize_visium_coordinates(
        visium_raw_xy,
        expected_pitch_um=visium_pitch_um
    )

    visium_layer = (
        layer_to_ordinal(
            visium[layer_col],
            name="Visium layer"
        )
        .to_numpy()
    )

    print(
        f"  Visium raw NN spacing : "
        f"{raw_nn:.4f}"
    )

    print(
        f"  Visium coordinate scale -> µm: "
        f"{visium_scale:.6f}"
    )

    # -------------------------
    # remove missing
    # -------------------------

    xenium_good = (
        np.isfinite(
            xenium_xy
        ).all(axis=1)
        & np.isfinite(
            xenium_layer
        )
    )

    visium_good = (
        np.isfinite(
            visium_xy
        ).all(axis=1)
        & np.isfinite(
            visium_layer
        )
    )

    xenium_xy = (
        xenium_xy[
            xenium_good
        ]
    )

    xenium_layer = (
        xenium_layer[
            xenium_good
        ]
    )

    visium_xy = (
        visium_xy[
            visium_good
        ]
    )

    visium_layer = (
        visium_layer[
            visium_good
        ]
    )

    print(
        f"  Xenium cells used : "
        f"{len(xenium_xy):,}"
    )

    print(
        f"  Visium spots used : "
        f"{len(visium_xy):,}"
    )

    # -------------------------
    # anatomical depth axes
    # -------------------------

    (
        xenium_center,
        xenium_depth,
        xenium_tangent
    ) = fit_depth_axis(
        xenium_xy,
        xenium_layer
    )

    (
        visium_center,
        visium_depth,
        visium_tangent
    ) = fit_depth_axis(
        visium_xy,
        visium_layer
    )

    rotation_deg = (
        rotation_angle_deg(
            xenium_depth,
            visium_depth
        )
    )

    # -------------------------
    # test both mirror states
    # -------------------------

    candidate_false = (
        initial_affine(
            xenium_xy,
            xenium_center,
            xenium_depth,
            xenium_tangent,
            visium_center,
            visium_depth,
            visium_tangent,
            mirror=False
        )
    )

    candidate_true = (
        initial_affine(
            xenium_xy,
            xenium_center,
            xenium_depth,
            xenium_tangent,
            visium_center,
            visium_depth,
            visium_tangent,
            mirror=True
        )
    )

    rho_false, bins_false = (
        layer_concordance(
            candidate_false,
            xenium_layer,
            visium_xy,
            visium_layer
        )
    )

    rho_true, bins_true = (
        layer_concordance(
            candidate_true,
            xenium_layer,
            visium_xy,
            visium_layer
        )
    )

    iou_false = outline_iou(
        candidate_false,
        visium_xy
    )

    iou_true = outline_iou(
        candidate_true,
        visium_xy
    )

    # -------------------------
    # choose best mirror
    # -------------------------

    if (
        np.isfinite(rho_false)
        and np.isfinite(rho_true)
    ):

        if rho_true > rho_false:
            chosen_mirror = True
        else:
            chosen_mirror = False

    elif np.isfinite(rho_true):

        chosen_mirror = True

    else:

        chosen_mirror = False

    print()
    print(
        f"  Rotation        : "
        f"{rotation_deg:+.2f}°"
    )

    print(
        f"  mirror=False    : "
        f"rho={rho_false:.4f}, "
        f"IoU={iou_false:.4f}, "
        f"bins={bins_false}"
    )

    print(
        f"  mirror=True     : "
        f"rho={rho_true:.4f}, "
        f"IoU={iou_true:.4f}, "
        f"bins={bins_true}"
    )

    print(
        f"  CHOSEN          : "
        f"mirror={chosen_mirror}"
    )

    # -------------------------
    # output
    # -------------------------

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    output_png = (
        output_dir
        / f"{br}_initial_affine_qc.png"
    )

    output_pdf = (
        output_dir
        / f"{br}_initial_affine_qc.pdf"
    )

    draw_qc(
        br=br,
        visium_xy=visium_xy,
        visium_layer=visium_layer,
        xenium_false=candidate_false,
        xenium_true=candidate_true,
        xenium_layer=xenium_layer,
        rho_false=rho_false,
        rho_true=rho_true,
        iou_false=iou_false,
        iou_true=iou_true,
        bins_false=bins_false,
        bins_true=bins_true,
        chosen_mirror=chosen_mirror,
        rotation_deg=rotation_deg,
        output_png=output_png,
        output_pdf=output_pdf
    )

    return {
        "BrNum": br,
        "status": "PASS",

        "aligned_csv": str(
            aligned_file
        ),

        "visium_csv": str(
            visium_file
        ),

        "n_xenium_cells": len(
            xenium_xy
        ),

        "n_visium_spots": len(
            visium_xy
        ),

        "visium_x_col": x_col,
        "visium_y_col": y_col,
        "visium_layer_col": layer_col,

        "visium_raw_nn_spacing": raw_nn,
        "visium_to_um_scale": visium_scale,

        "rotation_deg": rotation_deg,

        "rho_mirror_false": rho_false,
        "rho_mirror_true": rho_true,

        "iou_mirror_false": iou_false,
        "iou_mirror_true": iou_true,

        "shared_bins_mirror_false": bins_false,
        "shared_bins_mirror_true": bins_true,

        "chosen_mirror": chosen_mirror,

        "figure_png": str(
            output_png
        ),

        "figure_pdf": str(
            output_pdf
        )
    }


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_INPUT_ROOT,
        help=(
            "Root directory to recursively search "
            "for Br*_aligned_coordinates.csv"
        )
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Output directory"
    )

    parser.add_argument(
        "--visium-pitch-um",
        type=float,
        default=100.0,
        help=(
            "Expected Visium spot centre-to-centre "
            "distance in microns"
        )
    )

    parser.add_argument(
        "--br",
        type=str,
        default=None,
        help=(
            "Optional single donor, e.g. Br1113"
        )
    )

    args = parser.parse_args()

    root = args.root.resolve()
    output_dir = args.out.resolve()

    print("=" * 80)
    print("STALIGN INITIAL AFFINE QC")
    print("=" * 80)

    print(
        f"Search root : {root}"
    )

    print(
        f"Output      : {output_dir}"
    )

    samples = discover_samples(
        root
    )

    if args.br is not None:

        samples = [
            x
            for x in samples
            if x[0] == args.br
        ]

    print(
        f"Samples found: {len(samples)}"
    )

    if not samples:

        raise RuntimeError(
            "No valid sample pairs found."
        )

    results = []

    for (
        br,
        aligned_file,
        visium_file
    ) in samples:

        try:

            result = process_sample(
                br=br,
                aligned_file=aligned_file,
                visium_file=visium_file,
                output_dir=output_dir,
                visium_pitch_um=args.visium_pitch_um
            )

        except Exception as exc:

            print()
            print(
                f"ERROR: {br}"
            )

            print(
                str(exc)
            )

            traceback.print_exc()

            result = {
                "BrNum": br,
                "status": "FAIL",
                "error": str(exc)
            }

        results.append(
            result
        )

    summary = pd.DataFrame(
        results
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    summary_file = (
        output_dir
        / "initial_affine_qc_summary.csv"
    )

    summary.to_csv(
        summary_file,
        index=False
    )

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    summary_cols = [
        c
        for c in [
            "BrNum",
            "status",
            "rotation_deg",
            "rho_mirror_false",
            "rho_mirror_true",
            "iou_mirror_false",
            "iou_mirror_true",
            "chosen_mirror"
        ]
        if c in summary.columns
    ]

    print(
        summary[
            summary_cols
        ].to_string(
            index=False
        )
    )

    print()
    print(
        f"Summary CSV: {summary_file}"
    )

    print(
        f"QC figures : {output_dir}"
    )


if __name__ == "__main__":
    main()
