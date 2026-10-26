#!/usr/bin/env python3

"""
02_make_final_aligned_qc.py

Advisor-style post-hoc QC of exported Xenium -> Visium STalign results.

IMPORTANT
---------
This script DOES NOT rerun STalign/LDDMM.

It uses the already exported:
    x_aligned
    y_aligned
    match_prob
    layer
    spot_dist_um
    inside_visium
    usable

and the corresponding:
    BrXXXX_visium_coordinates.csv

It reproduces the advisor's Cell-8 style evaluation as far as possible
from the exported CSV files:

    match_prob > 0.85
    nearest Visium spot
    within actual spot radius (~32.54 um)
    valid Xenium layer
    valid Visium domain

Metrics:
    exact layer agreement
    +/- 1 layer agreement
    Spearman layer concordance
    permutation baseline
    GM/WM coarse agreement
    within-GM layer agreement
    within-WM layer agreement
    per-layer landing/bias
    geometric/export QC

Also creates:
    final alignment figure
    layer landing figure
    per-layer CSV
    cohort summary CSV

Layer ordering follows advisor notebook exactly, ZERO-BASED:

    0  L1/M   SpD07
    1  L2/3   SpD06
    2  L3/4   SpD02
    3  L5     SpD05
    4  L6     SpD03
    5  WM     SpD04
    6  WMtz   SpD01

Author: GeneBridge
"""

from pathlib import Path
import argparse
import traceback

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from scipy.stats import spearmanr


# ============================================================
# PROJECT PATHS
# ============================================================

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
    / "stalign_final_aligned_qc_advisor_style"
)


# ============================================================
# PARAMETERS
# ============================================================

MATCH_THRESH = 0.85

# Advisor v4 value.
# Derived there from SpaceRanger spot_diameter_fullres.
DEFAULT_SPOT_RADIUS_UM = 32.54

SPOT_PITCH_UM = 100.0

N_PERMUTATIONS = 200


# ============================================================
# ADVISOR LAYER ORDER -- ZERO BASED
# ============================================================

LAYER_NAMES = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WM",
    "WMtz",
]

SPD_NAMES = [
    "SpD07-L1",
    "SpD06-L2/3",
    "SpD02-L3/4",
    "SpD05-L5",
    "SpD03-L6",
    "SpD04-WM",
    "SpD01-WMtz",
]


def normalize_label(x):
    if pd.isna(x):
        return None

    return (
        str(x)
        .strip()
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
    )


LABEL_MAP = {
    # ----------------------------------
    # L1/M = 0
    # ----------------------------------
    "l1/m": 0,
    "l1m": 0,
    "l1": 0,
    "layer1": 0,
    "spd07": 0,
    "spd07l1": 0,
    "spd07l1/m": 0,

    # ----------------------------------
    # L2/3 = 1
    # ----------------------------------
    "l2/3": 1,
    "l23": 1,
    "layer2/3": 1,
    "layer23": 1,
    "spd06": 1,
    "spd06l2/3": 1,
    "spd06l23": 1,

    # ----------------------------------
    # L3/4 = 2
    # ----------------------------------
    "l3/4": 2,
    "l34": 2,
    "layer3/4": 2,
    "layer34": 2,
    "spd02": 2,
    "spd02l3/4": 2,
    "spd02l34": 2,

    # ----------------------------------
    # L5 = 3
    # ----------------------------------
    "l5": 3,
    "layer5": 3,
    "spd05": 3,
    "spd05l5": 3,

    # ----------------------------------
    # L6 = 4
    # ----------------------------------
    "l6": 4,
    "layer6": 4,
    "spd03": 4,
    "spd03l6": 4,

    # ----------------------------------
    # WM = 5
    # ----------------------------------
    "wm": 5,
    "whitematter": 5,
    "spd04": 5,
    "spd04wm": 5,

    # ----------------------------------
    # WMtz = 6
    # ----------------------------------
    "wmtz": 6,
    "wmtransition": 6,
    "whitemattertransition": 6,
    "spd01": 6,
    "spd01wmtz": 6,
}


def layer_to_ordinal(series):
    """
    Convert Xenium layer or Visium domain to advisor-style
    zero-based ordinal 0..6.
    """

    result = pd.Series(
        np.nan,
        index=series.index,
        dtype=float,
    )

    # Allow zero-based numeric values already present
    numeric = pd.to_numeric(
        series,
        errors="coerce",
    )

    numeric_mask = numeric.between(
        0,
        6,
    )

    result.loc[numeric_mask] = (
        numeric.loc[numeric_mask]
    )

    for idx, value in series.items():

        if pd.notna(result.loc[idx]):
            continue

        key = normalize_label(value)

        if key in LABEL_MAP:
            result.loc[idx] = LABEL_MAP[key]

    return result


def parse_bool(series):

    if series.dtype == bool:
        return series.to_numpy()

    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map(
            {
                "true": True,
                "false": False,
                "1": True,
                "0": False,
                "yes": True,
                "no": False,
            }
        )
        .fillna(False)
        .astype(bool)
        .to_numpy()
    )


def safe_spearman(x, y):

    if len(x) < 3:
        return np.nan

    if (
        len(np.unique(x)) < 2
        or len(np.unique(y)) < 2
    ):
        return np.nan

    return float(
        spearmanr(
            x,
            y,
            nan_policy="omit",
        ).statistic
    )


# ============================================================
# DISCOVERY
# ============================================================

def discover_samples(root, wanted_br=None):

    samples = []

    aligned_files = sorted(
        root.rglob(
            "Br*_aligned_coordinates.csv"
        )
    )

    for aligned_file in aligned_files:

        br = (
            aligned_file.name
            .replace(
                "_aligned_coordinates.csv",
                "",
            )
        )

        if (
            wanted_br is not None
            and br != wanted_br
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
                visium_file,
            )
        )

    return samples


# ============================================================
# PERMUTATION BASELINE
# ============================================================

def permutation_baseline(
    xo,
    vo,
    n_perm=N_PERMUTATIONS,
):

    rng = np.random.default_rng(0)

    baseline = np.array(
        [
            (
                np.abs(
                    rng.permutation(xo)
                    - vo
                )
                <= 1
            ).mean()
            for _ in range(n_perm)
        ]
    )

    mean = float(
        baseline.mean()
    )

    sd = float(
        baseline.std()
    )

    observed = float(
        (
            np.abs(xo - vo)
            <= 1
        ).mean()
    )

    if sd > 0:
        lift = (
            observed - mean
        ) / sd
    else:
        lift = np.nan

    return (
        mean,
        sd,
        float(lift),
    )


# ============================================================
# FINE COMPARTMENT METRICS
# ============================================================

def fine_metrics(x, v):

    if len(x) == 0:

        return {
            "n": 0,
            "exact": np.nan,
            "pm1": np.nan,
            "rho": np.nan,
        }

    return {
        "n": len(x),

        "exact": float(
            (x == v).mean()
        ),

        "pm1": float(
            (
                np.abs(x - v)
                <= 1
            ).mean()
        ),

        "rho": safe_spearman(
            x,
            v,
        ),
    }


# ============================================================
# PER-LAYER METRICS
# ============================================================

def make_layer_table(
    br,
    xo,
    vo,
):

    rows = []

    for layer_idx, name in enumerate(
        LAYER_NAMES
    ):

        sel = (
            xo == layer_idx
        )

        n = int(
            sel.sum()
        )

        if n == 0:

            rows.append(
                {
                    "BrNum": br,
                    "layer_index": layer_idx,
                    "layer": name,
                    "n": 0,
                    "exact": np.nan,
                    "pm1": np.nan,
                    "mean_landed": np.nan,
                    "sd_landed": np.nan,
                    "bias": np.nan,
                }
            )

            continue

        landed = (
            vo[sel]
        )

        exact = float(
            (
                landed
                == layer_idx
            ).mean()
        )

        pm1 = float(
            (
                np.abs(
                    landed
                    - layer_idx
                )
                <= 1
            ).mean()
        )

        mean_landed = float(
            np.mean(
                landed
            )
        )

        sd_landed = float(
            np.std(
                landed
            )
        )

        bias = float(
            mean_landed
            - layer_idx
        )

        rows.append(
            {
                "BrNum": br,
                "layer_index": layer_idx,
                "layer": name,
                "n": n,
                "exact": exact,
                "pm1": pm1,
                "mean_landed": mean_landed,
                "sd_landed": sd_landed,
                "bias": bias,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# FIGURE 1
# Advisor-style final alignment plot
# ============================================================

def plot_final_alignment(
    br,
    xen_xy,
    xen_ord,
    vis_xy,
    vis_ord,
    match_prob,
    matched,
    score_mask,
    result,
    outdir,
):

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 6),
    )

    fig.suptitle(
        (
            f"{br} — final Xenium → Visium "
            f"STalign coordinates"
        ),
        fontsize=14,
        fontweight="bold",
    )

    valid_x = (
        np.isfinite(
            xen_xy
        ).all(axis=1)
    )

    valid_v = (
        np.isfinite(
            vis_xy
        ).all(axis=1)
        & np.isfinite(
            vis_ord
        )
    )

    # --------------------------------------------------------
    # Panel 1: continuous STalign match probability
    # --------------------------------------------------------

    p1 = axes[0].scatter(
        xen_xy[
            valid_x,
            0
        ] / 1000,
        xen_xy[
            valid_x,
            1
        ] / 1000,
        c=match_prob[
            valid_x
        ],
        cmap="RdYlGn",
        vmin=0,
        vmax=1,
        s=0.8,
        linewidths=0,
        rasterized=True,
    )

    axes[0].set_title(
        (
            "STalign match probability\n"
            f"{result['matched_pct']:.1f}% > "
            f"{MATCH_THRESH}"
        )
    )

    fig.colorbar(
        p1,
        ax=axes[0],
        fraction=0.046,
    )

    # --------------------------------------------------------
    # Panel 2: matched Xenium over Visium
    # --------------------------------------------------------

    axes[1].scatter(
        vis_xy[
            valid_v,
            0
        ] / 1000,
        vis_xy[
            valid_v,
            1
        ] / 1000,
        c=vis_ord[
            valid_v
        ],
        cmap="viridis",
        vmin=0,
        vmax=6,
        s=9,
        alpha=0.45,
        linewidths=0,
        rasterized=True,
    )

    matched_valid = (
        valid_x
        & matched
        & np.isfinite(
            xen_ord
        )
    )

    axes[1].scatter(
        xen_xy[
            matched_valid,
            0
        ] / 1000,
        xen_xy[
            matched_valid,
            1
        ] / 1000,
        c=xen_ord[
            matched_valid
        ],
        cmap="turbo",
        vmin=0,
        vmax=6,
        s=0.8,
        alpha=0.45,
        linewidths=0,
        rasterized=True,
    )

    axes[1].set_title(
        (
            "Matched Xenium over Visium\n"
            f"match_prob > {MATCH_THRESH}"
        )
    )

    # --------------------------------------------------------
    # Panel 3: Visium reference
    # --------------------------------------------------------

    axes[2].scatter(
        vis_xy[
            valid_v,
            0
        ] / 1000,
        vis_xy[
            valid_v,
            1
        ] / 1000,
        c=vis_ord[
            valid_v
        ],
        cmap="viridis",
        vmin=0,
        vmax=6,
        s=10,
        linewidths=0,
        rasterized=True,
    )

    axes[2].set_title(
        "Visium spatial-domain reference"
    )

    # --------------------------------------------------------
    # shared ranges
    # --------------------------------------------------------

    combined = np.vstack(
        [
            xen_xy[
                valid_x
            ],
            vis_xy[
                valid_v
            ],
        ]
    )

    mins = np.nanmin(
        combined,
        axis=0,
    )

    maxs = np.nanmax(
        combined,
        axis=0,
    )

    pad = (
        maxs - mins
    ) * 0.03

    for ax in axes:

        ax.set_xlim(
            (
                mins[0]
                - pad[0]
            ) / 1000,
            (
                maxs[0]
                + pad[0]
            ) / 1000,
        )

        ax.set_ylim(
            (
                mins[1]
                - pad[1]
            ) / 1000,
            (
                maxs[1]
                + pad[1]
            ) / 1000,
        )

        ax.invert_yaxis()

        ax.set_aspect(
            "equal"
        )

        ax.set_xlabel(
            "mm"
        )

        ax.set_ylabel(
            "mm"
        )

    footer = (
        f"scored in physical Visium spots: "
        f"n={result['n_scored']:,} "
        f"({result['scored_pct']:.1f}%); "
        f"exact={100*result['exact']:.1f}%; "
        f"±1={100*result['pm1']:.1f}%; "
        f"Spearman={result['rho']:.3f}; "
        f"GM/WM={100*result['gmwm_accuracy']:.1f}%"
    )

    fig.text(
        0.5,
        0.015,
        footer,
        ha="center",
        fontsize=8,
    )

    plt.tight_layout(
        rect=[
            0,
            0.055,
            1,
            0.94,
        ]
    )

    png = (
        outdir
        / f"{br}_final_alignment_advisor_style.png"
    )

    pdf = (
        outdir
        / f"{br}_final_alignment_advisor_style.pdf"
    )

    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    plt.close(fig)

    return png


# ============================================================
# FIGURE 2
# Where each Xenium layer lands
# ============================================================

def plot_layer_landing(
    br,
    xen_xy,
    xen_ord,
    vis_xy,
    vis_ord,
    score_mask,
    nearest_idx,
    layer_table,
    outdir,
):

    fig, axes = plt.subplots(
        2,
        4,
        figsize=(18, 10),
    )

    axes = axes.ravel()

    valid_vis = (
        np.isfinite(
            vis_xy
        ).all(axis=1)
        & np.isfinite(
            vis_ord
        )
    )

    for layer_idx, name in enumerate(
        LAYER_NAMES
    ):

        ax = axes[
            layer_idx
        ]

        # all Visium spots in grey
        ax.scatter(
            vis_xy[
                valid_vis,
                0
            ] / 1000,
            vis_xy[
                valid_vis,
                1
            ] / 1000,
            s=6,
            c="lightgrey",
            alpha=0.35,
            linewidths=0,
            rasterized=True,
        )

        # Visium spots of this same layer
        target_layer = (
            valid_vis
            & (
                vis_ord
                == layer_idx
            )
        )

        ax.scatter(
            vis_xy[
                target_layer,
                0
            ] / 1000,
            vis_xy[
                target_layer,
                1
            ] / 1000,
            s=18,
            facecolors="none",
            edgecolors="black",
            linewidths=0.7,
            rasterized=True,
        )

        # scored Xenium cells of this layer
        cells = (
            score_mask
            & (
                xen_ord
                == layer_idx
            )
        )

        ax.scatter(
            xen_xy[
                cells,
                0
            ] / 1000,
            xen_xy[
                cells,
                1
            ] / 1000,
            s=1.2,
            alpha=0.55,
            rasterized=True,
        )

        row = (
            layer_table
            .iloc[layer_idx]
        )

        ax.set_title(
            (
                f"{name}   n={int(row['n']):,}\n"
                f"exact={100*row['exact']:.1f}%  "
                f"±1={100*row['pm1']:.1f}%  "
                f"bias={row['bias']:+.2f}"
            )
            if row["n"] > 0
            else f"{name}   n=0"
        )

        ax.invert_yaxis()
        ax.set_aspect(
            "equal"
        )

    # final panel = explanation
    axes[7].axis(
        "off"
    )

    axes[7].text(
        0.02,
        0.98,
        (
            "coloured dots = aligned Xenium cells\n"
            "black rings = Visium spots of same layer\n"
            "grey = all Visium spots\n\n"
            "Cells close to same-layer rings indicate\n"
            "laminar correspondence.\n\n"
            "Systematic displacement indicates a\n"
            "laminar offset."
        ),
        va="top",
        fontsize=10,
    )

    fig.suptitle(
        (
            "Where each Xenium layer lands "
            "on the Visium section — "
            f"{br}"
        ),
        fontsize=13,
        fontweight="bold",
    )

    plt.tight_layout(
        rect=[
            0,
            0,
            1,
            0.96,
        ]
    )

    png = (
        outdir
        / f"{br}_layer_landing.png"
    )

    pdf = (
        outdir
        / f"{br}_layer_landing.pdf"
    )

    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    plt.close(fig)


# ============================================================
# ONE SAMPLE
# ============================================================

def process_sample(
    br,
    aligned_file,
    visium_file,
    outdir,
    spot_radius_um,
):

    print()
    print("=" * 90)
    print(br)
    print("=" * 90)

    print(
        f"Aligned : {aligned_file}"
    )

    print(
        f"Visium  : {visium_file}"
    )

    xen = pd.read_csv(
        aligned_file
    )

    vis = pd.read_csv(
        visium_file
    )

    # ========================================================
    # REQUIRED COLUMNS
    # ========================================================

    required_xen = [
        "x_aligned",
        "y_aligned",
        "match_prob",
        "layer",
    ]

    required_vis = [
        "x",
        "y",
        "domain",
    ]

    missing_x = [
        c
        for c in required_xen
        if c not in xen.columns
    ]

    missing_v = [
        c
        for c in required_vis
        if c not in vis.columns
    ]

    if missing_x:
        raise ValueError(
            f"Missing Xenium columns: {missing_x}"
        )

    if missing_v:
        raise ValueError(
            f"Missing Visium columns: {missing_v}"
        )

    # ========================================================
    # FINAL COORDINATES
    # ========================================================

    xen_xy = np.c_[
        pd.to_numeric(
            xen["x_aligned"],
            errors="coerce",
        ),
        pd.to_numeric(
            xen["y_aligned"],
            errors="coerce",
        ),
    ]

    vis_xy = np.c_[
        pd.to_numeric(
            vis["x"],
            errors="coerce",
        ),
        pd.to_numeric(
            vis["y"],
            errors="coerce",
        ),
    ]

    match_prob = (
        pd.to_numeric(
            xen["match_prob"],
            errors="coerce",
        )
        .to_numpy()
    )

    xen_ord = (
        layer_to_ordinal(
            xen["layer"]
        )
        .to_numpy()
    )

    vis_ord = (
        layer_to_ordinal(
            vis["domain"]
        )
        .to_numpy()
    )

    # ========================================================
    # COORDINATE VALIDITY
    # ========================================================

    valid_x_coord = (
        np.isfinite(
            xen_xy
        ).all(axis=1)
    )

    valid_v_coord = (
        np.isfinite(
            vis_xy
        ).all(axis=1)
    )

    valid_vis_rows = (
        valid_v_coord
    )

    if valid_vis_rows.sum() == 0:
        raise ValueError(
            "No valid Visium coordinates."
        )

    # ========================================================
    # NEAREST VISIUM SPOT
    # Exactly analogous to advisor Cell 8 KD-tree step.
    # ========================================================

    vis_tree_xy = (
        vis_xy[
            valid_vis_rows
        ]
    )

    vis_tree_index = np.where(
        valid_vis_rows
    )[0]

    tree = cKDTree(
        vis_tree_xy
    )

    nearest_dist = np.full(
        len(xen),
        np.nan,
    )

    nearest_idx = np.full(
        len(xen),
        -1,
        dtype=int,
    )

    query_rows = np.where(
        valid_x_coord
    )[0]

    dtmp, itmp = tree.query(
        xen_xy[
            query_rows
        ],
        k=1,
    )

    nearest_dist[
        query_rows
    ] = dtmp

    nearest_idx[
        query_rows
    ] = (
        vis_tree_index[
            itmp
        ]
    )

    # ========================================================
    # FILTERS
    #
    # Advisor:
    #
    # keep    = wm > 0.85
    # near8   = distance <= SPOT_R_UM
    # has_lab = Xenium layer not null
    #
    # sel8 = near8 & keep & has_lab & valid Visium domain
    # ========================================================

    matched = (
        np.isfinite(
            match_prob
        )
        & (
            match_prob
            > MATCH_THRESH
        )
    )

    in_spot = (
        np.isfinite(
            nearest_dist
        )
        & (
            nearest_dist
            <= spot_radius_um
        )
    )

    has_xen_layer = (
        np.isfinite(
            xen_ord
        )
    )

    target_layer = np.full(
        len(xen),
        np.nan,
    )

    has_nearest = (
        nearest_idx
        >= 0
    )

    target_layer[
        has_nearest
    ] = (
        vis_ord[
            nearest_idx[
                has_nearest
            ]
        ]
    )

    has_vis_layer = (
        np.isfinite(
            target_layer
        )
    )

    score_mask = (
        valid_x_coord
        & matched
        & in_spot
        & has_xen_layer
        & has_vis_layer
    )

    if score_mask.sum() == 0:

        raise RuntimeError(
            "No cells survived advisor-style scoring filter."
        )

    xo = (
        xen_ord[
            score_mask
        ]
        .astype(int)
    )

    vo = (
        target_layer[
            score_mask
        ]
        .astype(int)
    )

    dd = (
        nearest_dist[
            score_mask
        ]
    )

    # ========================================================
    # OVERALL METRICS
    # ========================================================

    exact = float(
        (
            xo == vo
        ).mean()
    )

    pm1 = float(
        (
            np.abs(
                xo - vo
            )
            <= 1
        ).mean()
    )

    rho = safe_spearman(
        xo,
        vo,
    )

    (
        perm_mean,
        perm_sd,
        perm_lift,
    ) = permutation_baseline(
        xo,
        vo,
    )

    # ========================================================
    # COARSE GM / WM
    #
    # Advisor ZERO-BASED definition:
    #
    #   0 = L1/M -- excluded
    #   1..4 = GM = L2/3 ... L6
    #   5..6 = WM = WM, WMtz
    #
    # The following intentionally reproduces the advisor logic:
    # target "not GM" counts as the WM side for the coarse metric.
    # ========================================================

    GM_I = [
        1,
        2,
        3,
        4,
    ]

    WM_I = [
        5,
        6,
    ]

    is_gm_x = np.isin(
        xo,
        GM_I,
    )

    is_wm_x = np.isin(
        xo,
        WM_I,
    )

    is_gm_v = np.isin(
        vo,
        GM_I,
    )

    is_wm_v = np.isin(
        vo,
        WM_I,
    )

    coarse_sel = (
        is_gm_x
        | is_wm_x
    )

    gmwm_accuracy = float(
        (
            is_gm_x[
                coarse_sel
            ]
            == is_gm_v[
                coarse_sel
            ]
        ).mean()
    )

    gm_to_gm = float(
        is_gm_v[
            is_gm_x
        ].mean()
    ) if is_gm_x.any() else np.nan

    # exact advisor interpretation:
    # source WM -> target not-GM
    wm_to_wm_advisor = float(
        (
            ~is_gm_v[
                is_wm_x
            ]
        ).mean()
    ) if is_wm_x.any() else np.nan

    # stricter additional check:
    # source WM -> target specifically WM/WMtz
    wm_to_wm_strict = float(
        is_wm_v[
            is_wm_x
        ].mean()
    ) if is_wm_x.any() else np.nan

    # ========================================================
    # FINE METRICS GIVEN CORRECT COMPARTMENT
    # ========================================================

    gm_correct = (
        is_gm_x
        & is_gm_v
    )

    # advisor logic
    wm_correct = (
        is_wm_x
        & ~is_gm_v
    )

    fine_gm = fine_metrics(
        xo[
            gm_correct
        ],
        vo[
            gm_correct
        ],
    )

    fine_wm = fine_metrics(
        xo[
            wm_correct
        ],
        vo[
            wm_correct
        ],
    )

    # ========================================================
    # PER-LAYER
    # ========================================================

    layer_table = make_layer_table(
        br,
        xo,
        vo,
    )

    layer_csv = (
        outdir
        / f"{br}_layer_landing_metrics.csv"
    )

    layer_table.to_csv(
        layer_csv,
        index=False,
    )

    # GM biases = L2/3 through L6
    gm_biases = (
        layer_table.loc[
            layer_table[
                "layer_index"
            ].isin(
                GM_I
            ),
            "bias",
        ]
        .dropna()
        .to_numpy()
    )

    gm_bias_mean = (
        float(
            np.mean(
                gm_biases
            )
        )
        if len(gm_biases)
        else np.nan
    )

    gm_bias_same_sign = (
        bool(
            np.all(
                gm_biases > 0
            )
            or np.all(
                gm_biases < 0
            )
        )
        if len(gm_biases)
        else False
    )

    # ========================================================
    # EXPORTED GEOMETRY METRICS
    # These are retained from our previous cohort QC.
    # ========================================================

    inside_n = np.nan
    inside_pct = np.nan

    if "inside_visium" in xen.columns:

        inside = parse_bool(
            xen[
                "inside_visium"
            ]
        )

        inside_n = int(
            inside.sum()
        )

        inside_pct = float(
            100
            * inside.mean()
        )

    usable_n = np.nan
    usable_pct = np.nan

    if "usable" in xen.columns:

        usable = parse_bool(
            xen[
                "usable"
            ]
        )

        usable_n = int(
            usable.sum()
        )

        usable_pct = float(
            100
            * usable.mean()
        )

    assigned_pct = np.nan

    if "visium_spot" in xen.columns:

        spot = (
            xen[
                "visium_spot"
            ]
            .astype(str)
            .str.strip()
        )

        assigned = (
            xen[
                "visium_spot"
            ]
            .notna()
            & (
                spot != ""
            )
            & (
                spot.str.lower()
                != "nan"
            )
        )

        assigned_pct = float(
            100
            * assigned.mean()
        )

    # ========================================================
    # EXPORTED vs RECOMPUTED DISTANCE CHECK
    # ========================================================

    exported_distance_mae = np.nan

    if "spot_dist_um" in xen.columns:

        exported_dist = (
            pd.to_numeric(
                xen[
                    "spot_dist_um"
                ],
                errors="coerce",
            )
            .to_numpy()
        )

        compare = (
            np.isfinite(
                exported_dist
            )
            & np.isfinite(
                nearest_dist
            )
        )

        if compare.any():

            exported_distance_mae = float(
                np.median(
                    np.abs(
                        exported_dist[
                            compare
                        ]
                        - nearest_dist[
                            compare
                        ]
                    )
                )
            )

    # ========================================================
    # SELF-CHECK
    # Based on advisor notebook.
    # ========================================================

    all_layers_present = (
        len(
            np.unique(
                xo
            )
        )
        == 7
    )

    scored_fraction = float(
        score_mask.mean()
    )

    checks = {
        "enough_cells_gt_5000":
            len(xo) > 5000,

        "pm1_above_baseline_5sd":
            (
                pm1
                > (
                    perm_mean
                    + 5
                    * perm_sd
                )
            ),

        "gmwm_binary_gt_70pct":
            gmwm_accuracy > 0.70,

        "spearman_gt_0_3":
            (
                np.isfinite(
                    rho
                )
                and rho > 0.30
            ),

        "all_7_layers_represented":
            all_layers_present,

        "in_spot_scored_fraction_gt_20pct":
            scored_fraction > 0.20,

        "all_visium_domains_mapped":
            (
                np.isfinite(
                    vis_ord
                ).sum()
                == len(
                    vis_ord
                )
            ),

        "match_prob_is_numeric":
            np.issubdtype(
                match_prob.dtype,
                np.number,
            ),
    }

    all_pass = bool(
        all(
            checks.values()
        )
    )

    # ========================================================
    # RESULT RECORD
    # ========================================================

    result = {
        "BrNum": br,
        "status": "PASS",

        "n_cells": len(
            xen
        ),

        "n_visium_spots": len(
            vis
        ),

        "match_threshold": MATCH_THRESH,
        "spot_radius_um": spot_radius_um,

        "n_matched": int(
            matched.sum()
        ),

        "matched_pct": float(
            100
            * matched.mean()
        ),

        "n_within_spot_radius": int(
            in_spot.sum()
        ),

        "within_spot_radius_pct": float(
            100
            * in_spot.mean()
        ),

        "n_has_xenium_layer": int(
            has_xen_layer.sum()
        ),

        "n_scored": int(
            score_mask.sum()
        ),

        "scored_pct": float(
            100
            * scored_fraction
        ),

        "median_distance_matched_um": float(
            np.nanmedian(
                nearest_dist[
                    matched
                ]
            )
        ),

        "median_distance_scored_um": float(
            np.nanmedian(
                dd
            )
        ),

        "exact": exact,
        "pm1": pm1,
        "rho": rho,

        "permutation_pm1_mean": perm_mean,
        "permutation_pm1_sd": perm_sd,
        "permutation_lift_sigma": perm_lift,

        "gmwm_n": int(
            coarse_sel.sum()
        ),

        "gmwm_accuracy": gmwm_accuracy,

        "gm_to_gm": gm_to_gm,

        "wm_to_wm_advisor": wm_to_wm_advisor,
        "wm_to_wm_strict": wm_to_wm_strict,

        "fine_gm_n": fine_gm[
            "n"
        ],

        "fine_gm_exact": fine_gm[
            "exact"
        ],

        "fine_gm_pm1": fine_gm[
            "pm1"
        ],

        "fine_gm_rho": fine_gm[
            "rho"
        ],

        "fine_wm_n": fine_wm[
            "n"
        ],

        "fine_wm_exact": fine_wm[
            "exact"
        ],

        "fine_wm_pm1": fine_wm[
            "pm1"
        ],

        "fine_wm_rho": fine_wm[
            "rho"
        ],

        "gm_bias_mean": gm_bias_mean,

        "gm_bias_all_same_sign":
            gm_bias_same_sign,

        "inside_visium_n": inside_n,
        "inside_visium_pct": inside_pct,

        "usable_n": usable_n,
        "usable_pct": usable_pct,

        "assigned_to_spot_pct":
            assigned_pct,

        "exported_vs_recomputed_dist_median_abs_diff_um":
            exported_distance_mae,

        "self_check_all_pass":
            all_pass,
    }

    # ========================================================
    # TERMINAL REPORT
    # ========================================================

    print(
        f"cells                         : "
        f"{len(xen):,}"
    )

    print(
        f"  matched (match_prob > "
        f"{MATCH_THRESH})       : "
        f"{matched.sum():6d} "
        f"({matched.mean()*100:5.1f}%)"
    )

    print(
        f"  within {spot_radius_um:.2f} um "
        f"of a spot      : "
        f"{in_spot.sum():6d} "
        f"({in_spot.mean()*100:5.1f}%)"
    )

    print(
        f"  has a Xenium layer            : "
        f"{has_xen_layer.sum():6d} "
        f"({has_xen_layer.mean()*100:5.1f}%)"
    )

    print(
        f"  scored here (all above)       : "
        f"{score_mask.sum():6d} "
        f"({score_mask.mean()*100:5.1f}%)"
    )

    print(
        f"  median nearest distance "
        f"(matched) : "
        f"{result['median_distance_matched_um']:.1f} um"
    )

    spot_area_fraction = (
        np.pi
        * spot_radius_um ** 2
        / (
            SPOT_PITCH_UM
            * SPOT_PITCH_UM
            * np.sqrt(3)
            / 2
        )
    )

    print(
        f"  spot discs cover ~"
        f"{spot_area_fraction*100:.0f}% "
        f"of tissue area"
    )

    print()
    print(
        f"--- overall "
        f"(n = {len(xo):,}) ---"
    )

    print(
        f"  exact {exact*100:.1f}%"
        f"  |  +/-1 {pm1*100:.1f}%"
        f"  |  Spearman {rho:.3f}"
    )

    print(
        f"  +/-1 permutation baseline "
        f"{perm_mean*100:.1f}%"
        f"  ->  {perm_lift:.1f} sigma"
    )

    print()
    print(
        f"--- coarse: GM/WM binary "
        f"(n = {coarse_sel.sum():,}, "
        f"L1/M excluded) ---"
    )

    print(
        f"  landed on correct side : "
        f"{gmwm_accuracy*100:.1f}%"
    )

    print(
        f"  GM cells -> GM spot    : "
        f"n={is_gm_x.sum():6d}  "
        f"{gm_to_gm*100:.1f}%"
    )

    print(
        f"  WM cells -> WM side    : "
        f"n={is_wm_x.sum():6d}  "
        f"{wm_to_wm_advisor*100:.1f}%"
    )

    print(
        f"  WM cells -> actual "
        f"WM/WMtz spot (strict): "
        f"{wm_to_wm_strict*100:.1f}%"
    )

    print()
    print(
        "--- fine: laminar position "
        "given correct compartment ---"
    )

    print(
        f"  within GM "
        f"n={fine_gm['n']:6d}  "
        f"exact {fine_gm['exact']*100:5.1f}%  "
        f"+/-1 {fine_gm['pm1']*100:5.1f}%  "
        f"rho {fine_gm['rho']:+.3f}"
    )

    print(
        f"  within WM "
        f"n={fine_wm['n']:6d}  "
        f"exact {fine_wm['exact']*100:5.1f}%  "
        f"+/-1 {fine_wm['pm1']*100:5.1f}%  "
        f"rho {fine_wm['rho']:+.3f}"
    )

    print()
    print(
        layer_table[
            [
                "layer",
                "n",
                "exact",
                "pm1",
                "mean_landed",
                "sd_landed",
                "bias",
            ]
        ].to_string(
            index=False,
            formatters={
                "exact":
                    lambda x: (
                        f"{x*100:.1f}%"
                        if pd.notna(x)
                        else "NA"
                    ),

                "pm1":
                    lambda x: (
                        f"{x*100:.1f}%"
                        if pd.notna(x)
                        else "NA"
                    ),

                "mean_landed":
                    lambda x: (
                        f"{x:.2f}"
                        if pd.notna(x)
                        else "NA"
                    ),

                "sd_landed":
                    lambda x: (
                        f"{x:.2f}"
                        if pd.notna(x)
                        else "NA"
                    ),

                "bias":
                    lambda x: (
                        f"{x:+.2f}"
                        if pd.notna(x)
                        else "NA"
                    ),
            },
        )
    )

    print()
    print(
        "--- self-check ---"
    )

    for key, value in checks.items():

        print(
            f"  "
            f"[{'PASS ' if value else 'CHECK'}] "
            f"{key}"
        )

    print(
        "ALL PASS"
        if all_pass
        else ">>> some CHECK <<<"
    )

    # ========================================================
    # FIGURES
    # ========================================================

    png = plot_final_alignment(
        br=br,
        xen_xy=xen_xy,
        xen_ord=xen_ord,
        vis_xy=vis_xy,
        vis_ord=vis_ord,
        match_prob=match_prob,
        matched=matched,
        score_mask=score_mask,
        result=result,
        outdir=outdir,
    )

    plot_layer_landing(
        br=br,
        xen_xy=xen_xy,
        xen_ord=xen_ord,
        vis_xy=vis_xy,
        vis_ord=vis_ord,
        score_mask=score_mask,
        nearest_idx=nearest_idx,
        layer_table=layer_table,
        outdir=outdir,
    )

    result[
        "figure"
    ] = str(
        png
    )

    return (
        result,
        layer_table,
    )


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
    )

    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
    )

    parser.add_argument(
        "--br",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--spot-radius-um",
        type=float,
        default=DEFAULT_SPOT_RADIUS_UM,
        help=(
            "Physical Visium spot radius. "
            "Advisor v4 reference = 32.54 um."
        ),
    )

    args = parser.parse_args()

    root = (
        args.root.resolve()
    )

    outdir = (
        args.out.resolve()
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    samples = discover_samples(
        root,
        args.br,
    )

    print("=" * 90)
    print(
        "STALIGN FINAL ALIGNMENT QC "
        "— ADVISOR STYLE"
    )
    print("=" * 90)

    print(
        f"Root             : {root}"
    )

    print(
        f"Output           : {outdir}"
    )

    print(
        f"Match threshold  : {MATCH_THRESH}"
    )

    print(
        f"Spot radius      : "
        f"{args.spot_radius_um:.2f} um"
    )

    print(
        f"Samples found    : "
        f"{len(samples)}"
    )

    if not samples:

        raise RuntimeError(
            "No valid sample pairs found."
        )

    results = []
    all_layer_tables = []

    for (
        br,
        aligned_file,
        visium_file,
    ) in samples:

        try:

            (
                result,
                layer_table,
            ) = process_sample(
                br=br,
                aligned_file=aligned_file,
                visium_file=visium_file,
                outdir=outdir,
                spot_radius_um=args.spot_radius_um,
            )

            all_layer_tables.append(
                layer_table
            )

        except Exception as exc:

            print()
            print(
                f"ERROR {br}: {exc}"
            )

            traceback.print_exc()

            result = {
                "BrNum": br,
                "status": "FAIL",
                "error": str(exc),
            }

        results.append(
            result
        )

    # ========================================================
    # COHORT SUMMARY
    # ========================================================

    summary = pd.DataFrame(
        results
    )

    summary_file = (
        outdir
        / "final_aligned_qc_summary.csv"
    )

    summary.to_csv(
        summary_file,
        index=False,
    )

    if all_layer_tables:

        cohort_layers = pd.concat(
            all_layer_tables,
            ignore_index=True,
        )

        cohort_layers.to_csv(
            (
                outdir
                / "final_aligned_qc_all_layer_metrics.csv"
            ),
            index=False,
        )

    print()
    print("=" * 90)
    print(
        "COHORT SUMMARY"
    )
    print("=" * 90)

    show_cols = [
        c
        for c in [
            "BrNum",
            "status",
            "matched_pct",
            "within_spot_radius_pct",
            "scored_pct",
            "exact",
            "pm1",
            "rho",
            "gmwm_accuracy",
            "fine_gm_rho",
            "usable_pct",
            "self_check_all_pass",
        ]
        if c in summary.columns
    ]

    print(
        summary[
            show_cols
        ].to_string(
            index=False
        )
    )

    print()
    print(
        f"Summary: {summary_file}"
    )

    print(
        f"Figures: {outdir}"
    )


if __name__ == "__main__":
    main()
