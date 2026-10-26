#!/usr/bin/env python3

"""
02_make_final_aligned_qc.py

Post-hoc cohort QC of final Xenium -> Visium STalign coordinates.

IMPORTANT
---------
This script DOES NOT rerun STalign/LDDMM.

It evaluates the already-exported final coordinates:
    x_aligned
    y_aligned
    match_prob

The alignment figure follows the advisor notebook plotting style:
    Panel 1: Xenium match probability
             RdYlGn, vmin=0, vmax=1, s=0.6

    Panel 2: Visium domains + matched Xenium layers
             Visium: viridis, s=9, alpha=0.45
             Xenium: autumn, s=0.6, alpha=0.35,
                     vmin=1, vmax=7

    Panel 3: Visium spatial-domain reference
             viridis, s=9

Cell-to-spot scoring follows the advisor Cell 8 logic:
    matched = match_prob > 0.85
    near    = distance to nearest Visium spot <= spot radius
    scored  = matched & near & annotated Xenium layer
              & valid nearest Visium domain

Metrics:
    exact layer agreement
    +/-1 layer agreement
    Spearman rho
    +/-1 permutation baseline
    GM/WM boundary accuracy
    GM -> GM
    WM -> WM
    within-GM exact / +/-1 / rho
    within-WM exact / +/-1 / rho
    per-layer recall
    per-layer +/-1
    mean landed domain
    systematic bias
    confusion matrices

Default spot radius = 32.57 um, based on the exported STalign
cohort metadata. Override using --spot-radius-um if needed.

Layer ordering follows the advisor notebook:

    0 L1/M   <-> SpD07-L1
    1 L2/3   <-> SpD06-L2/3
    2 L3/4   <-> SpD02-L3/4
    3 L5     <-> SpD05-L5
    4 L6     <-> SpD03-L6
    5 WM     <-> SpD04-WM
    6 WMtz   <-> SpD01-WMtz

Plotting uses 1..7 only for matching the advisor's colour scale.
Metrics use 0..6 exactly like advisor Cell 8.
"""

from pathlib import Path
import argparse
import traceback
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

from scipy.spatial import cKDTree
from scipy.stats import spearmanr


# ============================================================
# PATHS
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
    / "stalign_final_aligned_qc"
)


# ============================================================
# ADVISOR PARAMETERS
# ============================================================

MATCH_THRESH = 0.85
SPOT_PITCH_UM = 100.0
DEFAULT_SPOT_RADIUS_UM = 32.57
N_PERM = 200

LAY = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WM",
    "WMtz",
]

SPD_ORDER = [
    "SpD07-L1",
    "SpD06-L2/3",
    "SpD02-L3/4",
    "SpD05-L5",
    "SpD03-L6",
    "SpD04-WM",
    "SpD01-WMtz",
]


# ============================================================
# LABEL MAPPING
# ============================================================

def clean_label(x):

    if pd.isna(x):
        return None

    return re.sub(
        r"[^a-z0-9]",
        "",
        str(x).lower()
    )


# 1..7 for plotting
LABEL_TO_ORD1 = {

    # L1/M
    "l1": 1,
    "l1m": 1,
    "spd07": 1,
    "spd07l1": 1,

    # L2/3
    "l23": 2,
    "spd06": 2,
    "spd06l23": 2,

    # L3/4
    "l34": 3,
    "spd02": 3,
    "spd02l34": 3,

    # L5
    "l5": 4,
    "spd05": 4,
    "spd05l5": 4,

    # L6
    "l6": 5,
    "spd03": 5,
    "spd03l6": 5,

    # WM
    "wm": 6,
    "spd04": 6,
    "spd04wm": 6,

    # WM transition zone
    "wmtz": 7,
    "spd01": 7,
    "spd01wmtz": 7,
}


def labels_to_ord1(series):

    out = np.full(
        len(series),
        np.nan,
        dtype=float
    )

    for j, val in enumerate(series):

        if pd.isna(val):
            continue

        # Numeric 1..7
        try:
            f = float(val)

            if (
                np.isfinite(f)
                and f >= 1
                and f <= 7
                and float(f).is_integer()
            ):
                out[j] = int(f)
                continue

        except Exception:
            pass

        key = clean_label(val)

        if key in LABEL_TO_ORD1:
            out[j] = LABEL_TO_ORD1[key]

    return out


# ============================================================
# UTILITIES
# ============================================================

def parse_bool(series):

    if series.dtype == bool:
        return series.to_numpy()

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
        .to_numpy()
    )


def safe_mean(x):

    x = np.asarray(x)

    if len(x) == 0:
        return np.nan

    return float(
        np.mean(x)
    )


def safe_rho(a, b):

    a = np.asarray(a)
    b = np.asarray(b)

    if len(a) < 3:
        return np.nan

    if (
        len(np.unique(a)) < 2
        or len(np.unique(b)) < 2
    ):
        return np.nan

    return float(
        spearmanr(
            a,
            b
        ).statistic
    )


def get_spot_radius(
    xen,
    vis,
    fallback
):

    candidates = [
        "spot_radius_um",
        "spot_r_um",
        "SPOT_R_UM",
    ]

    for df in (vis, xen):

        for col in candidates:

            if col not in df.columns:
                continue

            v = pd.to_numeric(
                df[col],
                errors="coerce"
            )

            v = v[
                np.isfinite(v)
                & (v > 0)
            ]

            if len(v):

                return (
                    float(np.median(v)),
                    f"column:{col}"
                )

    return (
        float(fallback),
        "fallback"
    )



def soften_cmap(name="RdYlGn", blend=0.45):
    """
    Blend a colormap with white to create a softer pastel version.
    blend=0 -> original colormap
    blend=1 -> fully white
    """
    base = plt.get_cmap(name, 256)
    cols = base(np.linspace(0, 1, 256))

    white = np.array([1.0, 1.0, 1.0, 1.0])
    cols = (1 - blend) * cols + blend * white
    cols[:, 3] = 1.0

    return mcolors.LinearSegmentedColormap.from_list(
        f"soft_{name}",
        cols
    )


# ============================================================
# ALIGNMENT FIGURE
# ============================================================

def make_alignment_figure(
    br,
    xen_xy,
    match_prob,
    matched,
    xen_ord1,
    vis_xy,
    vis_ord1,
    outdir
):

    fig, ax = plt.subplots(
        1,
        3,
        figsize=(18, 5)
    )

    # --------------------------------------------------------
    # PANEL 1
    # Exact advisor style:
    # Xenium match probability
    # --------------------------------------------------------

    good_x = (
        np.isfinite(xen_xy).all(axis=1)
        & np.isfinite(match_prob)
    )

    soft_prob_cmap = soften_cmap(
        "RdYlGn",
        blend=0.45
    )

    s = ax[0].scatter(
        xen_xy[good_x, 0] / 1000,
        xen_xy[good_x, 1] / 1000,
        s=0.6,
        c=match_prob[good_x],
        cmap=soft_prob_cmap,
        vmin=0,
        vmax=1,
        linewidths=0,
        rasterized=True
    )

    ax[0].set_title(
        f"match probability "
        f"({matched.mean()*100:.0f}% > "
        f"{MATCH_THRESH})"
    )

    plt.colorbar(
        s,
        ax=ax[0],
        fraction=.046
    )

    # --------------------------------------------------------
    # PANEL 2
    # Visium underlay + matched Xenium
    # --------------------------------------------------------

    good_v = (
        np.isfinite(vis_xy).all(axis=1)
        & np.isfinite(vis_ord1)
    )

    ax[1].scatter(
        vis_xy[good_v, 0] / 1000,
        vis_xy[good_v, 1] / 1000,
        s=9,
        c=vis_ord1[good_v],
        cmap="viridis",
        alpha=.45,
        linewidths=0,
        rasterized=True
    )

    good_overlay = (
        matched
        & np.isfinite(xen_xy).all(axis=1)
        & np.isfinite(xen_ord1)
    )

    ax[1].scatter(
        xen_xy[good_overlay, 0] / 1000,
        xen_xy[good_overlay, 1] / 1000,
        s=.6,
        c=xen_ord1[good_overlay],
        cmap="autumn",
        vmin=1,
        vmax=7,
        alpha=.35,
        linewidths=0,
        rasterized=True
    )

    ax[1].set_title(
        "aligned Xenium (matched) on Visium spd"
    )

    # --------------------------------------------------------
    # PANEL 3
    # Visium reference
    # --------------------------------------------------------

    ax[2].scatter(
        vis_xy[good_v, 0] / 1000,
        vis_xy[good_v, 1] / 1000,
        s=9,
        c=vis_ord1[good_v],
        cmap="viridis",
        linewidths=0,
        rasterized=True
    )

    ax[2].set_title(
        "Visium spd reference"
    )

    for a in ax:

        a.set_aspect("equal")
        a.invert_yaxis()

        a.set_xlabel("mm")
        a.set_ylabel("mm")

    fig.suptitle(
        f"{br} — Xenium-to-Visium alignment after LDDMM",
        fontsize=13,
        fontweight="bold",
        y=1.01
    )

    plt.tight_layout()

    png = (
        outdir
        / f"{br}_final_alignment_qc.png"
    )

    pdf = (
        outdir
        / f"{br}_final_alignment_qc.pdf"
    )

    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf,
        bbox_inches="tight"
    )

    plt.close(fig)

    return png, pdf


# ============================================================
# CELL-TO-SPOT QC
# ============================================================

def calculate_cell_spot_qc(
    br,
    xen,
    vis,
    xen_xy,
    vis_xy,
    match_prob,
    xen_ord1,
    vis_ord1,
    spot_radius_um
):

    n = len(xen)

    # --------------------------------------------------------
    # Nearest Visium spot exactly from final coordinates
    # --------------------------------------------------------

    vis_coord_ok = (
        np.isfinite(vis_xy).all(axis=1)
    )

    if vis_coord_ok.sum() == 0:
        raise ValueError(
            f"{br}: no valid Visium coordinates"
        )

    vis_xy_tree = (
        vis_xy[
            vis_coord_ok
        ]
    )

    vis_ord1_tree = (
        vis_ord1[
            vis_coord_ok
        ]
    )

    tree = cKDTree(
        vis_xy_tree
    )

    coord_ok = (
        np.isfinite(xen_xy).all(axis=1)
    )

    nearest_dist = np.full(
        n,
        np.inf,
        dtype=float
    )

    nearest_idx = np.full(
        n,
        -1,
        dtype=int
    )

    d, j = tree.query(
        xen_xy[
            coord_ok
        ],
        k=1
    )

    nearest_dist[
        coord_ok
    ] = d

    nearest_idx[
        coord_ok
    ] = j

    nearest_vis_ord1 = np.full(
        n,
        np.nan,
        dtype=float
    )

    nearest_vis_ord1[
        coord_ok
    ] = (
        vis_ord1_tree[j]
    )

    # --------------------------------------------------------
    # Advisor Cell 8 definitions
    # --------------------------------------------------------

    matched = (
        np.isfinite(match_prob)
        & (
            match_prob
            > MATCH_THRESH
        )
    )

    has_lab = (
        np.isfinite(
            xen_ord1
        )
    )

    near = (
        nearest_dist
        <= spot_radius_um
    )

    valid_target = (
        np.isfinite(
            nearest_vis_ord1
        )
    )

    scored = (
        matched
        & near
        & has_lab
        & valid_target
        & coord_ok
    )

    if scored.sum() == 0:

        raise ValueError(
            f"{br}: no cells survived "
            "matched + in-spot + layer filter"
        )

    # Advisor Cell 8 uses zero-based ordinals
    xo = (
        xen_ord1[
            scored
        ].astype(int)
        - 1
    )

    vo = (
        nearest_vis_ord1[
            scored
        ].astype(int)
        - 1
    )

    dd = (
        nearest_dist[
            scored
        ]
    )

    # --------------------------------------------------------
    # Overall
    # --------------------------------------------------------

    exact = float(
        (xo == vo).mean()
    )

    adj = float(
        (
            np.abs(
                xo - vo
            ) <= 1
        ).mean()
    )

    rho = safe_rho(
        xo,
        vo
    )

    rng = np.random.default_rng(
        0
    )

    base = np.array(
        [
            (
                np.abs(
                    rng.permutation(xo)
                    - vo
                ) <= 1
            ).mean()
            for _ in range(
                N_PERM
            )
        ],
        dtype=float
    )

    baseline = float(
        base.mean()
    )

    baseline_sd = float(
        base.std()
    )

    if baseline_sd > 0:

        lift_sigma = float(
            (
                adj
                - baseline
            )
            / baseline_sd
        )

    else:

        lift_sigma = np.nan

    # --------------------------------------------------------
    # Coarse GM/WM
    #
    # EXACT advisor Cell 8 zero-based definition:
    #
    # 0 = L1/M excluded
    # 1..4 = GM
    # 5..6 = WM / WMtz
    # --------------------------------------------------------

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
        GM_I
    )

    is_wm_x = np.isin(
        xo,
        WM_I
    )

    is_gm_v = np.isin(
        vo,
        GM_I
    )

    coarse_sel = (
        is_gm_x
        | is_wm_x
    )

    gmwm_acc = (
        safe_mean(
            (
                is_gm_x
                == is_gm_v
            )[
                coarse_sel
            ]
        )
    )

    gm_to_gm = (
        safe_mean(
            is_gm_v[
                is_gm_x
            ]
        )
    )

    wm_to_wm = (
        safe_mean(
            (
                ~is_gm_v[
                    is_wm_x
                ]
            )
        )
    )

    # --------------------------------------------------------
    # Fine: within correct compartment
    # --------------------------------------------------------

    within = {}

    for (
        name,
        mx,
        right
    ) in [

        (
            "gm",
            is_gm_x,
            is_gm_v
        ),

        (
            "wm",
            is_wm_x,
            ~is_gm_v
        ),
    ]:

        mm = (
            mx
            & right
        )

        nn = int(
            mm.sum()
        )

        if nn > 0:

            ex = float(
                (
                    xo[mm]
                    == vo[mm]
                ).mean()
            )

            ad = float(
                (
                    np.abs(
                        xo[mm]
                        - vo[mm]
                    ) <= 1
                ).mean()
            )

            rr = safe_rho(
                xo[mm],
                vo[mm]
            )

        else:

            ex = np.nan
            ad = np.nan
            rr = np.nan

        within[name] = {
            "n": nn,
            "exact": ex,
            "adj": ad,
            "rho": rr,
        }

    # --------------------------------------------------------
    # Per-layer metrics
    # --------------------------------------------------------

    layer_rows = []

    mu = np.full(
        7,
        np.nan
    )

    sd = np.full(
        7,
        np.nan
    )

    for i, layer in enumerate(
        LAY
    ):

        m = (
            xo == i
        )

        n_i = int(
            m.sum()
        )

        if n_i > 0:

            exact_i = float(
                (
                    vo[m]
                    == i
                ).mean()
            )

            adj_i = float(
                (
                    np.abs(
                        vo[m]
                        - i
                    ) <= 1
                ).mean()
            )

            mean_i = float(
                vo[m].mean()
            )

            sd_i = float(
                vo[m].std()
            )

            bias_i = (
                mean_i
                - i
            )

            mu[i] = mean_i
            sd[i] = sd_i

        else:

            exact_i = np.nan
            adj_i = np.nan
            mean_i = np.nan
            sd_i = np.nan
            bias_i = np.nan

        layer_rows.append(
            {
                "BrNum": br,
                "layer": layer,
                "layer_ordinal_0based": i,
                "n": n_i,
                "exact": exact_i,
                "plusminus1": adj_i,
                "mean_landed_visium_spd": mean_i,
                "sd_landed_visium_spd": sd_i,
                "bias": bias_i,
            }
        )

    gm_bias = (
        mu[1:5]
        - np.arange(
            1,
            5
        )
    )

    ok_bias = (
        np.isfinite(
            gm_bias
        )
    )

    mean_gm_bias = (
        float(
            np.nanmean(
                gm_bias
            )
        )
        if ok_bias.any()
        else np.nan
    )

    if (
        ok_bias.any()
        and np.isfinite(
            mean_gm_bias
        )
        and mean_gm_bias != 0
    ):

        gm_bias_same_sign = bool(
            np.all(
                np.sign(
                    gm_bias[
                        ok_bias
                    ]
                )
                ==
                np.sign(
                    mean_gm_bias
                )
            )
        )

    else:

        gm_bias_same_sign = False

    # --------------------------------------------------------
    # Confusion matrices
    # --------------------------------------------------------

    C = np.zeros(
        (
            7,
            7
        ),
        dtype=float
    )

    np.add.at(
        C,
        (
            xo,
            vo
        ),
        1
    )

    row_den = (
        C.sum(
            axis=1,
            keepdims=True
        )
    )

    col_den = (
        C.sum(
            axis=0,
            keepdims=True
        )
    )

    Crow = np.divide(
        C,
        row_den,
        out=np.zeros_like(C),
        where=row_den > 0
    )

    Ccol = np.divide(
        C,
        col_den,
        out=np.zeros_like(C),
        where=col_den > 0
    )

    # --------------------------------------------------------
    # Extra exported-QC fields
    # --------------------------------------------------------

    inside_pct = np.nan
    usable_pct = np.nan

    if "inside_visium" in xen.columns:

        inside_bool = parse_bool(
            xen[
                "inside_visium"
            ]
        )

        inside_pct = (
            100
            * inside_bool.mean()
        )

    if "usable" in xen.columns:

        usable_bool = parse_bool(
            xen[
                "usable"
            ]
        )

        usable_pct = (
            100
            * usable_bool.mean()
        )

    assigned_pct = np.nan

    if "visium_spot" in xen.columns:

        v = (
            xen[
                "visium_spot"
            ]
            .astype(str)
            .str.strip()
        )

        assigned = (
            xen[
                "visium_spot"
            ].notna()
            & (v != "")
            & (v != "nan")
        )

        assigned_pct = (
            100
            * assigned.mean()
        )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary = {

        "BrNum": br,
        "status": "PASS",

        "n_cells": n,

        "n_matched": int(
            matched.sum()
        ),

        "matched_pct": (
            100
            * matched.mean()
        ),

        "spot_radius_um": spot_radius_um,

        "n_within_spot_radius": int(
            near.sum()
        ),

        "within_spot_radius_pct": (
            100
            * near.mean()
        ),

        "n_has_layer": int(
            has_lab.sum()
        ),

        "has_layer_pct": (
            100
            * has_lab.mean()
        ),

        "n_scored": int(
            scored.sum()
        ),

        "scored_pct": (
            100
            * scored.mean()
        ),

        "median_dist_matched_um": (
            float(
                np.median(
                    nearest_dist[
                        matched
                        & coord_ok
                    ]
                )
            )
            if (
                matched
                & coord_ok
            ).any()
            else np.nan
        ),

        "median_dist_scored_um": (
            float(
                np.median(
                    dd
                )
            )
        ),

        "inside_visium_pct": inside_pct,
        "usable_pct": usable_pct,
        "assigned_to_spot_pct": assigned_pct,

        "exact_layer_agreement": exact,
        "plusminus1_layer_agreement": adj,
        "layer_spearman": rho,

        "plusminus1_permutation_baseline": baseline,
        "permutation_baseline_sd": baseline_sd,
        "plusminus1_lift_sigma": lift_sigma,

        "gmwm_accuracy": gmwm_acc,
        "gm_to_gm": gm_to_gm,
        "wm_to_wm": wm_to_wm,

        "within_gm_n": within["gm"]["n"],
        "within_gm_exact": within["gm"]["exact"],
        "within_gm_plusminus1": within["gm"]["adj"],
        "within_gm_spearman": within["gm"]["rho"],

        "within_wm_n": within["wm"]["n"],
        "within_wm_exact": within["wm"]["exact"],
        "within_wm_plusminus1": within["wm"]["adj"],
        "within_wm_spearman": within["wm"]["rho"],

        "mean_gm_bias": mean_gm_bias,
        "gm_bias_same_sign": gm_bias_same_sign,
    }

    return (
        summary,
        layer_rows,
        C,
        Crow,
        Ccol,
        mu,
        sd,
        matched,
        nearest_dist
    )


# ============================================================
# ADVISOR-STYLE LAYER CONCORDANCE FIGURE
# ============================================================

def make_layer_qc_figure(
    br,
    summary,
    Crow,
    Ccol,
    mu,
    sd,
    outdir
):

    fig, ax = plt.subplots(
        1,
        3,
        figsize=(19, 5.5)
    )

    for a, M, ttl in [

        (
            ax[0],
            Crow,
            (
                "row-normalized (recall)\n"
                f"exact "
                f"{summary['exact_layer_agreement']*100:.0f}%, "
                f"+/-1 "
                f"{summary['plusminus1_layer_agreement']*100:.0f}%, "
                f"{summary['plusminus1_lift_sigma']:.0f} sigma"
            )
        ),

        (
            ax[1],
            Ccol,
            "column-normalized (precision)"
        ),
    ]:

        im = a.imshow(
            M,
            cmap="magma",
            vmin=0,
            vmax=1
        )

        a.set_xticks(
            range(7)
        )

        a.set_xticklabels(
            LAY,
            rotation=45,
            fontsize=7
        )

        a.set_yticks(
            range(7)
        )

        a.set_yticklabels(
            LAY,
            fontsize=7
        )

        a.set_xlabel(
            "landed Visium spd"
        )

        a.set_ylabel(
            "Xenium Lieber layer"
        )

        a.set_title(
            ttl,
            fontsize=9
        )

        a.plot(
            [-.5, 6.5],
            [-.5, 6.5],
            c="#00E5FF",
            lw=.7,
            alpha=.6
        )

        a.axvline(
            4.5,
            c="#00E5FF",
            lw=.7,
            ls="--",
            alpha=.5
        )

        a.axhline(
            4.5,
            c="#00E5FF",
            lw=.7,
            ls="--",
            alpha=.5
        )

        plt.colorbar(
            im,
            ax=a,
            fraction=.046,
            pad=.02
        )

    # systematic offset panel
    ax[2].plot(
        [0, 6],
        [0, 6],
        c="#bbb",
        lw=1,
        ls="--",
        label="perfect"
    )

    ax[2].errorbar(
        range(7),
        mu,
        yerr=sd,
        fmt="-o",
        ms=5,
        lw=1.2,
        c="#0072B2",
        capsize=3,
        label="observed"
    )

    ax[2].axvspan(
        .5,
        4.5,
        color="#009E73",
        alpha=.07
    )

    ax[2].text(
        2.5,
        5.6,
        "GM",
        ha="center",
        fontsize=9,
        color="#009E73"
    )

    ax[2].set_xticks(
        range(7)
    )

    ax[2].set_xticklabels(
        LAY,
        rotation=45,
        fontsize=7
    )

    ax[2].set_xlabel(
        "Xenium Lieber layer"
    )

    ax[2].set_ylabel(
        "mean landed Visium spd"
    )

    ax[2].set_title(
        "systematic offset vs noise",
        fontsize=9
    )

    ax[2].legend(
        fontsize=8
    )

    fig.suptitle(
        f"{br} — Xenium Lieber layer vs landed Visium spd",
        fontsize=13,
        fontweight="bold"
    )

    plt.tight_layout()

    png = (
        outdir
        / f"{br}_layer_concordance_qc.png"
    )

    pdf = (
        outdir
        / f"{br}_layer_concordance_qc.pdf"
    )

    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        pdf,
        bbox_inches="tight"
    )

    plt.close(fig)

    return png, pdf


# ============================================================
# ONE SAMPLE
# ============================================================

def process_sample(
    br,
    aligned_file,
    visium_file,
    outdir,
    fallback_spot_radius
):

    print()
    print("=" * 88)
    print(br)
    print("=" * 88)

    print(
        "Aligned:",
        aligned_file
    )

    print(
        "Visium :",
        visium_file
    )

    xen = pd.read_csv(
        aligned_file
    )

    vis = pd.read_csv(
        visium_file
    )

    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

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

    for col in required_xen:

        if col not in xen.columns:

            raise ValueError(
                f"{br}: missing Xenium "
                f"column '{col}'"
            )

    for col in required_vis:

        if col not in vis.columns:

            raise ValueError(
                f"{br}: missing Visium "
                f"column '{col}'"
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
    ].astype(float)

    vis_xy = np.c_[

        pd.to_numeric(
            vis["x"],
            errors="coerce"
        ),

        pd.to_numeric(
            vis["y"],
            errors="coerce"
        )
    ].astype(float)

    match_prob = pd.to_numeric(
        xen["match_prob"],
        errors="coerce"
    ).to_numpy(
        dtype=float
    )

    xen_ord1 = labels_to_ord1(
        xen["layer"]
    )

    vis_ord1 = labels_to_ord1(
        vis["domain"]
    )

    print(
        "Xenium mapped layers : "
        f"{np.isfinite(xen_ord1).sum():,}/"
        f"{len(xen_ord1):,}"
    )

    print(
        "Visium mapped domains: "
        f"{np.isfinite(vis_ord1).sum():,}/"
        f"{len(vis_ord1):,}"
    )

    spot_radius_um, radius_source = (
        get_spot_radius(
            xen,
            vis,
            fallback_spot_radius
        )
    )

    print(
        f"Spot radius used      : "
        f"{spot_radius_um:.2f} um "
        f"({radius_source})"
    )

    # --------------------------------------------------------
    # Cell-to-spot QC
    # --------------------------------------------------------

    (
        summary,
        layer_rows,
        C,
        Crow,
        Ccol,
        mu,
        sd,
        matched,
        nearest_dist
    ) = calculate_cell_spot_qc(
        br=br,
        xen=xen,
        vis=vis,
        xen_xy=xen_xy,
        vis_xy=vis_xy,
        match_prob=match_prob,
        xen_ord1=xen_ord1,
        vis_ord1=vis_ord1,
        spot_radius_um=spot_radius_um
    )

    # --------------------------------------------------------
    # Alignment figure
    # --------------------------------------------------------

    align_png, align_pdf = (
        make_alignment_figure(
            br=br,
            xen_xy=xen_xy,
            match_prob=match_prob,
            matched=matched,
            xen_ord1=xen_ord1,
            vis_xy=vis_xy,
            vis_ord1=vis_ord1,
            outdir=outdir
        )
    )

    # --------------------------------------------------------
    # Layer QC figure
    # --------------------------------------------------------

    layer_png, layer_pdf = (
        make_layer_qc_figure(
            br=br,
            summary=summary,
            Crow=Crow,
            Ccol=Ccol,
            mu=mu,
            sd=sd,
            outdir=outdir
        )
    )

    summary[
        "alignment_figure_png"
    ] = str(
        align_png
    )

    summary[
        "layer_qc_figure_png"
    ] = str(
        layer_png
    )

    # --------------------------------------------------------
    # Terminal report
    # --------------------------------------------------------

    print()
    print("--- FILTERING ---")

    print(
        f"cells                         : "
        f"{summary['n_cells']:,}"
    )

    print(
        f"matched (prob > {MATCH_THRESH}) "
        f"        : "
        f"{summary['n_matched']:6d} "
        f"({summary['matched_pct']:5.1f}%)"
    )

    print(
        f"within {spot_radius_um:.2f} um "
        f"of a spot  : "
        f"{summary['n_within_spot_radius']:6d} "
        f"({summary['within_spot_radius_pct']:5.1f}%)"
    )

    print(
        f"has Lieber layer               : "
        f"{summary['n_has_layer']:6d} "
        f"({summary['has_layer_pct']:5.1f}%)"
    )

    print(
        f"scored (all above)             : "
        f"{summary['n_scored']:6d} "
        f"({summary['scored_pct']:5.1f}%)"
    )

    print(
        f"median nearest distance, "
        f"matched : "
        f"{summary['median_dist_matched_um']:.2f} um"
    )

    disc_cover = (
        np.pi
        * spot_radius_um**2
        /
        (
            SPOT_PITCH_UM**2
            * np.sqrt(3)
            / 2
        )
        * 100
    )

    print(
        f"spot discs cover approximately "
        f"{disc_cover:.0f}% of tissue area"
    )

    print()
    print(
        f"--- OVERALL "
        f"(n={summary['n_scored']:,}) ---"
    )

    print(
        f"exact "
        f"{summary['exact_layer_agreement']*100:.1f}%"
        f" | +/-1 "
        f"{summary['plusminus1_layer_agreement']*100:.1f}%"
        f" | Spearman "
        f"{summary['layer_spearman']:.3f}"
    )

    print(
        f"+/-1 permutation baseline "
        f"{summary['plusminus1_permutation_baseline']*100:.1f}%"
        f" -> "
        f"{summary['plusminus1_lift_sigma']:.1f} sigma"
    )

    print()
    print(
        "--- COARSE GM/WM "
        "(L1/M excluded) ---"
    )

    print(
        f"correct side : "
        f"{summary['gmwm_accuracy']*100:.1f}%"
    )

    print(
        f"GM -> GM     : "
        f"{summary['gm_to_gm']*100:.1f}%"
    )

    print(
        f"WM -> WM     : "
        f"{summary['wm_to_wm']*100:.1f}%"
    )

    print()
    print(
        "--- FINE: CORRECT COMPARTMENT ---"
    )

    print(
        f"within GM n={summary['within_gm_n']:6d}"
        f" | exact "
        f"{summary['within_gm_exact']*100:.1f}%"
        f" | +/-1 "
        f"{summary['within_gm_plusminus1']*100:.1f}%"
        f" | rho "
        f"{summary['within_gm_spearman']:.3f}"
    )

    print(
        f"within WM n={summary['within_wm_n']:6d}"
        f" | exact "
        f"{summary['within_wm_exact']*100:.1f}%"
        f" | +/-1 "
        f"{summary['within_wm_plusminus1']*100:.1f}%"
        f" | rho "
        f"{summary['within_wm_spearman']:.3f}"
    )

    print()
    print(
        f"mean GM bias : "
        f"{summary['mean_gm_bias']:+.2f}"
        f" | all same sign: "
        f"{summary['gm_bias_same_sign']}"
    )

    return (
        summary,
        layer_rows
    )


# ============================================================
# MAIN
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

    parser.add_argument(
        "--spot-radius-um",
        type=float,
        default=DEFAULT_SPOT_RADIUS_UM,
        help=(
            "Fallback Visium spot radius "
            "if not stored in exported CSV "
            f"(default {DEFAULT_SPOT_RADIUS_UM} um)"
        )
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
                f"WARNING {br}: "
                f"missing "
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

    print("=" * 88)
    print("FINAL STALIGN COHORT QC")
    print("=" * 88)

    print(
        "Root        :",
        root
    )

    print(
        "Output      :",
        outdir
    )

    print(
        "Match cutoff:",
        MATCH_THRESH
    )

    print(
        "Fallback spot radius:",
        f"{args.spot_radius_um:.2f} um"
    )

    print(
        "Samples found:",
        len(samples)
    )

    summaries = []
    per_layer_rows = []

    for (
        br,
        aligned_file,
        visium_file
    ) in samples:

        try:

            summary, layer_rows = (
                process_sample(
                    br=br,
                    aligned_file=aligned_file,
                    visium_file=visium_file,
                    outdir=outdir,
                    fallback_spot_radius=(
                        args.spot_radius_um
                    )
                )
            )

            summaries.append(
                summary
            )

            per_layer_rows.extend(
                layer_rows
            )

        except Exception as e:

            traceback.print_exc()

            summaries.append(
                {
                    "BrNum": br,
                    "status": "FAIL",
                    "error": str(e),
                }
            )

    # --------------------------------------------------------
    # Cohort summary
    # --------------------------------------------------------

    summary_df = pd.DataFrame(
        summaries
    )

    summary_file = (
        outdir
        / "final_aligned_qc_summary.csv"
    )

    summary_df.to_csv(
        summary_file,
        index=False
    )

    per_layer_df = pd.DataFrame(
        per_layer_rows
    )

    per_layer_file = (
        outdir
        / "final_aligned_qc_per_layer.csv"
    )

    per_layer_df.to_csv(
        per_layer_file,
        index=False
    )

    print()
    print("=" * 88)
    print("FINAL ALIGNMENT SUMMARY")
    print("=" * 88)

    show_cols = [
        "BrNum",
        "status",
        "matched_pct",
        "scored_pct",
        "exact_layer_agreement",
        "plusminus1_layer_agreement",
        "layer_spearman",
        "gmwm_accuracy",
        "within_gm_spearman",
        "mean_gm_bias",
    ]

    show_cols = [
        c
        for c in show_cols
        if c in summary_df.columns
    ]

    print(
        summary_df[
            show_cols
        ].to_string(
            index=False
        )
    )

    print()
    print(
        "Summary CSV :",
        summary_file
    )

    print(
        "Per-layer CSV:",
        per_layer_file
    )

    print(
        "Figures      :",
        outdir
    )


if __name__ == "__main__":
    main()
