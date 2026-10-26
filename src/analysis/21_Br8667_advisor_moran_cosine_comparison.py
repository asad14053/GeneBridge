#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr, pearsonr, wilcoxon


# =============================================================================
# CONFIG
# =============================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BR = "Br8667"

BASE = (
    ROOT
    / "results/analysis/three_source_validation"
    / BR
)

METRICS = (
    BASE
    / f"{BR}_three_source_gene_metrics.csv"
)

OUT = (
    BASE
    / "advisor_moran_cosine_validation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# COLUMNS
# =============================================================================

MORAN_A = "moran_A_Visium"
MORAN_B = "moran_B_measured"
MORAN_C = "moran_C_OOF"

COS_AB = "cosine_A_Visium_vs_B_measured"
COS_AC = "cosine_A_Visium_vs_C_OOF"
COS_BC = "cosine_B_measured_vs_C_OOF"

PATTERN_A = "patterned_A_Visium"
EVALUABLE = "evaluable_all3"


# =============================================================================
# FIXED COLORS
# =============================================================================

# Moran distributions
COLOR_VISIUM = "#4c78a8"
COLOR_MEASURED = "#f58518"
COLOR_OOF = "#54a24b"

# Moran scatter classification
COLOR_PATTERNED = "#f28e2b"
COLOR_NONPATTERNED = "#9ecae1"

# Cosine pair colors
COLOR_AB = "#4c78a8"   # Visium <-> Measured
COLOR_AC = "#f58518"   # Visium <-> OOF
COLOR_BC = "#54a24b"   # Measured <-> OOF


# =============================================================================
# HELPERS
# =============================================================================

def as_bool(s):

    if pd.api.types.is_bool_dtype(s):

        return s.to_numpy(
            bool
        )

    return (
        s.astype(str)
        .str.lower()
        .isin([
            "true",
            "1",
            "yes",
        ])
        .to_numpy(
            bool
        )
    )


def padded_limits(
    values,
    pad_fraction=0.06
):

    values = np.asarray(
        values,
        dtype=float
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    if len(values) == 0:

        return (
            -0.05,
            1.0
        )


    lo = float(
        np.min(
            values
        )
    )

    hi = float(
        np.max(
            values
        )
    )


    span = hi - lo

    if span <= 0:

        span = 1.0


    pad = (
        span
        *
        pad_fraction
    )


    return (
        lo - pad,
        hi + pad
    )


def get_stats(
    x,
    y
):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )


    valid = (
        np.isfinite(
            x
        )
        &
        np.isfinite(
            y
        )
    )


    x = x[
        valid
    ]

    y = y[
        valid
    ]


    rho = float(
        spearmanr(
            x,
            y
        ).statistic
    )


    pear = float(
        pearsonr(
            x,
            y
        ).statistic
    )


    slope, intercept = np.polyfit(
        x,
        y,
        1
    )


    return {
        "x":
            x,

        "y":
            y,

        "n":
            len(
                x
            ),

        "spearman":
            rho,

        "pearson":
            pear,

        "slope":
            float(
                slope
            ),

        "intercept":
            float(
                intercept
            ),

        "median_delta":
            float(
                np.median(
                    y - x
                )
            ),
    }


def plot_moran_scatter(
    ax,
    x,
    y,
    pattern,
    xlim,
    ylim,
    xlabel,
    ylabel,
    title
):

    stats = get_stats(
        x,
        y
    )


    valid = (
        np.isfinite(
            x
        )
        &
        np.isfinite(
            y
        )
    )


    x = np.asarray(
        x,
        dtype=float
    )[
        valid
    ]

    y = np.asarray(
        y,
        dtype=float
    )[
        valid
    ]

    p = np.asarray(
        pattern,
        dtype=bool
    )[
        valid
    ]


    # -------------------------------------------------------------------------
    # NON-PATTERNED
    # -------------------------------------------------------------------------

    ax.scatter(
        x[
            ~p
        ],
        y[
            ~p
        ],
        s=18,
        color=COLOR_NONPATTERNED,
        alpha=0.65,
        linewidths=0,
        label="not patterned",
        rasterized=True
    )


    # -------------------------------------------------------------------------
    # PATTERNED
    # -------------------------------------------------------------------------

    ax.scatter(
        x[
            p
        ],
        y[
            p
        ],
        s=20,
        color=COLOR_PATTERNED,
        alpha=0.78,
        linewidths=0,
        label="patterned",
        rasterized=True
    )


    # -------------------------------------------------------------------------
    # Y = X
    # -------------------------------------------------------------------------

    ref_lo = min(
        xlim[0],
        ylim[0]
    )

    ref_hi = max(
        xlim[1],
        ylim[1]
    )


    ax.plot(
        [
            ref_lo,
            ref_hi
        ],
        [
            ref_lo,
            ref_hi
        ],
        linestyle="--",
        linewidth=1.0,
        color="#377eb8"
    )


    # -------------------------------------------------------------------------
    # LINEAR FIT
    # -------------------------------------------------------------------------

    xx = np.linspace(
        xlim[0],
        xlim[1],
        200
    )


    yy = (
        stats[
            "intercept"
        ]
        +
        stats[
            "slope"
        ]
        *
        xx
    )


    ax.plot(
        xx,
        yy,
        color="#ff7f0e",
        linewidth=1.3
    )


    ax.set_xlim(
        xlim
    )

    ax.set_ylim(
        ylim
    )


    ax.set_xlabel(
        xlabel,
        fontsize=9
    )

    ax.set_ylabel(
        ylabel,
        fontsize=9
    )


    ax.set_title(
        title,
        loc="left",
        fontsize=11,
        fontweight="bold"
    )


    ax.tick_params(
        labelsize=8
    )


    ax.legend(
        frameon=False,
        fontsize=7,
        loc="lower right"
    )


    ax.text(
        0.03,
        0.97,
        (
            f"genes n={stats['n']}\n"
            f"Spearman ρ={stats['spearman']:.3f}\n"
            f"Pearson r={stats['pearson']:.3f}\n"
            f"slope={stats['slope']:.3f}\n"
            f"median ΔI={stats['median_delta']:+.3f}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
        bbox=dict(
            facecolor="white",
            edgecolor="0.8",
            alpha=0.88,
            boxstyle="round,pad=0.2"
        )
    )


    return stats


def horizontal_distribution(
    ax,
    arrays,
    labels,
    colors,
    bins,
    title,
    ylabel
):

    # -------------------------------------------------------------------------
    # SAME BINS FOR ALL DISTRIBUTIONS
    # -------------------------------------------------------------------------

    clean = []

    for values in arrays:

        values = np.asarray(
            values,
            dtype=float
        )

        values = values[
            np.isfinite(
                values
            )
        ]

        clean.append(
            values
        )


    histograms = []

    edges = bins


    for values in clean:

        h, _ = np.histogram(
            values,
            bins=edges
        )

        histograms.append(
            h
        )


    centers = (
        edges[:-1]
        +
        edges[1:]
    ) / 2


    bin_width = (
        edges[1]
        -
        edges[0]
    )


    # Split each bin vertically among 3 series.
    bar_height = (
        bin_width
        *
        0.26
    )


    offsets = [
        +bar_height,
        0,
        -bar_height,
    ]


    for h, label, color, offset in zip(
        histograms,
        labels,
        colors,
        offsets
    ):

        ax.barh(
            centers
            +
            offset,
            h,
            height=bar_height,
            alpha=0.65,
            color=color,
            label=label
        )


    medians = [
        float(
            np.median(
                x
            )
        )
        for x
        in clean
    ]


    for med, color in zip(
        medians,
        colors
    ):

        ax.axhline(
            med,
            color=color,
            linestyle="--",
            linewidth=1.2
        )


    ax.set_xlabel(
        "gene count",
        fontsize=9
    )

    ax.set_ylabel(
        ylabel,
        fontsize=9
    )


    ax.set_title(
        title,
        loc="left",
        fontsize=11,
        fontweight="bold"
    )


    ax.tick_params(
        labelsize=8
    )


    ax.legend(
        frameon=False,
        fontsize=7
    )


    text = []

    for label, med, values in zip(
        labels,
        medians,
        clean
    ):

        text.append(
            f"{label}: n={len(values)}, median={med:.3f}"
        )


    ax.text(
        0.03,
        0.97,
        "\n".join(
            text
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.75
        )
    )


    return {
        label: {
            "n":
                len(
                    values
                ),

            "median":
                med,
        }

        for label, med, values
        in zip(
            labels,
            medians,
            clean
        )
    }


def paired_wilcoxon(
    x,
    y
):

    x = np.asarray(
        x,
        dtype=float
    )

    y = np.asarray(
        y,
        dtype=float
    )


    valid = (
        np.isfinite(
            x
        )
        &
        np.isfinite(
            y
        )
    )


    x = x[
        valid
    ]

    y = y[
        valid
    ]


    try:

        p = float(
            wilcoxon(
                x,
                y
            ).pvalue
        )

    except ValueError:

        p = np.nan


    return {
        "n":
            len(
                x
            ),

        "median_x":
            float(
                np.median(
                    x
                )
            ),

        "median_y":
            float(
                np.median(
                    y
                )
            ),

        "median_delta_y_minus_x":
            float(
                np.median(
                    y - x
                )
            ),

        "wilcoxon_p":
            p,
    }


# =============================================================================
# LOAD
# =============================================================================

if not METRICS.exists():

    raise FileNotFoundError(
        METRICS
    )


df = pd.read_csv(
    METRICS
)


required = [
    "gene",
    MORAN_A,
    MORAN_B,
    MORAN_C,
    COS_AB,
    COS_AC,
    COS_BC,
    PATTERN_A,
    EVALUABLE,
]


missing = [
    c
    for c
    in required
    if c not in df.columns
]


if missing:

    raise RuntimeError(
        f"Missing required columns: {missing}"
    )


for c in [
    MORAN_A,
    MORAN_B,
    MORAN_C,
    COS_AB,
    COS_AC,
    COS_BC,
]:

    df[
        c
    ] = pd.to_numeric(
        df[
            c
        ],
        errors="coerce"
    )


evaluable = as_bool(
    df[
        EVALUABLE
    ]
)


patterned = as_bool(
    df[
        PATTERN_A
    ]
)


# =============================================================================
# MORAN DATA
# =============================================================================

moran_valid = (
    evaluable
    &
    np.isfinite(
        df[
            MORAN_A
        ].to_numpy(float)
    )
    &
    np.isfinite(
        df[
            MORAN_B
        ].to_numpy(float)
    )
    &
    np.isfinite(
        df[
            MORAN_C
        ].to_numpy(float)
    )
)


A = df.loc[
    moran_valid,
    MORAN_A
].to_numpy(float)


B = df.loc[
    moran_valid,
    MORAN_B
].to_numpy(float)


C = df.loc[
    moran_valid,
    MORAN_C
].to_numpy(float)


P = patterned[
    moran_valid
]


print("=" * 100)
print("BR8667 MORAN + RAW COSINE ADVISOR FIGURE")
print("=" * 100)

print(
    "Moran genes:",
    len(
        A
    )
)

print(
    "Visium patterned:",
    int(
        P.sum()
    )
)

print(
    "Visium non-patterned:",
    int(
        (
            ~P
        ).sum()
    )
)


# =============================================================================
# MATCHED AXES FOR d/e
# =============================================================================

de_xlim = padded_limits(
    A
)


de_ylim = padded_limits(
    np.concatenate([
        B,
        C
    ])
)


# =============================================================================
# COSINE DATA
# =============================================================================

cos_valid = (
    evaluable
    &
    np.isfinite(
        df[
            COS_AB
        ].to_numpy(float)
    )
    &
    np.isfinite(
        df[
            COS_AC
        ].to_numpy(float)
    )
    &
    np.isfinite(
        df[
            COS_BC
        ].to_numpy(float)
    )
)


AB = df[
    COS_AB
].to_numpy(float)

AC = df[
    COS_AC
].to_numpy(float)

BC = df[
    COS_BC
].to_numpy(float)


pat_mask = (
    cos_valid
    &
    patterned
)


non_mask = (
    cos_valid
    &
    ~patterned
)


# =============================================================================
# FIGURE LAYOUT
#
# d | e | f
# g spans full width
# h | i
# =============================================================================

fig = plt.figure(
    figsize=(16, 14)
)


gs = fig.add_gridspec(
    nrows=3,
    ncols=6,
    height_ratios=[
        1.0,
        0.72,
        0.82
    ],
    left=0.065,
    right=0.975,
    top=0.91,
    bottom=0.075,
    wspace=0.55,
    hspace=0.45
)


# Moran pairwise scatterplots
ax_d = fig.add_subplot(
    gs[0, 0:2]
)

ax_e = fig.add_subplot(
    gs[0, 2:4]
)

ax_f = fig.add_subplot(
    gs[0, 4:6]
)


# Moran distributions
ax_g = fig.add_subplot(
    gs[1, 1:5]
)


# Cosine distributions
ax_h = fig.add_subplot(
    gs[2, 0:3]
)

ax_i = fig.add_subplot(
    gs[2, 3:6]
)


# =============================================================================
# d — VISIUM vs MEASURED MORAN
# =============================================================================

stats_d = plot_moran_scatter(
    ax=ax_d,
    x=A,
    y=B,
    pattern=P,
    xlim=de_xlim,
    ylim=de_ylim,
    xlabel="Visium Moran's I",
    ylabel="Measured Xenium Moran's I",
    title="d   Visium vs measured Xenium"
)


# =============================================================================
# e — VISIUM vs OOF MORAN
# =============================================================================

stats_e = plot_moran_scatter(
    ax=ax_e,
    x=A,
    y=C,
    pattern=P,
    xlim=de_xlim,
    ylim=de_ylim,
    xlabel="Visium Moran's I",
    ylabel="ENVI OOF Moran's I",
    title="e   Visium vs ENVI OOF"
)


# =============================================================================
# f — MEASURED XENIUM vs OOF MORAN
#
# Direct test of whether ENVI OOF preserves the Moran's I measured
# experimentally by Xenium.
# =============================================================================

f_lim = padded_limits(
    np.concatenate([
        B,
        C
    ])
)


stats_f = plot_moran_scatter(
    ax=ax_f,
    x=B,
    y=C,
    pattern=P,
    xlim=f_lim,
    ylim=f_lim,
    xlabel="Measured Xenium Moran's I",
    ylabel="ENVI OOF Moran's I",
    title="f   Measured Xenium vs ENVI OOF"
)


ax_f.set_aspect(
    "equal",
    adjustable="box"
)


# =============================================================================
# g — MORAN'S-I DISTRIBUTIONS
#
# x = gene count
# y = Moran's I
#
# Three SOURCE distributions:
#     Visium
#     Measured Xenium
#     ENVI OOF
# =============================================================================

moran_min = min(
    np.nanmin(A),
    np.nanmin(B),
    np.nanmin(C)
)

moran_max = max(
    np.nanmax(A),
    np.nanmax(B),
    np.nanmax(C)
)


moran_bins = np.linspace(
    moran_min,
    moran_max,
    26
)


stats_g = horizontal_distribution(
    ax=ax_g,
    arrays=[
        A,
        B,
        C
    ],
    labels=[
        "Visium",
        "Measured Xenium",
        "ENVI OOF"
    ],
    colors=[
        COLOR_VISIUM,
        COLOR_MEASURED,
        COLOR_OOF
    ],
    bins=moran_bins,
    title="g   Moran's-I distributions",
    ylabel="Moran's I"
)


# =============================================================================
# h — VISIUM-PATTERNED RAW-COSINE DISTRIBUTIONS
#
# x = gene count
# y = raw cosine similarity
#
# NOTE:
# Cosine is pairwise, so the three curves are:
#     Visium <-> Measured
#     Visium <-> OOF
#     Measured <-> OOF
# =============================================================================

cos_bins = np.linspace(
    0,
    1,
    26
)


stats_h = horizontal_distribution(
    ax=ax_h,
    arrays=[
        AB[pat_mask],
        AC[pat_mask],
        BC[pat_mask]
    ],
    labels=[
        "Visium ↔ Measured",
        "Visium ↔ OOF",
        "Measured ↔ OOF"
    ],
    colors=[
        COLOR_AB,
        COLOR_AC,
        COLOR_BC
    ],
    bins=cos_bins,
    title=(
        "h   Visium-patterned genes — "
        "raw cosine distributions"
    ),
    ylabel="raw cosine similarity"
)


ax_h.set_ylim(
    0,
    1
)


# =============================================================================
# i — VISIUM NON-PATTERNED RAW-COSINE DISTRIBUTIONS
# =============================================================================

stats_i = horizontal_distribution(
    ax=ax_i,
    arrays=[
        AB[non_mask],
        AC[non_mask],
        BC[non_mask]
    ],
    labels=[
        "Visium ↔ Measured",
        "Visium ↔ OOF",
        "Measured ↔ OOF"
    ],
    colors=[
        COLOR_AB,
        COLOR_AC,
        COLOR_BC
    ],
    bins=cos_bins,
    title=(
        "i   Visium non-patterned genes — "
        "raw cosine distributions"
    ),
    ylabel="raw cosine similarity"
)


ax_i.set_ylim(
    0,
    1
)


# =============================================================================
# PAIRED COSINE STATISTICS
# =============================================================================

pat_AB_AC = paired_wilcoxon(
    AB[
        pat_mask
    ],
    AC[
        pat_mask
    ]
)


pat_AB_BC = paired_wilcoxon(
    AB[
        pat_mask
    ],
    BC[
        pat_mask
    ]
)


non_AB_AC = paired_wilcoxon(
    AB[
        non_mask
    ],
    AC[
        non_mask
    ]
)


non_AB_BC = paired_wilcoxon(
    AB[
        non_mask
    ],
    BC[
        non_mask
    ]
)


# =============================================================================
# TITLE
# =============================================================================

fig.suptitle(
    (
        "Br8667: spatial correspondence of measured Xenium "
        "and ENVI OOF predictions"
    ),
    fontsize=15,
    fontweight="bold",
    y=0.96
)


# ---- commented out overlapping annotation ----
# fig.text(
#     0.33,
#     0.91,
#     (
#         "d/e use the same 300-gene universe "
#         "and identical axis ranges"
#     ),
#     ha="center",
#     fontsize=8,
#     fontweight="bold"
# )


# ---- commented out overlapping annotation ----
# fig.text(
#     0.81,
#     0.91,
#     (
#         "f/g/h use identical bins within each panel"
#     ),
#     ha="center",
#     fontsize=8,
#     fontweight="bold"
# )


# =============================================================================
# CAPTION
# =============================================================================

caption = (
    "Panels d and e compare gene-level Moran's I for measured Xenium and ENVI OOF "
    "against the same Visium reference using identical axis limits. "
    "Panel f compares the complete Moran's-I distributions of Visium, measured Xenium, "
    "and ENVI OOF with gene count on the x-axis and Moran's I on the y-axis. "
    "Panels g and h show raw-cosine distributions for Visium-patterned and Visium "
    "non-patterned genes, respectively. Three pairwise correspondence measures are shown: "
    "Visium↔Measured, Visium↔OOF, and Measured↔OOF. "
    "All three curves/bars within each distribution panel use identical bins."
)


fig.text(
    0.50,
    0.025,
    caption,
    ha="center",
    va="bottom",
    fontsize=7,
    wrap=True
)


# =============================================================================
# SAVE SUMMARY
# =============================================================================

summary_rows = [
    {
        "comparison":
            "Moran_Visium_vs_Measured",

        "n":
            stats_d[
                "n"
            ],

        "spearman":
            stats_d[
                "spearman"
            ],

        "pearson":
            stats_d[
                "pearson"
            ],

        "slope":
            stats_d[
                "slope"
            ],

        "median_delta":
            stats_d[
                "median_delta"
            ],
    },

    {
        "comparison":
            "Moran_Visium_vs_OOF",

        "n":
            stats_e[
                "n"
            ],

        "spearman":
            stats_e[
                "spearman"
            ],

        "pearson":
            stats_e[
                "pearson"
            ],

        "slope":
            stats_e[
                "slope"
            ],

        "median_delta":
            stats_e[
                "median_delta"
            ],
    },


    {
        "comparison":
            "Moran_Measured_vs_OOF",

        "n":
            stats_f[
                "n"
            ],

        "spearman":
            stats_f[
                "spearman"
            ],

        "pearson":
            stats_f[
                "pearson"
            ],

        "slope":
            stats_f[
                "slope"
            ],

        "median_delta":
            stats_f[
                "median_delta"
            ],
    },

    {
        "comparison":
            "Patterned_cosine_VisiumMeasured_vs_VisiumOOF",

        **pat_AB_AC,
    },

    {
        "comparison":
            "Patterned_cosine_VisiumMeasured_vs_MeasuredOOF",

        **pat_AB_BC,
    },

    {
        "comparison":
            "Nonpatterned_cosine_VisiumMeasured_vs_VisiumOOF",

        **non_AB_AC,
    },

    {
        "comparison":
            "Nonpatterned_cosine_VisiumMeasured_vs_MeasuredOOF",

        **non_AB_BC,
    },
]


summary = pd.DataFrame(
    summary_rows
)


SUMMARY = (
    OUT
    / f"{BR}_advisor_moran_cosine_summary.csv"
)


summary.to_csv(
    SUMMARY,
    index=False
)


# =============================================================================
# SAVE FULL GENE TABLE
# =============================================================================

gene_table = df[
    [
        "gene",
        MORAN_A,
        MORAN_B,
        MORAN_C,
        COS_AB,
        COS_AC,
        COS_BC,
        PATTERN_A,
        EVALUABLE,
    ]
].copy()


GENE_CSV = (
    OUT
    / f"{BR}_advisor_moran_cosine_gene_metrics.csv"
)


gene_table.to_csv(
    GENE_CSV,
    index=False
)


# =============================================================================
# SAVE FIGURE
# =============================================================================

PNG = (
    OUT
    / f"fig_{BR}_advisor_moran_cosine_validation.png"
)

PDF = (
    OUT
    / f"fig_{BR}_advisor_moran_cosine_validation.pdf"
)


fig.savefig(
    PNG,
    dpi=400,
    bbox_inches="tight"
)


fig.savefig(
    PDF,
    bbox_inches="tight"
)


plt.close(
    fig
)


# =============================================================================
# DONE
# =============================================================================

print()
print("=" * 100)
print("COMPLETE")
print("=" * 100)

print()
print(
    summary.to_string(
        index=False
    )
)

print()
print(
    "PNG:",
    PNG
)

print(
    "PDF:",
    PDF
)

print(
    "Summary:",
    SUMMARY
)

print(
    "Gene table:",
    GENE_CSV
)
