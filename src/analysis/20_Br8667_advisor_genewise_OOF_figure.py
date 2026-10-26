#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy import sparse
from scipy.stats import spearmanr, pearsonr


# =============================================================================
# CONFIG
# =============================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BR = "Br8667"
MATCH_THRESH = 0.85


BASE = (
    ROOT
    / "results/analysis/three_source_validation"
    / BR
)

METRICS = (
    BASE
    / f"{BR}_three_source_gene_metrics.csv"
)

SPOTS = (
    BASE
    / f"{BR}_three_source_spot_matrices.npz"
)

OOF = (
    ROOT
    / "outputs/imputation_full/ENVI_N24"
    / BR
    / "ex5_1/envi/oof_300/oof_predicted_300genes.h5ad"
)

ALIGN = (
    ROOT
    / "data/metadata/stalign_results_update"
    / f"{BR}_STalign_result"
    / f"{BR}_aligned_coordinates.csv"
)

OUT = (
    BASE
    / "advisor_genewise_validation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# FIXED DISPLAY GENES
#
# Chosen once from A-vs-C and frozen for all comparisons.
# =============================================================================

DISPLAY_GENES = [
    (
        "MT-ND1",
        "A–C selected: high patterned"
    ),
    (
        "RASGRF2",
        "A–C selected: median patterned"
    ),
    (
        "NXPH4",
        "A–C selected: low patterned"
    ),
    (
        "NPIPA1",
        "A–C selected: non-patterned control"
    ),
]


# =============================================================================
# HELPERS
# =============================================================================

def dense_float(X):

    if sparse.issparse(X):
        X = X.toarray()

    return np.asarray(
        X,
        dtype=np.float32
    )


def cell_cpm(X):

    X = np.asarray(
        X,
        dtype=np.float32
    )

    lib = np.asarray(
        X.sum(
            axis=1
        ),
        dtype=np.float64
    ).ravel()

    lib[
        lib <= 0
    ] = 1.0

    return (
        X
        /
        lib[:, None]
        *
        1e6
    )


def shared_map_vmax(
    x,
    y
):

    z = np.concatenate([
        np.asarray(
            x,
            dtype=float
        ),
        np.asarray(
            y,
            dtype=float
        )
    ])

    z = z[
        np.isfinite(
            z
        )
    ]

    if len(z) == 0:
        return 1.0

    vmax = float(
        np.nanpercentile(
            z,
            99
        )
    )

    if (
        not np.isfinite(vmax)
        or
        vmax <= 0
    ):
        vmax = 1.0

    return vmax


def plot_spatial(
    ax,
    xy,
    values,
    vmin,
    vmax,
    point_size
):

    sc = ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=values,
        s=point_size,
        cmap="Reds",
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        rasterized=True
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_linewidth(
            0.6
        )

    return sc


def padded_limit(
    values,
    lower_zero=True,
    pad_fraction=0.04
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
            0,
            1
        )

    lo = float(
        np.nanmin(
            values
        )
    )

    hi = float(
        np.nanmax(
            values
        )
    )

    if lower_zero:
        lo = min(
            0.0,
            lo
        )

    span = (
        hi - lo
    )

    if span <= 0:
        span = 1.0

    pad = (
        span
        *
        pad_fraction
    )

    if lower_zero:
        lo = 0.0
    else:
        lo -= pad

    hi += pad

    return (
        lo,
        hi
    )


# =============================================================================
# LOAD THREE-SOURCE PSEUDOSPOT RESULTS
# =============================================================================

print("=" * 100)
print("BR8667 ADVISOR GENE-WISE OOF FIGURE")
print("=" * 100)


metrics = pd.read_csv(
    METRICS
)


z = np.load(
    SPOTS,
    allow_pickle=True
)


genes = np.asarray(
    z[
        "genes"
    ],
    dtype=str
)


xy_spot = np.asarray(
    z[
        "xy"
    ],
    dtype=float
)


A_cpm = np.asarray(
    z[
        "visium_cpm"
    ],
    dtype=float
)


B_cpm = np.asarray(
    z[
        "measured_xenium_cpm"
    ],
    dtype=float
)


C_cpm = np.asarray(
    z[
        "oof_envi_cpm"
    ],
    dtype=float
)


A_log = np.log10(
    A_cpm + 1
)

B_log = np.log10(
    B_cpm + 1
)

C_log = np.log10(
    C_cpm + 1
)


spot_gene_idx = {
    gene: i
    for i, gene
    in enumerate(
        genes
    )
}


print(
    "Pseudospots:",
    len(
        xy_spot
    )
)

print(
    "Genes:",
    len(
        genes
    )
)


# =============================================================================
# LOAD OOF CELL-LEVEL DATA
#
# ground_truth_counts = measured Xenium
# count_scale         = OOF ENVI prediction
# =============================================================================

o = ad.read_h5ad(
    OOF
)


if "ground_truth_counts" not in o.layers:
    raise RuntimeError(
        "ground_truth_counts missing from OOF object."
    )


if "count_scale" not in o.layers:
    raise RuntimeError(
        "count_scale missing from OOF object."
    )


oof_genes = np.asarray(
    o.var_names,
    dtype=str
)


oof_gene_idx = {
    gene: i
    for i, gene
    in enumerate(
        oof_genes
    )
}


for gene, _ in DISPLAY_GENES:

    if gene not in spot_gene_idx:
        raise RuntimeError(
            f"{gene} absent from pseudospot matrices."
        )

    if gene not in oof_gene_idx:
        raise RuntimeError(
            f"{gene} absent from OOF object."
        )


B_cell = dense_float(
    o.layers[
        "ground_truth_counts"
    ]
)


C_cell = dense_float(
    o.layers[
        "count_scale"
    ]
)


# Same 300-gene CPM normalization for observed and OOF cell maps.
B_cell_log = np.log10(
    cell_cpm(
        B_cell
    )
    +
    1
)


C_cell_log = np.log10(
    cell_cpm(
        C_cell
    )
    +
    1
)


# =============================================================================
# ALIGN OOF CELLS TO STALIGN
#
# Spatial maps b/c use all confidently aligned cells.
# Scatter panels use the previously generated physical pseudospots.
# =============================================================================

aa = pd.read_csv(
    ALIGN,
    low_memory=False
)


st_lookup = {
    str(v): i
    for i, v
    in enumerate(
        aa[
            "cell_id"
        ].astype(str)
    )
}


oof_ids = (
    o.obs[
        "original_cell_id"
    ]
    .astype(str)
    .to_numpy()
)


st_row = np.array([
    st_lookup.get(
        cid,
        -1
    )
    for cid
    in oof_ids
])


has_st = (
    st_row >= 0
)


ii = np.flatnonzero(
    has_st
)


match_prob = np.full(
    o.n_obs,
    np.nan
)


x_aligned = np.full(
    o.n_obs,
    np.nan
)


y_aligned = np.full(
    o.n_obs,
    np.nan
)


sub = aa.iloc[
    st_row[
        ii
    ]
]


match_prob[
    ii
] = pd.to_numeric(
    sub[
        "match_prob"
    ],
    errors="coerce"
).to_numpy()


x_aligned[
    ii
] = pd.to_numeric(
    sub[
        "x_aligned"
    ],
    errors="coerce"
).to_numpy()


y_aligned[
    ii
] = pd.to_numeric(
    sub[
        "y_aligned"
    ],
    errors="coerce"
).to_numpy()


map_use = (
    has_st
    &
    np.isfinite(
        match_prob
    )
    &
    (
        match_prob
        >
        MATCH_THRESH
    )
    &
    np.isfinite(
        x_aligned
    )
    &
    np.isfinite(
        y_aligned
    )
)


map_rows = np.flatnonzero(
    map_use
)


cell_xy = np.column_stack([
    x_aligned[
        map_rows
    ],
    y_aligned[
        map_rows
    ]
])


print(
    "Confidently aligned cells for maps:",
    len(
        map_rows
    )
)


# =============================================================================
# FIGURE
# =============================================================================

fig = plt.figure(
    figsize=(
        22,
        13
    )
)


gs = fig.add_gridspec(
    nrows=4,
    ncols=6,
    width_ratios=[
        1.0,
        1.0,
        1.0,
        1.13,
        1.13,
        1.13
    ],
    left=0.08,
    right=0.985,
    top=0.90,
    bottom=0.08,
    wspace=0.38,
    hspace=0.42
)


axes = [
    [
        fig.add_subplot(
            gs[
                r,
                c
            ]
        )
        for c in range(
            6
        )
    ]
    for r in range(
        4
    )
]


# =============================================================================
# ROWS
# =============================================================================

summary_rows = []


for row_i, (
    gene,
    row_label
) in enumerate(
    DISPLAY_GENES
):

    print()
    print(
        "-" * 100
    )

    print(
        gene,
        "|",
        row_label
    )


    j_spot = spot_gene_idx[
        gene
    ]

    j_cell = oof_gene_idx[
        gene
    ]


    # -------------------------------------------------------------------------
    # DATA
    # -------------------------------------------------------------------------

    visium = A_log[
        :,
        j_spot
    ]


    measured_spot = B_log[
        :,
        j_spot
    ]


    oof_spot = C_log[
        :,
        j_spot
    ]


    measured_cell = B_cell_log[
        map_rows,
        j_cell
    ]


    oof_cell = C_cell_log[
        map_rows,
        j_cell
    ]


    # -------------------------------------------------------------------------
    # a — VISIUM
    # -------------------------------------------------------------------------

    ax = axes[
        row_i
    ][
        0
    ]


    vis_vmax = float(
        np.nanpercentile(
            visium,
            99
        )
    )


    if (
        not np.isfinite(
            vis_vmax
        )
        or
        vis_vmax <= 0
    ):
        vis_vmax = 1.0


    plot_spatial(
        ax,
        xy_spot,
        visium,
        vmin=0,
        vmax=vis_vmax,
        point_size=7
    )


    ax.text(
        -0.12,
        0.50,
        (
            f"{gene}\n"
            f"{row_label}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=8
    )


    # -------------------------------------------------------------------------
    # b/c — OBSERVED / OOF WITH SAME COLOR SCALE
    # -------------------------------------------------------------------------

    shared_vmax = shared_map_vmax(
        measured_cell,
        oof_cell
    )


    ax_b = axes[
        row_i
    ][
        1
    ]


    ax_c = axes[
        row_i
    ][
        2
    ]


    plot_spatial(
        ax_b,
        cell_xy,
        measured_cell,
        vmin=0,
        vmax=shared_vmax,
        point_size=0.9
    )


    plot_spatial(
        ax_c,
        cell_xy,
        oof_cell,
        vmin=0,
        vmax=shared_vmax,
        point_size=0.9
    )


    # Shared gene-specific colorbar for b and c.
    norm = Normalize(
        vmin=0,
        vmax=shared_vmax
    )


    sm = ScalarMappable(
        norm=norm,
        cmap="Reds"
    )

    sm.set_array([])


    cbar = fig.colorbar(
        sm,
        ax=[
            ax_b,
            ax_c
        ],
        orientation="horizontal",
        fraction=0.035,
        pad=0.035,
        aspect=35
    )


    cbar.ax.tick_params(
        labelsize=6
    )


    cbar.set_label(
        (
            f"{gene} log10(CPM + 1) — "
            "shared observed/OOF scale"
        ),
        fontsize=6
    )


    # -------------------------------------------------------------------------
    # d/e — MATCHED AXES
    #
    # d: Visium vs measured
    # e: Visium vs OOF
    # -------------------------------------------------------------------------

    # Same x-axis for d/e because both use same Visium values.
    xlim = padded_limit(
        visium,
        lower_zero=True
    )


    # Same y-axis across observed and OOF.
    ylim = padded_limit(
        np.concatenate([
            measured_spot,
            oof_spot
        ]),
        lower_zero=True
    )


    # d
    ax_d = axes[
        row_i
    ][
        3
    ]


    ax_d.scatter(
        visium,
        measured_spot,
        s=6,
        alpha=0.35,
        c="0.25",
        linewidths=0,
        rasterized=True
    )


    ax_d.set_xlim(
        xlim
    )

    ax_d.set_ylim(
        ylim
    )


    ax_d.set_xlabel(
        "Visium\nlog10(CPM + 1)",
        fontsize=7
    )

    ax_d.set_ylabel(
        "Xenium observed\nlog10(CPM + 1)",
        fontsize=7
    )


    # e
    ax_e = axes[
        row_i
    ][
        4
    ]


    ax_e.scatter(
        visium,
        oof_spot,
        s=6,
        alpha=0.35,
        c="0.25",
        linewidths=0,
        rasterized=True
    )


    ax_e.set_xlim(
        xlim
    )

    ax_e.set_ylim(
        ylim
    )


    ax_e.set_xlabel(
        "Visium\nlog10(CPM + 1)",
        fontsize=7
    )

    ax_e.set_ylabel(
        "ENVI OOF\nlog10(CPM + 1)",
        fontsize=7
    )


    # Find metrics row.
    mr = metrics.loc[
        metrics[
            "gene"
        ].astype(str)
        ==
        gene
    ]


    if len(
        mr
    ) != 1:

        raise RuntimeError(
            f"Expected exactly one metrics row for {gene}"
        )


    mr = mr.iloc[
        0
    ]


    cos_ab = float(
        mr[
            "cosine_A_Visium_vs_B_measured"
        ]
    )

    sp_ab = float(
        mr[
            "spearman_A_Visium_vs_B_measured"
        ]
    )

    cos_ac = float(
        mr[
            "cosine_A_Visium_vs_C_OOF"
        ]
    )

    sp_ac = float(
        mr[
            "spearman_A_Visium_vs_C_OOF"
        ]
    )


    ax_d.text(
        0.03,
        0.97,
        (
            f"raw cosine {cos_ab:.2f}\n"
            f"Spearman {sp_ab:+.2f}"
        ),
        transform=ax_d.transAxes,
        va="top",
        fontsize=7,
        bbox=dict(
            facecolor="white",
            edgecolor="0.7",
            alpha=0.85,
            boxstyle="round,pad=0.2"
        )
    )


    ax_e.text(
        0.03,
        0.97,
        (
            f"raw cosine {cos_ac:.2f}\n"
            f"Spearman {sp_ac:+.2f}"
        ),
        transform=ax_e.transAxes,
        va="top",
        fontsize=7,
        bbox=dict(
            facecolor="white",
            edgecolor="0.7",
            alpha=0.85,
            boxstyle="round,pad=0.2"
        )
    )


    # -------------------------------------------------------------------------
    # f — DIRECT XENIUM OBSERVED vs OOF
    #
    # Equal x/y range and y=x line.
    # -------------------------------------------------------------------------

    ax_f = axes[
        row_i
    ][
        5
    ]


    direct_values = np.concatenate([
        measured_spot,
        oof_spot
    ])


    direct_lim = padded_limit(
        direct_values,
        lower_zero=True
    )


    ax_f.scatter(
        measured_spot,
        oof_spot,
        s=7,
        alpha=0.40,
        linewidths=0,
        rasterized=True
    )


    ax_f.plot(
        direct_lim,
        direct_lim,
        linestyle="--",
        linewidth=1.0,
        color="0.35"
    )


    ax_f.set_xlim(
        direct_lim
    )

    ax_f.set_ylim(
        direct_lim
    )

    ax_f.set_aspect(
        "equal",
        adjustable="box"
    )


    ax_f.set_xlabel(
        "Xenium observed\nlog10(CPM + 1)",
        fontsize=7
    )

    ax_f.set_ylabel(
        "ENVI OOF\nlog10(CPM + 1)",
        fontsize=7
    )


    valid_direct = (
        np.isfinite(
            measured_spot
        )
        &
        np.isfinite(
            oof_spot
        )
    )


    x_direct = measured_spot[
        valid_direct
    ]

    y_direct = oof_spot[
        valid_direct
    ]


    direct_spearman = spearmanr(
        x_direct,
        y_direct
    ).statistic


    direct_pearson = pearsonr(
        x_direct,
        y_direct
    ).statistic


    direct_rmse = float(
        np.sqrt(
            np.mean(
                (
                    x_direct
                    -
                    y_direct
                )
                **
                2
            )
        )
    )


    ax_f.text(
        0.03,
        0.97,
        (
            f"Spearman {direct_spearman:+.2f}\n"
            f"Pearson {direct_pearson:+.2f}\n"
            f"RMSE {direct_rmse:.2f}"
        ),
        transform=ax_f.transAxes,
        va="top",
        fontsize=7,
        bbox=dict(
            facecolor="white",
            edgecolor="0.7",
            alpha=0.85,
            boxstyle="round,pad=0.2"
        )
    )


    # Tick consistency.
    for ax in [
        ax_d,
        ax_e,
        ax_f
    ]:

        ax.tick_params(
            labelsize=6
        )


    summary_rows.append({
        "gene":
            gene,

        "display_role":
            row_label,

        "map_shared_vmin":
            0.0,

        "map_shared_vmax":
            shared_vmax,

        "DE_xmin":
            xlim[0],

        "DE_xmax":
            xlim[1],

        "DE_ymin":
            ylim[0],

        "DE_ymax":
            ylim[1],

        "AB_raw_cosine":
            cos_ab,

        "AB_spearman":
            sp_ab,

        "AC_raw_cosine":
            cos_ac,

        "AC_spearman":
            sp_ac,

        "BC_direct_spearman":
            direct_spearman,

        "BC_direct_pearson":
            direct_pearson,

        "BC_direct_RMSE":
            direct_rmse,
    })


# =============================================================================
# HEADINGS
# =============================================================================

column_titles = [
    "a   Visium",
    "b   Xenium observed",
    "c   ENVI OOF imputed",
    "d   Visium vs observed",
    "e   Visium vs OOF",
    "f   Observed vs OOF",
]


for col_i, title in enumerate(
    column_titles
):

    axes[
        0
    ][
        col_i
    ].set_title(
        title,
        fontsize=10,
        fontweight="bold",
        loc="left"
    )


# =============================================================================
# TITLE + ANNOTATIONS
# =============================================================================

fig.suptitle(
    (
        "Br8667: gene-wise validation of ENVI out-of-fold predictions "
        "against measured Xenium expression"
    ),
    fontsize=15,
    fontweight="bold",
    y=0.97
)


fig.text(
    0.345,
    0.925,
    "Same gene-specific color range for observed and OOF",
    ha="center",
    fontsize=8,
    fontweight="bold"
)


fig.text(
    0.685,
    0.925,
    "Matched x- and y-axis ranges within each gene",
    ha="center",
    fontsize=8,
    fontweight="bold"
)


fig.text(
    0.925,
    0.925,
    "Equal x/y range + y=x",
    ha="center",
    fontsize=8,
    fontweight="bold"
)


# =============================================================================
# CAPTION
# =============================================================================

caption = (
    "Four fixed genes from the Br8667 300-gene OOF benchmark are shown. "
    "Panel a shows paired Visium expression. Panels b and c show measured Xenium "
    "and ENVI OOF expression in confidently aligned Xenium cells; for each gene, "
    "the observed and OOF maps use exactly the same color scale, while different "
    "genes retain gene-specific scales. Panels d and e compare Visium with measured "
    "and OOF pseudospot expression across the same occupied Visium spots and use "
    "identical x- and y-axis limits within each gene. Panel f directly compares "
    "measured Xenium and ENVI OOF pseudospot expression with equal x/y limits and "
    "a y=x reference line."
)


fig.text(
    0.5,
    0.025,
    caption,
    ha="center",
    va="bottom",
    fontsize=7,
    wrap=True
)


# =============================================================================
# SAVE
# =============================================================================

PNG = (
    OUT
    / f"fig_{BR}_advisor_genewise_OOF_validation.png"
)

PDF = (
    OUT
    / f"fig_{BR}_advisor_genewise_OOF_validation.pdf"
)

CSV = (
    OUT
    / f"{BR}_advisor_genewise_OOF_validation_metrics.csv"
)


pd.DataFrame(
    summary_rows
).to_csv(
    CSV,
    index=False
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


if getattr(
    o,
    "file",
    None
) is not None:

    o.file.close()


print()
print("=" * 100)
print("COMPLETE")
print("=" * 100)

print(
    "PNG:",
    PNG
)

print(
    "PDF:",
    PDF
)

print(
    "Metrics:",
    CSV
)
