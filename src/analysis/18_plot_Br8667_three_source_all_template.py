#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    / "PI_style_pairwise_final"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


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


def as_bool(s):

    if pd.api.types.is_bool_dtype(s):
        return s.to_numpy(bool)

    return (
        s.astype(str)
        .str.lower()
        .isin([
            "true",
            "1",
            "yes"
        ])
        .to_numpy(bool)
    )


def cell_cpm(X):

    lib = np.asarray(
        X.sum(axis=1),
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


def display_normalize(x):
    """
    Map expression to [0,1] using the 1st-99th percentile.

    DISPLAY ONLY.
    Metrics are never computed from this transformed vector.
    """

    x = np.asarray(
        x,
        dtype=float
    )

    out = np.zeros_like(
        x,
        dtype=float
    )

    finite = np.isfinite(x)

    if not finite.any():
        return out

    vals = x[
        finite
    ]

    lo = np.nanpercentile(
        vals,
        1
    )

    hi = np.nanpercentile(
        vals,
        99
    )

    if (
        not np.isfinite(lo)
        or not np.isfinite(hi)
        or hi <= lo
    ):
        hi = lo + 1e-12

    out[
        finite
    ] = np.clip(
        (
            x[
                finite
            ]
            -
            lo
        )
        /
        (
            hi
            -
            lo
        ),
        0,
        1
    )

    return out


def select_display_genes(
    df,
    cosine_col,
    pattern_col
):

    evaluable = as_bool(
        df[
            "evaluable_all3"
        ]
    )

    patterned_flag = as_bool(
        df[
            pattern_col
        ]
    )

    cosine = pd.to_numeric(
        df[
            cosine_col
        ],
        errors="coerce"
    ).to_numpy(float)

    valid = (
        evaluable
        &
        np.isfinite(
            cosine
        )
    )

    pat = df.loc[
        valid
        &
        patterned_flag
    ].copy()

    non = df.loc[
        valid
        &
        ~patterned_flag
    ].copy()

    pat[
        cosine_col
    ] = pd.to_numeric(
        pat[
            cosine_col
        ],
        errors="coerce"
    )

    non[
        cosine_col
    ] = pd.to_numeric(
        non[
            cosine_col
        ],
        errors="coerce"
    )

    if len(pat) < 3:
        raise RuntimeError(
            f"Too few patterned genes: {len(pat)}"
        )

    if len(non) < 1:
        raise RuntimeError(
            "No non-patterned genes."
        )

    pat = pat.sort_values(
        cosine_col,
        ascending=True
    )

    non = non.sort_values(
        cosine_col,
        ascending=True
    )

    # highest patterned
    high = pat.iloc[
        -1
    ]

    # closest gene to patterned median
    pat_median = float(
        pat[
            cosine_col
        ].median()
    )

    mid_idx = (
        pat[
            cosine_col
        ]
        .sub(
            pat_median
        )
        .abs()
        .idxmin()
    )

    median_gene = pat.loc[
        mid_idx
    ]

    # lowest patterned
    low = pat.iloc[
        0
    ]

    # non-patterned median control
    non_median = float(
        non[
            cosine_col
        ].median()
    )

    ctrl_idx = (
        non[
            cosine_col
        ]
        .sub(
            non_median
        )
        .abs()
        .idxmin()
    )

    control = non.loc[
        ctrl_idx
    ]

    return [
        (
            high,
            "patterned, highest raw cosine"
        ),
        (
            median_gene,
            "patterned, median raw cosine"
        ),
        (
            low,
            "patterned, lowest raw cosine"
        ),
        (
            control,
            "not patterned control"
        ),
    ]


def format_map_axis(
    ax,
    xmin,
    xmax,
    ymin,
    ymax
):

    ax.set_xlim(
        xmin,
        xmax
    )

    # same orientation for all sources
    ax.set_ylim(
        ymax,
        ymin
    )

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.7)
        spine.set_color("0.25")


def plot_spot_map(
    ax,
    values,
    xmin,
    xmax,
    ymin,
    ymax
):

    disp = display_normalize(
        values
    )

    ax.scatter(
        xy_spot[:, 0],
        xy_spot[:, 1],
        c=disp,
        cmap="Reds",
        vmin=0,
        vmax=1,
        s=8,
        linewidths=0,
        rasterized=True
    )

    format_map_axis(
        ax,
        xmin,
        xmax,
        ymin,
        ymax
    )


def plot_cell_map(
    ax,
    values,
    xmin,
    xmax,
    ymin,
    ymax
):

    disp = display_normalize(
        values
    )

    ax.scatter(
        cell_xy[:, 0],
        cell_xy[:, 1],
        c=disp,
        cmap="Reds",
        vmin=0,
        vmax=1,
        s=1.0,
        linewidths=0,
        rasterized=True
    )

    format_map_axis(
        ax,
        xmin,
        xmax,
        ymin,
        ymax
    )


def add_panel_letter(
    ax,
    letter,
    x=-0.14,
    y=1.04
):

    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="bottom"
    )


def plot_cosine_hist(
    ax,
    values,
    selected,
    cosine_col,
    patterned,
    letter,
    title
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

    bins = np.linspace(
        0,
        1,
        26
    )

    ax.hist(
        values,
        bins=bins,
        orientation="horizontal",
        color="0.78",
        edgecolor="white",
        linewidth=0.6
    )

    med = float(
        np.nanmedian(
            values
        )
    )

    ax.axhline(
        med,
        color="black",
        linestyle="--",
        linewidth=0.8
    )

    ax.set_ylim(
        0,
        1
    )

    ax.set_ylabel(
        "raw cosine similarity",
        fontsize=8
    )

    ax.set_xlabel(
        "gene count",
        fontsize=8
    )

    ax.tick_params(
        labelsize=7
    )

    ax.set_title(
        f"{letter}   {title} (n={len(values):,})",
        fontsize=9,
        fontweight="bold",
        loc="left"
    )

    ax.text(
        0.03,
        0.94,
        f"median={med:.3f}",
        transform=ax.transAxes,
        fontsize=7,
        va="top"
    )

    xmax = ax.get_xlim()[1]

    for row, desc in selected:

        is_patterned = (
            desc.startswith(
                "patterned"
            )
        )

        if (
            patterned
            and
            not is_patterned
        ):
            continue

        if (
            not patterned
            and
            is_patterned
        ):
            continue

        value = float(
            row[
                cosine_col
            ]
        )

        gene = str(
            row[
                "gene"
            ]
        )

        ax.axhline(
            value,
            color="red",
            linewidth=1.2
        )

        ax.text(
            xmax * 1.02,
            value,
            gene,
            fontsize=7,
            fontstyle="italic",
            va="center",
            ha="left",
            clip_on=False
        )


# =============================================================================
# LOAD METRICS + PSEUDOSPOT MATRICES
# =============================================================================

for f in [
    METRICS,
    SPOTS,
    OOF,
    ALIGN
]:

    if not f.exists():
        raise FileNotFoundError(
            f
        )


res = pd.read_csv(
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


if not (
    A_cpm.shape
    ==
    B_cpm.shape
    ==
    C_cpm.shape
):
    raise RuntimeError(
        "Spot matrices do not have identical dimensions."
    )


if A_cpm.shape[1] != len(
    genes
):
    raise RuntimeError(
        "Gene dimension mismatch."
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


gene_to_idx = {
    g: i
    for i, g
    in enumerate(
        genes
    )
}


print("=" * 100)
print("BR8667 THREE-SOURCE FIGURE GENERATOR")
print("=" * 100)

print(
    "Genes:",
    len(genes)
)

print(
    "Occupied spots:",
    len(xy_spot)
)


# =============================================================================
# LOAD CELL-LEVEL B/C + STALIGN COORDINATES
# =============================================================================

o = ad.read_h5ad(
    OOF
)

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


o_ids = (
    o.obs[
        "original_cell_id"
    ]
    .astype(str)
    .to_numpy()
)


st_row = np.array([
    st_lookup.get(
        x,
        -1
    )
    for x
    in o_ids
])


has_st = (
    st_row >= 0
)


idx = np.flatnonzero(
    has_st
)


match_prob = np.full(
    o.n_obs,
    np.nan
)

vis_spot = np.full(
    o.n_obs,
    "",
    dtype=object
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
        idx
    ]
]


match_prob[
    idx
] = pd.to_numeric(
    sub[
        "match_prob"
    ],
    errors="coerce"
).to_numpy()


vis_spot[
    idx
] = (
    sub[
        "visium_spot"
    ]
    .fillna("")
    .astype(str)
    .to_numpy()
)


x_aligned[
    idx
] = pd.to_numeric(
    sub[
        "x_aligned"
    ],
    errors="coerce"
).to_numpy()


y_aligned[
    idx
] = pd.to_numeric(
    sub[
        "y_aligned"
    ],
    errors="coerce"
).to_numpy()


use = (
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
    (
        vis_spot != ""
    )
    &
    (
        vis_spot != "nan"
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


cell_rows = np.flatnonzero(
    use
)


cell_xy = np.column_stack([
    x_aligned[
        cell_rows
    ],
    y_aligned[
        cell_rows
    ]
])


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


print(
    "Aligned cells inside physical Visium spots:",
    len(
        cell_rows
    )
)


# =============================================================================
# COMMON MAP EXTENT
# =============================================================================

all_x = np.concatenate([
    xy_spot[:, 0],
    cell_xy[:, 0]
])

all_y = np.concatenate([
    xy_spot[:, 1],
    cell_xy[:, 1]
])


xmin = float(
    np.nanmin(
        all_x
    )
)

xmax = float(
    np.nanmax(
        all_x
    )
)

ymin = float(
    np.nanmin(
        all_y
    )
)

ymax = float(
    np.nanmax(
        all_y
    )
)


xpad = (
    xmax
    -
    xmin
) * 0.02

ypad = (
    ymax
    -
    ymin
) * 0.02


xmin -= xpad
xmax += xpad
ymin -= ypad
ymax += ypad


# =============================================================================
# FIXED DISPLAY GENES FOR THE 300-GENE THREE-SOURCE BENCHMARK
#
# Selected ONCE from the A-vs-C external-validation comparison and then
# frozen across A-vs-B, A-vs-C, and B-vs-C.
#
# The four display genes from the previous full-transcriptome figure
# (MT-CO1, EEF1G, KCNQ1OT1, NMT1) are not present in this 300-gene
# OOF benchmark and therefore cannot be used here.
# =============================================================================

FIXED_DISPLAY_GENES = [
    ("MT-ND1",  "patterned, highest raw cosine"),
    ("RASGRF2", "patterned, median raw cosine"),
    ("NXPH4",   "patterned, lowest raw cosine"),
    ("NPIPA1",  "not patterned control"),
]


def get_fixed_display_genes(df):

    available = set(
        df["gene"].astype(str)
    )

    out = []

    for gene, description in FIXED_DISPLAY_GENES:

        if gene not in available:
            raise RuntimeError(
                f"Fixed display gene absent from 300-gene benchmark: {gene}"
            )

        hit = df.loc[
            df["gene"].astype(str) == gene
        ]

        if len(hit) != 1:
            raise RuntimeError(
                f"Expected one row for {gene}; found {len(hit)}"
            )

        out.append(
            (
                hit.iloc[0],
                description
            )
        )

    return out


# =============================================================================
# THREE PAIRS
# =============================================================================

PAIRS = [
    {
        "code":
            "A_vs_B",

        "title":
            "Visium vs Measured Xenium",

        "source1":
            "Visium",

        "source2":
            "Measured Xenium",

        "map1_type":
            "spot",

        "map2_type":
            "cell_B",

        "spot1":
            A_log,

        "spot2":
            B_log,

        "cos":
            "cosine_A_Visium_vs_B_measured",

        "spearman":
            "spearman_A_Visium_vs_B_measured",

        "moran1":
            "moran_A_Visium",

        "moran2":
            "moran_B_measured",

        "pattern":
            "patterned_A_Visium",

        "pattern_source":
            "Visium",
    },

    {
        "code":
            "A_vs_C",

        "title":
            "Visium vs ENVI OOF",

        "source1":
            "Visium",

        "source2":
            "ENVI OOF",

        "map1_type":
            "spot",

        "map2_type":
            "cell_C",

        "spot1":
            A_log,

        "spot2":
            C_log,

        "cos":
            "cosine_A_Visium_vs_C_OOF",

        "spearman":
            "spearman_A_Visium_vs_C_OOF",

        "moran1":
            "moran_A_Visium",

        "moran2":
            "moran_C_OOF",

        "pattern":
            "patterned_A_Visium",

        "pattern_source":
            "Visium",
    },

    {
        "code":
            "B_vs_C",

        "title":
            "Measured Xenium vs ENVI OOF",

        "source1":
            "Measured Xenium",

        "source2":
            "ENVI OOF",

        "map1_type":
            "cell_B",

        "map2_type":
            "cell_C",

        "spot1":
            B_log,

        "spot2":
            C_log,

        "cos":
            "cosine_B_measured_vs_C_OOF",

        "spearman":
            "spearman_B_measured_vs_C_OOF",

        "moran1":
            "moran_B_measured",

        "moran2":
            "moran_C_OOF",

        "pattern":
            "patterned_B_measured",

        "pattern_source":
            "Measured Xenium",
    },
]


# =============================================================================
# GENERATE ALL THREE
# =============================================================================

for spec in PAIRS:

    print()
    print("=" * 100)
    print(
        spec[
            "code"
        ],
        ":",
        spec[
            "title"
        ]
    )
    print("=" * 100)


    # Same four 300-gene display genes in all three pairwise figures.
    selected = get_fixed_display_genes(
        res
    )


    selected_records = []

    for row_label, (
        row,
        description
    ) in enumerate(
        selected
    ):

        print(
            f"{description:36s}",
            f"{row['gene']:15s}",
            f"cos={float(row[spec['cos']]):.3f}",
            f"sp={float(row[spec['spearman']]):+.3f}"
        )

        selected_records.append({
            "row":
                row_label + 1,

            "description":
                description,

            "gene":
                row[
                    "gene"
                ],

            "raw_cosine":
                float(
                    row[
                        spec[
                            "cos"
                        ]
                    ]
                ),

            "spearman":
                float(
                    row[
                        spec[
                            "spearman"
                        ]
                    ]
                ),
        })


    pd.DataFrame(
        selected_records
    ).to_csv(
        OUT
        / f"{BR}_{spec['code']}_selected_display_genes.csv",
        index=False
    )


    # =========================================================================
    # SAME GEOMETRY AS UPLOADED TEMPLATE:
    #
    # row 1: a | b | c | d
    # row 2: a | b | c | d
    # row 3: a | b | c | e
    # row 4: a | b | c | f
    # =========================================================================

    fig = plt.figure(
        figsize=(
            14.5,
            13.5
        )
    )


    gs = fig.add_gridspec(
        nrows=4,
        ncols=4,
        width_ratios=[
            1.0,
            1.0,
            1.10,
            1.35
        ],
        height_ratios=[
            1,
            1,
            1,
            1
        ],
        left=0.17,
        right=0.94,
        top=0.89,
        bottom=0.12,
        wspace=0.38,
        hspace=0.30
    )


    axes_a = [
        fig.add_subplot(
            gs[
                i,
                0
            ]
        )
        for i in range(4)
    ]


    axes_b = [
        fig.add_subplot(
            gs[
                i,
                1
            ]
        )
        for i in range(4)
    ]


    axes_c = [
        fig.add_subplot(
            gs[
                i,
                2
            ]
        )
        for i in range(4)
    ]


    # d occupies top TWO rows
    axd = fig.add_subplot(
        gs[
            0:2,
            3
        ]
    )

    # e = row 3
    axe = fig.add_subplot(
        gs[
            2,
            3
        ]
    )

    # f = row 4
    axf = fig.add_subplot(
        gs[
            3,
            3
        ]
    )


    # =========================================================================
    # PANELS a / b / c
    # =========================================================================

    for i, (
        r,
        description
    ) in enumerate(
        selected
    ):

        gene = str(
            r[
                "gene"
            ]
        )

        j = gene_to_idx[
            gene
        ]


        # ---------------------------------------------------------------------
        # a
        # ---------------------------------------------------------------------

        if spec[
            "map1_type"
        ] == "spot":

            plot_spot_map(
                axes_a[
                    i
                ],
                spec[
                    "spot1"
                ][
                    :,
                    j
                ],
                xmin,
                xmax,
                ymin,
                ymax
            )

        elif spec[
            "map1_type"
        ] == "cell_B":

            plot_cell_map(
                axes_a[
                    i
                ],
                B_cell_log[
                    cell_rows,
                    j
                ],
                xmin,
                xmax,
                ymin,
                ymax
            )

        else:

            raise RuntimeError(
                spec[
                    "map1_type"
                ]
            )


        axes_a[
            i
        ].text(
            -0.11,
            0.50,
            f"{gene}\n{description}",
            transform=axes_a[
                i
            ].transAxes,
            ha="right",
            va="center",
            fontsize=8
        )


        # ---------------------------------------------------------------------
        # b
        # ---------------------------------------------------------------------

        if spec[
            "map2_type"
        ] == "cell_B":

            plot_cell_map(
                axes_b[
                    i
                ],
                B_cell_log[
                    cell_rows,
                    j
                ],
                xmin,
                xmax,
                ymin,
                ymax
            )

        elif spec[
            "map2_type"
        ] == "cell_C":

            plot_cell_map(
                axes_b[
                    i
                ],
                C_cell_log[
                    cell_rows,
                    j
                ],
                xmin,
                xmax,
                ymin,
                ymax
            )

        else:

            raise RuntimeError(
                spec[
                    "map2_type"
                ]
            )


        # ---------------------------------------------------------------------
        # c
        # ---------------------------------------------------------------------

        x = spec[
            "spot1"
        ][
            :,
            j
        ]

        y = spec[
            "spot2"
        ][
            :,
            j
        ]


        axes_c[
            i
        ].scatter(
            x,
            y,
            s=6,
            c="0.30",
            alpha=0.40,
            linewidths=0,
            rasterized=True
        )


        axes_c[
            i
        ].set_xlabel(
            f"{spec['source1']}\nlog10(CPM + 1)",
            fontsize=7
        )

        axes_c[
            i
        ].set_ylabel(
            f"{spec['source2']}\nlog10(CPM + 1)",
            fontsize=7
        )


        axes_c[
            i
        ].tick_params(
            labelsize=7
        )


        axes_c[
            i
        ].text(
            0.03,
            0.97,
            (
                f"raw cosine "
                f"{float(r[spec['cos']]):.2f}\n"
                f"Spearman "
                f"{float(r[spec['spearman']]):+.2f}"
            ),
            transform=axes_c[
                i
            ].transAxes,
            ha="left",
            va="top",
            fontsize=7,
            bbox=dict(
                facecolor="white",
                edgecolor="0.7",
                alpha=0.85,
                boxstyle="round,pad=0.2"
            )
        )


    # headings
    axes_a[
        0
    ].set_title(
        spec[
            "source1"
        ],
        fontsize=10,
        fontweight="bold",
        loc="left"
    )


    axes_b[
        0
    ].set_title(
        (
            spec[
                "source2"
            ]
            +
            (
                " aligned"
                if spec[
                    "source2"
                ] != "Visium"
                else
                ""
            )
        ),
        fontsize=10,
        fontweight="bold",
        loc="left"
    )


    axes_c[
        0
    ].set_title(
        "Normalized expression\ncorrespondence",
        fontsize=10,
        fontweight="bold",
        loc="left"
    )


    add_panel_letter(
        axes_a[
            0
        ],
        "a"
    )

    add_panel_letter(
        axes_b[
            0
        ],
        "b"
    )

    add_panel_letter(
        axes_c[
            0
        ],
        "c"
    )


    # =========================================================================
    # PANEL d — MORAN'S I: SOURCE 1 vs SOURCE 2
    # =========================================================================

    mx = pd.to_numeric(
        res[
            spec[
                "moran1"
            ]
        ],
        errors="coerce"
    ).to_numpy(float)

    my = pd.to_numeric(
        res[
            spec[
                "moran2"
            ]
        ],
        errors="coerce"
    ).to_numpy(float)


    evaluable = as_bool(
        res[
            "evaluable_all3"
        ]
    )

    pattern = as_bool(
        res[
            spec[
                "pattern"
            ]
        ]
    )


    valid = (
        evaluable
        &
        np.isfinite(
            mx
        )
        &
        np.isfinite(
            my
        )
    )


    xv = mx[
        valid
    ]

    yv = my[
        valid
    ]

    pv = pattern[
        valid
    ]


    # non-patterned first
    axd.scatter(
        xv[
            ~pv
        ],
        yv[
            ~pv
        ],
        s=14,
        alpha=0.50,
        linewidths=0,
        c="#9ecae1",
        label="not patterned",
        rasterized=True
    )


    # patterned
    axd.scatter(
        xv[
            pv
        ],
        yv[
            pv
        ],
        s=16,
        alpha=0.72,
        linewidths=0,
        c="#f28e2b",
        label="patterned",
        rasterized=True
    )


    lo = float(
        min(
            np.nanmin(
                xv
            ),
            np.nanmin(
                yv
            )
        )
    )

    hi = float(
        max(
            np.nanmax(
                xv
            ),
            np.nanmax(
                yv
            )
        )
    )


    pad = max(
        0.02,
        (
            hi - lo
        )
        *
        0.05
    )

    lo -= pad
    hi += pad


    # y = x reference
    axd.plot(
        [
            lo,
            hi
        ],
        [
            lo,
            hi
        ],
        linestyle="--",
        linewidth=1.0,
        color="#377eb8"
    )


    # linear fit
    slope, intercept = np.polyfit(
        xv,
        yv,
        1
    )


    xx = np.linspace(
        lo,
        hi,
        200
    )


    axd.plot(
        xx,
        intercept
        +
        slope
        *
        xx,
        color="#ff7f0e",
        linewidth=1.2
    )


    rho = spearmanr(
        xv,
        yv
    ).statistic

    pear = pearsonr(
        xv,
        yv
    ).statistic


    delta = yv - xv


    axd.set_xlim(
        lo,
        hi
    )

    axd.set_ylim(
        lo,
        hi
    )

    axd.set_aspect(
        "equal",
        adjustable="box"
    )


    axd.set_xlabel(
        f"{spec['source1']} Moran's I",
        fontsize=8
    )

    axd.set_ylabel(
        f"{spec['source2']} Moran's I",
        fontsize=8
    )

    axd.tick_params(
        labelsize=7
    )


    axd.set_title(
        "d   Moran's-I correspondence",
        fontsize=9,
        fontweight="bold",
        loc="left"
    )


    axd.text(
        0.03,
        0.97,
        (
            f"genes n={len(xv)}\n"
            f"Spearman ρ={rho:.3f}\n"
            f"Pearson r={pear:.3f}\n"
            f"slope={slope:.3f}\n"
            f"median ΔI="
            f"{np.nanmedian(delta):+.3f}"
        ),
        transform=axd.transAxes,
        va="top",
        fontsize=7,
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.75
        )
    )


    axd.legend(
        frameon=False,
        fontsize=6,
        loc="lower right"
    )


    # =========================================================================
    # PANEL e — PATTERNED RAW COSINE
    # =========================================================================

    cosine = pd.to_numeric(
        res[
            spec[
                "cos"
            ]
        ],
        errors="coerce"
    ).to_numpy(float)


    pat_values = cosine[
        evaluable
        &
        pattern
        &
        np.isfinite(
            cosine
        )
    ]


    plot_cosine_hist(
        axe,
        pat_values,
        selected,
        spec[
            "cos"
        ],
        patterned=True,
        letter="e",
        title="spatially patterned"
    )


    # =========================================================================
    # PANEL f — NON-PATTERNED RAW COSINE
    # =========================================================================

    non_values = cosine[
        evaluable
        &
        ~pattern
        &
        np.isfinite(
            cosine
        )
    ]


    plot_cosine_hist(
        axf,
        non_values,
        selected,
        spec[
            "cos"
        ],
        patterned=False,
        letter="f",
        title="not patterned"
    )


    # =========================================================================
    # TITLE
    # =========================================================================

    fig.suptitle(
        (
            f"{BR}: "
            f"{spec['title']}"
        ),
        fontsize=14,
        fontweight="bold",
        y=0.965
    )


    # =========================================================================
    # DISPLAY COLOR STRIP
    # =========================================================================

    cax = fig.add_axes([
        0.19,
        0.073,
        0.18,
        0.009
    ])


    gradient = np.linspace(
        0,
        1,
        256
    )[
        None,
        :
    ]


    cax.imshow(
        gradient,
        aspect="auto",
        cmap="Reds",
        origin="lower"
    )

    cax.set_xticks([
        0,
        255
    ])

    cax.set_xticklabels([
        "low",
        "high"
    ])

    cax.tick_params(
        axis="x",
        labelsize=6,
        length=0
    )

    cax.set_yticks([])

    cax.set_title(
        "normalized expression (1st–99th pct; display only)",
        fontsize=6
    )


    # =========================================================================
    # CAPTION
    # =========================================================================

    caption = (
        f"Panels a–b show {spec['source1']} and {spec['source2']} spatial expression. "
        f"Spatial maps are independently normalized to the 1st–99th percentile for display only. "
        f"Panel c compares log10(CPM + 1) expression across the same "
        f"{len(xy_spot):,} occupied Visium spots. "
        f"Panel d compares Moran's I for the same 300 genes using the same "
        f"6-nearest-neighbor spatial graph for both sources. "
        f"Spatial patterning is defined independently from {spec['pattern_source']} "
        f"Moran's I with BH q < 0.05 and I > 0. "
        f"Panels e–f show raw cosine-similarity distributions for patterned and "
        f"non-patterned genes, respectively."
    )


    fig.text(
        0.50,
        0.025,
        caption,
        ha="center",
        va="bottom",
        fontsize=6.5,
        wrap=True
    )


    # =========================================================================
    # SAVE
    # =========================================================================

    PNG = (
        OUT
        / (
            f"fig_{BR}_"
            f"{spec['code']}_"
            "PI_template_raw_cosine.png"
        )
    )


    PDF = (
        OUT
        / (
            f"fig_{BR}_"
            f"{spec['code']}_"
            "PI_template_raw_cosine.pdf"
        )
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


    print(
        "Saved PNG:",
        PNG
    )

    print(
        "Saved PDF:",
        PDF
    )


# =============================================================================
# DONE
# =============================================================================

o.file.close()


print()
print("=" * 100)
print("ALL THREE BR8667 FIGURES COMPLETE")
print("=" * 100)

for f in sorted(
    OUT.glob(
        "fig_Br8667_*_PI_template_raw_cosine.png"
    )
):

    print(
        f.name
    )
