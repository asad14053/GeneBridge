#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.stats import spearmanr, pearsonr


# ============================================================
# CONFIG
# ============================================================

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
    / "PI_style_pairwise_v2"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def dense_float(X):

    if sparse.issparse(X):
        X = X.toarray()

    return np.asarray(
        X,
        dtype=np.float32
    )


def cell_cpm(X):

    lib = np.asarray(
        X.sum(axis=1),
        dtype=float
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


# ============================================================
# LOAD SPOT-LEVEL RESULTS
# ============================================================

res = pd.read_csv(
    METRICS
)

z = np.load(
    SPOTS,
    allow_pickle=True
)

genes = np.asarray(
    z["genes"],
    dtype=str
)

xy_spot = np.asarray(
    z["xy"],
    dtype=float
)

A_cpm = np.asarray(
    z["visium_cpm"],
    dtype=float
)

B_cpm = np.asarray(
    z["measured_xenium_cpm"],
    dtype=float
)

C_cpm = np.asarray(
    z["oof_envi_cpm"],
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


gene_to_idx = {
    g: i
    for i, g
    in enumerate(
        genes
    )
}


# ============================================================
# LOAD CELL-LEVEL XENIUM / OOF FOR MAPS
# ============================================================

o = ad.read_h5ad(
    OOF
)

aa = pd.read_csv(
    ALIGN,
    low_memory=False
)


lookup = {
    str(v): i
    for i, v
    in enumerate(
        aa[
            "cell_id"
        ].astype(str)
    )
}


okey = (
    o.obs[
        "original_cell_id"
    ]
    .astype(str)
    .to_numpy()
)


st_row = np.array([
    lookup.get(
        x,
        -1
    )
    for x
    in okey
])


has_st = (
    st_row >= 0
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


ii = np.flatnonzero(
    has_st
)

sub = aa.iloc[
    st_row[ii]
]


match_prob[ii] = pd.to_numeric(
    sub[
        "match_prob"
    ],
    errors="coerce"
).to_numpy()


vis_spot[ii] = (
    sub[
        "visium_spot"
    ]
    .fillna("")
    .astype(str)
    .to_numpy()
)


x_aligned[ii] = pd.to_numeric(
    sub[
        "x_aligned"
    ],
    errors="coerce"
).to_numpy()


y_aligned[ii] = pd.to_numeric(
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
        match_prob > 0.85
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


Bcell = dense_float(
    o.layers[
        "ground_truth_counts"
    ]
)

Ccell = dense_float(
    o.layers[
        "count_scale"
    ]
)


Bcell_log = np.log10(
    cell_cpm(
        Bcell
    )
    +
    1
)

Ccell_log = np.log10(
    cell_cpm(
        Ccell
    )
    +
    1
)


# ============================================================
# PAIRS
# ============================================================

PAIRS = [
    {
        "code":
            "A_vs_B",

        "name1":
            "Visium",

        "name2":
            "Measured Xenium",

        "spot1":
            A_log,

        "spot2":
            B_log,

        "map1":
            "visium",

        "map2":
            "measured",

        "cos":
            "cosine_A_Visium_vs_B_measured",

        "sp":
            "spearman_A_Visium_vs_B_measured",

        "moran1":
            "moran_A_Visium",

        "moran2":
            "moran_B_measured",

        "pattern":
            "patterned_A_Visium",
    },

    {
        "code":
            "A_vs_C",

        "name1":
            "Visium",

        "name2":
            "ENVI OOF",

        "spot1":
            A_log,

        "spot2":
            C_log,

        "map1":
            "visium",

        "map2":
            "oof",

        "cos":
            "cosine_A_Visium_vs_C_OOF",

        "sp":
            "spearman_A_Visium_vs_C_OOF",

        "moran1":
            "moran_A_Visium",

        "moran2":
            "moran_C_OOF",

        "pattern":
            "patterned_A_Visium",
    },

    {
        "code":
            "B_vs_C",

        "name1":
            "Measured Xenium",

        "name2":
            "ENVI OOF",

        "spot1":
            B_log,

        "spot2":
            C_log,

        "map1":
            "measured",

        "map2":
            "oof",

        "cos":
            "cosine_B_measured_vs_C_OOF",

        "sp":
            "spearman_B_measured_vs_C_OOF",

        "moran1":
            "moran_B_measured",

        "moran2":
            "moran_C_OOF",

        "pattern":
            "patterned_B_measured",
    },
]


# ============================================================
# DISPLAY-GENE SELECTION
# ============================================================

def choose_genes(spec):

    valid = (
        res[
            "evaluable_all3"
        ].astype(bool)
        &
        np.isfinite(
            pd.to_numeric(
                res[
                    spec[
                        "cos"
                    ]
                ],
                errors="coerce"
            )
        )
    )


    pat = (
        res.loc[
            valid
            &
            res[
                spec[
                    "pattern"
                ]
            ].astype(bool)
        ]
        .sort_values(
            spec[
                "cos"
            ],
            ascending=False
        )
        .copy()
    )


    non = (
        res.loc[
            valid
            &
            ~res[
                spec[
                    "pattern"
                ]
            ].astype(bool)
        ]
        .sort_values(
            spec[
                "cos"
            ],
            ascending=False
        )
        .copy()
    )


    if len(pat) < 3:
        raise RuntimeError(
            "Too few patterned genes."
        )

    high = pat.iloc[0]

    median_value = np.nanmedian(
        pat[
            spec[
                "cos"
            ]
        ]
    )

    med_idx = (
        pat[
            spec[
                "cos"
            ]
        ]
        .sub(
            median_value
        )
        .abs()
        .idxmin()
    )

    median_gene = pat.loc[
        med_idx
    ]

    low = pat.iloc[-1]


    non_median = np.nanmedian(
        non[
            spec[
                "cos"
            ]
        ]
    )

    non_idx = (
        non[
            spec[
                "cos"
            ]
        ]
        .sub(
            non_median
        )
        .abs()
        .idxmin()
    )

    control = non.loc[
        non_idx
    ]


    return [
        (
            high,
            "patterned, highest cosine"
        ),
        (
            median_gene,
            "patterned, median cosine"
        ),
        (
            low,
            "patterned, lowest cosine"
        ),
        (
            control,
            "not patterned control"
        ),
    ]


# ============================================================
# MAP FUNCTIONS
# ============================================================

def spot_map(
    ax,
    gene,
    M
):

    j = gene_to_idx[
        gene
    ]

    values = M[
        :,
        j
    ]

    vmax = np.nanpercentile(
        values,
        99
    )

    ax.scatter(
        xy_spot[:, 0],
        xy_spot[:, 1],
        c=values,
        s=7,
        linewidths=0,
        rasterized=True,
        vmin=0,
        vmax=vmax
    )

    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.axis(
        "off"
    )


def cell_map(
    ax,
    gene,
    source
):

    j = gene_to_idx[
        gene
    ]

    if source == "measured":

        values = Bcell_log[
            cell_rows,
            j
        ]

    elif source == "oof":

        values = Ccell_log[
            cell_rows,
            j
        ]

    else:

        raise ValueError(
            source
        )


    vmax = np.nanpercentile(
        values,
        99
    )

    ax.scatter(
        cell_xy[:, 0],
        cell_xy[:, 1],
        c=values,
        s=1.0,
        linewidths=0,
        rasterized=True,
        vmin=0,
        vmax=vmax
    )

    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.axis(
        "off"
    )


# ============================================================
# PLOT
# ============================================================

for spec in PAIRS:

    selected = choose_genes(
        spec
    )


    print()
    print("=" * 100)
    print(
        spec[
            "name1"
        ],
        "vs",
        spec[
            "name2"
        ]
    )
    print("=" * 100)


    for r, label in selected:

        print(
            label,
            r[
                "gene"
            ],
            "cos=",
            f"{r[spec['cos']]:.3f}",
            "sp=",
            f"{r[spec['sp']]:+.3f}"
        )


    # ========================================================
    # SAME TEMPLATE AS PREVIOUS PI FIGURE
    # ========================================================

    fig = plt.figure(
        figsize=(
            14.5,
            14.0
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
        left=0.16,
        right=0.97,
        top=0.90,
        bottom=0.10,
        wspace=0.38,
        hspace=0.30
    )


    # ========================================================
    # a/b/c — FOUR GENE ROWS
    # ========================================================

    for row, (
        r,
        row_label
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


        # ----------------------------------------------------
        # PANEL a
        # ----------------------------------------------------

        axa = fig.add_subplot(
            gs[
                row,
                0
            ]
        )

        if spec[
            "map1"
        ] == "visium":

            spot_map(
                axa,
                gene,
                spec[
                    "spot1"
                ]
            )

        else:

            cell_map(
                axa,
                gene,
                spec[
                    "map1"
                ]
            )


        if row == 0:

            axa.set_title(
                spec[
                    "name1"
                ],
                loc="left",
                fontsize=11,
                fontweight="bold"
            )


        axa.text(
            -0.10,
            0.50,
            f"{gene}\n{row_label}",
            transform=axa.transAxes,
            ha="right",
            va="center",
            fontsize=9
        )


        # ----------------------------------------------------
        # PANEL b
        # ----------------------------------------------------

        axb = fig.add_subplot(
            gs[
                row,
                1
            ]
        )

        if spec[
            "map2"
        ] == "visium":

            spot_map(
                axb,
                gene,
                spec[
                    "spot2"
                ]
            )

        else:

            cell_map(
                axb,
                gene,
                spec[
                    "map2"
                ]
            )


        if row == 0:

            axb.set_title(
                spec[
                    "name2"
                ],
                loc="left",
                fontsize=11,
                fontweight="bold"
            )


        # ----------------------------------------------------
        # PANEL c
        # ----------------------------------------------------

        axc = fig.add_subplot(
            gs[
                row,
                2
            ]
        )


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


        axc.scatter(
            x,
            y,
            s=7,
            alpha=0.35,
            linewidths=0,
            rasterized=True
        )


        if row == 0:

            axc.set_title(
                "Normalized expression\ncorrespondence",
                loc="left",
                fontsize=11,
                fontweight="bold"
            )


        axc.set_xlabel(
            f"{spec['name1']} log10(CPM + 1)",
            fontsize=8
        )

        axc.set_ylabel(
            f"{spec['name2']} log10(CPM + 1)",
            fontsize=8
        )


        axc.text(
            0.03,
            0.97,
            (
                f"raw cosine "
                f"{r[spec['cos']]:.2f}\n"
                f"Spearman "
                f"{r[spec['sp']]:+.2f}"
            ),
            transform=axc.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            bbox=dict(
                facecolor="white",
                alpha=0.7,
                edgecolor="0.7",
                linewidth=0.5
            )
        )


    # ========================================================
    # PANEL d — MORAN'S I A vs B
    # SAME RIGHT-SIDE LARGE PANEL POSITION AS OLD FIGURE
    # ========================================================

    axd = fig.add_subplot(
        gs[
            0:3,
            3
        ]
    )


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


    valid = (
        res[
            "evaluable_all3"
        ].astype(bool).to_numpy()
        &
        np.isfinite(
            mx
        )
        &
        np.isfinite(
            my
        )
    )


    pattern = (
        res[
            spec[
                "pattern"
            ]
        ]
        .astype(bool)
        .to_numpy()
    )


    # non-patterned first
    axd.scatter(
        mx[
            valid
            &
            ~pattern
        ],
        my[
            valid
            &
            ~pattern
        ],
        s=16,
        alpha=0.35,
        linewidths=0,
        label="not patterned"
    )


    # patterned
    axd.scatter(
        mx[
            valid
            &
            pattern
        ],
        my[
            valid
            &
            pattern
        ],
        s=18,
        alpha=0.65,
        linewidths=0,
        label="patterned"
    )


    xv = mx[
        valid
    ]

    yv = my[
        valid
    ]


    rho = spearmanr(
        xv,
        yv
    ).statistic


    pear = pearsonr(
        xv,
        yv
    ).statistic


    slope, intercept = np.polyfit(
        xv,
        yv,
        1
    )


    lo = min(
        np.nanmin(
            xv
        ),
        np.nanmin(
            yv
        )
    )

    hi = max(
        np.nanmax(
            xv
        ),
        np.nanmax(
            yv
        )
    )

    pad = max(
        0.02,
        (
            hi - lo
        )
        *
        0.06
    )

    lo -= pad
    hi += pad


    axd.plot(
        [
            lo,
            hi
        ],
        [
            lo,
            hi
        ],
        "--",
        lw=1
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
        slope * xx,
        lw=1.2
    )


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
        f"{spec['name1']} Moran's I"
    )

    axd.set_ylabel(
        f"{spec['name2']} Moran's I"
    )


    axd.set_title(
        "d   Moran's-I correspondence",
        loc="left",
        fontsize=11,
        fontweight="bold"
    )


    axd.text(
        0.03,
        0.97,
        (
            f"genes n={valid.sum()}\n"
            f"Spearman ρ={rho:.3f}\n"
            f"Pearson r={pear:.3f}\n"
            f"slope={slope:.3f}\n"
            f"median ΔI="
            f"{np.nanmedian(yv-xv):+.3f}"
        ),
        transform=axd.transAxes,
        ha="left",
        va="top",
        fontsize=8
    )


    axd.legend(
        frameon=False,
        fontsize=7,
        loc="lower right"
    )


    # ========================================================
    # PANEL e — RAW COSINE DISTRIBUTIONS
    #
    # Same bottom-right location as old "not patterned" panel.
    # Both patterned/non-patterned shown here because panel d
    # is now reserved for Moran's I.
    # ========================================================

    axe = fig.add_subplot(
        gs[
            3,
            3
        ]
    )


    patterned_cos = (
        res.loc[
            res[
                "evaluable_all3"
            ].astype(bool)
            &
            res[
                spec[
                    "pattern"
                ]
            ].astype(bool),
            spec[
                "cos"
            ]
        ]
        .dropna()
        .to_numpy(float)
    )


    nonpatterned_cos = (
        res.loc[
            res[
                "evaluable_all3"
            ].astype(bool)
            &
            ~res[
                spec[
                    "pattern"
                ]
            ].astype(bool),
            spec[
                "cos"
            ]
        ]
        .dropna()
        .to_numpy(float)
    )


    bins = np.linspace(
        0,
        1,
        21
    )


    # horizontal histograms; same orientation as prior template
    hp, edges = np.histogram(
        patterned_cos,
        bins=bins
    )

    hn, _ = np.histogram(
        nonpatterned_cos,
        bins=bins
    )


    centers = (
        edges[:-1]
        +
        edges[1:]
    ) / 2


    height = (
        edges[1]
        -
        edges[0]
    ) * 0.40


    axe.barh(
        centers
        +
        height / 2,
        hp,
        height=height,
        alpha=0.65,
        label=(
            f"patterned "
            f"(n={len(patterned_cos)})"
        )
    )


    axe.barh(
        centers
        -
        height / 2,
        hn,
        height=height,
        alpha=0.45,
        label=(
            f"not patterned "
            f"(n={len(nonpatterned_cos)})"
        )
    )


    med_pat = np.nanmedian(
        patterned_cos
    )

    med_non = np.nanmedian(
        nonpatterned_cos
    )


    axe.axhline(
        med_pat,
        lw=1.0,
        ls="--"
    )

    axe.axhline(
        med_non,
        lw=1.0,
        ls=":"
    )


    axe.set_ylim(
        0,
        1
    )

    axe.set_xlabel(
        "gene count"
    )

    axe.set_ylabel(
        "raw cosine similarity"
    )


    axe.set_title(
        "e   raw cosine distribution",
        loc="left",
        fontsize=11,
        fontweight="bold"
    )


    axe.text(
        0.03,
        0.97,
        (
            f"patterned median "
            f"{med_pat:.2f}\n"
            f"non-patterned median "
            f"{med_non:.2f}"
        ),
        transform=axe.transAxes,
        ha="left",
        va="top",
        fontsize=7
    )


    axe.legend(
        frameon=False,
        fontsize=6,
        loc="lower right"
    )


    # ========================================================
    # PANEL LETTERS
    # ========================================================

    fig.text(
        0.155,
        0.915,
        "a",
        fontsize=12,
        fontweight="bold"
    )

    fig.text(
        0.367,
        0.915,
        "b",
        fontsize=12,
        fontweight="bold"
    )

    fig.text(
        0.577,
        0.915,
        "c",
        fontsize=12,
        fontweight="bold"
    )


    # ========================================================
    # TITLE
    # ========================================================

    fig.suptitle(
        (
            f"External spatial validation: "
            f"{spec['name1']} vs {spec['name2']} "
            f"({BR})"
        ),
        fontsize=14,
        fontweight="bold",
        y=0.965
    )


    # ========================================================
    # CAPTION — SAME BOTTOM STYLE
    # ========================================================

    pattern_source = (
        "Visium"
        if spec[
            "pattern"
        ] == "patterned_A_Visium"
        else
        "measured Xenium"
    )


    caption = (
        f"Panels a–c show four representative genes: "
        f"high-, median-, and low-raw-cosine spatially patterned genes "
        f"plus one non-patterned control. "
        f"Expression correspondence is calculated across the same "
        f"{len(xy_spot):,} occupied Visium spots after CPM normalization. "
        f"Panel d compares Moran's I for the same 300 genes using the "
        f"same 6-nearest-neighbor spatial graph for both sources. "
        f"Spatial patterning was defined independently from "
        f"{pattern_source} Moran's I with BH q < 0.05 and I > 0. "
        f"Panel e shows raw cosine distributions for patterned and "
        f"non-patterned genes."
    )


    fig.text(
        0.50,
        0.035,
        caption,
        ha="center",
        va="bottom",
        fontsize=7,
        wrap=True
    )


    # ========================================================
    # SAVE
    # ========================================================

    PNG = (
        OUT
        / f"fig_{BR}_{spec['code']}_PI_template_raw_cosine.png"
    )

    PDF = (
        OUT
        / f"fig_{BR}_{spec['code']}_PI_template_raw_cosine.pdf"
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
        "Saved:",
        PNG
    )


o.file.close()


print()
print("=" * 100)
print("ALL THREE PI-TEMPLATE FIGURES COMPLETE")
print("=" * 100)
