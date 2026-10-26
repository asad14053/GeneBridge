#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    / "Br8667_three_source_gene_metrics.csv"
)

SPOTS = (
    BASE
    / "Br8667_three_source_spot_matrices.npz"
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
    / "PI_style_pairwise"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# LOAD GENE METRICS + SPOT MATRICES
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

A = np.asarray(
    z["visium_cpm"],
    dtype=float
)

B = np.asarray(
    z["measured_xenium_cpm"],
    dtype=float
)

C = np.asarray(
    z["oof_envi_cpm"],
    dtype=float
)


# plotting expression
Alog = np.log10(
    A + 1
)

Blog = np.log10(
    B + 1
)

Clog = np.log10(
    C + 1
)


gene_to_idx = {
    g: i
    for i, g
    in enumerate(genes)
}


# ============================================================
# LOAD CELL-LEVEL B/C FOR SPATIAL MAPS
# ============================================================

o = ad.read_h5ad(
    OOF
)

aa = pd.read_csv(
    ALIGN,
    low_memory=False
)


# map OOF original cell ID -> STalign cell ID
lookup = {
    str(x): i
    for i, x
    in enumerate(
        aa["cell_id"].astype(str)
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


has = (
    st_row >= 0
)

idx = np.flatnonzero(
    has
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
    st_row[idx]
]


match_prob[idx] = pd.to_numeric(
    sub["match_prob"],
    errors="coerce"
).to_numpy()


vis_spot[idx] = (
    sub["visium_spot"]
    .fillna("")
    .astype(str)
    .to_numpy()
)


x_aligned[idx] = pd.to_numeric(
    sub["x_aligned"],
    errors="coerce"
).to_numpy()


y_aligned[idx] = pd.to_numeric(
    sub["y_aligned"],
    errors="coerce"
).to_numpy()


use = (
    has
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


# B = measured ground truth
# C = ENVI OOF prediction
# AnnData layers may be scipy sparse matrices, so explicitly
# materialize them before converting to numpy.

def dense_float(X):
    if hasattr(X, "toarray"):
        X = X.toarray()

    return np.asarray(
        X,
        dtype=np.float32
    )


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

print(
    "Cell-level B measured:",
    Bcell.shape,
    Bcell.dtype
)

print(
    "Cell-level C OOF:",
    Ccell.shape,
    Ccell.dtype
)


# cell-level CPM over same 300 genes
def cell_cpm(X):

    lib = X.sum(
        axis=1
    )

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


cell_xy = np.column_stack([
    x_aligned[
        cell_rows
    ],
    y_aligned[
        cell_rows
    ]
])


# ============================================================
# PAIR DEFINITIONS
# ============================================================

pairs = [
    {
        "code": "A_vs_B",
        "source1": "Visium",
        "source2": "Measured Xenium",
        "X": Alog,
        "Y": Blog,
        "cos_col":
            "cosine_A_Visium_vs_B_measured",
        "sp_col":
            "spearman_A_Visium_vs_B_measured",
        "moran_x":
            "moran_A_Visium",
        "moran_y":
            "moran_B_measured",
        "pattern":
            "patterned_A_Visium",
        "map1":
            "spot_A",
        "map2":
            "cell_B",
    },

    {
        "code": "A_vs_C",
        "source1": "Visium",
        "source2": "ENVI OOF",
        "X": Alog,
        "Y": Clog,
        "cos_col":
            "cosine_A_Visium_vs_C_OOF",
        "sp_col":
            "spearman_A_Visium_vs_C_OOF",
        "moran_x":
            "moran_A_Visium",
        "moran_y":
            "moran_C_OOF",
        "pattern":
            "patterned_A_Visium",
        "map1":
            "spot_A",
        "map2":
            "cell_C",
    },

    {
        "code": "B_vs_C",
        "source1": "Measured Xenium",
        "source2": "ENVI OOF",
        "X": Blog,
        "Y": Clog,
        "cos_col":
            "cosine_B_measured_vs_C_OOF",
        "sp_col":
            "spearman_B_measured_vs_C_OOF",
        "moran_x":
            "moran_B_measured",
        "moran_y":
            "moran_C_OOF",
        "pattern":
            "patterned_B_measured",
        "map1":
            "cell_B",
        "map2":
            "cell_C",
    },
]


# ============================================================
# HELPERS
# ============================================================

def choose_genes(
    df,
    cosine_col,
    pattern_col
):

    ok = (
        df["evaluable_all3"].astype(bool)
        &
        np.isfinite(
            df[
                cosine_col
            ]
        )
    )

    patterned = (
        df.loc[
            ok
            &
            df[
                pattern_col
            ].astype(bool)
        ]
        .sort_values(
            cosine_col,
            ascending=False
        )
        .copy()
    )

    nonpatterned = (
        df.loc[
            ok
            &
            ~df[
                pattern_col
            ].astype(bool)
        ]
        .sort_values(
            cosine_col,
            ascending=False
        )
        .copy()
    )

    if len(patterned) < 3:
        raise RuntimeError(
            "Too few patterned genes."
        )

    high = patterned.iloc[0]

    med = patterned.iloc[
        len(patterned) // 2
    ]

    low = patterned.iloc[-1]

    if len(nonpatterned):

        ctrl = nonpatterned.iloc[
            len(nonpatterned) // 2
        ]

    else:

        raise RuntimeError(
            "No non-patterned control genes."
        )

    return [
        (
            high,
            "patterned, highest cosine"
        ),
        (
            med,
            "patterned, median cosine"
        ),
        (
            low,
            "patterned, lowest cosine"
        ),
        (
            ctrl,
            "non-patterned control"
        ),
    ]


def plot_spot_map(
    ax,
    gene,
    M,
    title
):

    j = gene_to_idx[
        gene
    ]

    val = M[
        :,
        j
    ]

    ax.scatter(
        xy_spot[:, 0],
        xy_spot[:, 1],
        c=val,
        s=8,
        linewidths=0,
        rasterized=True
    )

    ax.set_title(
        title,
        fontsize=9
    )

    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.axis(
        "off"
    )


def plot_cell_map(
    ax,
    gene,
    source,
    title
):

    j = gene_to_idx[
        gene
    ]

    if source == "B":

        val = Bcell_log[
            cell_rows,
            j
        ]

    else:

        val = Ccell_log[
            cell_rows,
            j
        ]

    # robust display scaling only
    vmax = np.nanpercentile(
        val,
        99
    )

    ax.scatter(
        cell_xy[:, 0],
        cell_xy[:, 1],
        c=val,
        s=1.0,
        linewidths=0,
        rasterized=True,
        vmin=0,
        vmax=vmax
    )

    ax.set_title(
        title,
        fontsize=9
    )

    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.axis(
        "off"
    )


# ============================================================
# GENERATE ONE FIGURE PER PAIR
# ============================================================

for spec in pairs:

    selected = choose_genes(
        res,
        spec[
            "cos_col"
        ],
        spec[
            "pattern"
        ]
    )

    print()
    print("=" * 100)
    print(
        spec["code"],
        spec["source1"],
        "vs",
        spec["source2"]
    )
    print("=" * 100)

    for r, label in selected:

        print(
            label,
            ":",
            r["gene"],
            f"cos={r[spec['cos_col']]:.3f}",
            f"sp={r[spec['sp_col']]:+.3f}"
        )


    # --------------------------------------------------------
    # Layout
    #
    # cols 0-2 = four selected gene rows
    # col 3    = Moran scatter, spans all rows
    # col 4    = patterned raw-cosine distribution
    # col 5    = non-patterned raw-cosine distribution
    # --------------------------------------------------------

    fig = plt.figure(
        figsize=(
            20,
            11
        )
    )

    gs = fig.add_gridspec(
        4,
        6,
        width_ratios=[
            1.0,
            1.0,
            1.15,
            1.35,
            0.9,
            0.9,
        ],
        wspace=0.42,
        hspace=0.38
    )


    # ========================================================
    # a / b / c — selected genes
    # ========================================================

    for i, (
        r,
        label
    ) in enumerate(
        selected
    ):

        gene = str(
            r["gene"]
        )

        j = gene_to_idx[
            gene
        ]


        # ----------------------------------------------------
        # panel a
        # ----------------------------------------------------

        ax = fig.add_subplot(
            gs[
                i,
                0
            ]
        )

        if spec[
            "map1"
        ] == "spot_A":

            plot_spot_map(
                ax,
                gene,
                Alog,
                (
                    spec[
                        "source1"
                    ]
                    if i == 0
                    else ""
                )
            )

        elif spec[
            "map1"
        ] == "cell_B":

            plot_cell_map(
                ax,
                gene,
                "B",
                (
                    spec[
                        "source1"
                    ]
                    if i == 0
                    else ""
                )
            )


        ax.text(
            -0.03,
            0.5,
            f"{gene}\n{label}",
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=9
        )


        # ----------------------------------------------------
        # panel b
        # ----------------------------------------------------

        ax = fig.add_subplot(
            gs[
                i,
                1
            ]
        )

        if spec[
            "map2"
        ] == "cell_B":

            plot_cell_map(
                ax,
                gene,
                "B",
                (
                    spec[
                        "source2"
                    ]
                    if i == 0
                    else ""
                )
            )

        elif spec[
            "map2"
        ] == "cell_C":

            plot_cell_map(
                ax,
                gene,
                "C",
                (
                    spec[
                        "source2"
                    ]
                    if i == 0
                    else ""
                )
            )


        # ----------------------------------------------------
        # panel c
        # ----------------------------------------------------

        ax = fig.add_subplot(
            gs[
                i,
                2
            ]
        )

        x = spec[
            "X"
        ][
            :,
            j
        ]

        y = spec[
            "Y"
        ][
            :,
            j
        ]

        ax.scatter(
            x,
            y,
            s=8,
            alpha=0.35,
            linewidths=0,
            rasterized=True
        )

        ax.set_xlabel(
            spec[
                "source1"
            ]
            +
            "\nlog10(CPM + 1)",
            fontsize=8
        )

        ax.set_ylabel(
            spec[
                "source2"
            ]
            +
            "\nlog10(CPM + 1)",
            fontsize=8
        )

        if i == 0:

            ax.set_title(
                "normalized expression correspondence",
                fontsize=10
            )


        ax.text(
            0.03,
            0.97,
            (
                f"raw cosine "
                f"{r[spec['cos_col']]:.2f}\n"
                f"Spearman "
                f"{r[spec['sp_col']]:+.2f}"
            ),
            transform=ax.transAxes,
            va="top",
            fontsize=8
        )


    # ========================================================
    # d — MORAN I GENE-LEVEL CORRESPONDENCE
    # ========================================================

    axd = fig.add_subplot(
        gs[
            :,
            3
        ]
    )

    mx = pd.to_numeric(
        res[
            spec[
                "moran_x"
            ]
        ],
        errors="coerce"
    ).to_numpy(float)

    my = pd.to_numeric(
        res[
            spec[
                "moran_y"
            ]
        ],
        errors="coerce"
    ).to_numpy(float)


    mok = (
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


    pat = (
        res[
            spec[
                "pattern"
            ]
        ]
        .astype(bool)
        .to_numpy()
    )


    # non-patterned
    axd.scatter(
        mx[
            mok
            &
            ~pat
        ],
        my[
            mok
            &
            ~pat
        ],
        s=18,
        alpha=0.28,
        linewidths=0,
        label="not patterned"
    )


    # patterned
    axd.scatter(
        mx[
            mok
            &
            pat
        ],
        my[
            mok
            &
            pat
        ],
        s=22,
        alpha=0.65,
        linewidths=0,
        label="patterned"
    )


    xv = mx[
        mok
    ]

    yv = my[
        mok
    ]


    lo = float(
        min(
            np.nanmin(xv),
            np.nanmin(yv)
        )
    )

    hi = float(
        max(
            np.nanmax(xv),
            np.nanmax(yv)
        )
    )

    pad = max(
        0.02,
        0.05
        *
        (
            hi - lo
        )
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
        linestyle="--",
        linewidth=1
    )


    rho, _ = spearmanr(
        xv,
        yv
    )

    pr, _ = pearsonr(
        xv,
        yv
    )

    slope, intercept = np.polyfit(
        xv,
        yv,
        1
    )


    xx = np.linspace(
        lo,
        hi,
        100
    )

    axd.plot(
        xx,
        intercept
        +
        slope
        *
        xx,
        linewidth=1.2
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
        f"{spec['source1']} Moran's I"
    )

    axd.set_ylabel(
        f"{spec['source2']} Moran's I"
    )

    axd.set_title(
        "d   Moran's-I correspondence",
        loc="left",
        fontweight="bold"
    )


    axd.text(
        0.03,
        0.97,
        (
            f"genes n={len(xv)}\n"
            f"Spearman ρ={rho:.3f}\n"
            f"Pearson r={pr:.3f}\n"
            f"slope={slope:.3f}\n"
            f"median ΔI="
            f"{np.median(yv-xv):+.3f}"
        ),
        transform=axd.transAxes,
        va="top",
        fontsize=9
    )


    axd.legend(
        frameon=False,
        fontsize=8,
        loc="lower right"
    )


    # ========================================================
    # e — PATTERNED RAW COSINE
    # ========================================================

    axe = fig.add_subplot(
        gs[
            :2,
            4
        ]
    )

    pvals = (
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
                "cos_col"
            ]
        ]
        .dropna()
        .to_numpy(float)
    )


    axe.hist(
        pvals,
        bins=np.linspace(
            0,
            1,
            26
        ),
        orientation="horizontal"
    )

    medp = np.median(
        pvals
    )

    axe.axhline(
        medp,
        linestyle="--",
        linewidth=1
    )

    axe.set_ylim(
        0,
        1
    )

    axe.set_xlabel(
        "genes"
    )

    axe.set_ylabel(
        "raw cosine similarity"
    )

    axe.set_title(
        "e   spatially patterned",
        loc="left",
        fontweight="bold"
    )

    axe.text(
        0.03,
        0.95,
        (
            f"n={len(pvals)}\n"
            f"median={medp:.3f}"
        ),
        transform=axe.transAxes,
        va="top",
        fontsize=9
    )


    # selected patterned genes
    for r, label in selected[:3]:

        c = float(
            r[
                spec[
                    "cos_col"
                ]
            ]
        )

        axe.axhline(
            c,
            linewidth=0.8,
            alpha=0.55
        )


    # ========================================================
    # f — NON-PATTERNED RAW COSINE
    # ========================================================

    axf = fig.add_subplot(
        gs[
            2:,
            4
        ]
    )

    npvals = (
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
                "cos_col"
            ]
        ]
        .dropna()
        .to_numpy(float)
    )


    axf.hist(
        npvals,
        bins=np.linspace(
            0,
            1,
            26
        ),
        orientation="horizontal"
    )

    medn = np.median(
        npvals
    )

    axf.axhline(
        medn,
        linestyle="--",
        linewidth=1
    )

    axf.set_ylim(
        0,
        1
    )

    axf.set_xlabel(
        "genes"
    )

    axf.set_ylabel(
        "raw cosine similarity"
    )

    axf.set_title(
        "f   non-patterned",
        loc="left",
        fontweight="bold"
    )

    axf.text(
        0.03,
        0.95,
        (
            f"n={len(npvals)}\n"
            f"median={medn:.3f}"
        ),
        transform=axf.transAxes,
        va="top",
        fontsize=9
    )


    # non-patterned selected control
    r = selected[3][0]

    axf.axhline(
        float(
            r[
                spec[
                    "cos_col"
                ]
            ]
        ),
        linewidth=0.8,
        alpha=0.55
    )


    # ========================================================
    # SUMMARY COLUMN
    # ========================================================

    axs = fig.add_subplot(
        gs[
            :,
            5
        ]
    )

    axs.axis(
        "off"
    )

    overall_cos = (
        res.loc[
            res[
                "evaluable_all3"
            ].astype(bool),
            spec[
                "cos_col"
            ]
        ]
        .dropna()
        .to_numpy(float)
    )

    overall_sp = (
        res.loc[
            res[
                "evaluable_all3"
            ].astype(bool),
            spec[
                "sp_col"
            ]
        ]
        .dropna()
        .to_numpy(float)
    )


    summary_text = (
        f"{spec['source1']}\n"
        f"vs\n"
        f"{spec['source2']}\n\n"
        f"Genes: {len(overall_cos)}\n\n"
        f"Median raw cosine\n"
        f"{np.median(overall_cos):.3f}\n\n"
        f"Median Spearman\n"
        f"{np.median(overall_sp):+.3f}\n\n"
        f"Patterned genes\n"
        f"{int(pat[mok].sum())}\n\n"
        f"Pattern definition:\n"
        f"{'Visium' if spec['pattern']=='patterned_A_Visium' else 'measured Xenium'}"
    )


    axs.text(
        0.05,
        0.95,
        summary_text,
        va="top",
        fontsize=11
    )


    # ========================================================
    # PANEL LETTERS
    # ========================================================

    fig.text(
        0.01,
        0.97,
        "a",
        fontsize=15,
        fontweight="bold"
    )

    fig.text(
        0.17,
        0.97,
        "b",
        fontsize=15,
        fontweight="bold"
    )

    fig.text(
        0.33,
        0.97,
        "c",
        fontsize=15,
        fontweight="bold"
    )


    fig.suptitle(
        (
            f"{BR}: "
            f"{spec['source1']} vs {spec['source2']}"
        ),
        fontsize=15,
        fontweight="bold",
        y=0.995
    )


    # ========================================================
    # SAVE
    # ========================================================

    PNG = (
        OUT
        / (
            f"fig_{BR}_"
            f"{spec['code']}_"
            f"PI_style_raw_cosine.png"
        )
    )

    PDF = (
        OUT
        / (
            f"fig_{BR}_"
            f"{spec['code']}_"
            f"PI_style_raw_cosine.pdf"
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
        "Saved:",
        PNG
    )


o.file.close()

print()
print("=" * 100)
print("ALL THREE FIGURES COMPLETE")
print("=" * 100)
