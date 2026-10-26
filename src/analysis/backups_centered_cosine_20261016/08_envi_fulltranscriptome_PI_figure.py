#!/usr/bin/env python3

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.io import mmread


# ============================================================
# CONFIG
# ============================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BR = os.environ.get(
    "BR",
    "Br8667"
)

MATCH_THRESH = 0.85

ENVI_FILE = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/envi"
    / BR
    / f"spatial_data_xenium_{BR}_ENVI_full_transcriptome.h5ad"
)

STDIR = (
    ROOT
    / "data/metadata/stalign_results_update"
    / f"{BR}_STalign_result"
)

ALIGN_FILE = (
    STDIR
    / f"{BR}_aligned_coordinates.csv"
)

VCOORD_FILE = (
    STDIR
    / f"{BR}_visium_coordinates.csv"
)

METRICS_FILE = (
    ROOT
    / "results/analysis/envi_fulltranscriptome_crosstech"
    / BR
    / f"{BR}_ENVI_fulltranscriptome_gene_metrics.csv.gz"
)

VRAW = (
    ROOT
    / "data/raw/visium"
    / BR
)

OUT = (
    ROOT
    / "results/analysis/envi_fulltranscriptome_crosstech"
    / BR
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def pick_one(path, pattern):

    hits = sorted(
        path.glob(pattern)
    )

    if len(hits) != 1:
        raise RuntimeError(
            f"{BR}: expected one {pattern} under {path}, "
            f"found {len(hits)}:\n"
            +
            "\n".join(
                str(x)
                for x in hits
            )
        )

    return hits[0]


def dense(x):

    if sparse.issparse(x):
        return x.toarray()

    return np.asarray(x)


def n01(v):

    v = np.asarray(
        v,
        dtype=float
    ).ravel()

    finite = np.isfinite(v)

    out = np.zeros(
        len(v),
        dtype=float
    )

    if finite.sum() == 0:
        return out

    a, b = np.percentile(
        v[finite],
        [1, 99]
    )

    out[finite] = np.clip(
        (
            v[finite]
            -
            a
        )
        /
        (
            b
            -
            a
            +
            1e-9
        ),
        0,
        1
    )

    return out


def as_bool(x):

    if pd.api.types.is_bool_dtype(x):
        return x.to_numpy(bool)

    return (
        x.astype(str)
        .str.strip()
        .str.lower()
        .isin([
            "true",
            "1",
            "t",
            "yes"
        ])
        .to_numpy(bool)
    )


# ============================================================
# CHECK FILES
# ============================================================

for f in [
    ENVI_FILE,
    ALIGN_FILE,
    VCOORD_FILE,
    METRICS_FILE,
]:

    if not f.exists():
        raise FileNotFoundError(f)


VMTX = pick_one(
    VRAW,
    "*_matrix.mtx.gz"
)

VBAR = pick_one(
    VRAW,
    "*_barcodes.tsv.gz"
)

VFEAT = pick_one(
    VRAW,
    "*_features.tsv.gz"
)


# ============================================================
# LOAD GENE METRICS
# ============================================================

res = pd.read_csv(
    METRICS_FILE
)

for c in [
    "measured_in_xenium",
    "evaluable",
    "visium_highconf",
    "spatially_patterned",
]:

    res[c] = as_bool(
        res[c]
    )


# ============================================================
# IMPUTED-ONLY GENE UNIVERSE
# ============================================================

imp = (
    (~res["measured_in_xenium"])
    &
    res["evaluable"]
)

imp_pat = (
    imp
    &
    res["spatially_patterned"]
)

imp_non = (
    imp
    &
    (~res["spatially_patterned"])
)


# Display examples require stronger Visium expression.
display_ok = (
    imp
    &
    res["visium_highconf"]
)


patd = (
    res[
        display_ok
        &
        res["spatially_patterned"]
    ]
    .sort_values(
        "cosine",
        ascending=False
    )
    .reset_index(drop=True)
)


nonp = (
    res[
        display_ok
        &
        (~res["spatially_patterned"])
    ]
    .sort_values(
        "cosine",
        ascending=False
    )
    .reset_index(drop=True)
)


if len(patd) < 3:
    raise RuntimeError(
        f"{BR}: not enough depth-qualified patterned "
        "imputed genes for plotting."
    )

if len(nonp) < 1:
    raise RuntimeError(
        f"{BR}: no depth-qualified non-patterned "
        "imputed control available."
    )


G = [
    patd.iloc[0]["gene"],
    patd.iloc[len(patd) // 2]["gene"],
    patd.iloc[-1]["gene"],
    nonp.iloc[len(nonp) // 2]["gene"],
]


GLAB = [
    "imputed, patterned, highest cosine",
    "imputed, patterned, median cosine",
    "imputed, patterned, lowest cosine",
    "imputed, not patterned (control)",
]


print("=" * 80)
print(f"ENVI PI-STYLE FIGURE: {BR}")
print("=" * 80)

print(
    "imputed evaluable:",
    int(imp.sum())
)

print(
    "imputed patterned:",
    int(imp_pat.sum())
)

print(
    "imputed non-patterned:",
    int(imp_non.sum())
)

print()
print("DISPLAY GENES")

for nm, g in zip(
    GLAB,
    G
):

    r = (
        res
        .set_index("gene")
        .loc[g]
    )

    print(
        f"{nm:40s} "
        f"{g:15s} "
        f"cos={r.cosine:.3f} "
        f"sp={r.spearman:+.3f} "
        f"Moran={r.moran_I:.3f} "
        f"q={r.moran_q:.2g}"
    )


# ============================================================
# LOAD ENVI + STALIGN
# ============================================================

e = ad.read_h5ad(
    ENVI_FILE,
    backed="r"
)

aa = pd.read_csv(
    ALIGN_FILE,
    low_memory=False
)

vcoord = pd.read_csv(
    VCOORD_FILE
)


# ============================================================
# MAP ENVI CELLS TO STALIGN
# ============================================================

if (
    "cell_id" not in e.obs.columns
    or
    "cell_id" not in aa.columns
):
    raise RuntimeError(
        "cell_id missing from ENVI or STalign."
    )


lookup = {
    str(k): i
    for i, k
    in enumerate(
        aa["cell_id"]
        .astype(str)
    )
}


st_row = np.array([
    lookup.get(
        str(k),
        -1
    )
    for k in e.obs[
        "cell_id"
    ].astype(str)
])


has_st = (
    st_row >= 0
)


idx = np.flatnonzero(
    has_st
)


sub = aa.iloc[
    st_row[idx]
]


match_prob = np.full(
    e.n_obs,
    np.nan
)

x_aligned = np.full(
    e.n_obs,
    np.nan
)

y_aligned = np.full(
    e.n_obs,
    np.nan
)

vis_spot = np.full(
    e.n_obs,
    "",
    dtype=object
)


match_prob[idx] = pd.to_numeric(
    sub["match_prob"],
    errors="coerce"
).to_numpy()


x_aligned[idx] = pd.to_numeric(
    sub["x_aligned"],
    errors="coerce"
).to_numpy()


y_aligned[idx] = pd.to_numeric(
    sub["y_aligned"],
    errors="coerce"
).to_numpy()


vis_spot[idx] = (
    sub["visium_spot"]
    .fillna("")
    .astype(str)
    .to_numpy()
)


# ============================================================
# PANEL-b CELL UNIVERSE:
# all high-confidence aligned Xenium cells
# ============================================================

matched = (
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


matched_rows = np.flatnonzero(
    matched
)


# ============================================================
# PSEUDOSPOT UNIVERSE:
# matched cells whose centroids fall inside a physical Visium spot
# ============================================================

assigned = (
    matched
    &
    (
        vis_spot != ""
    )
    &
    (
        vis_spot != "nan"
    )
)


assigned_rows = np.flatnonzero(
    assigned
)


print()
print(
    "matched ENVI cells:",
    len(matched_rows)
)

print(
    "inside Visium spots:",
    len(assigned_rows)
)


# ============================================================
# TARGET VISIUM SPOTS
# ============================================================

target_barcodes = (
    vcoord["barcode"]
    .astype(str)
    .tolist()
)


barcode_to_spot = {
    b: i
    for i, b
    in enumerate(
        target_barcodes
    )
}


spot_idx = np.array([
    barcode_to_spot.get(
        str(x),
        -1
    )
    for x in vis_spot[
        assigned_rows
    ]
])


if (
    spot_idx < 0
).any():

    raise RuntimeError(
        f"{BR}: some exported visium_spot values "
        "are absent from vcoord."
    )


n_spots = len(
    target_barcodes
)


ncell = np.bincount(
    spot_idx,
    minlength=n_spots
)


occupied = (
    ncell > 0
)


occ_idx = np.flatnonzero(
    occupied
)


print(
    "occupied pseudospots:",
    int(occupied.sum()),
    "/",
    n_spots
)


# ============================================================
# LOAD RAW VISIUM
# ============================================================

vbar = pd.read_csv(
    VBAR,
    sep="\t",
    header=None,
    compression="gzip"
)


vfeat = pd.read_csv(
    VFEAT,
    sep="\t",
    header=None,
    compression="gzip"
)


VM = mmread(
    VMTX
).tocsr()


vbarcodes = (
    vbar.iloc[:, 0]
    .astype(str)
    .to_numpy()
)


vgenes = (
    vfeat.iloc[:, 1]
    .astype(str)
    .to_numpy()
)


if VM.shape == (
    len(vgenes),
    len(vbarcodes)
):

    VM = VM.T.tocsr()


elif VM.shape != (
    len(vbarcodes),
    len(vgenes)
):

    raise RuntimeError(
        f"Unexpected Visium matrix shape: "
        f"{VM.shape}"
    )


raw_bar_lookup = {
    b: i
    for i, b
    in enumerate(
        vbarcodes
    )
}


missing = [
    b
    for b in target_barcodes
    if b not in raw_bar_lookup
]


if missing:
    raise RuntimeError(
        f"{BR}: {len(missing)} STalign Visium "
        "barcodes absent from raw matrix."
    )


vrows = np.array([
    raw_bar_lookup[b]
    for b in target_barcodes
])


VT = VM[
    vrows
].tocsr()


# First occurrence of each Visium gene symbol.
vgene_lookup = {}

for j, g in enumerate(
    vgenes
):

    if g not in vgene_lookup:
        vgene_lookup[g] = j


# ============================================================
# GENE INDEXING
# ============================================================

egenes = pd.Index(
    e.var_names.astype(str)
)


eval_genes = (
    res.loc[
        res["evaluable"],
        "gene"
    ]
    .astype(str)
    .tolist()
)


# Keep only evaluable genes actually available in both matrices.
eval_genes = [
    g
    for g in eval_genes
    if (
        g in egenes
        and
        g in vgene_lookup
    )
]


eval_eidx = np.array([
    egenes.get_loc(g)
    for g in eval_genes
])


eval_vidx = np.array([
    vgene_lookup[g]
    for g in eval_genes
])


selected_eidx = np.array([
    egenes.get_loc(g)
    for g in G
])


selected_vidx = np.array([
    vgene_lookup[g]
    for g in G
])


print()
print(
    "evaluable shared genes used for CPM:",
    len(eval_genes)
)


# ============================================================
# LOAD ASSIGNED ENVI MATRIX ONCE
#
# Br8667 is ~15.5K cells x 35K genes (~2 GB float32).
# This keeps the plotting job fast enough.
# ============================================================

print()
print(
    "Loading assigned-cell ENVI count_scale..."
)


E_assigned = e.layers[
    "count_scale"
][
    assigned_rows,
    :
]


E_assigned = np.asarray(
    E_assigned,
    dtype=np.float32
)


print(
    "assigned ENVI matrix:",
    E_assigned.shape,
    f"{E_assigned.nbytes / 1024**3:.2f} GB"
)


# ============================================================
# PSEUDOSPOT LIBRARY SIZE
# over the SAME evaluable shared-gene universe used in Task 2
# ============================================================

print()
print(
    "Computing ENVI pseudospot library sizes..."
)


cell_lib = np.zeros(
    len(assigned_rows),
    dtype=np.float64
)


CHUNK = 1024


for s in range(
    0,
    len(eval_eidx),
    CHUNK
):

    t = min(
        s + CHUNK,
        len(eval_eidx)
    )

    cell_lib += (
        E_assigned[
            :,
            eval_eidx[s:t]
        ]
        .sum(axis=1)
    )


P_lib = np.bincount(
    spot_idx,
    weights=cell_lib,
    minlength=n_spots
).astype(float)


P_lib[
    P_lib == 0
] = 1.0


# ============================================================
# VISIUM LIBRARY SIZE
# same evaluable shared-gene universe
# ============================================================

print(
    "Computing Visium library sizes..."
)


V_lib = np.zeros(
    n_spots,
    dtype=np.float64
)


for s in range(
    0,
    len(eval_vidx),
    CHUNK
):

    t = min(
        s + CHUNK,
        len(eval_vidx)
    )

    Vchunk = VT[
        :,
        eval_vidx[s:t]
    ]

    V_lib += np.asarray(
        Vchunk.sum(axis=1)
    ).ravel()


V_lib[
    V_lib == 0
] = 1.0


# ============================================================
# SELECTED-GENE PSEUDOSPOTS
# ============================================================

E_sel = E_assigned[
    :,
    selected_eidx
].astype(
    np.float64,
    copy=False
)


P_sel = np.zeros(
    (
        n_spots,
        len(G)
    ),
    dtype=np.float64
)


np.add.at(
    P_sel,
    spot_idx,
    E_sel
)


V_sel = (
    VT[
        :,
        selected_vidx
    ]
    .toarray()
    .astype(float)
)


Pcpm = (
    P_sel[
        occupied
    ]
    /
    P_lib[
        occupied,
        None
    ]
    *
    1e6
)


Vcpm = (
    V_sel[
        occupied
    ]
    /
    V_lib[
        occupied,
        None
    ]
    *
    1e6
)


xy_u = (
    vcoord.loc[
        occupied,
        [
            "x",
            "y"
        ]
    ]
    .to_numpy(float)
)


# ============================================================
# ENVI CELL VALUES FOR PANEL b
# all matched cells
# ============================================================

print(
    "Loading selected genes for aligned-cell maps..."
)


E_matched = np.zeros(
    (
        len(matched_rows),
        len(G)
    ),
    dtype=np.float32
)


for j, gi in enumerate(
    selected_eidx
):

    x = e.layers[
        "count_scale"
    ][
        matched_rows,
        gi
    ]

    E_matched[
        :,
        j
    ] = np.asarray(
        x
    ).ravel()


# No longer need 2 GB assigned matrix.
del E_assigned


# ============================================================
# FIGURE
# ============================================================

RED = plt.cm.Reds


nrow = 4


fig = plt.figure(
    figsize=(
        16,
        14.4
    )
)


gs = fig.add_gridspec(
    nrow,
    4,
    width_ratios=[
        1,
        1,
        1.05,
        1.25
    ],
    hspace=0.30,
    wspace=0.36
)


res_idx = (
    res
    .set_index("gene")
)


# ============================================================
# PANELS a / b / c
# ============================================================

for i, (
    g,
    nm
) in enumerate(
    zip(
        G,
        GLAB
    )
):

    r = res_idx.loc[
        g
    ]


    # --------------------------------------------------------
    # a — REAL VISIUM
    # --------------------------------------------------------

    ax = fig.add_subplot(
        gs[
            i,
            0
        ]
    )


    ax.scatter(
        xy_u[:, 0] / 1000,
        xy_u[:, 1] / 1000,
        c=n01(
            Vcpm[
                :,
                i
            ]
        ),
        cmap=RED,
        s=9,
        linewidths=0
    )


    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])


    ax.set_ylabel(
        f"{g}\n{nm}",
        fontsize=9.3,
        rotation=0,
        labelpad=70,
        va="center",
        ha="center"
    )


    if i == 0:

        ax.set_title(
            "a   Visium",
            fontsize=11,
            loc="left",
            fontweight="bold"
        )


    # --------------------------------------------------------
    # b — ALIGNED ENVI-IMPUTED CELL EXPRESSION
    # --------------------------------------------------------

    ax = fig.add_subplot(
        gs[
            i,
            1
        ]
    )


    ax.scatter(
        x_aligned[
            matched_rows
        ] / 1000,
        y_aligned[
            matched_rows
        ] / 1000,
        c=n01(
            E_matched[
                :,
                i
            ]
        ),
        cmap=RED,
        s=0.5,
        linewidths=0
    )


    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])


    if i == 0:

        ax.set_title(
            "b   ENVI imputed aligned",
            fontsize=11,
            loc="left",
            fontweight="bold"
        )


    # --------------------------------------------------------
    # c — PSEUDOSPOT vs VISIUM
    # --------------------------------------------------------

    ax = fig.add_subplot(
        gs[
            i,
            2
        ]
    )


    ax.scatter(
        Pcpm[
            :,
            i
        ],
        Vcpm[
            :,
            i
        ],
        s=9,
        c="#444444",
        alpha=0.45,
        linewidths=0
    )


    ax.set_xlabel(
        "ENVI pseudospot (CPM)",
        fontsize=9
    )

    ax.set_ylabel(
        "Visium spot (CPM)",
        fontsize=9
    )


    ax.ticklabel_format(
        style="sci",
        scilimits=(
            0,
            0
        ),
        axis="both"
    )


    ax.text(
        0.04,
        0.93,
        f"cosine {r.cosine:.2f}\n"
        f"Spearman {r.spearman:+.2f}",
        transform=ax.transAxes,
        fontsize=8.5,
        va="top",
        bbox=dict(
            boxstyle="round,pad=0.3",
            fc="white",
            ec="#999999",
            lw=0.5
        )
    )


    if i == 0:

        ax.set_title(
            "c   Normalized expression\n"
            "     correspondence",
            fontsize=10,
            loc="left",
            fontweight="bold"
        )


# ============================================================
# PANELS d / e
#
# Distribution uses ALL evaluable imputed-only genes,
# not just the high-confidence display subset.
# ============================================================

edges = np.linspace(
    0,
    1,
    26
)


dist_specs = [
    (
        imp_pat,
        "spatially patterned",
        0.55,
        slice(
            0,
            3
        )
    ),
    (
        imp_non,
        "not patterned",
        0.06,
        slice(
            3,
            4
        )
    ),
]


for k, (
    msk,
    lbl,
    ref,
    rows
) in enumerate(
    dist_specs
):

    axd = fig.add_subplot(
        gs[
            rows,
            3
        ]
    )


    vals = (
        res.loc[
            msk,
            "cosine"
        ]
        .to_numpy(float)
    )


    vals = vals[
        np.isfinite(vals)
    ]


    hist = np.histogram(
        vals,
        bins=edges
    )[0]


    centres = (
        edges[:-1]
        +
        edges[1:]
    ) / 2


    axd.barh(
        centres,
        hist,
        height=(
            edges[1]
            -
            edges[0]
        ) * 0.9,
        color="#c8c8c8",
        edgecolor="white",
        linewidth=0.3
    )


    # Clifton reference
    axd.axhline(
        ref,
        color="#0072B2",
        lw=1.4,
        ls="--"
    )


    ours = np.nanmedian(
        vals
    )


    axd.axhline(
        ours,
        color="#333333",
        lw=1.4
    )


    # Stable axis-relative labels:
    # reference right; ours left.
    axd.text(
        0.98,
        ref + 0.02,
        f"Clifton {ref:.2f}",
        transform=axd.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8,
        color="#0072B2"
    )


    axd.text(
        0.03,
        ours + 0.015,
        f"ours {ours:.2f}",
        transform=axd.get_yaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=8,
        color="#333333",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.80,
            pad=1
        )
    )


    # Red lines = four displayed genes
    for g in G:

        r = res_idx.loc[
            g
        ]


        belongs = (
            bool(
                r[
                    "spatially_patterned"
                ]
            )
            if k == 0
            else
            not bool(
                r[
                    "spatially_patterned"
                ]
            )
        )


        if not belongs:
            continue


        c = float(
            r[
                "cosine"
            ]
        )


        axd.axhline(
            c,
            color="#E8000B",
            lw=2
        )


        axd.text(
            1.02,
            c,
            g,
            transform=axd.get_yaxis_transform(),
            ha="left",
            va="center",
            fontsize=8.5,
            style="italic",
            clip_on=False
        )


    axd.set_ylim(
        0,
        1
    )


    axd.set_ylabel(
        "cosine similarity"
    )


    axd.set_xlabel(
        "gene count"
    )


    axd.set_title(
        f"{'d' if k == 0 else 'e'}   "
        f"{lbl} (n={len(vals):,})",
        fontweight="bold",
        loc="left",
        fontsize=9.5
    )


# ============================================================
# COLORBAR
# ============================================================

cax = fig.add_axes(
    [
        0.13,
        0.030,
        0.20,
        0.010
    ]
)


plt.colorbar(
    plt.cm.ScalarMappable(
        cmap=RED
    ),
    cax=cax,
    orientation="horizontal"
)


cax.set_xticks(
    [
        0,
        1
    ]
)


cax.set_xticklabels([
    "0",
    "high"
])


cax.set_title(
    "normalized expression (1st–99th pct)",
    fontsize=8
)


# ============================================================
# TITLE / CAPTION
# ============================================================

fig.suptitle(
    f"External spatial validation of ENVI-imputed expression "
    f"({BR})",
    fontsize=13,
    fontweight="bold",
    y=0.995
)


fig.text(
    0.5,
    0.002,
    f"ENVI full-transcriptome predictions were evaluated against the paired Visium section. "
    f"Panel b shows ENVI-imputed expression in aligned Xenium cells passing match probability "
    f"> {MATCH_THRESH}; panels c–e use only cells whose aligned centroids fall inside physical "
    f"Visium spots. These cells were summed into ENVI pseudospots and both platforms were "
    f"CPM-normalized over the evaluable shared-gene universe. "
    f"Rows show high-, median-, and low-cosine spatially patterned imputed-only genes plus "
    f"one non-patterned imputed-only control. Spatial patterning was defined independently "
    f"from Visium Moran's I with BH q < 0.05. Panels d/e show cosine distributions for all "
    f"evaluable imputed-only patterned/non-patterned genes.",
    ha="center",
    va="bottom",
    fontsize=8.3,
    color="#444444",
    wrap=True
)


# ============================================================
# SAVE
# ============================================================

PNG = (
    OUT
    / f"fig_ENVI_fulltranscriptome_{BR}_PI_style.png"
)

PDF = (
    OUT
    / f"fig_ENVI_fulltranscriptome_{BR}_PI_style.pdf"
)


fig.savefig(
    PNG,
    dpi=600,
    bbox_inches="tight"
)


fig.savefig(
    PDF,
    bbox_inches="tight"
)


plt.close(fig)


# ============================================================
# SAVE DISPLAY GENE TABLE
# ============================================================

selected = (
    res_idx.loc[
        G
    ]
    .reset_index()
)


selected.insert(
    1,
    "display_row",
    GLAB
)


SEL = (
    OUT
    / f"{BR}_ENVI_selected_display_genes.csv"
)


selected.to_csv(
    SEL,
    index=False
)


e.file.close()


print()
print("=" * 80)
print("OUTPUT")
print("=" * 80)

print(
    "PNG      :",
    PNG
)

print(
    "PDF      :",
    PDF
)

print(
    "selection:",
    SEL
)
