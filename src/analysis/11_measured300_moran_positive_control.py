#!/usr/bin/env python3

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.io import mmread
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, pearsonr


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
MORAN_K = int(os.environ.get("MORAN_K", "6"))
CHUNK = 512


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

VRAW = (
    ROOT
    / "data/raw/visium"
    / BR
)

OUTDIR = (
    ROOT
    / "results/analysis/measured300_moran_visium_vs_xenium"
    / BR
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# HELPERS
# ============================================================

def as_bool(s):

    if pd.api.types.is_bool_dtype(s):
        return s.to_numpy(bool)

    return (
        s.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "t", "yes"])
        .to_numpy(bool)
    )


def pick_one(path, pattern):

    hits = sorted(path.glob(pattern))

    if len(hits) != 1:
        raise RuntimeError(
            f"{BR}: expected exactly one {pattern} "
            f"under {path}; found {len(hits)}"
        )

    return hits[0]


def build_knn_weights(xy, k=6):

    xy = np.asarray(xy, dtype=float)

    n = len(xy)

    if n < 3:
        raise RuntimeError("Too few spots.")

    k = min(k, n - 1)

    tree = cKDTree(xy)

    _, nn = tree.query(
        xy,
        k=k + 1
    )

    nn = nn[:, 1:]

    rows = np.repeat(
        np.arange(n),
        k
    )

    cols = nn.reshape(-1)

    data = np.full(
        len(rows),
        1.0 / k,
        dtype=float
    )

    return sparse.csr_matrix(
        (
            data,
            (rows, cols)
        ),
        shape=(n, n)
    )


def moran_matrix(X, W):

    X = np.asarray(
        X,
        dtype=np.float64
    )

    n = X.shape[0]

    Z = X - X.mean(
        axis=0,
        keepdims=True
    )

    den = np.sum(
        Z * Z,
        axis=0
    )

    WZ = W @ Z

    num = np.sum(
        Z * WZ,
        axis=0
    )

    S0 = float(W.sum())

    out = np.full(
        X.shape[1],
        np.nan
    )

    ok = den > 0

    out[ok] = (
        n / S0
    ) * (
        num[ok]
        /
        den[ok]
    )

    return out


# ============================================================
# INPUT CHECK
# ============================================================

for f in [
    ENVI_FILE,
    ALIGN_FILE,
    VCOORD_FILE,
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
# LOAD ENVI
# ============================================================

e = ad.read_h5ad(
    ENVI_FILE,
    backed="r"
)

print("=" * 90)
print(f"MEASURED-300 MORAN POSITIVE CONTROL: {BR}")
print("=" * 90)

print(
    "ENVI shape:",
    e.shape
)


if "measured_in_xenium" not in e.var.columns:
    raise RuntimeError(
        "ENVI var missing measured_in_xenium."
    )


measured_flag = as_bool(
    e.var[
        "measured_in_xenium"
    ]
)


egenes = pd.Index(
    e.var_names.astype(str)
)


measured_genes_all = (
    egenes[
        measured_flag
    ]
    .astype(str)
    .tolist()
)


print(
    "measured Xenium genes in ENVI:",
    len(measured_genes_all)
)


# ============================================================
# STALIGN
# ============================================================

aa = pd.read_csv(
    ALIGN_FILE,
    low_memory=False
)

vcoord = pd.read_csv(
    VCOORD_FILE
)


if "cell_id" not in e.obs.columns:
    raise RuntimeError(
        "ENVI obs missing cell_id."
    )


lookup = {
    str(x): i
    for i, x in enumerate(
        aa["cell_id"].astype(str)
    )
}


st_row = np.array([
    lookup.get(
        str(x),
        -1
    )
    for x in e.obs[
        "cell_id"
    ].astype(str)
])


has_st = st_row >= 0

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

vis_spot = np.full(
    e.n_obs,
    "",
    dtype=object
)


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


assigned = (
    has_st
    &
    np.isfinite(match_prob)
    &
    (match_prob > MATCH_THRESH)
    &
    (vis_spot != "")
    &
    (vis_spot != "nan")
)


assigned_rows = np.flatnonzero(
    assigned
)


print(
    "cells inside physical Visium spots:",
    len(assigned_rows)
)


# ============================================================
# SPOT ASSIGNMENT
# ============================================================

target_barcodes = (
    vcoord["barcode"]
    .astype(str)
    .tolist()
)


barcode_to_spot = {
    b: i
    for i, b in enumerate(
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


if (spot_idx < 0).any():

    raise RuntimeError(
        "Some assigned cells map to unknown Visium spots."
    )


n_spots = len(
    target_barcodes
)


ncell = np.bincount(
    spot_idx,
    minlength=n_spots
)


occupied = ncell > 0


print(
    "occupied spots:",
    int(occupied.sum()),
    "/",
    n_spots
)


# ============================================================
# SAME GRAPH AS IMPUTED ANALYSIS
# ============================================================

xy = (
    vcoord.loc[
        occupied,
        ["x", "y"]
    ]
    .to_numpy(float)
)


W = build_knn_weights(
    xy,
    MORAN_K
)


print(
    "Moran kNN:",
    MORAN_K
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
        f"Unexpected Visium matrix shape {VM.shape}"
    )


raw_bar_lookup = {
    b: i
    for i, b in enumerate(
        vbarcodes
    )
}


vrows = np.array([
    raw_bar_lookup[b]
    for b in target_barcodes
])


VT = VM[
    vrows
].tocsr()


vgene_lookup = {}

for j, g in enumerate(
    vgenes
):

    if g not in vgene_lookup:
        vgene_lookup[g] = j


# ============================================================
# FULL SHARED TRANSCRIPTOME
#
# IMPORTANT:
# Use exactly the same normalization denominator as the
# full-imputed Moran analysis:
#
# all ENVI genes shared with paired Visium.
# ============================================================

shared_genes = [
    g
    for g in egenes.astype(str)
    if g in vgene_lookup
]


measured_genes = [
    g
    for g in measured_genes_all
    if g in vgene_lookup
]


print()
print("=" * 90)
print("GENE UNIVERSE")
print("=" * 90)

print(
    "all ENVI genes shared with Visium:",
    len(shared_genes)
)

print(
    "measured genes total:",
    len(measured_genes_all)
)

print(
    "measured genes shared with Visium:",
    len(measured_genes)
)


# ============================================================
# LOAD SAME ASSIGNED ENVI CELL MATRIX
# ============================================================

print()
print(
    "Loading assigned-cell count_scale matrix..."
)


E_assigned = e.layers[
    "count_scale"
][
    assigned_rows,
    :
]


if sparse.issparse(E_assigned):
    E_assigned = E_assigned.toarray()


E_assigned = np.asarray(
    E_assigned,
    dtype=np.float32
)


print(
    "matrix:",
    E_assigned.shape
)

print(
    "memory:",
    f"{E_assigned.nbytes / 1024**3:.2f} GB"
)


# ============================================================
# SAME FULL-TRANSCRIPTOME CPM DENOMINATOR
# ============================================================

shared_eidx = np.array([
    egenes.get_loc(g)
    for g in shared_genes
])


shared_vidx = np.array([
    vgene_lookup[g]
    for g in shared_genes
])


P_lib = np.zeros(
    n_spots,
    dtype=np.float64
)

V_lib = np.zeros(
    n_spots,
    dtype=np.float64
)


print()
print(
    "Computing full-shared-transcriptome library sizes..."
)


for s in range(
    0,
    len(shared_genes),
    CHUNK
):

    t = min(
        s + CHUNK,
        len(shared_genes)
    )


    Echunk = E_assigned[
        :,
        shared_eidx[s:t]
    ].astype(
        np.float64,
        copy=False
    )


    Pchunk = np.zeros(
        (
            n_spots,
            t - s
        ),
        dtype=np.float64
    )


    np.add.at(
        Pchunk,
        spot_idx,
        Echunk
    )


    P_lib += np.sum(
        Pchunk,
        axis=1
    )


    Vchunk = VT[
        :,
        shared_vidx[s:t]
    ]


    V_lib += np.asarray(
        Vchunk.sum(axis=1)
    ).ravel()


P_lib[
    P_lib <= 0
] = 1.0

V_lib[
    V_lib <= 0
] = 1.0


# ============================================================
# MEASURED 300 EXPRESSION
# ============================================================

m_eidx = np.array([
    egenes.get_loc(g)
    for g in measured_genes
])


m_vidx = np.array([
    vgene_lookup[g]
    for g in measured_genes
])


E_m = E_assigned[
    :,
    m_eidx
].astype(
    np.float64,
    copy=False
)


P_m = np.zeros(
    (
        n_spots,
        len(measured_genes)
    ),
    dtype=np.float64
)


np.add.at(
    P_m,
    spot_idx,
    E_m
)


V_m = (
    VT[
        :,
        m_vidx
    ]
    .toarray()
    .astype(float)
)


# same occupied spots
P_m = P_m[
    occupied,
    :
]

V_m = V_m[
    occupied,
    :
]


# same normalization
Pcpm = (
    P_m
    /
    P_lib[
        occupied,
        None
    ]
    *
    1e6
)


Vcpm = (
    V_m
    /
    V_lib[
        occupied,
        None
    ]
    *
    1e6
)


# same transformation
Pexpr = np.log1p(
    Pcpm
)

Vexpr = np.log1p(
    Vcpm
)


# ============================================================
# MORAN'S I
# ============================================================

moran_xenium = moran_matrix(
    Pexpr,
    W
)

moran_visium = moran_matrix(
    Vexpr,
    W
)


out = pd.DataFrame({
    "gene":
        measured_genes,

    "moran_visium":
        moran_visium,

    "moran_xenium_measured":
        moran_xenium,
})


out[
    "delta_moran_xenium_minus_visium"
] = (
    out[
        "moran_xenium_measured"
    ]
    -
    out[
        "moran_visium"
    ]
)


OUT_GENE = (
    OUTDIR
    / f"{BR}_measured300_Visium_vs_Xenium_MoranI.csv"
)


out.to_csv(
    OUT_GENE,
    index=False
)


# ============================================================
# VALID MORAN IN BOTH
# ============================================================

ok = (
    np.isfinite(
        out["moran_visium"]
    )
    &
    np.isfinite(
        out["moran_xenium_measured"]
    )
)


x = (
    out.loc[
        ok,
        "moran_visium"
    ]
    .to_numpy(float)
)

y = (
    out.loc[
        ok,
        "moran_xenium_measured"
    ]
    .to_numpy(float)
)


if len(x) >= 3:

    rho, rho_p = spearmanr(
        x,
        y
    )

    r, r_p = pearsonr(
        x,
        y
    )

    slope, intercept = np.polyfit(
        x,
        y,
        1
    )

else:

    rho = rho_p = np.nan
    r = r_p = np.nan
    slope = intercept = np.nan


delta = y - x


summary = pd.DataFrame([
    {
        "BrNum":
            BR,

        "n_assigned_cells":
            len(assigned_rows),

        "n_occupied_spots":
            int(occupied.sum()),

        "n_measured_total":
            len(measured_genes_all),

        "n_measured_shared":
            len(measured_genes),

        "n_measured_valid_moran_both":
            len(x),

        "moran_knn_k":
            MORAN_K,

        "spearman_moran":
            rho,

        "spearman_p":
            rho_p,

        "pearson_moran":
            r,

        "pearson_p":
            r_p,

        "regression_slope":
            slope,

        "regression_intercept":
            intercept,

        "median_visium_moran":
            np.nanmedian(x),

        "median_xenium_measured_moran":
            np.nanmedian(y),

        "median_delta_moran":
            np.nanmedian(delta),

        "fraction_xenium_gt_visium":
            np.mean(y > x),
    }
])


OUT_SUMMARY = (
    OUTDIR
    / f"{BR}_measured300_MoranI_summary.csv"
)


summary.to_csv(
    OUT_SUMMARY,
    index=False
)


# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.3, 7.1)
)


ax.scatter(
    x,
    y,
    s=25,
    alpha=0.65,
    linewidths=0
)


allv = np.concatenate([
    x,
    y
])


lo = float(
    np.nanmin(allv)
)

hi = float(
    np.nanmax(allv)
)


pad = max(
    0.03,
    0.05 * (hi - lo)
)

lo -= pad
hi += pad


# ideal line
ax.plot(
    [lo, hi],
    [lo, hi],
    linestyle="--",
    linewidth=1.2,
    label="y = x"
)


# regression
xx = np.linspace(
    lo,
    hi,
    200
)


ax.plot(
    xx,
    intercept + slope * xx,
    linewidth=1.3,
    label="Gene-level fit"
)


ax.set_xlim(
    lo,
    hi
)

ax.set_ylim(
    lo,
    hi
)

ax.set_aspect(
    "equal",
    adjustable="box"
)


ax.set_xlabel(
    "Visium Moran's I",
    fontsize=12
)

ax.set_ylabel(
    "Measured Xenium pseudospot Moran's I",
    fontsize=12
)


ax.set_title(
    f"{BR}: measured-gene spatial-autocorrelation control",
    fontsize=13,
    fontweight="bold"
)


ax.text(
    0.03,
    0.97,
    f"Measured genes shared: {len(measured_genes):,}\n"
    f"Valid Moran's I in both: {len(x):,}\n"
    f"Spearman ρ={rho:.3f}\n"
    f"Pearson r={r:.3f}\n"
    f"median ΔI={np.nanmedian(delta):+.3f}\n"
    f"slope={slope:.3f}",
    transform=ax.transAxes,
    va="top",
    fontsize=10,
    bbox=dict(
        boxstyle="round,pad=0.35",
        facecolor="white",
        edgecolor="0.7",
        alpha=0.92
    )
)


ax.legend(
    frameon=False,
    loc="lower right"
)


fig.tight_layout()


PNG = (
    OUTDIR
    / f"fig_{BR}_measured300_Visium_vs_Xenium_MoranI.png"
)

PDF = (
    OUTDIR
    / f"fig_{BR}_measured300_Visium_vs_Xenium_MoranI.pdf"
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

e.file.close()


print()
print("=" * 90)
print("MEASURED-300 POSITIVE CONTROL RESULT")
print("=" * 90)

print(
    "Measured shared       :",
    len(measured_genes)
)

print(
    "Valid Moran both      :",
    len(x)
)

print(
    f"Spearman Moran        : {rho:.4f}"
)

print(
    f"Pearson Moran         : {r:.4f}"
)

print(
    f"Median Visium I       : {np.nanmedian(x):.4f}"
)

print(
    f"Median measured X I   : {np.nanmedian(y):.4f}"
)

print(
    f"Median delta I        : {np.nanmedian(delta):+.4f}"
)

print(
    f"Regression slope      : {slope:.4f}"
)

print()
print("Gene table :", OUT_GENE)
print("Summary    :", OUT_SUMMARY)
print("PNG        :", PNG)
print("PDF        :", PDF)
