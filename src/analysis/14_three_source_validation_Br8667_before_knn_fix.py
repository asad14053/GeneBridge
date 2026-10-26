#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from scipy import sparse
from scipy.io import mmread
from scipy.spatial import Delaunay
from scipy.stats import norm, rankdata


# ============================================================
# CONFIG
# ============================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BR = "Br8667"

MATCH_THRESH = 0.85


OOF = (
    ROOT
    / "outputs/imputation_full/ENVI_N24"
    / BR
    / "ex5_1/envi/oof_300/oof_predicted_300genes.h5ad"
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

OUT = (
    ROOT
    / "results/analysis/three_source_validation"
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
            f"Expected one {pattern} in {path}; "
            f"found {len(hits)}"
        )

    return hits[0]


def dense(x):

    if sparse.issparse(x):
        return x.toarray()

    return np.asarray(x)


def bh(p):

    p = np.asarray(
        p,
        dtype=float
    )

    q = np.full(
        len(p),
        np.nan
    )

    ok = np.isfinite(p)

    pv = p[ok]

    if len(pv) == 0:
        return q

    order = np.argsort(pv)

    ranked = pv[order]

    adj = (
        ranked
        *
        len(ranked)
        /
        np.arange(
            1,
            len(ranked) + 1
        )
    )

    adj = np.minimum.accumulate(
        adj[::-1]
    )[::-1]

    adj = np.minimum(
        adj,
        1
    )

    tmp = np.empty_like(adj)

    tmp[order] = adj

    q[ok] = tmp

    return q


def cosine_cols(A, B):

    num = np.sum(
        A * B,
        axis=0
    )

    den = (
        np.sqrt(
            np.sum(
                A * A,
                axis=0
            )
        )
        *
        np.sqrt(
            np.sum(
                B * B,
                axis=0
            )
        )
    )

    out = np.full(
        A.shape[1],
        np.nan
    )

    ok = den > 0

    out[ok] = (
        num[ok]
        /
        den[ok]
    )

    return out


def spearman_cols(A, B):

    RA = rankdata(
        A,
        axis=0
    )

    RB = rankdata(
        B,
        axis=0
    )

    RA = (
        RA
        -
        RA.mean(
            axis=0,
            keepdims=True
        )
    )

    RB = (
        RB
        -
        RB.mean(
            axis=0,
            keepdims=True
        )
    )

    num = np.sum(
        RA * RB,
        axis=0
    )

    den = (
        np.sqrt(
            np.sum(
                RA**2,
                axis=0
            )
        )
        *
        np.sqrt(
            np.sum(
                RB**2,
                axis=0
            )
        )
    )

    out = np.full(
        A.shape[1],
        np.nan
    )

    ok = den > 0

    out[ok] = (
        num[ok]
        /
        den[ok]
    )

    return out


def build_spatial_W(xy):

    tri = Delaunay(xy)

    edges = set()

    for t in tri.simplices:

        for a, b in [
            (t[0], t[1]),
            (t[1], t[2]),
            (t[0], t[2]),
        ]:

            if a > b:
                a, b = b, a

            edges.add(
                (a, b)
            )

    edges = list(edges)

    dist = np.array([
        np.linalg.norm(
            xy[a] - xy[b]
        )
        for a, b in edges
    ])

    pitch = np.median(dist)

    max_edge = (
        1.5
        *
        pitch
    )

    keep = [
        (a, b)
        for (a, b), d
        in zip(
            edges,
            dist
        )
        if d <= max_edge
    ]

    rows = []
    cols = []

    for a, b in keep:

        rows.extend(
            [a, b]
        )

        cols.extend(
            [b, a]
        )

    data = np.ones(
        len(rows),
        dtype=float
    )

    W = sparse.csr_matrix(
        (
            data,
            (
                rows,
                cols
            )
        ),
        shape=(
            len(xy),
            len(xy)
        )
    )

    return W, pitch


def moran_chunk(X, W):

    X = np.asarray(
        X,
        dtype=float
    )

    n = X.shape[0]

    z = (
        X
        -
        np.mean(
            X,
            axis=0,
            keepdims=True
        )
    )

    zz = np.sum(
        z * z,
        axis=0
    )

    Wz = W @ z

    num = np.sum(
        z * Wz,
        axis=0
    )

    S0 = float(
        W.sum()
    )

    I = np.full(
        X.shape[1],
        np.nan
    )

    ok = zz > 0

    I[ok] = (
        n
        /
        S0
    ) * (
        num[ok]
        /
        zz[ok]
    )

    # same analytical positive-autocorrelation p-value
    EI = (
        -1.0
        /
        (n - 1)
    )

    Wplus = W + W.T

    S1 = (
        0.5
        *
        float(
            Wplus.multiply(
                Wplus
            ).sum()
        )
    )

    rsum = np.asarray(
        W.sum(axis=1)
    ).ravel()

    csum = np.asarray(
        W.sum(axis=0)
    ).ravel()

    S2 = float(
        np.sum(
            (
                rsum
                +
                csum
            ) ** 2
        )
    )

    b2 = np.full(
        X.shape[1],
        np.nan
    )

    b2[ok] = (
        n
        *
        np.sum(
            z[:, ok] ** 4,
            axis=0
        )
        /
        (
            zz[ok] ** 2
        )
    )

    A = (
        n
        *
        (
            (
                n**2
                -
                3*n
                +
                3
            )
            *
            S1
            -
            n * S2
            +
            3 * S0**2
        )
    )

    B = (
        b2
        *
        (
            (
                n**2
                -
                n
            )
            *
            S1
            -
            2*n*S2
            +
            6*S0**2
        )
    )

    VI = np.full(
        X.shape[1],
        np.nan
    )

    VI[ok] = (
        (
            A
            -
            B[ok]
        )
        /
        (
            (n - 1)
            *
            (n - 2)
            *
            (n - 3)
            *
            S0**2
        )
        -
        EI**2
    )

    p = np.full(
        X.shape[1],
        np.nan
    )

    vok = (
        ok
        &
        np.isfinite(VI)
        &
        (VI > 0)
    )

    zscore = (
        I[vok]
        -
        EI
    ) / np.sqrt(
        VI[vok]
    )

    p[vok] = (
        1
        -
        norm.cdf(
            zscore
        )
    )

    return I, p


# ============================================================
# INPUT FILES
# ============================================================

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
# LOAD OOF
# ============================================================

o = ad.read_h5ad(
    OOF
)

print("=" * 100)
print("BR8667 THREE-SOURCE VALIDATION")
print("=" * 100)

print(
    "OOF:",
    o.shape
)

if o.n_vars != 300:
    raise RuntimeError(
        f"Expected 300 OOF genes, got {o.n_vars}"
    )


genes = np.asarray(
    o.var_names,
    dtype=str
)


# ============================================================
# B = MEASURED XENIUM
# C = ENVI OOF PREDICTION
# ============================================================

B_cell = dense(
    o.layers[
        "ground_truth_counts"
    ]
).astype(
    np.float64
)

C_cell = dense(
    o.layers[
        "count_scale"
    ]
).astype(
    np.float64
)

print(
    "B measured:",
    B_cell.shape,
    "min/max",
    B_cell.min(),
    B_cell.max()
)

print(
    "C OOF:",
    C_cell.shape,
    "min/max",
    C_cell.min(),
    C_cell.max()
)


# ============================================================
# STALIGN CELL MATCH
# ============================================================

aa = pd.read_csv(
    ALIGN_FILE,
    low_memory=False
)

vcoord = pd.read_csv(
    VCOORD_FILE
)


# audit showed this is full 66148 overlap
okey = (
    o.obs[
        "original_cell_id"
    ]
    .astype(str)
    .to_numpy()
)

akey = (
    aa[
        "cell_id"
    ]
    .astype(str)
)

lookup = {
    k: i
    for i, k
    in enumerate(
        akey
    )
}


st_row = np.array([
    lookup.get(
        k,
        -1
    )
    for k
    in okey
])


has_st = (
    st_row >= 0
)


print(
    "OOF cells matching STalign:",
    has_st.sum(),
    "/",
    o.n_obs
)


idx = np.flatnonzero(
    has_st
)

sub = aa.iloc[
    st_row[idx]
]


match_prob = np.full(
    o.n_obs,
    np.nan
)

vis_spot = np.full(
    o.n_obs,
    "",
    dtype=object
)


match_prob[idx] = pd.to_numeric(
    sub[
        "match_prob"
    ],
    errors="coerce"
).to_numpy()


vis_spot[idx] = (
    sub[
        "visium_spot"
    ]
    .fillna("")
    .astype(str)
    .to_numpy()
)


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
)


cell_rows = np.flatnonzero(
    use
)


print(
    "OOF cells used inside Visium spots:",
    len(cell_rows)
)


# ============================================================
# VISIUM SPOT ASSIGNMENT
# ============================================================

target_barcodes = (
    vcoord[
        "barcode"
    ]
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
    for x
    in vis_spot[
        cell_rows
    ]
])


if (
    spot_idx < 0
).any():

    raise RuntimeError(
        "Unknown Visium spot assignment."
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
    "occupied spots:",
    occupied.sum(),
    "/",
    n_spots
)


# ============================================================
# AGGREGATE B AND C USING IDENTICAL CELLS
# ============================================================

Aagg = sparse.csr_matrix(
    (
        np.ones(
            len(cell_rows),
            dtype=np.float64
        ),
        (
            spot_idx,
            np.arange(
                len(cell_rows)
            )
        )
    ),
    shape=(
        n_spots,
        len(cell_rows)
    )
)


B_spot = (
    Aagg
    @
    B_cell[
        cell_rows,
        :
    ]
)

C_spot = (
    Aagg
    @
    C_cell[
        cell_rows,
        :
    ]
)


B_spot = np.asarray(
    B_spot[
        occ_idx,
        :
    ],
    dtype=np.float64
)

C_spot = np.asarray(
    C_spot[
        occ_idx,
        :
    ],
    dtype=np.float64
)


# ============================================================
# LOAD A = REAL VISIUM
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


raw_bar_lookup = {
    b: i
    for i, b
    in enumerate(
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


missing = [
    g
    for g in genes
    if g not in vgene_lookup
]


if missing:
    raise RuntimeError(
        f"{len(missing)} OOF genes absent from Visium: "
        f"{missing[:20]}"
    )


vidx = np.array([
    vgene_lookup[g]
    for g in genes
])


A_spot = (
    VT[
        occ_idx
    ][
        :,
        vidx
    ]
    .toarray()
    .astype(
        np.float64
    )
)


print(
    "A Visium:",
    A_spot.shape
)

print(
    "B measured:",
    B_spot.shape
)

print(
    "C OOF:",
    C_spot.shape
)


# ============================================================
# EVALUABILITY
# ============================================================

A_total = A_spot.sum(
    axis=0
)

B_total = B_spot.sum(
    axis=0
)

C_total = C_spot.sum(
    axis=0
)


gene_ok = (
    (A_total > 0)
    &
    (B_total > 0)
    &
    (C_total > 0)
)


print(
    "genes evaluable in all 3:",
    gene_ok.sum(),
    "/ 300"
)


# ============================================================
# NORMALIZATION
#
# Since all sources contain the SAME 300 genes here,
# normalize each spot over these 300 genes.
# ============================================================

def cpm(X):

    lib = X[
        :,
        gene_ok
    ].sum(
        axis=1
    )

    lib = np.asarray(
        lib,
        dtype=float
    )

    lib[
        lib <= 0
    ] = 1.0

    return (
        X
        /
        lib[
            :,
            None
        ]
        *
        1e6
    )


A_cpm = cpm(
    A_spot
)

B_cpm = cpm(
    B_spot
)

C_cpm = cpm(
    C_spot
)


# raw cosine uses same definition as existing workflow
A_log = np.log10(
    A_cpm + 1
)

B_log = np.log10(
    B_cpm + 1
)

C_log = np.log10(
    C_cpm + 1
)


# ============================================================
# PAIRWISE EXPRESSION CORRESPONDENCE
# ============================================================

cos_AB = cosine_cols(
    A_log,
    B_log
)

cos_AC = cosine_cols(
    A_log,
    C_log
)

cos_BC = cosine_cols(
    B_log,
    C_log
)


sp_AB = spearman_cols(
    A_log,
    B_log
)

sp_AC = spearman_cols(
    A_log,
    C_log
)

sp_BC = spearman_cols(
    B_log,
    C_log
)


# ============================================================
# SAME SPATIAL GRAPH FOR ALL THREE
# ============================================================

xy = (
    vcoord.loc[
        occupied,
        [
            "x",
            "y"
        ]
    ]
    .to_numpy(
        float
    )
)


W, pitch = build_spatial_W(
    xy
)


print(
    "spatial graph spots:",
    len(xy)
)

print(
    "estimated pitch:",
    pitch
)


# ============================================================
# MORAN'S I — SAME SPOTS + SAME W
# ============================================================

moran_A, p_A = moran_chunk(
    A_cpm,
    W
)

moran_B, p_B = moran_chunk(
    B_cpm,
    W
)

moran_C, p_C = moran_chunk(
    C_cpm,
    W
)


q_A = bh(
    p_A
)

q_B = bh(
    p_B
)

q_C = bh(
    p_C
)


pat_A = (
    gene_ok
    &
    np.isfinite(
        q_A
    )
    &
    (
        q_A < 0.05
    )
    &
    (
        moran_A > 0
    )
)


pat_B = (
    gene_ok
    &
    np.isfinite(
        q_B
    )
    &
    (
        q_B < 0.05
    )
    &
    (
        moran_B > 0
    )
)


pat_C = (
    gene_ok
    &
    np.isfinite(
        q_C
    )
    &
    (
        q_C < 0.05
    )
    &
    (
        moran_C > 0
    )
)


# ============================================================
# OUTPUT GENE TABLE
# ============================================================

res = pd.DataFrame({
    "gene":
        genes,

    "evaluable_all3":
        gene_ok,

    # expression correspondence
    "cosine_A_Visium_vs_B_measured":
        cos_AB,

    "spearman_A_Visium_vs_B_measured":
        sp_AB,

    "cosine_A_Visium_vs_C_OOF":
        cos_AC,

    "spearman_A_Visium_vs_C_OOF":
        sp_AC,

    "cosine_B_measured_vs_C_OOF":
        cos_BC,

    "spearman_B_measured_vs_C_OOF":
        sp_BC,

    # Moran
    "moran_A_Visium":
        moran_A,

    "moran_B_measured":
        moran_B,

    "moran_C_OOF":
        moran_C,

    "moran_q_A_Visium":
        q_A,

    "moran_q_B_measured":
        q_B,

    "moran_q_C_OOF":
        q_C,

    # pattern calls
    "patterned_A_Visium":
        pat_A,

    "patterned_B_measured":
        pat_B,

    "patterned_C_OOF":
        pat_C,
})


OUT_GENE = (
    OUT
    / "Br8667_three_source_gene_metrics.csv"
)


res.to_csv(
    OUT_GENE,
    index=False
)


# ============================================================
# SAVE SPOT-LEVEL MATRICES FOR FIGURE
# ============================================================

np.savez_compressed(
    OUT
    / "Br8667_three_source_spot_matrices.npz",

    genes=genes,

    barcodes=np.asarray(
        target_barcodes,
        dtype=object
    )[
        occupied
    ],

    xy=xy,

    visium_cpm=A_cpm,

    measured_xenium_cpm=B_cpm,

    oof_envi_cpm=C_cpm,
)


# ============================================================
# SUMMARY
# ============================================================

def med(x):

    x = np.asarray(
        x,
        dtype=float
    )

    x = x[
        np.isfinite(x)
    ]

    if len(x) == 0:
        return np.nan

    return float(
        np.median(x)
    )


summary = pd.DataFrame([
    {
        "BrNum":
            BR,

        "n_oof_cells":
            o.n_obs,

        "n_cells_inside_visium_spots":
            len(cell_rows),

        "n_occupied_visium_spots":
            int(
                occupied.sum()
            ),

        "n_genes":
            len(genes),

        "n_evaluable_all3":
            int(
                gene_ok.sum()
            ),

        "n_patterned_visium":
            int(
                pat_A.sum()
            ),

        "n_patterned_measured_xenium":
            int(
                pat_B.sum()
            ),

        "n_patterned_oof":
            int(
                pat_C.sum()
            ),

        "median_cosine_AB":
            med(
                cos_AB[
                    gene_ok
                ]
            ),

        "median_cosine_AC":
            med(
                cos_AC[
                    gene_ok
                ]
            ),

        "median_cosine_BC":
            med(
                cos_BC[
                    gene_ok
                ]
            ),

        "median_spearman_AB":
            med(
                sp_AB[
                    gene_ok
                ]
            ),

        "median_spearman_AC":
            med(
                sp_AC[
                    gene_ok
                ]
            ),

        "median_spearman_BC":
            med(
                sp_BC[
                    gene_ok
                ]
            ),

        "median_moran_A":
            med(
                moran_A[
                    gene_ok
                ]
            ),

        "median_moran_B":
            med(
                moran_B[
                    gene_ok
                ]
            ),

        "median_moran_C":
            med(
                moran_C[
                    gene_ok
                ]
            ),
    }
])


OUT_SUMMARY = (
    OUT
    / "Br8667_three_source_summary.csv"
)


summary.to_csv(
    OUT_SUMMARY,
    index=False
)


print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    summary.T.to_string(
        header=False
    )
)

print()
print("Gene metrics :", OUT_GENE)
print("Summary      :", OUT_SUMMARY)
print(
    "Spot data    :",
    OUT
    / "Br8667_three_source_spot_matrices.npz"
)

o.file.close()
