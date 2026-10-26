#!/usr/bin/env python3

import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

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

BR = os.environ.get(
    "BR",
    "Br8667"
)

MATCH_THRESH = 0.85

# Parallel settings.
# SLURM_CPUS_PER_TASK will be 16 in production.
N_JOBS = int(
    os.environ.get(
        "SLURM_CPUS_PER_TASK",
        os.environ.get("N_JOBS", "16")
    )
)

# Gene chunk sizes.
AGG_CHUNK = 1024
METRIC_CHUNK = 512

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
            f"{BR}: expected exactly one {pattern} "
            f"under {path}; found {len(hits)}:\n"
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

    order = np.argsort(
        pv
    )

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

    tmp = np.empty_like(
        adj
    )

    tmp[order] = adj

    q[ok] = tmp

    return q


def build_spatial_W(xy):

    tri = Delaunay(
        xy
    )

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

    # Estimate typical Visium pitch from triangulation edges.
    edges = list(edges)

    dist = np.array([
        np.linalg.norm(
            xy[a] - xy[b]
        )
        for a, b in edges
    ])

    pitch = np.median(
        dist
    )

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

    """
    X:
        spots x genes
        Visium CPM values.
    """

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

    EI = (
        -1.0
        /
        (n - 1)
    )

    Wplus = (
        W
        +
        W.T
    )

    S1 = 0.5 * float(
        Wplus.multiply(
            Wplus
        ).sum()
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
            n
            *
            S2
            +
            3
            *
            S0**2
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
            2*n
            *
            S2
            +
            6
            *
            S0**2
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

    # positive spatial autocorrelation
    p[vok] = (
        1
        -
        norm.cdf(
            zscore
        )
    )

    return I, p


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


# ============================================================
# INPUTS
# ============================================================

print(
    "=" * 80
)

print(
    f"ENVI FULL-TRANSCRIPTOME CROSS-TECH: {BR}"
)

print(
    "=" * 80
)

for f in [
    ENVI_FILE,
    ALIGN_FILE,
    VCOORD_FILE
]:

    if not f.exists():
        raise FileNotFoundError(
            f
        )


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


print(
    "ENVI   :",
    ENVI_FILE
)

print(
    "STalign:",
    STDIR
)

print(
    "Visium :",
    VMTX
)


# ============================================================
# LOAD ENVI
# ============================================================

e = ad.read_h5ad(
    ENVI_FILE,
    backed="r"
)

print()
print(
    "ENVI shape:",
    e.shape
)

if (
    "count_scale"
    not in e.layers.keys()
):

    raise RuntimeError(
        "ENVI count_scale layer missing."
    )

if (
    "measured_in_xenium"
    not in e.var.columns
):

    raise RuntimeError(
        "measured_in_xenium missing."
    )


# ============================================================
# LOAD STALIGN
# ============================================================

aa = pd.read_csv(
    ALIGN_FILE
)

vcoord = pd.read_csv(
    VCOORD_FILE
)

print(
    "STalign cells:",
    len(aa)
)

print(
    "STalign Visium spots:",
    len(vcoord)
)


# ============================================================
# MATCH ENVI CELLS TO STALIGN CELLS
# ============================================================

e_candidates = {
    "obs_names":
        pd.Index(
            e.obs_names.astype(str)
        )
}

for c in [
    "cell_id",
    "original_cell_id",
    "donor_cell_id"
]:

    if c in e.obs.columns:

        e_candidates[c] = pd.Index(
            e.obs[c]
            .astype(str)
        )


a_candidates = {}

for c in [
    "cell_id",
    "lieber_key"
]:

    if c in aa.columns:

        a_candidates[c] = pd.Index(
            aa[c]
            .astype(str)
        )


best = None

for ename, ei in e_candidates.items():

    eset = set(
        ei
    )

    for aname, ai in a_candidates.items():

        n_overlap = sum(
            x in eset
            for x in ai
        )

        rec = (
            n_overlap,
            ename,
            aname
        )

        if (
            best is None
            or rec[0] > best[0]
        ):
            best = rec


if (
    best is None
    or best[0] == 0
):

    raise RuntimeError(
        f"{BR}: could not match ENVI "
        "cells to STalign."
    )


n_overlap, ename, aname = best

print()
print(
    "Best cell-ID mapping:",
    f"ENVI {ename} <-> STalign {aname}"
)

print(
    "overlap:",
    n_overlap,
    "/",
    e.n_obs
)


ekey = np.asarray(
    e_candidates[ename],
    dtype=str
)

akey = aa[
    aname
].astype(str)


if akey.duplicated().any():

    raise RuntimeError(
        f"{BR}: duplicate STalign "
        f"{aname} identifiers."
    )


a_lookup = {
    k: i
    for i, k
    in enumerate(
        akey
    )
}


st_row = np.array([
    a_lookup.get(
        k,
        -1
    )
    for k
    in ekey
])

has_st = (
    st_row >= 0
)


# ============================================================
# MATCHED + ASSIGNED CELLS
# ============================================================

match_prob = np.full(
    e.n_obs,
    np.nan
)

vis_spot = np.full(
    e.n_obs,
    "",
    dtype=object
)


idx = np.flatnonzero(
    has_st
)

match_prob[idx] = (
    pd.to_numeric(
        aa.iloc[
            st_row[idx]
        ]["match_prob"],
        errors="coerce"
    )
    .to_numpy()
)

vis_spot[idx] = (
    aa.iloc[
        st_row[idx]
    ]["visium_spot"]
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

print()
print(
    "ENVI cells assigned to "
    "Visium spots:",
    len(cell_rows)
)


# ============================================================
# VISIUM BARCODE ORDER
# ============================================================

if "barcode" not in vcoord.columns:

    raise RuntimeError(
        "barcode column missing from "
        "STalign Visium coordinates."
    )


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
    for x
    in vis_spot[
        cell_rows
    ]
])


bad_spot = (
    spot_idx < 0
)

if bad_spot.any():

    raise RuntimeError(
        f"{BR}: {bad_spot.sum()} "
        "assigned cells have Visium "
        "barcodes absent from vcoord."
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
    occupied.sum(),
    "/",
    n_spots,
    f"({100*occupied.mean():.1f}%)"
)


# Sparse cell -> spot aggregation matrix
A = sparse.csr_matrix(
    (
        np.ones(
            len(cell_rows),
            dtype=np.float32
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


# ============================================================
# RAW VISIUM
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


# Ensure spots x genes.
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
        f"Unexpected Visium matrix "
        f"shape {VM.shape}"
    )


raw_bar_lookup = {
    b: i
    for i, b
    in enumerate(
        vbarcodes
    )
}


missing_barcodes = [
    b
    for b in target_barcodes
    if b not in raw_bar_lookup
]

if missing_barcodes:

    raise RuntimeError(
        f"{BR}: {len(missing_barcodes)} "
        "STalign target barcodes absent "
        "from raw Visium."
    )


vrows = np.array([
    raw_bar_lookup[b]
    for b
    in target_barcodes
])


VT = VM[
    vrows
].tocsr()


print()
print(
    "Raw Visium:",
    VT.shape
)


# ============================================================
# GENE MATCH
# ============================================================

# Preserve first raw Visium feature if symbols repeat.
vgene_lookup = {}

for j, g in enumerate(
    vgenes
):

    if g not in vgene_lookup:
        vgene_lookup[g] = j


egenes = np.asarray(
    e.var_names,
    dtype=str
)


shared_mask = np.array([
    g in vgene_lookup
    for g in egenes
])


shared_eidx = np.flatnonzero(
    shared_mask
)

shared_genes = egenes[
    shared_eidx
]

shared_vidx = np.array([
    vgene_lookup[g]
    for g in shared_genes
])


measured_flag = (
    e.var["measured_in_xenium"]
    .astype(bool)
    .to_numpy()[
        shared_eidx
    ]
)


print()
print(
    "ENVI genes:",
    e.n_vars
)

print(
    "shared with Visium:",
    len(
        shared_genes
    )
)

print(
    "shared measured:",
    measured_flag.sum()
)

print(
    "shared imputed-only:",
    (
        ~measured_flag
    ).sum()
)


# ============================================================
# FAST PARALLEL IMPLEMENTATION
#
# Old implementation repeatedly read the backed ENVI h5ad
# three times in small chunks.
#
# New implementation:
#   1. Read matched cells from count_scale ONCE.
#   2. Keep only genes shared with Visium.
#   3. Aggregate cells -> pseudospots once.
#   4. Materialize corresponding Visium matrix once.
#   5. Parallelize expensive per-gene metrics across CPUs.
# ============================================================

ng = len(shared_genes)

print()
print("=" * 80)
print("FAST PARALLEL MODE")
print("=" * 80)
print("CPUs                 :", N_JOBS)
print("matched ENVI cells   :", len(cell_rows))
print("shared genes         :", ng)
print("aggregation chunk    :", AGG_CHUNK)
print("metric chunk         :", METRIC_CHUNK)


# ============================================================
# LOAD MATCHED ENVI CELLS ONCE
# ============================================================

print()
print(
    "Loading count_scale for matched cells once..."
)

# First select rows from backed HDF5.
# This is ~2 GB for Br8667 rather than loading all layers/cells.
E_used = e.layers[
    "count_scale"
][
    cell_rows,
    :
]

E_used = np.asarray(
    E_used,
    dtype=np.float32
)

# Retain only genes shared with raw Visium.
E_used = E_used[
    :,
    shared_eidx
]

print(
    "ENVI matched matrix:",
    E_used.shape,
    E_used.dtype
)

print(
    "ENVI matrix memory  :",
    f"{E_used.nbytes / 1024**3:.2f} GB"
)


# ============================================================
# AGGREGATE ENVI -> VISIUM PSEUDOSPOTS IN PARALLEL
# ============================================================

print()
print(
    "Aggregating ENVI cells into Visium pseudospots..."
)

P_full = np.empty(
    (
        occupied.sum(),
        ng
    ),
    dtype=np.float32
)


agg_bounds = [
    (
        start,
        min(
            start + AGG_CHUNK,
            ng
        )
    )
    for start in range(
        0,
        ng,
        AGG_CHUNK
    )
]


def aggregate_chunk(bounds):

    s0, t0 = bounds

    out = A @ E_used[
        :,
        s0:t0
    ]

    out = np.asarray(
        out[
            occ_idx,
            :
        ],
        dtype=np.float32
    )

    return (
        s0,
        t0,
        out
    )


done = 0

with ThreadPoolExecutor(
    max_workers=N_JOBS
) as ex:

    futures = {
        ex.submit(
            aggregate_chunk,
            b
        ): b
        for b in agg_bounds
    }

    for fut in as_completed(
        futures
    ):

        s0, t0, out = fut.result()

        P_full[
            :,
            s0:t0
        ] = out

        done += (
            t0 - s0
        )

        print(
            f"  aggregation: "
            f"{done:5d}/{ng}",
            end="\r",
            flush=True
        )

print()

print(
    "Pseudospot matrix:",
    P_full.shape
)

print(
    "Pseudospot memory:",
    f"{P_full.nbytes / 1024**3:.2f} GB"
)


# ENVI cell-level matrix is no longer required.
del E_used


# ============================================================
# MATERIALIZE MATCHED VISIUM MATRIX ONCE
# ============================================================

print()
print(
    "Loading corresponding raw Visium matrix once..."
)

V_full = VT[
    occ_idx
][
    :,
    shared_vidx
]

V_full = V_full.toarray().astype(
    np.float32,
    copy=False
)

print(
    "Visium matrix:",
    V_full.shape
)

print(
    "Visium memory:",
    f"{V_full.nbytes / 1024**3:.2f} GB"
)


# ============================================================
# EXPRESSION ELIGIBILITY — VECTORIZED
# ============================================================

print()
print(
    "Computing expression eligibility..."
)

p_total = np.asarray(
    P_full.sum(
        axis=0
    )
).ravel()

v_total = np.asarray(
    V_full.sum(
        axis=0
    )
).ravel()

v_detect = np.asarray(
    (
        V_full > 0
    ).sum(
        axis=0
    )
).ravel()

v_max = np.asarray(
    V_full.max(
        axis=0
    )
).ravel()


gene_ok = (
    (p_total >= 1)
    &
    (v_total >= 1)
    &
    (v_detect > 1)
)


# Independent Visium-only display/high-confidence criterion.
visium_highconf = (
    (
        v_detect
        /
        occupied.sum()
    )
    >= 0.20
) & (
    v_max >= 10
)


print(
    "basic evaluable genes:",
    int(
        gene_ok.sum()
    )
)

print(
    "high-confidence Visium-expressed:",
    int(
        (
            gene_ok
            &
            visium_highconf
        ).sum()
    )
)


# ============================================================
# LIBRARY SIZES
# ============================================================

print()
print(
    "Computing CPM library sizes..."
)

P_lib = np.asarray(
    P_full[
        :,
        gene_ok
    ].sum(
        axis=1
    ),
    dtype=np.float64
).ravel()

V_lib = np.asarray(
    V_full[
        :,
        gene_ok
    ].sum(
        axis=1
    ),
    dtype=np.float64
).ravel()


P_lib[
    P_lib == 0
] = 1.0

V_lib[
    V_lib == 0
] = 1.0


# ============================================================
# SPATIAL GRAPH
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

print()
print(
    "Visium graph spots:",
    len(xy)
)

print(
    "estimated pitch:",
    f"{pitch:.2f} um"
)


# ============================================================
# PARALLEL CORRESPONDENCE + MORAN
# ============================================================

cosine = np.full(
    ng,
    np.nan,
    dtype=float
)

spearman = np.full(
    ng,
    np.nan,
    dtype=float
)

moran = np.full(
    ng,
    np.nan,
    dtype=float
)

moran_p = np.full(
    ng,
    np.nan,
    dtype=float
)


metric_bounds = [
    (
        start,
        min(
            start + METRIC_CHUNK,
            ng
        )
    )
    for start in range(
        0,
        ng,
        METRIC_CHUNK
    )
]


def metric_chunk(bounds):

    s0, t0 = bounds

    local_ok = gene_ok[
        s0:t0
    ]

    pos = np.flatnonzero(
        local_ok
    )

    if len(pos) == 0:

        return (
            s0,
            t0,
            np.array([], dtype=int),
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )


    P = P_full[
        :,
        s0:t0
    ][
        :,
        local_ok
    ].astype(
        np.float64,
        copy=False
    )

    V = V_full[
        :,
        s0:t0
    ][
        :,
        local_ok
    ].astype(
        np.float64,
        copy=False
    )


    # Same full-transcriptome CPM denominator for every chunk.
    Pcpm = (
        P
        /
        P_lib[
            :,
            None
        ]
        *
        1e6
    )

    Vcpm = (
        V
        /
        V_lib[
            :,
            None
        ]
        *
        1e6
    )


    # Same correspondence definition as previous script.
    Plog = np.log10(
        Pcpm + 1
    )

    Vlog = np.log10(
        Vcpm + 1
    )


    cs = cosine_cols(
        Plog,
        Vlog
    )

    sp = spearman_cols(
        Plog,
        Vlog
    )

    mi, mp = moran_chunk(
        Vcpm,
        W
    )


    return (
        s0,
        t0,
        pos,
        cs,
        sp,
        mi,
        mp
    )


print()
print(
    f"Computing correspondence + Moran "
    f"using {N_JOBS} threads..."
)

done = 0

with ThreadPoolExecutor(
    max_workers=N_JOBS
) as ex:

    futures = {
        ex.submit(
            metric_chunk,
            b
        ): b
        for b in metric_bounds
    }

    for fut in as_completed(
        futures
    ):

        (
            s0,
            t0,
            pos,
            cs,
            sp,
            mi,
            mp
        ) = fut.result()


        if len(pos):

            gp = (
                s0
                +
                pos
            )

            cosine[
                gp
            ] = cs

            spearman[
                gp
            ] = sp

            moran[
                gp
            ] = mi

            moran_p[
                gp
            ] = mp


        done += (
            t0 - s0
        )

        print(
            f"  metrics: "
            f"{done:5d}/{ng}",
            end="\r",
            flush=True
        )

print()


# ============================================================
# BH ACROSS EVALUABLE GENES
# ============================================================

moran_q = np.full(
    ng,
    np.nan
)

moran_q[
    gene_ok
] = bh(
    moran_p[
        gene_ok
    ]
)


patterned = (
    gene_ok
    &
    np.isfinite(
        moran_q
    )
    &
    (
        moran_q < 0.05
    )
    &
    (
        moran > 0
    )
)


# ============================================================
# OUTPUT TABLE
# ============================================================

res = pd.DataFrame({
    "gene":
        shared_genes,

    "measured_in_xenium":
        measured_flag,

    "expression_source":
        np.where(
            measured_flag,
            "measured_xenium",
            "envi_imputed"
        ),

    "visium_total":
        v_total,

    "visium_detected_spots":
        v_detect,

    "visium_detection_fraction":
        (
            v_detect
            /
            occupied.sum()
        ),

    "visium_max_count":
        v_max,

    "envi_pseudospot_total":
        p_total,

    "evaluable":
        gene_ok,

    "visium_highconf":
        visium_highconf,

    "moran_I":
        moran,

    "moran_p":
        moran_p,

    "moran_q":
        moran_q,

    "spatially_patterned":
        patterned,

    "cosine":
        cosine,

    "spearman":
        spearman,
})


GENE_OUT = (
    OUT
    / f"{BR}_ENVI_fulltranscriptome_gene_metrics.csv.gz"
)

res.to_csv(
    GENE_OUT,
    index=False,
    compression="gzip"
)


# ============================================================
# SUMMARY
# ============================================================

def med(x):

    x = np.asarray(
        x,
        float
    )

    x = x[
        np.isfinite(
            x
        )
    ]

    if len(x) == 0:
        return np.nan

    return float(
        np.median(
            x
        )
    )


meas = (
    res["measured_in_xenium"]
    &
    res["evaluable"]
)

imp = (
    ~res["measured_in_xenium"]
    &
    res["evaluable"]
)

meas_pat = (
    meas
    &
    res["spatially_patterned"]
)

imp_pat = (
    imp
    &
    res["spatially_patterned"]
)


summary = pd.DataFrame([
    {
        "BrNum":
            BR,

        "n_envi_cells":
            e.n_obs,

        "n_stalign_cells":
            len(aa),

        "n_cells_used":
            len(cell_rows),

        "n_visium_spots":
            n_spots,

        "n_occupied_pseudospots":
            int(
                occupied.sum()
            ),

        "pseudospot_coverage":
            float(
                occupied.mean()
            ),

        "n_envi_genes":
            e.n_vars,

        "n_shared_genes":
            ng,

        "n_evaluable_genes":
            int(
                gene_ok.sum()
            ),

        "n_measured_evaluable":
            int(
                meas.sum()
            ),

        "n_imputed_evaluable":
            int(
                imp.sum()
            ),

        "n_patterned_total":
            int(
                patterned.sum()
            ),

        "n_patterned_measured":
            int(
                meas_pat.sum()
            ),

        "n_patterned_imputed":
            int(
                imp_pat.sum()
            ),

        "median_cosine_measured":
            med(
                res.loc[
                    meas,
                    "cosine"
                ]
            ),

        "median_cosine_imputed":
            med(
                res.loc[
                    imp,
                    "cosine"
                ]
            ),

        "median_cosine_patterned_measured":
            med(
                res.loc[
                    meas_pat,
                    "cosine"
                ]
            ),

        "median_cosine_patterned_imputed":
            med(
                res.loc[
                    imp_pat,
                    "cosine"
                ]
            ),

        "median_spearman_patterned_measured":
            med(
                res.loc[
                    meas_pat,
                    "spearman"
                ]
            ),

        "median_spearman_patterned_imputed":
            med(
                res.loc[
                    imp_pat,
                    "spearman"
                ]
            ),

        "n_highconf_imputed":
            int(
                (
                    imp
                    &
                    res["visium_highconf"]
                ).sum()
            ),

        "median_cosine_highconf_imputed":
            med(
                res.loc[
                    imp
                    &
                    res["visium_highconf"],
                    "cosine"
                ]
            ),
    }
])


SUM_OUT = (
    OUT
    / f"{BR}_ENVI_fulltranscriptome_summary.csv"
)

summary.to_csv(
    SUM_OUT,
    index=False
)


print()
print(
    "=" * 80
)

print(
    "RESULTS"
)

print(
    "=" * 80
)

print(
    summary.T.to_string(
        header=False
    )
)

print()
print(
    "Gene metrics:",
    GENE_OUT
)

print(
    "Summary     :",
    SUM_OUT
)


e.file.close()
