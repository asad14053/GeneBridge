#!/usr/bin/env python3

import os
import re
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.spatial import Delaunay
from scipy.stats import norm, spearmanr, mannwhitneyu


# ============================================================
# CONFIG
# ============================================================

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BR = "Br8667"

MATCH_THRESH = 0.85
SPOT_PITCH = 100.0

# Prefer the unfiltered measured Xenium object.
# Fall back to the imputation target if necessary.
X_CANDIDATES = [
    ROOT / f"data/processed/xenium/raw_by_donor/{BR}.h5ad",
    ROOT / f"data/processed/imputation_full/ex1_ntc/targets/"
           f"spatial_data_xenium_{BR}_ex1_ntc.h5ad",
]

VFILE = (
    ROOT
    / f"data/processed/imputation_beta/{BR}/"
      f"spatial_data_huuki_visium_{BR}.h5ad"
)

STDIR = (
    ROOT
    / "data/metadata/stalign_results_update"
    / f"{BR}_STalign_result"
)

ACSV = STDIR / f"{BR}_aligned_coordinates.csv"
VCSV = STDIR / f"{BR}_visium_coordinates.csv"

OUT = (
    ROOT
    / "results/analysis/crosstech_measured300"
    / BR
)

OUT.mkdir(parents=True, exist_ok=True)


# ============================================================
# HELPERS
# ============================================================

def get_count_matrix(adata):
    """
    Prefer raw counts stored in layers['counts'].
    Fall back to X only when counts layer is unavailable.
    """
    if "counts" in adata.layers:
        return adata.layers["counts"], "layers['counts']"
    return adata.X, "X"


def dense(M):
    if sparse.issparse(M):
        return M.toarray()
    return np.asarray(M)


def dense_col(M, j):
    z = M[:, j]
    if sparse.issparse(z):
        z = z.toarray()
    return np.asarray(z).ravel()


def cpm(M):
    M = np.asarray(M, dtype=float)

    lib = M.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0

    return M / lib * 1e6


def log_cpm(M):
    return np.log10(cpm(M) + 1.0)


def cosine_columns(A, B):
    num = np.sum(A * B, axis=0)

    den = (
        np.linalg.norm(A, axis=0)
        *
        np.linalg.norm(B, axis=0)
    )

    out = np.full(A.shape[1], np.nan, dtype=float)

    good = den > 0
    out[good] = num[good] / den[good]

    return out


def bh(p):
    """Benjamini-Hochberg, matching advisor notebook."""
    p = np.asarray(p, dtype=float)
    n = len(p)

    order = np.argsort(p)

    q = np.empty(n, dtype=float)

    tmp = (
        p[order]
        * n
        / np.arange(1, n + 1)
    )

    tmp = np.minimum.accumulate(
        tmp[::-1]
    )[::-1]

    q[order] = tmp

    return np.clip(q, 0, 1)


def spatial_weights(xy, filter_dist):
    """
    MERINGUE-style:
    Delaunay graph, then remove edges > filter_dist.
    """

    tri = Delaunay(xy)

    edges = set()

    for s in tri.simplices:
        for a in range(3):
            for b in range(a + 1, 3):

                i, j = sorted(
                    (
                        int(s[a]),
                        int(s[b])
                    )
                )

                edges.add((i, j))

    E = np.asarray(
        sorted(edges),
        dtype=int
    )

    dist = np.linalg.norm(
        xy[E[:, 0]]
        -
        xy[E[:, 1]],
        axis=1
    )

    E = E[
        dist <= filter_dist
    ]

    n = len(xy)

    rows = np.r_[
        E[:, 0],
        E[:, 1]
    ]

    cols = np.r_[
        E[:, 1],
        E[:, 0]
    ]

    vals = np.ones(
        len(rows),
        dtype=float
    )

    W = sparse.csr_matrix(
        (vals, (rows, cols)),
        shape=(n, n)
    )

    return W


def moran_test(x, W, S0, S1, S2):
    """
    Same randomization-assumption Moran test used
    in the advisor notebook.
    """

    x = np.asarray(x, dtype=float)

    n = len(x)

    z = x - x.mean()

    zz = np.sum(z ** 2)

    if zz == 0:
        return 0.0, 1.0

    I = (
        (n / S0)
        * np.dot(
            z,
            W.dot(z)
        )
        / zz
    )

    EI = -1.0 / (n - 1)

    b2 = (
        n
        * np.sum(z ** 4)
        / zz ** 2
    )

    A = n * (
        (n**2 - 3*n + 3) * S1
        - n * S2
        + 3 * S0**2
    )

    B = b2 * (
        (n**2 - n) * S1
        - 2*n*S2
        + 6*S0**2
    )

    VI = (
        (A - B)
        /
        (
            (n - 1)
            * (n - 2)
            * (n - 3)
            * S0**2
        )
        -
        EI**2
    )

    if VI <= 0:
        return I, 1.0

    zscore = (
        I - EI
    ) / np.sqrt(VI)

    # one-sided positive spatial autocorrelation
    p = float(
        1 - norm.cdf(zscore)
    )

    return I, p


def n01(v):
    """1st–99th percentile scaling for visualization only."""

    v = np.asarray(v, dtype=float).ravel()

    lo, hi = np.percentile(
        v,
        [1, 99]
    )

    return np.clip(
        (v - lo)
        /
        (hi - lo + 1e-9),
        0,
        1
    )


# ============================================================
# LOAD STALIGN EXPORT
# ============================================================

print("=" * 80)
print("STALIGN")
print("=" * 80)

align = pd.read_csv(
    ACSV,
    low_memory=False
)

vcoord = pd.read_csv(
    VCSV,
    low_memory=False
)

align["cell_id"] = (
    align["cell_id"]
    .astype(str)
)

align["lieber_key"] = (
    align["lieber_key"]
    .astype(str)
)

vcoord["barcode"] = (
    vcoord["barcode"]
    .astype(str)
)

print("aligned cells :", len(align))
print("Visium spots  :", len(vcoord))


# ============================================================
# LOAD XENIUM — CHOOSE FILE WITH BEST ID OVERLAP
# ============================================================

print()
print("=" * 80)
print("XENIUM SOURCE SELECTION")
print("=" * 80)

best = None

for path in X_CANDIDATES:

    if not path.exists():
        continue

    xx = ad.read_h5ad(path)

    candidate_ids = {
        "obs_names":
            pd.Series(
                xx.obs_names.astype(str),
                index=np.arange(xx.n_obs)
            )
    }

    for col in [
        "cell_id",
        "original_cell_id",
        "donor_cell_id"
    ]:
        if col in xx.obs.columns:
            candidate_ids[col] = (
                xx.obs[col]
                .astype(str)
                .reset_index(drop=True)
            )

    for src_name, ids in candidate_ids.items():

        ids_set = set(ids)

        for align_col in [
            "cell_id",
            "lieber_key"
        ]:

            score = len(
                ids_set
                &
                set(
                    align[
                        align_col
                    ].astype(str)
                )
            )

            print(
                f"{path.name:55s} "
                f"{src_name:18s} -> "
                f"{align_col:12s}: "
                f"{score}"
            )

            rec = (
                score,
                path,
                xx,
                src_name,
                align_col,
                ids
            )

            if (
                best is None
                or score > best[0]
            ):
                best = rec


if best is None:
    raise RuntimeError(
        "No Xenium H5AD found."
    )


(
    x_overlap,
    XFILE,
    x,
    X_ID_SOURCE,
    X_ALIGN_COL,
    x_ids,
) = best


print()
print("SELECTED XENIUM:")
print(" file       :", XFILE)
print(" shape      :", x.shape)
print(" ID source  :", X_ID_SOURCE)
print(" align col  :", X_ALIGN_COL)
print(" overlap    :", x_overlap)


if x_overlap < 0.90 * x.n_obs:
    warnings.warn(
        "Xenium/STalign overlap <90%. "
        "Inspect identifiers."
    )


# unique ID -> row
x_id_table = pd.DataFrame({
    "id": x_ids.astype(str).values,
    "row": np.arange(x.n_obs)
})

x_id_table = (
    x_id_table
    .drop_duplicates(
        "id",
        keep="first"
    )
)

x_to_row = pd.Series(
    x_id_table["row"].values,
    index=x_id_table["id"].values
)


align["_xrow"] = (
    align[X_ALIGN_COL]
    .astype(str)
    .map(x_to_row)
)

has_x = align["_xrow"].notna()

print(
    "STalign cells with expression:",
    has_x.sum(),
    "/",
    len(align)
)


XMAT, x_matrix_source = (
    get_count_matrix(x)
)

print(
    "Xenium count source:",
    x_matrix_source
)


# ============================================================
# VISIUM — IDENTIFY CORRECT CAPTURE USING BARCODE + PIXEL GEOMETRY
#
# IMPORTANT:
# 10x barcodes can repeat across capture areas.
# Therefore barcode overlap alone is not sufficient.
# The correct capture must also reproduce the STalign
# full-resolution pixel coordinates x_px / y_px.
# ============================================================

print()
print("=" * 80)
print("VISIUM CAPTURE SELECTION")
print("=" * 80)

v = ad.read_h5ad(VFILE)

print("Visium h5ad:", VFILE)
print("shape       :", v.shape)

if "spatial" not in v.obsm:
    raise RuntimeError("Visium h5ad has no obsm['spatial'].")

if "sample_id" not in v.obs.columns:
    raise RuntimeError("Visium h5ad has no sample_id column.")

vxy_px = np.asarray(v.obsm["spatial"], dtype=float)

obs_tmp = pd.DataFrame({
    "vrow": np.arange(v.n_obs),
    "obs_name": v.obs_names.astype(str),
    "sample_id": v.obs["sample_id"].astype(str).values,
    "x_px_h5": vxy_px[:,0],
    "y_px_h5": vxy_px[:,1],
})

# Example:
# AAACAAGTATCTCCCA-1.24 -> AAACAAGTATCTCCCA-1
obs_tmp["barcode"] = (
    obs_tmp["obs_name"]
    .str.replace(r"\.\d+$", "", regex=True)
)

vc_geom = vcoord[
    ["barcode","x_px","y_px"]
].copy()

vc_geom["barcode"] = vc_geom["barcode"].astype(str)

scores = []

for sample, d in obs_tmp.groupby("sample_id"):

    m = d.merge(
        vc_geom,
        on="barcode",
        how="inner"
    )

    if len(m) == 0:
        continue

    dist = np.sqrt(
        (m["x_px_h5"] - m["x_px"])**2 +
        (m["y_px_h5"] - m["y_px"])**2
    )

    scores.append({
        "sample_id": sample,
        "n_sample_spots": len(d),
        "n_barcode_overlap": len(m),
        "median_px_error": float(np.median(dist)),
        "p95_px_error": float(np.percentile(dist,95)),
    })

score_df = (
    pd.DataFrame(scores)
    .sort_values(
        ["median_px_error","n_barcode_overlap"],
        ascending=[True,False]
    )
)

print(score_df.to_string(index=False))

if score_df.empty:
    raise RuntimeError("No Visium sample overlaps STalign barcodes.")

best = score_df.iloc[0]

V_SAMPLE = str(best["sample_id"])

print()
print("SELECTED VISIUM CAPTURE:", V_SAMPLE)
print("barcode overlap         :", int(best["n_barcode_overlap"]))
print("median pixel error      :", best["median_px_error"])
print("p95 pixel error         :", best["p95_px_error"])

# Correct capture should match STalign pixel geometry almost exactly.
if best["median_px_error"] > 2:
    raise RuntimeError(
        "Best Visium capture still has >2 px median geometry error. "
        "Do not continue."
    )

sub = obs_tmp[
    obs_tmp["sample_id"] == V_SAMPLE
].copy()

# Map using barcode WITHIN the now-validated capture.
sub = (
    sub
    .sort_values("vrow")
    .drop_duplicates("barcode")
)

v_to_row = pd.Series(
    sub["vrow"].values,
    index=sub["barcode"].values
)

spot_present = vcoord["barcode"].astype(str).isin(v_to_row.index)

vcoord = (
    vcoord.loc[spot_present]
    .reset_index(drop=True)
)

vrows = (
    vcoord["barcode"]
    .astype(str)
    .map(v_to_row)
    .astype(int)
    .to_numpy()
)

print("usable target spots      :", len(vcoord))

if len(vcoord) < 0.95 * 4800:
    raise RuntimeError(
        f"Only {len(vcoord)} / 4800 STalign target spots mapped."
    )

VMAT, v_matrix_source = get_count_matrix(v)

print("Visium count source      :", v_matrix_source)


# ============================================================
# SHARED GENES
# ============================================================

print()
print("=" * 80)
print("SHARED GENES")
print("=" * 80)

shared = sorted(
    set(
        x.var_names.astype(str)
    )
    &
    set(
        v.var_names.astype(str)
    )
)

print(
    "shared genes:",
    len(shared)
)


if len(shared) < 250:
    raise RuntimeError(
        f"Only {len(shared)} shared genes. "
        "Expected ~300."
    )


xi = x.var_names.get_indexer(
    shared
)

vi = v.var_names.get_indexer(
    shared
)


# ============================================================
# VISIUM COUNTS IN STALIGN SPOT ORDER
# ============================================================

V = dense(
    VMAT[vrows][:, vi]
).astype(float)


spot_to_i = pd.Series(
    np.arange(
        len(vcoord)
    ),
    index=vcoord[
        "barcode"
    ].astype(str)
)


# ============================================================
# XENIUM -> VISIUM PSEUDOSPOTS
#
# Use the exported final-STalign spot assignment.
#
# visium_spot is populated when the aligned cell centroid
# lies inside the physical Visium spot.
# ============================================================

aa = (
    align.loc[has_x]
    .copy()
)

aa["_xrow"] = (
    aa["_xrow"]
    .astype(int)
)


matched = (
    pd.to_numeric(
        aa["match_prob"],
        errors="coerce"
    )
    .fillna(0)
    .to_numpy()
    >
    MATCH_THRESH
)


assigned = (
    matched
    &
    aa["visium_spot"].notna().to_numpy()
    &
    aa["visium_spot"]
      .astype(str)
      .isin(
          spot_to_i.index
      )
      .to_numpy()
)


cell_use = (
    aa.loc[assigned]
    .copy()
)


xrows_use = (
    cell_use["_xrow"]
    .astype(int)
    .to_numpy()
)


sp_idx = (
    cell_use[
        "visium_spot"
    ]
    .astype(str)
    .map(
        spot_to_i
    )
    .astype(int)
    .to_numpy()
)


Xuse = dense(
    XMAT[
        xrows_use
    ][:, xi]
).astype(float)


P = np.zeros(
    (
        len(vcoord),
        len(shared)
    ),
    dtype=float
)


np.add.at(
    P,
    sp_idx,
    Xuse
)


ncell = np.bincount(
    sp_idx,
    minlength=len(vcoord)
)


occupied = (
    ncell >= 1
)


print()
print("=" * 80)
print("PSEUDOSPOTS")
print("=" * 80)

print(
    "cells with expression       :",
    len(aa)
)

print(
    f"matched > {MATCH_THRESH}          :",
    matched.sum()
)

print(
    "matched + inside spot       :",
    assigned.sum()
)

print(
    "occupied pseudospots        :",
    occupied.sum(),
    "/",
    len(vcoord),
    f"({occupied.mean()*100:.1f}%)"
)

if occupied.any():

    print(
        "cells/pseudospot median    :",
        f"{np.median(ncell[occupied]):.0f}"
    )

    print(
        "cells/pseudospot IQR       :",
        f"{np.percentile(ncell[occupied],25):.0f}"
        "–"
        f"{np.percentile(ncell[occupied],75):.0f}"
    )


if (
    "spot_dist_um"
    in cell_use.columns
):

    dd = pd.to_numeric(
        cell_use[
            "spot_dist_um"
        ],
        errors="coerce"
    )

    print(
        "assigned distance max (um) :",
        f"{dd.max():.2f}"
    )


# ============================================================
# GENE FILTER — SAME LOGIC AS ADVISOR CELL 9
# ============================================================

gene_ok = (
    (
        P[occupied]
        .sum(axis=0)
        >= 1
    )
    &
    (
        V[occupied]
        .sum(axis=0)
        >= 1
    )
    &
    (
        (
            V[occupied]
            > 0
        )
        .sum(axis=0)
        > 1
    )
)


genes = [
    g
    for g, good
    in zip(
        shared,
        gene_ok
    )
    if good
]


Pu = P[
    np.ix_(
        occupied,
        gene_ok
    )
]


Vu = V[
    np.ix_(
        occupied,
        gene_ok
    )
]


xy_u = (
    vcoord.loc[
        occupied,
        ["x", "y"]
    ]
    .to_numpy(dtype=float)
)


print()
print(
    "genes retained:",
    len(genes),
    "/",
    len(shared)
)


# ============================================================
# MORAN'S I — VISIUM
# ============================================================

print()
print("=" * 80)
print("VISIUM SPATIAL PATTERNING")
print("=" * 80)


W = spatial_weights(
    xy_u,
    filter_dist=1.5 * SPOT_PITCH
)


S0 = float(
    W.sum()
)


WT = (
    W
    +
    W.T
)


S1 = float(
    0.5
    *
    WT.multiply(WT).sum()
)


row_sum = np.asarray(
    W.sum(axis=1)
).ravel()

col_sum = np.asarray(
    W.sum(axis=0)
).ravel()


S2 = float(
    np.sum(
        (
            row_sum
            +
            col_sum
        ) ** 2
    )
)


print(
    "mean neighbours/spot:",
    f"{row_sum.mean():.2f}"
)


Vcpm = cpm(
    Vu
)


mor = np.array(
    [
        moran_test(
            Vcpm[:, j],
            W,
            S0,
            S1,
            S2
        )
        for j
        in range(
            len(genes)
        )
    ]
)


moran = mor[:, 0]
moran_p = mor[:, 1]

moran_q = bh(
    moran_p
)


patterned = (
    moran_q < 0.05
)


print(
    "spatially patterned:",
    patterned.sum(),
    "/",
    len(genes)
)


# ============================================================
# CROSS-TECHNOLOGY CORRESPONDENCE
# ============================================================

Pn = log_cpm(
    Pu
)

Vn = log_cpm(
    Vu
)


csim = cosine_columns(
    Pn,
    Vn
)


ssim = np.array(
    [
        (
            spearmanr(
                Pn[:, j],
                Vn[:, j]
            ).statistic
            if (
                np.std(
                    Pn[:, j]
                ) > 0
                and
                np.std(
                    Vn[:, j]
                ) > 0
            )
            else np.nan
        )
        for j
        in range(
            len(genes)
        )
    ]
)


good_pat = (
    patterned
    &
    np.isfinite(csim)
)

good_non = (
    (~patterned)
    &
    np.isfinite(csim)
)


if (
    good_pat.sum() > 0
    and
    good_non.sum() > 0
):

    mw = mannwhitneyu(
        csim[good_pat],
        csim[good_non],
        alternative="greater"
    )

    mw_p = mw.pvalue

else:

    mw_p = np.nan


print()
print(
    "patterned median cosine     :",
    f"{np.nanmedian(csim[patterned]):.3f}"
)

print(
    "non-patterned median cosine :",
    f"{np.nanmedian(csim[~patterned]):.3f}"
)

print(
    "patterned median Spearman   :",
    f"{np.nanmedian(ssim[patterned]):.3f}"
)

print(
    "Wilcoxon patterned > non    :",
    f"{mw_p:.3g}"
)


# ============================================================
# GENE METRICS
# ============================================================

metrics = pd.DataFrame({
    "gene": genes,
    "moran": moran,
    "moran_p": moran_p,
    "moran_q": moran_q,
    "patterned": patterned,
    "cosine": csim,
    "spearman": ssim,
    "visium_detection_fraction":
        (Vu > 0).mean(axis=0),
    "visium_max_count":
        Vu.max(axis=0),
})


metrics.to_csv(
    OUT
    /
    f"{BR}_measured300_gene_metrics.csv",
    index=False
)


# ============================================================
# SELECT HIGH / MEDIAN / LOW + NON-PATTERNED CONTROL
# ============================================================

sel = (
    metrics.copy()
)


sel = sel[
    ~sel["gene"].str.match(
        r"^(MT-|MTRNR|RP[LS])",
        case=True,
        na=False
    )
].copy()


NZ_MIN = 0.20
VMAX_MIN = 10


deep = (
    (
        sel[
            "visium_detection_fraction"
        ]
        >= NZ_MIN
    )
    &
    (
        sel[
            "visium_max_count"
        ]
        >= VMAX_MIN
    )
)


patd = (
    sel[
        sel["patterned"]
        &
        deep
    ]
    .sort_values(
        "cosine",
        ascending=False
    )
    .reset_index(drop=True)
)


nonp = (
    sel[
        (~sel["patterned"])
        &
        deep
    ]
    .sort_values(
        "cosine",
        ascending=False
    )
    .reset_index(drop=True)
)


if len(patd) < 3:
    raise RuntimeError(
        "Fewer than 3 eligible patterned genes."
    )


G = [
    patd.iloc[0]["gene"],
    patd.iloc[
        len(patd) // 2
    ]["gene"],
    patd.iloc[-1]["gene"],
]


GLAB = [
    "patterned, highest cosine",
    "patterned, median cosine",
    "patterned, lowest cosine",
]


if len(nonp) > 0:

    G.append(
        nonp.iloc[
            len(nonp) // 2
        ]["gene"]
    )

    GLAB.append(
        "not patterned (control)"
    )


print()
print("=" * 80)
print("DISPLAY GENES")
print("=" * 80)

print(
    f"depth requirement: "
    f"Visium detected >= {NZ_MIN*100:.0f}% "
    f"and max count >= {VMAX_MIN}"
)

print(
    "eligible patterned:",
    len(patd)
)

print(
    "eligible non-patterned:",
    len(nonp)
)


metrics_idx = (
    metrics
    .set_index("gene")
)


selected_rows = []

for label, g in zip(
    GLAB,
    G
):

    r = metrics_idx.loc[g]

    print(
        f"{label:28s} "
        f"{g:12s} "
        f"Moran={r['moran']:.3f} "
        f"q={r['moran_q']:.3g} "
        f"cos={r['cosine']:.3f} "
        f"rho={r['spearman']:+.3f}"
    )

    selected_rows.append({
        "class": label,
        "gene": g,
        "moran": r["moran"],
        "moran_q": r["moran_q"],
        "cosine": r["cosine"],
        "spearman": r["spearman"],
    })


pd.DataFrame(
    selected_rows
).to_csv(
    OUT
    /
    f"{BR}_selected_display_genes.csv",
    index=False
)


# ============================================================
# PREPARE CELL-LEVEL XENIUM CPM FOR PANEL b
# ============================================================

plot_cells = (
    aa.loc[matched]
    .copy()
)


plot_xrows = (
    plot_cells[
        "_xrow"
    ]
    .astype(int)
    .to_numpy()
)


plot_xy = (
    plot_cells[
        [
            "x_aligned",
            "y_aligned"
        ]
    ]
    .to_numpy(dtype=float)
)


Xplot = (
    XMAT[
        plot_xrows
    ]
)


xlib = np.asarray(
    Xplot.sum(axis=1)
).ravel()


xlib[
    xlib == 0
] = 1


# pseudospot CPM for panel c
Pcpm = cpm(
    Pu
)

Vcpm_all = cpm(
    Vu
)


# ============================================================
# FIGURE — SAME CONCEPT AS ATTACHMENT, WITHOUT SUBPLOT e
# ============================================================

RED = plt.cm.Reds

nrow = len(G)

fig = plt.figure(
    figsize=(
        16,
        3.6 * nrow
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


for i, (g, label) in enumerate(
    zip(
        G,
        GLAB
    )
):

    jg = genes.index(g)

    jx = (
        x.var_names
        .get_loc(g)
    )

    r = metrics_idx.loc[g]


    # --------------------------------------------------------
    # a — Visium
    # --------------------------------------------------------

    ax = fig.add_subplot(
        gs[i, 0]
    )

    ax.scatter(
        xy_u[:, 0] / 1000,
        xy_u[:, 1] / 1000,
        c=n01(
            Vcpm_all[:, jg]
        ),
        cmap=RED,
        s=9,
        linewidths=0,
        rasterized=True
    )

    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_ylabel(
        f"{g}\n{label}",
        fontsize=9.5,
        rotation=0,
        labelpad=60,
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
    # b — aligned Xenium
    # --------------------------------------------------------

    ax = fig.add_subplot(
        gs[i, 1]
    )


    xgene = (
        dense_col(
            Xplot,
            jx
        )
        /
        xlib
        *
        1e6
    )


    ax.scatter(
        plot_xy[:, 0] / 1000,
        plot_xy[:, 1] / 1000,
        c=n01(
            xgene
        ),
        cmap=RED,
        s=0.5,
        linewidths=0,
        rasterized=True
    )

    ax.set_aspect(
        "equal"
    )

    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])

    if i == 0:

        ax.set_title(
            "b   Xenium aligned",
            fontsize=11,
            loc="left",
            fontweight="bold"
        )


    # --------------------------------------------------------
    # c — normalized expression correspondence
    # --------------------------------------------------------

    ax = fig.add_subplot(
        gs[i, 2]
    )

    ax.scatter(
        Pcpm[:, jg],
        Vcpm_all[:, jg],
        s=9,
        c="#444444",
        alpha=0.45,
        linewidths=0,
        rasterized=True
    )

    ax.set_xlabel(
        "Xenium pseudospot (CPM)",
        fontsize=9
    )

    ax.set_ylabel(
        "Visium spot (CPM)",
        fontsize=9
    )

    ax.ticklabel_format(
        style="sci",
        scilimits=(0, 0),
        axis="both"
    )

    ax.text(
        0.04,
        0.93,
        (
            f"cosine {r['cosine']:.2f}\n"
            f"Spearman {r['spearman']:+.2f}"
        ),
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
# d — PATTERNED COSINE DISTRIBUTION ONLY
#
# No subplot e.
# ============================================================

rows_for_d = slice(
    0,
    min(
        3,
        nrow
    )
)


axd = fig.add_subplot(
    gs[
        rows_for_d,
        3
    ]
)


edges = np.linspace(
    0,
    1,
    26
)


vals = csim[
    patterned
    &
    np.isfinite(csim)
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


# Clifton cross-technology reference
CLIFTON = 0.55


axd.axhline(
    CLIFTON,
    color="#0072B2",
    lw=1.4,
    ls="--"
)


axd.text(
    0.98,
    CLIFTON + 0.02,
    f"Clifton {CLIFTON:.2f}",
    transform=axd.get_yaxis_transform(),
    ha="right",
    fontsize=8,
    color="#0072B2"
)


ours = np.nanmedian(
    csim[
        patterned
    ]
)


axd.axhline(
    ours,
    color="#333333",
    lw=1.4
)


axd.text(
    0.98,
    ours - 0.05,
    f"ours {ours:.2f}",
    transform=axd.get_yaxis_transform(),
    ha="right",
    fontsize=8
)


# mark only displayed PATTERNED genes
for g in G:

    rr = metrics_idx.loc[g]

    if not bool(
        rr["patterned"]
    ):
        continue

    c = rr["cosine"]

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
        va="center",
        fontsize=8.5,
        style="italic"
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
    f"d   spatially patterned (n={patterned.sum()})",
    fontweight="bold",
    loc="left",
    fontsize=9.5
)


# Bottom-right intentionally blank:
# subplot e REMOVED
if nrow >= 4:

    ax_blank = fig.add_subplot(
        gs[3, 3]
    )

    ax_blank.axis(
        "off"
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
    [0, 1]
)

cax.set_xticklabels(
    [
        "0",
        "high"
    ]
)

cax.set_title(
    "normalized expression (1st–99th pct)",
    fontsize=8
)


# ============================================================
# TITLE
# ============================================================

fig.suptitle(
    f"Cross-technology expression correspondence "
    f"(Xenium to Visium, {BR})",
    fontsize=13,
    fontweight="bold",
    y=1.0
)


fig.text(
    0.50,
    -0.025,
    (
        f"Final STalign Xenium cells with match probability > "
        f"{MATCH_THRESH} were assigned to their corresponding "
        f"Visium spot when the aligned cell centroid fell inside "
        f"the physical spot. Xenium counts were summed to form "
        f"one pseudospot per Visium spot. Both platforms were "
        f"CPM-normalized; cosine and Spearman were calculated "
        f"after log10(CPM+1). Spatial patterning was defined by "
        f"Visium Moran's I with BH-adjusted p < 0.05. "
        f"Rows show high-, median-, and low-cosine patterned "
        f"genes plus a non-patterned control. Subplot e is omitted."
    ),
    ha="center",
    va="top",
    fontsize=8.3,
    color="#444444",
    wrap=True
)


# ============================================================
# SAVE
# ============================================================

PNG = (
    OUT
    /
    f"fig_crosstech_{BR}_no_e.png"
)

PDF = (
    OUT
    /
    f"fig_crosstech_{BR}_no_e.pdf"
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

plt.close(
    fig
)


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame([
    {
        "BrNum": BR,
        "xenium_source": str(XFILE),
        "visium_sample": V_SAMPLE,
        "n_stalign_cells": len(align),
        "n_expression_cells": len(aa),
        "n_matched_cells": int(matched.sum()),
        "n_inside_spot": int(assigned.sum()),
        "n_visium_spots": len(vcoord),
        "n_occupied_pseudospots": int(occupied.sum()),
        "n_shared_genes": len(shared),
        "n_retained_genes": len(genes),
        "n_patterned_genes": int(patterned.sum()),
        "median_cosine_patterned":
            np.nanmedian(csim[patterned]),
        "median_cosine_nonpatterned":
            np.nanmedian(csim[~patterned]),
        "median_spearman_patterned":
            np.nanmedian(ssim[patterned]),
        "wilcoxon_p_patterned_gt_non":
            mw_p,
    }
])


summary.to_csv(
    OUT
    /
    f"{BR}_crosstech_summary.csv",
    index=False
)


print()
print("=" * 80)
print("DONE")
print("=" * 80)

print(
    "figure PNG :",
    PNG
)

print(
    "figure PDF :",
    PDF
)

print(
    "gene table :",
    OUT
    /
    f"{BR}_measured300_gene_metrics.csv"
)

print(
    "selection  :",
    OUT
    /
    f"{BR}_selected_display_genes.csv"
)

print(
    "summary    :",
    OUT
    /
    f"{BR}_crosstech_summary.csv"
)
