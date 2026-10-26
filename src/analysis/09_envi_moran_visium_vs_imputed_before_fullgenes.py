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

# Same graph for Visium and ENVI pseudospots
MORAN_K = int(
    os.environ.get(
        "MORAN_K",
        "6"
    )
)

GENE_CHUNK = 512


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

RESULT_DIR = (
    ROOT
    / "results/analysis/envi_fulltranscriptome_crosstech"
    / BR
)

METRICS_GZ = (
    RESULT_DIR
    / f"{BR}_ENVI_fulltranscriptome_gene_metrics.csv.gz"
)

METRICS_CSV = (
    RESULT_DIR
    / f"{BR}_ENVI_fulltranscriptome_gene_metrics.csv"
)

VRAW = (
    ROOT
    / "data/raw/visium"
    / BR
)

OUTDIR = (
    ROOT
    / "results/analysis/envi_moran_visium_vs_imputed"
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
        .isin([
            "true",
            "1",
            "t",
            "yes",
        ])
        .to_numpy(bool)
    )


def pick_one(path, pattern):

    hits = sorted(
        path.glob(pattern)
    )

    if len(hits) != 1:
        raise RuntimeError(
            f"{BR}: expected exactly one '{pattern}' "
            f"under {path}; found {len(hits)}"
        )

    return hits[0]


def build_knn_weights(xy, k=6):
    """
    Directed kNN, row-standardized.

    The important point for this analysis is that the SAME W
    is used for both Visium and ENVI pseudospot expression.
    """

    xy = np.asarray(
        xy,
        dtype=float
    )

    n = len(xy)

    if n < 3:
        raise RuntimeError(
            "Too few spots for Moran's I."
        )

    k = min(
        k,
        n - 1
    )

    tree = cKDTree(
        xy
    )

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
        dtype=np.float64
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
            n,
            n
        )
    )

    return W


def moran_matrix(X, W):
    """
    X: spots x genes

    Moran's I computed column-wise using one shared W.
    """

    X = np.asarray(
        X,
        dtype=np.float64
    )

    n = X.shape[0]

    Z = (
        X
        -
        np.mean(
            X,
            axis=0,
            keepdims=True
        )
    )

    denominator = np.sum(
        Z * Z,
        axis=0
    )

    WZ = W @ Z

    numerator = np.sum(
        Z * WZ,
        axis=0
    )

    S0 = float(
        W.sum()
    )

    out = np.full(
        X.shape[1],
        np.nan,
        dtype=float
    )

    ok = (
        denominator
        >
        0
    )

    out[ok] = (
        n
        /
        S0
    ) * (
        numerator[ok]
        /
        denominator[ok]
    )

    return out


# ============================================================
# CHECK INPUT
# ============================================================

for f in [
    ENVI_FILE,
    ALIGN_FILE,
    VCOORD_FILE,
]:

    if not f.exists():
        raise FileNotFoundError(f)


if METRICS_GZ.exists():

    METRICS_FILE = METRICS_GZ

elif METRICS_CSV.exists():

    METRICS_FILE = METRICS_CSV

else:

    raise FileNotFoundError(
        f"No Task-2 gene metrics file found for {BR}"
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


# ============================================================
# LOAD EXISTING TASK-2 METRICS
#
# We use this only to define the same evaluable gene universe.
# Moran's I itself is recomputed independently here for BOTH
# technologies using the same graph.
# ============================================================

metrics = pd.read_csv(
    METRICS_FILE
)

for c in [
    "measured_in_xenium",
    "evaluable",
    "spatially_patterned",
]:

    if c in metrics.columns:
        metrics[c] = as_bool(
            metrics[c]
        )


if "gene" not in metrics.columns:
    raise RuntimeError(
        "Expected a 'gene' column in Task-2 metrics."
    )


# ============================================================
# LOAD ENVI
# ============================================================

e = ad.read_h5ad(
    ENVI_FILE,
    backed="r"
)

print("=" * 90)
print(f"MORAN I CROSS-TECHNOLOGY ANALYSIS: {BR}")
print("=" * 90)

print(
    "ENVI shape:",
    e.shape
)


egenes = pd.Index(
    e.var_names.astype(str)
)


# ============================================================
# LOAD STALIGN
# ============================================================

aa = pd.read_csv(
    ALIGN_FILE,
    low_memory=False
)

vcoord = pd.read_csv(
    VCOORD_FILE
)


required_align = [
    "cell_id",
    "match_prob",
    "x_aligned",
    "y_aligned",
    "visium_spot",
]

missing = [
    c
    for c in required_align
    if c not in aa.columns
]

if missing:
    raise RuntimeError(
        f"Missing STalign columns: {missing}"
    )


if "cell_id" not in e.obs.columns:
    raise RuntimeError(
        "ENVI obs does not contain cell_id."
    )


# ============================================================
# MAP ENVI CELLS TO STALIGN
# ============================================================

lookup = {
    str(x): i
    for i, x in enumerate(
        aa["cell_id"]
        .astype(str)
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


assigned_rows = np.flatnonzero(
    assigned
)


print(
    "matched cells inside physical Visium spots:",
    len(assigned_rows)
)


# ============================================================
# TARGET VISIUM SPOTS
# ============================================================

if "barcode" not in vcoord.columns:
    raise RuntimeError(
        "Expected barcode column in Visium coordinates."
    )


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


if (
    spot_idx < 0
).any():

    raise RuntimeError(
        "Some STalign visium_spot values are absent "
        "from target Visium coordinates."
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
    int(
        occupied.sum()
    ),
    "/",
    n_spots
)


# ============================================================
# SAME COORDINATES FOR BOTH TECHNOLOGIES
# ============================================================

xy = (
    vcoord.loc[
        occupied,
        [
            "x",
            "y"
        ]
    ]
    .to_numpy(float)
)


print(
    "Moran graph spots:",
    len(xy)
)

print(
    "kNN:",
    MORAN_K
)


W = build_knn_weights(
    xy,
    k=MORAN_K
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
        f"Unexpected Visium matrix shape: {VM.shape}"
    )


raw_bar_lookup = {
    b: i
    for i, b in enumerate(
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
        f"{len(missing_barcodes)} target Visium barcodes "
        "absent from raw matrix."
    )


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
# SAME EVALUABLE UNIVERSE AS TASK 2
# ============================================================

eval_mask = metrics[
    "evaluable"
].to_numpy(bool)


eval_genes = (
    metrics.loc[
        eval_mask,
        "gene"
    ]
    .astype(str)
    .tolist()
)


eval_genes = [
    g
    for g in eval_genes
    if (
        g in egenes
        and
        g in vgene_lookup
    )
]


# Imputed-only genes are the actual genes used in the plot.
analysis_mask = (
    metrics["evaluable"].to_numpy(bool)
    &
    (
        ~metrics[
            "measured_in_xenium"
        ].to_numpy(bool)
    )
)


analysis_genes = (
    metrics.loc[
        analysis_mask,
        "gene"
    ]
    .astype(str)
    .tolist()
)


analysis_genes = [
    g
    for g in analysis_genes
    if (
        g in egenes
        and
        g in vgene_lookup
    )
]


print(
    "evaluable genes for CPM denominator:",
    len(eval_genes)
)

print(
    "evaluable imputed-only genes for Moran comparison:",
    len(analysis_genes)
)


# ============================================================
# LOAD ASSIGNED ENVI CELLS
#
# This is the same strategy used in the PI-style plotting
# workflow. It avoids repeated backed HDF5 reads.
# ============================================================

print()
print(
    "Loading assigned-cell ENVI count_scale matrix..."
)


E_assigned = e.layers[
    "count_scale"
][
    assigned_rows,
    :
]


if sparse.issparse(
    E_assigned
):

    E_assigned = (
        E_assigned
        .toarray()
    )


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
# FIRST PASS:
# LIBRARY SIZES USING SAME EVALUABLE GENE UNIVERSE
# ============================================================

P_lib = np.zeros(
    n_spots,
    dtype=np.float64
)

V_lib = np.zeros(
    n_spots,
    dtype=np.float64
)


eval_eidx = np.array([
    egenes.get_loc(g)
    for g in eval_genes
])

eval_vidx = np.array([
    vgene_lookup[g]
    for g in eval_genes
])


print()
print(
    "Computing matched pseudospot / Visium library sizes..."
)


for s in range(
    0,
    len(eval_genes),
    GENE_CHUNK
):

    t = min(
        s + GENE_CHUNK,
        len(eval_genes)
    )

    eidx = eval_eidx[
        s:t
    ]

    vidx = eval_vidx[
        s:t
    ]


    Echunk = E_assigned[
        :,
        eidx
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
        vidx
    ]


    V_lib += np.asarray(
        Vchunk.sum(
            axis=1
        )
    ).ravel()


P_lib[
    P_lib <= 0
] = 1.0

V_lib[
    V_lib <= 0
] = 1.0


# ============================================================
# SECOND PASS:
# GENE-WISE MORAN'S I
# ============================================================

analysis_eidx = np.array([
    egenes.get_loc(g)
    for g in analysis_genes
])

analysis_vidx = np.array([
    vgene_lookup[g]
    for g in analysis_genes
])


moran_vis = np.full(
    len(analysis_genes),
    np.nan
)

moran_envi = np.full(
    len(analysis_genes),
    np.nan
)


print()
print(
    "Computing gene-level Moran's I..."
)


for s in range(
    0,
    len(analysis_genes),
    GENE_CHUNK
):

    t = min(
        s + GENE_CHUNK,
        len(analysis_genes)
    )


    Echunk = E_assigned[
        :,
        analysis_eidx[
            s:t
        ]
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


    Vchunk = (
        VT[
            :,
            analysis_vidx[
                s:t
            ]
        ]
        .toarray()
        .astype(
            np.float64
        )
    )


    # Same spots.
    Pchunk = Pchunk[
        occupied,
        :
    ]

    Vchunk = Vchunk[
        occupied,
        :
    ]


    # Same type of normalization.
    Pcpm = (
        Pchunk
        /
        P_lib[
            occupied,
            None
        ]
        *
        1e6
    )


    Vcpm = (
        Vchunk
        /
        V_lib[
            occupied,
            None
        ]
        *
        1e6
    )


    # Same transformation.
    Pexpr = np.log1p(
        Pcpm
    )

    Vexpr = np.log1p(
        Vcpm
    )


    moran_envi[
        s:t
    ] = moran_matrix(
        Pexpr,
        W
    )


    moran_vis[
        s:t
    ] = moran_matrix(
        Vexpr,
        W
    )


    if (
        s == 0
        or
        t == len(
            analysis_genes
        )
        or
        t % 5000 < GENE_CHUNK
    ):

        print(
            f"  {t:,} / "
            f"{len(analysis_genes):,}"
        )


# ============================================================
# BUILD RESULT TABLE
# ============================================================

out = pd.DataFrame({
    "gene":
        analysis_genes,

    "moran_visium":
        moran_vis,

    "moran_envi_imputed":
        moran_envi,
})


out[
    "delta_moran_envi_minus_visium"
] = (
    out[
        "moran_envi_imputed"
    ]
    -
    out[
        "moran_visium"
    ]
)


# Bring previous Visium-defined patterned flag along
# ONLY for visualization.
mmeta = (
    metrics
    .set_index("gene")
)


out[
    "visium_patterned_task2"
] = [
    (
        bool(
            mmeta.loc[
                g,
                "spatially_patterned"
            ]
        )
        if (
            g in mmeta.index
            and
            "spatially_patterned"
            in mmeta.columns
        )
        else False
    )
    for g in analysis_genes
]


OUT_GENE = (
    OUTDIR
    / f"{BR}_Visium_vs_ENVI_MoranI_imputed_genes.csv.gz"
)


out.to_csv(
    OUT_GENE,
    index=False,
    compression="gzip"
)


# ============================================================
# CORRELATION BETWEEN THE TWO MORAN VALUES
# ============================================================

ok = (
    np.isfinite(
        out["moran_visium"]
    )
    &
    np.isfinite(
        out["moran_envi_imputed"]
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
        "moran_envi_imputed"
    ]
    .to_numpy(float)
)


rho, rho_p = spearmanr(
    x,
    y
)

r, r_p = pearsonr(
    x,
    y
)


delta = (
    y
    -
    x
)


# Regression is descriptive only.
if len(x) >= 2:

    slope, intercept = np.polyfit(
        x,
        y,
        1
    )

else:

    slope = np.nan
    intercept = np.nan


summary = pd.DataFrame([
    {
        "BrNum":
            BR,

        "n_assigned_cells":
            len(
                assigned_rows
            ),

        "n_occupied_spots":
            int(
                occupied.sum()
            ),

        "n_imputed_evaluable_genes":
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
            np.nanmedian(
                x
            ),

        "median_envi_moran":
            np.nanmedian(
                y
            ),

        "median_delta_moran":
            np.nanmedian(
                delta
            ),

        "fraction_envi_moran_gt_visium":
            np.mean(
                y > x
            ),
    }
])


OUT_SUMMARY = (
    OUTDIR
    / f"{BR}_Visium_vs_ENVI_MoranI_summary.csv"
)


summary.to_csv(
    OUT_SUMMARY,
    index=False
)


# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(
        7.4,
        7.1
    )
)


patterned = (
    out.loc[
        ok,
        "visium_patterned_task2"
    ]
    .to_numpy(bool)
)


# Non-patterned background
ax.scatter(
    x[
        ~patterned
    ],
    y[
        ~patterned
    ],
    s=9,
    alpha=0.20,
    linewidths=0,
    label="Not Visium-patterned"
)


# Visium-patterned genes
ax.scatter(
    x[
        patterned
    ],
    y[
        patterned
    ],
    s=12,
    alpha=0.65,
    linewidths=0,
    label="Visium-patterned"
)


# Same limits on x and y.
allv = np.concatenate([
    x,
    y
])


lo = float(
    np.nanmin(
        allv
    )
)

hi = float(
    np.nanmax(
        allv
    )
)


pad = max(
    0.03,
    0.04
    *
    (
        hi
        -
        lo
    )
)


lo -= pad
hi += pad


ax.plot(
    [
        lo,
        hi
    ],
    [
        lo,
        hi
    ],
    linestyle="--",
    linewidth=1.1,
    label="y = x"
)


# Descriptive regression line
xx = np.linspace(
    lo,
    hi,
    200
)


ax.plot(
    xx,
    intercept
    +
    slope
    *
    xx,
    linewidth=1.2,
    alpha=0.85,
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
    "ENVI-imputed Xenium pseudospot Moran's I",
    fontsize=12
)


ax.set_title(
    f"{BR}: preservation of spatial autocorrelation",
    fontsize=13,
    fontweight="bold"
)


ax.text(
    0.03,
    0.97,
    f"Imputed-only genes: n={len(x):,}\n"
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
    fontsize=9,
    loc="lower right"
)


fig.tight_layout()


PNG = (
    OUTDIR
    / f"fig_{BR}_Visium_vs_ENVI_MoranI.png"
)

PDF = (
    OUTDIR
    / f"fig_{BR}_Visium_vs_ENVI_MoranI.pdf"
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
print("RESULT")
print("=" * 90)

print(
    f"Genes               : {len(x):,}"
)

print(
    f"Spearman Moran      : {rho:.4f}"
)

print(
    f"Pearson Moran       : {r:.4f}"
)

print(
    f"Median Visium I     : {np.nanmedian(x):.4f}"
)

print(
    f"Median ENVI I       : {np.nanmedian(y):.4f}"
)

print(
    f"Median delta I      : {np.nanmedian(delta):+.4f}"
)

print(
    f"Regression slope    : {slope:.4f}"
)

print()
print(
    "Gene table:",
    OUT_GENE
)

print(
    "Summary   :",
    OUT_SUMMARY
)

print(
    "PNG       :",
    PNG
)

print(
    "PDF       :",
    PDF
)
