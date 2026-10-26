#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.spatial import cKDTree


# =============================================================================
# CONFIG
# =============================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BR = "Br8667"

MATCH_THRESH = 0.85
K = 6
CHUNK = 512


OOF = (
    ROOT
    / "outputs/imputation_full/ENVI_N24"
    / BR
    / "ex5_1/envi/oof_300/oof_predicted_300genes.h5ad"
)


FULL = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/envi"
    / BR
    / f"spatial_data_xenium_{BR}_ENVI_full_transcriptome.h5ad"
)


ALIGN = (
    ROOT
    / "data/metadata/stalign_results_update"
    / f"{BR}_STalign_result"
    / f"{BR}_aligned_coordinates.csv"
)


VIS_COORD = (
    ROOT
    / "data/metadata/stalign_results_update"
    / f"{BR}_STalign_result"
    / f"{BR}_visium_coordinates.csv"
)


OUT = (
    ROOT
    / "results/analysis/gene_variation"
    / BR
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


def get_cell_ids(adata):

    for c in [
        "original_cell_id",
        "cell_id",
        "donor_cell_id",
    ]:

        if c in adata.obs.columns:

            ids = (
                adata.obs[c]
                .astype(str)
                .to_numpy()
            )

            if len(np.unique(ids)) == len(ids):

                print(
                    "Using cell ID column:",
                    c
                )

                return ids

    ids = (
        adata.obs_names
        .astype(str)
        .to_numpy()
    )

    if len(np.unique(ids)) != len(ids):

        raise RuntimeError(
            "Could not find unique cell identifiers."
        )

    print(
        "Using obs_names as cell IDs."
    )

    return ids


def build_knn_W(
    xy,
    k=6
):

    n = len(
        xy
    )

    if n <= k:
        raise RuntimeError(
            f"Only {n} spots; cannot construct k={k} graph."
        )


    tree = cKDTree(
        xy
    )


    _, nbr = tree.query(
        xy,
        k=k + 1
    )


    nbr = nbr[
        :,
        1:
    ]


    rows = np.repeat(
        np.arange(
            n
        ),
        k
    )

    cols = nbr.reshape(
        -1
    )


    data = np.full(
        len(
            rows
        ),
        1.0 / k,
        dtype=np.float32
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


def moran_chunk(
    X,
    W
):

    """
    X:
        spots x genes, already transformed.
    """

    X = np.asarray(
        X,
        dtype=np.float32
    )


    Z = (
        X
        -
        X.mean(
            axis=0,
            keepdims=True
        )
    )


    denominator = np.sum(
        Z * Z,
        axis=0,
        dtype=np.float64
    )


    WZ = W @ Z


    numerator = np.sum(
        Z * WZ,
        axis=0,
        dtype=np.float64
    )


    # W is row standardized:
    # S0 = n, therefore n/S0 = 1.
    out = np.full(
        X.shape[1],
        np.nan,
        dtype=float
    )


    ok = (
        denominator > 0
    )


    out[
        ok
    ] = (
        numerator[
            ok
        ]
        /
        denominator[
            ok
        ]
    )


    return out


def calculate_metrics(
    spot_counts,
    library,
    genes,
    group,
    W
):

    """
    One common library denominator is supplied for BOTH groups.

    spot_counts:
        occupied spots x genes

    library:
        total predicted expression across
        300 OOF + 34,687 imputed-only genes.
    """

    spot_counts = np.asarray(
        spot_counts,
        dtype=np.float32
    )


    cpm = (
        spot_counts
        /
        library[
            :,
            None
        ]
        *
        1e6
    )


    logcpm = np.log1p(
        cpm
    ).astype(
        np.float32
    )


    mean_log = np.mean(
        logcpm,
        axis=0
    )


    sd_log = np.std(
        logcpm,
        axis=0,
        ddof=1
    )


    variance_log = (
        sd_log
        **
        2
    )


    q05 = np.quantile(
        logcpm,
        0.05,
        axis=0
    )


    q25 = np.quantile(
        logcpm,
        0.25,
        axis=0
    )


    q50 = np.quantile(
        logcpm,
        0.50,
        axis=0
    )


    q75 = np.quantile(
        logcpm,
        0.75,
        axis=0
    )


    q95 = np.quantile(
        logcpm,
        0.95,
        axis=0
    )


    iqr = (
        q75
        -
        q25
    )


    dynamic_range = (
        q95
        -
        q05
    )


    mean_cpm = np.mean(
        cpm,
        axis=0
    )


    sd_cpm = np.std(
        cpm,
        axis=0,
        ddof=1
    )


    cv = np.full(
        len(
            genes
        ),
        np.nan,
        dtype=float
    )


    nz = (
        mean_cpm > 0
    )


    cv[
        nz
    ] = (
        sd_cpm[
            nz
        ]
        /
        mean_cpm[
            nz
        ]
    )


    moran = moran_chunk(
        logcpm,
        W
    )


    out = pd.DataFrame({
        "gene":
            np.asarray(
                genes,
                dtype=str
            ),

        "group":
            group,

        "mean_log1p_CPM":
            mean_log,

        "sd_log1p_CPM":
            sd_log,

        "variance_log1p_CPM":
            variance_log,

        "q05_log1p_CPM":
            q05,

        "q25_log1p_CPM":
            q25,

        "median_log1p_CPM":
            q50,

        "q75_log1p_CPM":
            q75,

        "q95_log1p_CPM":
            q95,

        "IQR_log1p_CPM":
            iqr,

        "dynamic_range_95_5":
            dynamic_range,

        "mean_CPM":
            mean_cpm,

        "CV_CPM":
            cv,

        "moran_I":
            moran,
    })


    return out


GROUP_COLORS = {
    "Measured 300": "#1f77b4",
    "OOF 300": "#ff7f0e",
    "ENVI imputed-only 34K": "#2ca02c",
}


def ecdf(
    ax,
    x,
    label
):

    x = np.asarray(
        x,
        dtype=float
    )

    x = x[
        np.isfinite(
            x
        )
    ]

    x = np.sort(
        x
    )


    y = (
        np.arange(
            1,
            len(
                x
            ) + 1
        )
        /
        len(
            x
        )
    )


    ax.plot(
        x,
        y,
        linewidth=1.5,
        label=label,
        color=GROUP_COLORS[label]
    )


# =============================================================================
# INPUT CHECK
# =============================================================================

for f in [
    OOF,
    FULL,
    ALIGN,
]:

    if not f.exists():

        raise FileNotFoundError(
            f
        )


print("=" * 100)
print("BR8667: GENE-WISE VARIATION")
print("300 OOF vs 34K ENVI imputed-only")
print("=" * 100)


# =============================================================================
# LOAD DATA
# =============================================================================

print()
print("Loading OOF...")

oof = ad.read_h5ad(
    OOF
)

print(
    "OOF:",
    oof.shape
)


print()
print("Loading full ENVI transcriptome...")

full = ad.read_h5ad(
    FULL
)

print(
    "FULL:",
    full.shape
)


if "count_scale" not in oof.layers:

    raise RuntimeError(
        "OOF count_scale layer missing."
    )


if "ground_truth_counts" not in oof.layers:

    raise RuntimeError(
        "OOF ground_truth_counts layer missing."
    )


if "count_scale" not in full.layers:

    raise RuntimeError(
        "Full ENVI count_scale layer missing."
    )


oof_genes = np.asarray(
    oof.var_names,
    dtype=str
)


full_genes = np.asarray(
    full.var_names,
    dtype=str
)


if len(
    oof_genes
) != 300:

    raise RuntimeError(
        f"Expected 300 OOF genes; found {len(oof_genes)}"
    )


oof_set = set(
    oof_genes
)


imputed_mask = np.array([
    g not in oof_set
    for g
    in full_genes
])


imputed_idx = np.flatnonzero(
    imputed_mask
)


imputed_genes = full_genes[
    imputed_idx
]


print()
print(
    "OOF genes:",
    len(
        oof_genes
    )
)

print(
    "Full transcriptome genes:",
    len(
        full_genes
    )
)

print(
    "Imputed-only genes:",
    len(
        imputed_genes
    )
)


if len(
    imputed_genes
) != 34687:

    print(
        "WARNING: expected 34,687 imputed-only genes."
    )


# =============================================================================
# MATCH OOF CELLS TO FULL ENVI CELLS
# =============================================================================

print()
print("=" * 100)
print("MATCH OOF CELLS TO FULL ENVI")
print("=" * 100)


oof_ids = get_cell_ids(
    oof
)

full_ids = get_cell_ids(
    full
)


full_lookup = {
    cid: i
    for i, cid
    in enumerate(
        full_ids
    )
}


full_row_for_oof = np.array([
    full_lookup.get(
        cid,
        -1
    )
    for cid
    in oof_ids
])


has_full = (
    full_row_for_oof >= 0
)


print(
    "OOF cells with full-ENVI match:",
    f"{has_full.sum()}/{len(oof_ids)}"
)


# =============================================================================
# MATCH TO STALIGN + PHYSICAL VISIUM SPOTS
# =============================================================================

print()
print("=" * 100)
print("LOAD STALIGN ASSIGNMENTS")
print("=" * 100)


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


match_prob = np.full(
    len(
        oof_ids
    ),
    np.nan
)


spot_id = np.full(
    len(
        oof_ids
    ),
    "",
    dtype=object
)


x_aligned = np.full(
    len(
        oof_ids
    ),
    np.nan
)


y_aligned = np.full(
    len(
        oof_ids
    ),
    np.nan
)


ii = np.flatnonzero(
    has_st
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


spot_id[
    ii
] = (
    sub[
        "visium_spot"
    ]
    .fillna("")
    .astype(str)
    .to_numpy()
)


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


use = (
    has_full
    &
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
        spot_id != ""
    )
    &
    (
        spot_id != "nan"
    )
)


oof_rows = np.flatnonzero(
    use
)


full_rows = full_row_for_oof[
    oof_rows
]


used_spot_ids = spot_id[
    oof_rows
]


print(
    "Cells retained:",
    len(
        oof_rows
    )
)


# =============================================================================
# BUILD SPOT INDEX
# =============================================================================

spot_order = pd.Index(
    pd.unique(
        used_spot_ids
    )
)


spot_lookup = {
    s: i
    for i, s
    in enumerate(
        spot_order
    )
}


spot_index = np.array([
    spot_lookup[
        s
    ]
    for s
    in used_spot_ids
])


n_spots = len(
    spot_order
)


print(
    "Occupied spots:",
    n_spots
)


# cell -> spot aggregation matrix

A = sparse.csr_matrix(
    (
        np.ones(
            len(
                oof_rows
            ),
            dtype=np.float32
        ),
        (
            spot_index,
            np.arange(
                len(
                    oof_rows
                )
            )
        )
    ),
    shape=(
        n_spots,
        len(
            oof_rows
        )
    )
)


# =============================================================================
# SPOT COORDINATES
# =============================================================================

# Prefer actual Visium coordinates. If the coordinate file cannot be matched,
# fall back to the mean aligned Xenium coordinates within each physical spot.

spot_xy = None


if VIS_COORD.exists():

    try:

        vc = pd.read_csv(
            VIS_COORD,
            low_memory=False
        )


        id_candidates = [
            "barcode",
            "visium_spot",
            "spot",
            "spot_id",
        ]


        id_col = next(
            (
                c
                for c
                in id_candidates
                if c
                in vc.columns
            ),
            None
        )


        x_candidates = [
            "x",
            "x_coord",
            "pxl_col_in_fullres",
        ]


        y_candidates = [
            "y",
            "y_coord",
            "pxl_row_in_fullres",
        ]


        x_col = next(
            (
                c
                for c
                in x_candidates
                if c
                in vc.columns
            ),
            None
        )


        y_col = next(
            (
                c
                for c
                in y_candidates
                if c
                in vc.columns
            ),
            None
        )


        if (
            id_col is not None
            and
            x_col is not None
            and
            y_col is not None
        ):

            tmp = (
                vc[
                    [
                        id_col,
                        x_col,
                        y_col,
                    ]
                ]
                .copy()
            )


            tmp[
                id_col
            ] = tmp[
                id_col
            ].astype(str)


            tmp = tmp.set_index(
                id_col
            )


            if all(
                s in tmp.index
                for s
                in spot_order
            ):

                spot_xy = np.column_stack([
                    pd.to_numeric(
                        tmp.loc[
                            spot_order,
                            x_col
                        ],
                        errors="coerce"
                    ).to_numpy(),

                    pd.to_numeric(
                        tmp.loc[
                            spot_order,
                            y_col
                        ],
                        errors="coerce"
                    ).to_numpy(),
                ])


                print(
                    "Using true Visium coordinates."
                )

    except Exception as e:

        print(
            "Visium coordinate read warning:",
            repr(
                e
            )
        )


if spot_xy is None:

    print(
        "Using mean aligned Xenium coordinates per physical spot."
    )


    spot_xy = np.zeros(
        (
            n_spots,
            2
        ),
        dtype=float
    )


    used_x = x_aligned[
        oof_rows
    ]

    used_y = y_aligned[
        oof_rows
    ]


    for si in range(
        n_spots
    ):

        m = (
            spot_index
            ==
            si
        )


        spot_xy[
            si,
            0
        ] = np.nanmean(
            used_x[
                m
            ]
        )


        spot_xy[
            si,
            1
        ] = np.nanmean(
            used_y[
                m
            ]
        )


if not np.isfinite(
    spot_xy
).all():

    raise RuntimeError(
        "Non-finite spot coordinates."
    )


# =============================================================================
# SAME k=6 SPATIAL GRAPH
# =============================================================================

W = build_knn_W(
    spot_xy,
    k=K
)


print(
    "Spatial graph:",
    f"{K}-NN row-standardized"
)


# =============================================================================
# AGGREGATE 300 MEASURED XENIUM GENES
#
# ground_truth_counts is the true Xenium measurement for exactly the
# same OOF cell universe and the same 300 genes.
# =============================================================================

print()
print("=" * 100)
print("AGGREGATE 300 MEASURED XENIUM GENES")
print("=" * 100)


X_measured = oof.layers[
    "ground_truth_counts"
][
    oof_rows,
    :
]


MEASURED_spot = A @ X_measured


MEASURED_spot = dense_float(
    MEASURED_spot
)


print(
    "Measured pseudospot matrix:",
    MEASURED_spot.shape
)


# =============================================================================
# AGGREGATE 300 OOF GENES
# =============================================================================

print()
print("=" * 100)
print("AGGREGATE 300 OOF PREDICTIONS")
print("=" * 100)


X_oof = oof.layers[
    "count_scale"
][
    oof_rows,
    :
]


OOF_spot = A @ X_oof


OOF_spot = dense_float(
    OOF_spot
)


print(
    "OOF pseudospot matrix:",
    OOF_spot.shape
)


# =============================================================================
# AGGREGATE 34,687 IMPUTED-ONLY GENES
# =============================================================================

print()
print("=" * 100)
print("AGGREGATE FULL-TRANSCRIPTOME IMPUTED-ONLY GENES")
print("=" * 100)


n_imp = len(
    imputed_idx
)


IMP_spot = np.zeros(
    (
        n_spots,
        n_imp
    ),
    dtype=np.float32
)


for start in range(
    0,
    n_imp,
    CHUNK
):

    end = min(
        start + CHUNK,
        n_imp
    )


    idx_chunk = imputed_idx[
        start:end
    ]


    X = full.layers[
        "count_scale"
    ][
        full_rows,
        :
    ][
        :,
        idx_chunk
    ]


    agg = A @ X


    IMP_spot[
        :,
        start:end
    ] = dense_float(
        agg
    )


    print(
        f"Aggregated {start:,}-{end:,} / {n_imp:,}",
        flush=True
    )


print(
    "Imputed-only pseudospot matrix:",
    IMP_spot.shape
)


# =============================================================================
# ONE COMMON PREDICTED-TRANSCRIPTOME LIBRARY SIZE
#
# Use the SAME spot-specific size factor for measured-300, OOF-300,
# and imputed-only genes. This prevents differences between the three
# curves from being created by different normalization denominators.
#
# The denominator is based on:
#       300 OOF predictions + 34,687 imputed-only predictions
#       = 34,987 predicted genes.
#
# Measured-300 expression is normalized by this SAME external size
# factor solely to make its spatial-variation metrics directly
# comparable with OOF-300.
# =============================================================================

print()
print("=" * 100)
print("COMMON NORMALIZATION")
print("=" * 100)


library = (
    OOF_spot.sum(
        axis=1,
        dtype=np.float64
    )
    +
    IMP_spot.sum(
        axis=1,
        dtype=np.float64
    )
)


if np.any(
    library <= 0
):

    raise RuntimeError(
        "Found zero/negative pseudospot library."
    )


print(
    "Combined predicted transcriptome:",
    len(
        oof_genes
    )
    +
    len(
        imputed_genes
    ),
    "genes"
)


print(
    "Library median:",
    float(
        np.median(
            library
        )
    )
)


# =============================================================================
# 300 MEASURED METRICS
# =============================================================================

print()
print("=" * 100)
print("CALCULATE 300 MEASURED GENE-WISE VARIATION")
print("=" * 100)


metrics_measured = calculate_metrics(
    MEASURED_spot,
    library,
    oof_genes,
    "Measured_300",
    W
)


# =============================================================================
# 300 OOF METRICS
# =============================================================================

print()
print("=" * 100)
print("CALCULATE 300 OOF GENE-WISE VARIATION")
print("=" * 100)


metrics_oof = calculate_metrics(
    OOF_spot,
    library,
    oof_genes,
    "OOF_300",
    W
)


# =============================================================================
# 34K METRICS — CHUNKED
# =============================================================================

print()
print("=" * 100)
print("CALCULATE 34K IMPUTED-ONLY GENE-WISE VARIATION")
print("=" * 100)


parts = []


for start in range(
    0,
    n_imp,
    CHUNK
):

    end = min(
        start + CHUNK,
        n_imp
    )


    part = calculate_metrics(
        IMP_spot[
            :,
            start:end
        ],
        library,
        imputed_genes[
            start:end
        ],
        "ENVI_imputed_only_34K",
        W
    )


    parts.append(
        part
    )


    print(
        f"Metrics {start:,}-{end:,} / {n_imp:,}",
        flush=True
    )


metrics_imp = pd.concat(
    parts,
    ignore_index=True
)


metrics = pd.concat(
    [
        metrics_measured,
        metrics_oof,
        metrics_imp,
    ],
    ignore_index=True
)


# =============================================================================
# SAVE GENE-LEVEL TABLE
# =============================================================================

GENE_CSV = (
    OUT
    / f"{BR}_measured300_vs_OOF300_vs_imputed34K_gene_variation.csv"
)


metrics.to_csv(
    GENE_CSV,
    index=False
)


print(
    "Saved:",
    GENE_CSV
)


# =============================================================================
# SUMMARY TABLE
# =============================================================================

summary_rows = []


for group, d in metrics.groupby(
    "group",
    sort=False
):

    summary_rows.append({
        "group":
            group,

        "n_genes":
            len(
                d
            ),

        "median_SD_log1p_CPM":
            d[
                "sd_log1p_CPM"
            ].median(),

        "IQR_SD_log1p_CPM":
            (
                d[
                    "sd_log1p_CPM"
                ].quantile(
                    0.75
                )
                -
                d[
                    "sd_log1p_CPM"
                ].quantile(
                    0.25
                )
            ),

        "median_dynamic_range_95_5":
            d[
                "dynamic_range_95_5"
            ].median(),

        "median_IQR_log1p_CPM":
            d[
                "IQR_log1p_CPM"
            ].median(),

        "median_CV_CPM":
            d[
                "CV_CPM"
            ].median(),

        "median_Moran_I":
            d[
                "moran_I"
            ].median(),

        "fraction_Moran_gt_0_10":
            (
                d[
                    "moran_I"
                ]
                >
                0.10
            ).mean(),

        "fraction_Moran_gt_0_20":
            (
                d[
                    "moran_I"
                ]
                >
                0.20
            ).mean(),
    })


summary = pd.DataFrame(
    summary_rows
)


SUMMARY_CSV = (
    OUT
    / f"{BR}_measured300_vs_OOF300_vs_imputed34K_gene_variation_summary.csv"
)


summary.to_csv(
    SUMMARY_CSV,
    index=False
)


print()
print(summary.to_string(index=False))


# =============================================================================
# TOP VARIABLE GENES
# =============================================================================

for group, d in metrics.groupby(
    "group",
    sort=False
):

    top = (
        d.sort_values(
            "sd_log1p_CPM",
            ascending=False
        )
        .head(
            30
        )
    )


    fn = (
        OUT
        / f"{BR}_{group}_top30_variable_genes.csv"
    )


    top.to_csv(
        fn,
        index=False
    )


# =============================================================================
# FIGURE
# =============================================================================

m = metrics.loc[
    metrics[
        "group"
    ]
    ==
    "Measured_300"
].copy()


o = metrics.loc[
    metrics[
        "group"
    ]
    ==
    "OOF_300"
].copy()


i = metrics.loc[
    metrics[
        "group"
    ]
    ==
    "ENVI_imputed_only_34K"
].copy()


fig, axes = plt.subplots(
    2,
    2,
    figsize=(
        11,
        9
    )
)


# -----------------------------------------------------------------------------
# a. SD distribution
# -----------------------------------------------------------------------------

ax = axes[
    0,
    0
]


ecdf(
    ax,
    m[
        "sd_log1p_CPM"
    ],
    "Measured 300"
)


ecdf(
    ax,
    o[
        "sd_log1p_CPM"
    ],
    "OOF 300"
)


ecdf(
    ax,
    i[
        "sd_log1p_CPM"
    ],
    "ENVI imputed-only 34K"
)


ax.set_xlabel(
    "Gene-wise SD across spots\nlog1p(CPM)"
)

ax.set_ylabel(
    "Cumulative fraction of genes"
)

ax.set_title(
    "a   Spatial expression variation",
    loc="left",
    fontweight="bold"
)

ax.legend(
    frameon=False
)


# -----------------------------------------------------------------------------
# b. Dynamic range distribution
# -----------------------------------------------------------------------------

ax = axes[
    0,
    1
]


ecdf(
    ax,
    m[
        "dynamic_range_95_5"
    ],
    "Measured 300"
)


ecdf(
    ax,
    o[
        "dynamic_range_95_5"
    ],
    "OOF 300"
)


ecdf(
    ax,
    i[
        "dynamic_range_95_5"
    ],
    "ENVI imputed-only 34K"
)


ax.set_xlabel(
    "95th–5th percentile range\nlog1p(CPM)"
)

ax.set_ylabel(
    "Cumulative fraction of genes"
)

ax.set_title(
    "b   Spatial dynamic range",
    loc="left",
    fontweight="bold"
)

ax.legend(
    frameon=False
)


# -----------------------------------------------------------------------------
# c. Moran I distribution
# -----------------------------------------------------------------------------

ax = axes[
    1,
    0
]


ecdf(
    ax,
    m[
        "moran_I"
    ],
    "Measured 300"
)


ecdf(
    ax,
    o[
        "moran_I"
    ],
    "OOF 300"
)


ecdf(
    ax,
    i[
        "moran_I"
    ],
    "ENVI imputed-only 34K"
)


ax.axvline(
    0,
    linestyle="--",
    linewidth=0.8,
    color="0.5"
)


ax.set_xlabel(
    "Moran's I"
)

ax.set_ylabel(
    "Cumulative fraction of genes"
)

ax.set_title(
    "c   Spatial autocorrelation",
    loc="left",
    fontweight="bold"
)

ax.legend(
    frameon=False
)


# -----------------------------------------------------------------------------
# d. Mean versus SD
# -----------------------------------------------------------------------------

ax = axes[
    1,
    1
]


rng = np.random.default_rng(
    14053
)


if len(
    i
) > 5000:

    ii = rng.choice(
        len(
            i
        ),
        size=5000,
        replace=False
    )

    ii_df = i.iloc[
        ii
    ]

else:

    ii_df = i


ax.scatter(
    ii_df[
        "mean_log1p_CPM"
    ],
    ii_df[
        "sd_log1p_CPM"
    ],
    s=5,
    alpha=0.20,
    color=GROUP_COLORS["ENVI imputed-only 34K"],
    label="ENVI imputed-only 34K"
)


ax.scatter(
    m[
        "mean_log1p_CPM"
    ],
    m[
        "sd_log1p_CPM"
    ],
    s=14,
    alpha=0.75,
    color=GROUP_COLORS["Measured 300"],
    label="Measured 300"
)


ax.scatter(
    o[
        "mean_log1p_CPM"
    ],
    o[
        "sd_log1p_CPM"
    ],
    s=12,
    alpha=0.65,
    color=GROUP_COLORS["OOF 300"],
    label="OOF 300"
)


ax.set_xlabel(
    "Mean log1p(CPM)"
)

ax.set_ylabel(
    "SD log1p(CPM)"
)

ax.set_title(
    "d   Mean–variation relationship",
    loc="left",
    fontweight="bold"
)

ax.legend(
    frameon=False,
    fontsize=8
)


fig.suptitle(
    (
        "Br8667: gene-wise spatial variation of "
        "300 measured, 300 OOF, and 34K ENVI-imputed genes"
    ),
    fontweight="bold",
    fontsize=13
)


fig.text(
    0.5,
    0.015,
    (
        "Both gene sets use the same matched cells, the same physical Visium pseudospots, "
        "one combined 34,987-gene library-size denominator, and the same 6-NN spatial graph. "
        "Measured 300 are the Xenium ground-truth counts for the same cells; "
        "OOF 300 are held-out ENVI predictions of those same genes; "
        "the 34K group contains ENVI imputed-only genes."
    ),
    ha="center",
    fontsize=7
)


fig.tight_layout(
    rect=[
        0,
        0.05,
        1,
        0.95
    ]
)


PNG = (
    OUT
    / f"fig_{BR}_measured300_vs_OOF300_vs_imputed34K_gene_variation.png"
)

PDF = (
    OUT
    / f"fig_{BR}_measured300_vs_OOF300_vs_imputed34K_gene_variation.pdf"
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


print()
print("=" * 100)
print("COMPLETE")
print("=" * 100)

print(
    "Gene metrics:",
    GENE_CSV
)

print(
    "Summary:",
    SUMMARY_CSV
)

print(
    "Figure:",
    PNG
)


oof.file.close()
full.file.close()
