#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse, stats


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

META = (
    ROOT
    / "data/metadata/patient_xenium_visium_24_common_with_dx.csv"
)

EX1_DIR = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/targets"
)

EX2_DIR = (
    ROOT
    / "data/processed/imputation_full/ex2_scz/targets"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/before_imputation"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


def section(title):
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def resolve_col(df, candidates):
    lookup = {
        str(c).strip().lower(): c
        for c in df.columns
    }

    for candidate in candidates:
        key = candidate.lower()

        if key in lookup:
            return lookup[key]

    raise KeyError(
        f"Could not find any of {candidates}. "
        f"Columns={df.columns.tolist()}"
    )


def bh_fdr(pvalues):
    """
    Benjamini-Hochberg FDR correction.
    """

    p = np.asarray(
        pvalues,
        dtype=float,
    )

    n = len(p)

    order = np.argsort(p)

    ranked = p[order]

    adjusted = (
        ranked
        * n
        / np.arange(
            1,
            n + 1,
        )
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.minimum(
        adjusted,
        1.0,
    )

    result = np.empty(
        n,
        dtype=float,
    )

    result[order] = adjusted

    return result


# ============================================================
# Metadata
# ============================================================

section("LOAD DONOR METADATA")

meta = pd.read_csv(META)

patient_col = resolve_col(
    meta,
    ["patient_id", "BrNum"],
)

dx_col = resolve_col(
    meta,
    ["Dx", "diagnosis"],
)

age_col = resolve_col(
    meta,
    ["AGE", "age"],
)

sex_col = resolve_col(
    meta,
    ["SEX", "sex"],
)


meta = meta.rename(
    columns={
        patient_col: "patient_id",
        dx_col: "Dx",
        age_col: "Age",
        sex_col: "Sex",
    }
)


meta["patient_id"] = (
    meta["patient_id"]
    .astype(str)
    .str.strip()
)

meta["Dx"] = (
    meta["Dx"]
    .astype(str)
    .str.strip()
    .str.upper()
)

meta["Sex"] = (
    meta["Sex"]
    .astype(str)
    .str.strip()
    .str.upper()
)

meta["Age"] = pd.to_numeric(
    meta["Age"],
    errors="coerce",
)


# Br6432 intentionally excluded.
meta = meta.loc[
    meta["patient_id"] != "Br6432"
].copy()


meta = meta.loc[
    meta["Dx"].isin(
        ["NTC", "SCZ"]
    )
].copy()


print(
    meta[
        [
            "patient_id",
            "Dx",
            "Age",
            "Sex",
        ]
    ]
    .sort_values(
        ["Dx", "patient_id"]
    )
    .to_string(index=False)
)


print("\nDiagnosis counts:")
print(
    meta["Dx"].value_counts()
)


if (
    (meta["Dx"] == "NTC").sum()
    != 11
):
    raise RuntimeError(
        "Expected 11 NTC donors."
    )


if (
    (meta["Dx"] == "SCZ").sum()
    != 12
):
    raise RuntimeError(
        "Expected 12 SCZ donors."
    )


# ============================================================
# Locate target files
# ============================================================

section("LOCATE ORIGINAL XENIUM TARGET FILES")

records = []


for row in meta.itertuples(
    index=False
):

    donor = str(
        row.patient_id
    )

    dx = str(
        row.Dx
    )

    if dx == "NTC":

        path = (
            EX1_DIR
            / f"spatial_data_xenium_{donor}_ex1_ntc.h5ad"
        )

    else:

        path = (
            EX2_DIR
            / f"spatial_data_xenium_{donor}_ex2_scz.h5ad"
        )


    if not path.is_file():
        raise FileNotFoundError(
            path
        )


    records.append(
        {
            "patient_id":
                donor,

            "Dx":
                dx,

            "Age":
                float(
                    row.Age
                ),

            "Sex":
                str(
                    row.Sex
                ),

            "file":
                path,
        }
    )


print(
    f"Found {len(records)} / 23 target files."
)


# ============================================================
# Pseudobulk
# ============================================================

section("BUILD DONOR PSEUDOBULK")

gene_order = None

pseudobulk_counts = []

donor_rows = []


for i, record in enumerate(
    records,
    start=1,
):

    donor = record[
        "patient_id"
    ]

    print(
        f"[{i:02d}/23] {donor} | {record['Dx']}",
        flush=True,
    )


    a = ad.read_h5ad(
        record["file"],
        backed="r",
    )


    genes = (
        a.var_names
        .astype(str)
        .tolist()
    )


    if gene_order is None:

        gene_order = genes

    elif genes != gene_order:

        raise RuntimeError(
            f"{donor}: gene order differs."
        )


    if len(genes) != 300:
        raise RuntimeError(
            f"{donor}: expected 300 genes, "
            f"found {len(genes)}"
        )


    if "counts" in a.layers:

        matrix = (
            a.layers["counts"]
        )

    else:

        matrix = a.X


    # Sum cells → one donor-level profile.
    if sparse.issparse(
        matrix
    ):

        summed = np.asarray(
            matrix.sum(axis=0)
        ).ravel()

    else:

        # backed dense HDF5
        summed = np.asarray(
            matrix[:]
        ).sum(
            axis=0,
            dtype=np.float64,
        )


    summed = np.asarray(
        summed,
        dtype=np.float64,
    )


    if (
        not np.isfinite(
            summed
        ).all()
        or np.any(
            summed < 0
        )
    ):
        raise RuntimeError(
            f"{donor}: invalid pseudobulk counts."
        )


    pseudobulk_counts.append(
        summed
    )


    donor_rows.append(
        {
            "patient_id":
                donor,

            "Dx":
                record["Dx"],

            "Age":
                record["Age"],

            "Sex":
                record["Sex"],

            "n_cells":
                int(
                    a.n_obs
                ),

            "panel_library":
                float(
                    summed.sum()
                ),
        }
    )


    a.file.close()


counts = np.vstack(
    pseudobulk_counts
)


donors = pd.DataFrame(
    donor_rows
)


print("\nPseudobulk matrix:")
print(counts.shape)


# ============================================================
# Save raw donor pseudobulk
# ============================================================

raw_df = pd.DataFrame(
    counts,
    index=donors[
        "patient_id"
    ],
    columns=gene_order,
)


raw_df.to_csv(
    OUT
    / "original_xenium_pseudobulk_counts_23donors.csv"
)


donors.to_csv(
    OUT
    / "original_xenium_pseudobulk_donor_metadata.csv",
    index=False,
)


# ============================================================
# Panel CPM + log2
# ============================================================

section("NORMALIZE PSEUDOBULK")

library = counts.sum(
    axis=1,
    dtype=np.float64,
)


if np.any(
    library <= 0
):
    raise RuntimeError(
        "Zero pseudobulk library detected."
    )


cpm = (
    counts
    / library[:, None]
    * 1_000_000.0
)


log2cpm = np.log2(
    cpm + 1.0
)


pd.DataFrame(
    log2cpm,
    index=donors[
        "patient_id"
    ],
    columns=gene_order,
).to_csv(
    OUT
    / "original_xenium_pseudobulk_log2cpm_23donors.csv"
)


# ============================================================
# Design matrix:
#
# expression ~ Dx + Age + Sex
#
# NTC = reference
# SCZ coefficient = adjusted SCZ-vs-NTC log2FC
# ============================================================

section("BUILD DE DESIGN MATRIX")


dx_scz = (
    donors["Dx"]
    .eq("SCZ")
    .astype(float)
    .to_numpy()
)


age = (
    donors["Age"]
    .astype(float)
    .to_numpy()
)


# center age for numerical stability
age = (
    age
    - np.nanmean(age)
)


sex_m = (
    donors["Sex"]
    .eq("M")
    .astype(float)
    .to_numpy()
)


valid = (
    np.isfinite(age)
)


if not valid.all():

    print(
        "WARNING: dropping donors "
        "with missing Age."
    )

    donors = (
        donors.loc[
            valid
        ]
        .reset_index(
            drop=True
        )
    )

    log2cpm = (
        log2cpm[
            valid,
            :
        ]
    )

    dx_scz = dx_scz[
        valid
    ]

    age = age[
        valid
    ]

    sex_m = sex_m[
        valid
    ]


X = np.column_stack(
    [
        np.ones(
            len(donors)
        ),

        dx_scz,

        age,

        sex_m,
    ]
)


design_columns = [
    "Intercept",
    "Dx_SCZ",
    "Age_centered",
    "Sex_M",
]


design = pd.DataFrame(
    X,
    index=donors[
        "patient_id"
    ],
    columns=design_columns,
)


design.to_csv(
    OUT
    / "original_xenium_de_design_matrix.csv"
)


print(design)


# ============================================================
# Vectorized OLS across 300 genes
# ============================================================

section("SCZ VS NTC DIFFERENTIAL EXPRESSION")


n, p = X.shape

rank = np.linalg.matrix_rank(
    X
)


df_resid = (
    n - rank
)


if df_resid <= 0:
    raise RuntimeError(
        "No residual degrees of freedom."
    )


xtx_inv = np.linalg.pinv(
    X.T @ X
)


beta = (
    xtx_inv
    @ X.T
    @ log2cpm
)


fitted = (
    X
    @ beta
)


residual = (
    log2cpm
    - fitted
)


rss = np.sum(
    residual ** 2,
    axis=0,
)


sigma2 = (
    rss
    / df_resid
)


# Dx_SCZ is coefficient index 1.
log2fc = (
    beta[
        1,
        :
    ]
)


variance_log2fc = (
    sigma2
    * xtx_inv[
        1,
        1
    ]
)


se_log2fc = np.sqrt(
    variance_log2fc
)


t_stat = np.divide(
    log2fc,
    se_log2fc,
    out=np.zeros_like(
        log2fc
    ),
    where=(
        se_log2fc > 0
    ),
)


pvalue = (
    2.0
    * stats.t.sf(
        np.abs(
            t_stat
        ),
        df=df_resid,
    )
)


fdr = bh_fdr(
    pvalue
)


# ============================================================
# Group summaries
# ============================================================

ntc_mask = (
    donors["Dx"]
    .eq("NTC")
    .to_numpy()
)


scz_mask = (
    donors["Dx"]
    .eq("SCZ")
    .to_numpy()
)


mean_ntc = (
    log2cpm[
        ntc_mask,
        :
    ]
    .mean(
        axis=0
    )
)


mean_scz = (
    log2cpm[
        scz_mask,
        :
    ]
    .mean(
        axis=0
    )
)


result = pd.DataFrame(
    {
        "gene":
            gene_order,

        "gene_source":
            "measured_xenium",

        "mean_log2CPM_NTC":
            mean_ntc,

        "mean_log2CPM_SCZ":
            mean_scz,

        "log2FC":
            log2fc,

        "P-value":
            pvalue,

        "FDR":
            fdr,

        "t_statistic":
            t_stat,

        "n_NTC":
            int(
                ntc_mask.sum()
            ),

        "n_SCZ":
            int(
                scz_mask.sum()
            ),
    }
)


result[
    "direction"
] = np.where(
    result[
        "log2FC"
    ] > 0,
    "Up_in_SCZ",
    "Down_in_SCZ",
)


result[
    "significant_FDR_0.05"
] = (
    result[
        "FDR"
    ] < 0.05
)


result = result.sort_values(
    [
        "FDR",
        "P-value",
    ]
).reset_index(
    drop=True
)


# ============================================================
# Full table
# ============================================================

full_path = (
    OUT
    / "original_xenium_DE_SCZ_vs_NTC_full.csv"
)


result.to_csv(
    full_path,
    index=False,
)


# ============================================================
# Advisor-format table
# ============================================================

advisor = result[
    [
        "gene",
        "P-value",
        "FDR",
        "log2FC",
    ]
].copy()


advisor_path = (
    OUT
    / "original_xenium_DE_SCZ_vs_NTC_advisor_table.csv"
)


advisor.to_csv(
    advisor_path,
    index=False,
)


# ============================================================
# Significant genes
# ============================================================

significant = (
    result.loc[
        result[
            "FDR"
        ] < 0.05
    ]
    .copy()
)


significant.to_csv(
    OUT
    / "original_xenium_DE_SCZ_vs_NTC_FDR05.csv",
    index=False,
)


# ============================================================
# Summary
# ============================================================

section("DE SUMMARY")


print(
    "Biological replicates:",
    len(
        donors
    ),
)

print(
    "NTC donors:",
    int(
        ntc_mask.sum()
    ),
)

print(
    "SCZ donors:",
    int(
        scz_mask.sum()
    ),
)

print(
    "Genes tested:",
    len(
        result
    ),
)

print(
    "FDR < 0.05:",
    len(
        significant
    ),
)

print(
    "Up in SCZ:",
    int(
        (
            significant[
                "log2FC"
            ] > 0
        ).sum()
    ),
)

print(
    "Down in SCZ:",
    int(
        (
            significant[
                "log2FC"
            ] < 0
        ).sum()
    ),
)


print("\nTop 20:")
print(
    result[
        [
            "gene",
            "P-value",
            "FDR",
            "log2FC",
            "direction",
        ]
    ]
    .head(20)
    .to_string(
        index=False
    )
)


print("\nAdvisor table:")
print(advisor_path)

print("\nFull table:")
print(full_path)


section("DONE")
