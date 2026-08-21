#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = ROOT / "data/metadata/envi_production_23donors.tsv"
META = ROOT / "data/metadata/patient_xenium_visium_24_common_with_dx.csv"

OUT = ROOT / "outputs/imputation_full/DE/after_imputation_measured300"
OUT.mkdir(parents=True, exist_ok=True)

BASELINE_DE = (
    ROOT
    / "outputs/imputation_full/DE/before_imputation"
    / "original_xenium_DE_SCZ_vs_NTC_full.csv"
)


def dense(x):
    if sparse.issparse(x):
        return x.toarray()

    if hasattr(x, "toarray"):
        return x.toarray()

    return np.asarray(x)


def bh_fdr(p):
    p = np.asarray(p, dtype=float)

    order = np.argsort(p)
    ranked = p[order]
    n = len(p)

    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.minimum(q, 1.0)

    out = np.empty(n, dtype=float)
    out[order] = q

    return out


print("=" * 100)
print("AFTER-ENVI MEASURED-300 QC")
print("=" * 100)


# =============================================================================
# Manifest
# =============================================================================

manifest = pd.read_csv(MANIFEST, sep="\t")

required = {"donor", "Dx", "experiment", "target"}

if not required.issubset(manifest.columns):
    raise RuntimeError(
        f"Manifest columns: {manifest.columns.tolist()}\n"
        f"Required: {sorted(required)}"
    )

manifest["donor"] = manifest["donor"].astype(str).str.strip()
manifest["Dx"] = manifest["Dx"].astype(str).str.strip().str.upper()

assert len(manifest) == 23
assert manifest["donor"].nunique() == 23
assert (manifest["Dx"] == "NTC").sum() == 11
assert (manifest["Dx"] == "SCZ").sum() == 12

print("\nDonors:")
print(manifest[["donor", "Dx", "experiment"]].to_string(index=False))

print("\nNTC :", (manifest["Dx"] == "NTC").sum())
print("SCZ :", (manifest["Dx"] == "SCZ").sum())
print("Total:", len(manifest))


# =============================================================================
# Direct BEFORE-vs-AFTER comparison of the 300 measured genes
# =============================================================================

print("\n" + "=" * 100)
print("DIRECT CELL-LEVEL BEFORE vs AFTER COMPARISON")
print("=" * 100)

pb_rows = []
qc_rows = []

canonical_genes = None

CHUNK = 5000


for i, r in enumerate(manifest.itertuples(index=False), 1):

    donor = r.donor
    dx = r.Dx
    experiment = r.experiment

    before_path = Path(r.target)

    after_path = (
        ROOT
        / f"data/processed/imputation_full/{experiment}/envi/{donor}"
        / f"spatial_data_xenium_{donor}_ENVI_full_transcriptome.h5ad"
    )

    if not before_path.exists():
        raise FileNotFoundError(before_path)

    if not after_path.exists():
        raise FileNotFoundError(after_path)

    print(
        f"\n[{i:02d}/23] {donor} | {dx}",
        flush=True,
    )

    before = ad.read_h5ad(before_path, backed="r")
    after = ad.read_h5ad(after_path, backed="r")

    if before.n_obs != after.n_obs:
        raise RuntimeError(
            f"{donor}: cell number differs: "
            f"before={before.n_obs}, after={after.n_obs}"
        )

    # Production output should retain cell order.
    before_obs = before.obs_names.astype(str).to_numpy()
    after_obs = after.obs_names.astype(str).to_numpy()

    if not np.array_equal(before_obs, after_obs):
        raise RuntimeError(
            f"{donor}: cell names/order differ before vs after."
        )

    if "expression_source" not in after.var.columns:
        raise RuntimeError(
            f"{donor}: expression_source missing."
        )

    source = after.var["expression_source"].astype(str).to_numpy()

    measured_idx = np.flatnonzero(
        source == "measured_xenium"
    )

    if len(measured_idx) != 300:
        raise RuntimeError(
            f"{donor}: expected 300 measured genes, "
            f"found {len(measured_idx)}"
        )

    measured_genes = (
        after.var_names[
            measured_idx
        ]
        .astype(str)
        .tolist()
    )

    if canonical_genes is None:
        canonical_genes = measured_genes

    elif measured_genes != canonical_genes:
        raise RuntimeError(
            f"{donor}: measured-gene order differs across donors."
        )

    # Target contains exactly the original Xenium panel.
    before_gene_index = before.var_names.get_indexer(
        measured_genes
    )

    if np.any(before_gene_index < 0):
        missing = np.asarray(measured_genes)[
            before_gene_index < 0
        ]

        raise RuntimeError(
            f"{donor}: measured genes missing from target: "
            f"{missing[:20].tolist()}"
        )

    # Production runner uses original target counts.
    if "counts" in before.layers:
        before_matrix = before.layers["counts"]
        before_source = 'layers["counts"]'
    else:
        before_matrix = before.X
        before_source = "X"

    after_matrix = after.X

    pb_sum = np.zeros(
        300,
        dtype=np.float64,
    )

    n_diff = 0
    max_abs_diff = 0.0
    sum_abs_diff = 0.0
    n_values = 0

    for start in range(
        0,
        before.n_obs,
        CHUNK,
    ):

        stop = min(
            start + CHUNK,
            before.n_obs,
        )

        # BEFORE has only 300 genes, so load the entire 300-gene chunk,
        # then reorder in memory if needed.
        b = dense(
            before_matrix[
                start:stop,
                :
            ]
        )

        b = np.asarray(
            b,
            dtype=np.float64,
        )

        b = b[
            :,
            before_gene_index
        ]

        # AFTER has 34,987 genes; read only the 300 measured columns.
        a = dense(
            after_matrix[
                start:stop,
                measured_idx
            ]
        )

        a = np.asarray(
            a,
            dtype=np.float64,
        )

        if b.shape != a.shape:
            raise RuntimeError(
                f"{donor}: chunk shape mismatch "
                f"{b.shape} vs {a.shape}"
            )

        d = a - b

        n_diff += int(
            np.count_nonzero(d)
        )

        if d.size:
            max_abs_diff = max(
                max_abs_diff,
                float(
                    np.max(
                        np.abs(d)
                    )
                ),
            )

            sum_abs_diff += float(
                np.sum(
                    np.abs(d)
                )
            )

            n_values += d.size

        pb_sum += a.sum(
            axis=0,
            dtype=np.float64,
        )

    mean_abs_diff = (
        sum_abs_diff / n_values
        if n_values
        else 0.0
    )

    exact = (
        n_diff == 0
        and max_abs_diff == 0.0
    )

    print(
        f"  cells             : {after.n_obs:,}"
    )
    print(
        f"  before source     : {before_source}"
    )
    print(
        f"  after total genes : {after.n_vars:,}"
    )
    print(
        f"  measured genes    : {len(measured_idx)}"
    )
    print(
        f"  exact equality    : {exact}"
    )
    print(
        f"  differing values  : {n_diff:,}"
    )
    print(
        f"  max abs diff      : {max_abs_diff}"
    )

    qc_rows.append({
        "donor": donor,
        "Dx": dx,
        "n_cells": after.n_obs,
        "before_source": before_source,
        "after_total_genes": after.n_vars,
        "measured_genes": len(measured_idx),
        "exact_equal": exact,
        "n_different_values": n_diff,
        "max_abs_diff": max_abs_diff,
        "mean_abs_diff": mean_abs_diff,
    })

    pb_rows.append(
        pd.Series(
            pb_sum,
            index=measured_genes,
            name=donor,
        )
    )

    before.file.close()
    after.file.close()


qc = pd.DataFrame(qc_rows)

qc.to_csv(
    OUT / "before_vs_after_measured300_per_donor_qc.csv",
    index=False,
)


if not qc["exact_equal"].all():

    print("\nFAILED DONORS:")
    print(
        qc.loc[
            ~qc["exact_equal"]
        ].to_string(index=False)
    )

    raise RuntimeError(
        "FAIL: at least one donor has changed measured Xenium values."
    )


print("\n" + "=" * 100)
print("ALL 23 DONORS PASS CELL-LEVEL MEASURED-300 QC")
print("=" * 100)


# =============================================================================
# AFTER-ENVI donor pseudobulk
# =============================================================================

pb = pd.DataFrame(pb_rows)

# Match manifest donor order.
pb = pb.loc[
    manifest["donor"].tolist()
]

assert pb.shape == (23, 300)

pb.to_csv(
    OUT / "after_ENVI_measured300_pseudobulk_counts_23donors.csv"
)


# =============================================================================
# Metadata for DE
# =============================================================================

meta = pd.read_csv(META)

patient_col = (
    "patient_id"
    if "patient_id" in meta.columns
    else "BrNum"
)

age_col = (
    "AGE"
    if "AGE" in meta.columns
    else "Age"
)

sex_col = (
    "SEX"
    if "SEX" in meta.columns
    else "Sex"
)

meta[patient_col] = (
    meta[patient_col]
    .astype(str)
    .str.strip()
)

meta = meta.set_index(
    patient_col
)

donor_meta = manifest[
    ["donor", "Dx"]
].copy()

donor_meta["Age"] = [
    meta.loc[d, age_col]
    for d in donor_meta["donor"]
]

donor_meta["Sex"] = [
    meta.loc[d, sex_col]
    for d in donor_meta["donor"]
]

donor_meta["Age"] = pd.to_numeric(
    donor_meta["Age"],
    errors="coerce",
)

donor_meta["Sex"] = (
    donor_meta["Sex"]
    .astype(str)
    .str.strip()
    .str.upper()
)

if donor_meta["Age"].isna().any():
    raise RuntimeError(
        "Missing Age values:\n"
        + donor_meta.loc[
            donor_meta["Age"].isna()
        ].to_string(index=False)
    )

donor_meta.to_csv(
    OUT / "after_ENVI_measured300_donor_metadata.csv",
    index=False,
)


# =============================================================================
# Same normalization used for original 300-gene baseline
# =============================================================================

counts = pb.to_numpy(
    dtype=np.float64,
)

library = counts.sum(
    axis=1,
)

if np.any(library <= 0):
    raise RuntimeError(
        "Zero panel library detected."
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
    index=pb.index,
    columns=pb.columns,
).to_csv(
    OUT / "after_ENVI_measured300_log2CPM_23donors.csv"
)


# =============================================================================
# expression ~ Dx + Age + Sex
# =============================================================================

dx = (
    donor_meta["Dx"]
    .eq("SCZ")
    .astype(float)
    .to_numpy()
)

age = donor_meta["Age"].to_numpy(
    dtype=float
)

age = age - age.mean()

sex_m = (
    donor_meta["Sex"]
    .eq("M")
    .astype(float)
    .to_numpy()
)

X = np.column_stack([
    np.ones(len(donor_meta)),
    dx,
    age,
    sex_m,
])

design = pd.DataFrame(
    X,
    index=donor_meta["donor"],
    columns=[
        "Intercept",
        "Dx_SCZ",
        "Age_centered",
        "Sex_M",
    ],
)

design.to_csv(
    OUT / "after_ENVI_measured300_design_matrix.csv"
)


rank = np.linalg.matrix_rank(X)
df_resid = len(donor_meta) - rank

xtx_inv = np.linalg.pinv(
    X.T @ X
)

beta = (
    xtx_inv
    @ X.T
    @ log2cpm
)

fitted = X @ beta
resid = log2cpm - fitted

rss = np.sum(
    resid ** 2,
    axis=0,
)

sigma2 = rss / df_resid

log2fc = beta[1, :]

se = np.sqrt(
    sigma2
    * xtx_inv[1, 1]
)

tstat = np.divide(
    log2fc,
    se,
    out=np.zeros_like(log2fc),
    where=se > 0,
)

pvalue = (
    2.0
    * stats.t.sf(
        np.abs(tstat),
        df=df_resid,
    )
)

fdr = bh_fdr(pvalue)


ntc = (
    donor_meta["Dx"]
    .eq("NTC")
    .to_numpy()
)

scz = (
    donor_meta["Dx"]
    .eq("SCZ")
    .to_numpy()
)


result = pd.DataFrame({
    "gene": pb.columns,
    "gene_source": "measured_xenium",
    "mean_log2CPM_NTC": log2cpm[ntc].mean(axis=0),
    "mean_log2CPM_SCZ": log2cpm[scz].mean(axis=0),
    "log2FC": log2fc,
    "P-value": pvalue,
    "FDR": fdr,
    "t_statistic": tstat,
    "n_NTC": int(ntc.sum()),
    "n_SCZ": int(scz.sum()),
})

result["direction"] = np.where(
    result["log2FC"] > 0,
    "Higher_in_SCZ",
    np.where(
        result["log2FC"] < 0,
        "Higher_in_NTC",
        "No_change",
    ),
)

result = result.sort_values(
    ["FDR", "P-value"]
).reset_index(drop=True)


result.to_csv(
    OUT / "after_ENVI_measured300_DE_SCZ_vs_NTC_full.csv",
    index=False,
)

result[
    ["gene", "P-value", "FDR", "log2FC"]
].to_csv(
    OUT / "after_ENVI_measured300_DE_SCZ_vs_NTC_advisor_table.csv",
    index=False,
)


# =============================================================================
# Summary
# =============================================================================

nominal = result[
    result["P-value"] < 0.05
]

significant = result[
    result["FDR"] < 0.05
]


print("\n" + "=" * 100)
print("AFTER-IMPUTATION MEASURED-300 DE SUMMARY")
print("=" * 100)

print("\nAll 300 genes")
print(
    "Higher in SCZ :",
    int(
        (result["log2FC"] > 0).sum()
    ),
)
print(
    "Higher in NTC :",
    int(
        (result["log2FC"] < 0).sum()
    ),
)

print("\nNominal P < 0.05")
print("Total         :", len(nominal))
print(
    "Higher in SCZ :",
    int(
        (nominal["log2FC"] > 0).sum()
    ),
)
print(
    "Higher in NTC :",
    int(
        (nominal["log2FC"] < 0).sum()
    ),
)

print("\nFDR < 0.05")
print("Total         :", len(significant))
print(
    "Higher in SCZ :",
    int(
        (significant["log2FC"] > 0).sum()
    ),
)
print(
    "Higher in NTC :",
    int(
        (significant["log2FC"] < 0).sum()
    ),
)


# =============================================================================
# Compare against original baseline DE if available
# =============================================================================

if BASELINE_DE.exists():

    baseline = pd.read_csv(BASELINE_DE)

    cmp = baseline[
        ["gene", "log2FC", "P-value", "FDR"]
    ].merge(
        result[
            ["gene", "log2FC", "P-value", "FDR"]
        ],
        on="gene",
        suffixes=("_before", "_after"),
        validate="one_to_one",
    )

    for col in [
        "log2FC",
        "P-value",
        "FDR",
    ]:

        cmp[
            f"abs_diff_{col}"
        ] = np.abs(
            cmp[f"{col}_before"]
            - cmp[f"{col}_after"]
        )

    cmp.to_csv(
        OUT / "before_vs_after_measured300_DE_comparison.csv",
        index=False,
    )

    print("\nBEFORE vs AFTER DE STATISTICS")

    print(
        "Max |Δ log2FC| :",
        cmp["abs_diff_log2FC"].max(),
    )

    print(
        "Max |Δ P-value|:",
        cmp["abs_diff_P-value"].max(),
    )

    print(
        "Max |Δ FDR|    :",
        cmp["abs_diff_FDR"].max(),
    )

else:

    print(
        "\nBaseline DE table not found; "
        "direct measured-count QC still completed."
    )


print("\n" + "=" * 100)
print("PASS")
print("=" * 100)

print(
    "All 23 donors retain the original 300 measured Xenium genes "
    "unchanged after ENVI imputation."
)

print(
    "\nNext analysis: 34,687 ENVI-imputed genes, "
    "SCZ vs NTC."
)
