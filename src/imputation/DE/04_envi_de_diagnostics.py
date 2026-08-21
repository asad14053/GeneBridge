#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BASE = (
    ROOT
    / "outputs/imputation_full/DE/envi_transcriptome"
)

PB_FILE = (
    BASE
    / "ENVI_full34987_pseudobulk_countscale_23donors.csv.gz"
)

DX_FILE = (
    BASE
    / "ENVI_full34987_diagnosis.csv"
)

DE_FILE = (
    BASE
    / "ENVI_imputed34687_DE_SCZ_vs_NTC_full.csv"
)

OUT = BASE / "diagnostics"
OUT.mkdir(parents=True, exist_ok=True)


print("=" * 100)
print("ENVI DE DIAGNOSTICS")
print("=" * 100)


# =============================================================================
# Load data
# =============================================================================

pb = pd.read_csv(
    PB_FILE,
    index_col=0,
)

dx = pd.read_csv(
    DX_FILE
)

de = pd.read_csv(
    DE_FILE
)


dx["donor"] = (
    dx["donor"]
    .astype(str)
    .str.strip()
)

dx["Dx"] = (
    dx["Dx"]
    .astype(str)
    .str.strip()
    .str.upper()
)

dx = (
    dx.set_index("donor")
    .loc[pb.index]
)


print("Pseudobulk:", pb.shape)
print("DE genes  :", len(de))
print("\nDiagnosis:")
print(dx["Dx"].value_counts())


assert pb.shape == (23, 34987)
assert len(de) == 34687
assert (dx["Dx"] == "SCZ").sum() == 12
assert (dx["Dx"] == "NTC").sum() == 11


# =============================================================================
# Same normalization used in limma analysis
# =============================================================================

counts = pb.to_numpy(
    dtype=np.float64
)

library = counts.sum(axis=1)

cpm = (
    counts
    / library[:, None]
    * 1_000_000.0
)

log2cpm = np.log2(
    cpm + 1.0
)


log2cpm_df = pd.DataFrame(
    log2cpm,
    index=pb.index,
    columns=pb.columns,
)


# =============================================================================
# 1. P-value histogram
# =============================================================================

p = (
    pd.to_numeric(
        de["P.value"],
        errors="coerce"
    )
    .dropna()
    .to_numpy()
)

p = p[
    (p >= 0)
    & (p <= 1)
]


fig = plt.figure(
    figsize=(8, 6)
)

plt.hist(
    p,
    bins=40
)

expected_per_bin = len(p) / 40

plt.axhline(
    expected_per_bin,
    linestyle="--",
    label="Expected under uniform null"
)

plt.xlabel("Nominal P-value")
plt.ylabel("Number of genes")

plt.title(
    "ENVI-imputed genes: P-value distribution\n"
    "SCZ vs NTC"
)

plt.legend()

plt.tight_layout()

fig.savefig(
    OUT / "01_imputed34687_pvalue_histogram.png",
    dpi=300,
)

plt.close(fig)


# =============================================================================
# 2. QQ plot
# =============================================================================

observed = np.sort(
    np.clip(
        p,
        1e-300,
        1.0
    )
)

n = len(observed)

expected = (
    np.arange(
        1,
        n + 1
    )
    - 0.5
) / n


expected_log = -np.log10(
    expected
)

observed_log = -np.log10(
    observed
)


fig = plt.figure(
    figsize=(7, 7)
)

plt.scatter(
    expected_log,
    observed_log,
    s=8,
    alpha=0.7
)

lim = max(
    expected_log.max(),
    observed_log.max()
)

plt.plot(
    [0, lim],
    [0, lim],
    linestyle="--"
)

plt.xlabel(
    "Expected -log10(P)"
)

plt.ylabel(
    "Observed -log10(P)"
)

plt.title(
    "QQ plot: 34,687 ENVI-imputed genes"
)

plt.tight_layout()

fig.savefig(
    OUT / "02_imputed34687_QQ_plot.png",
    dpi=300,
)

plt.close(fig)


# =============================================================================
# 3. Donor PCA
#
# Standardize each gene across the 23 donors.
# PCA performed by SVD.
# =============================================================================

X = log2cpm.copy()

gene_sd = X.std(
    axis=0,
    ddof=1
)

keep = (
    np.isfinite(gene_sd)
    & (gene_sd > 0)
)

X = X[:, keep]


gene_mean = X.mean(
    axis=0
)

gene_sd = X.std(
    axis=0,
    ddof=1
)

Z = (
    X - gene_mean
) / gene_sd


U, S, VT = np.linalg.svd(
    Z,
    full_matrices=False
)

scores = U * S

variance = S ** 2

variance_ratio = (
    variance
    / variance.sum()
)


pca = pd.DataFrame({
    "donor": pb.index,
    "Dx": dx.loc[
        pb.index,
        "Dx"
    ].to_numpy(),
    "PC1": scores[:, 0],
    "PC2": scores[:, 1],
    "library_size": library,
})


pca.to_csv(
    OUT
    / "03_ENVI_full34987_donor_PCA_coordinates.csv",
    index=False,
)


fig = plt.figure(
    figsize=(8, 7)
)


for diagnosis in [
    "NTC",
    "SCZ"
]:

    subset = pca[
        pca["Dx"] == diagnosis
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        s=65,
        label=diagnosis,
    )


for _, row in pca.iterrows():

    plt.annotate(
        row["donor"],
        (
            row["PC1"],
            row["PC2"]
        ),
        xytext=(4, 3),
        textcoords="offset points",
        fontsize=7,
    )


plt.xlabel(
    f"PC1 ({variance_ratio[0] * 100:.1f}% variance)"
)

plt.ylabel(
    f"PC2 ({variance_ratio[1] * 100:.1f}% variance)"
)

plt.title(
    "Donor PCA — ENVI full transcriptome"
)

plt.legend()

plt.tight_layout()

fig.savefig(
    OUT
    / "03_ENVI_full34987_donor_PCA.png",
    dpi=300,
)

plt.close(fig)


# =============================================================================
# 4. Top-20 imputed genes — donor-level expression
# =============================================================================

de_sorted = de.sort_values(
    [
        "P.value",
        "FDR"
    ]
).reset_index(
    drop=True
)

top20 = (
    de_sorted[
        "gene"
    ]
    .head(20)
    .tolist()
)


missing = [
    gene
    for gene in top20
    if gene not in log2cpm_df.columns
]

if missing:

    raise RuntimeError(
        f"Top genes missing from pseudobulk: {missing}"
    )


long_rows = []

for gene in top20:

    for donor in log2cpm_df.index:

        long_rows.append({
            "gene": gene,
            "donor": donor,
            "Dx": dx.loc[
                donor,
                "Dx"
            ],
            "log2CPM": log2cpm_df.loc[
                donor,
                gene
            ],
        })


top_expr = pd.DataFrame(
    long_rows
)


top_expr.to_csv(
    OUT
    / "04_top20_imputed_genes_donor_expression.csv",
    index=False,
)


# Group-level summary
top_summary = (
    top_expr
    .groupby(
        [
            "gene",
            "Dx"
        ]
    )["log2CPM"]
    .agg(
        [
            "mean",
            "std",
            "median",
            "min",
            "max",
        ]
    )
    .reset_index()
)


top_summary.to_csv(
    OUT
    / "04_top20_imputed_genes_group_summary.csv",
    index=False,
)


# =============================================================================
# 5. Statistical diagnostic summary
# =============================================================================

n_genes = len(p)

n_p05 = int(
    (p < 0.05).sum()
)

n_p01 = int(
    (p < 0.01).sum()
)

n_p001 = int(
    (p < 0.001).sum()
)


expected_p05 = (
    0.05
    * n_genes
)

expected_p01 = (
    0.01
    * n_genes
)

expected_p001 = (
    0.001
    * n_genes
)


summary_lines = [
    "=" * 80,
    "ENVI IMPUTED-GENE DE DIAGNOSTICS",
    "=" * 80,
    "",
    f"Genes tested: {n_genes:,}",
    "",
    "Observed vs uniform-null expectation:",
    (
        f"P < 0.05  : observed={n_p05:,}, "
        f"expected≈{expected_p05:.1f}"
    ),
    (
        f"P < 0.01  : observed={n_p01:,}, "
        f"expected≈{expected_p01:.1f}"
    ),
    (
        f"P < 0.001 : observed={n_p001:,}, "
        f"expected≈{expected_p001:.1f}"
    ),
    "",
    (
        f"Smallest P-value: "
        f"{p.min():.8g}"
    ),
    "",
    (
        f"PC1 variance explained: "
        f"{variance_ratio[0] * 100:.2f}%"
    ),
    (
        f"PC2 variance explained: "
        f"{variance_ratio[1] * 100:.2f}%"
    ),
]


summary_text = "\n".join(
    summary_lines
)


print("\n" + summary_text)


with open(
    OUT / "diagnostic_summary.txt",
    "w"
) as f:

    f.write(
        summary_text
        + "\n"
    )


print("\nOutputs:")
for path in sorted(
    OUT.iterdir()
):

    print(
        " ",
        path.name
    )


print("\nSUCCESS")
