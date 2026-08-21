#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BASE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results"
)

PB_DIR = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/pseudobulk"
)

SIG_FILE = (
    BASE
    / "significant_candidates"
    / "ENVI_L6_Ex_all141_withinFDR05.csv"
)

META_FILE = (
    ROOT
    / "data/metadata/patient_xenium_visium_24_common_with_dx.csv"
)

GENE_INFO = (
    PB_DIR
    / "ENVI_full34987_gene_info.csv"
)

OUT = (
    BASE
    / "significant_candidates/L6_validation/robustness_qc"
)

OUT.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------
# Utility: read CSV even if .gz extension is misleading
# ------------------------------------------------------------------

def read_csv_safe(path):
    path = Path(path)

    with open(path, "rb") as f:
        magic = f.read(2)

    compression = "gzip" if magic == b"\x1f\x8b" else None

    return pd.read_csv(
        path,
        compression=compression
    )


# ------------------------------------------------------------------
# Locate combined DE file
# ------------------------------------------------------------------

de_candidates = [
    BASE / "ENVI_10celltypes_imputed34687_DE_combined.csv",
    BASE / "ENVI_10celltypes_imputed34687_DE_combined.csv.gz",
]

DE_FILE = next(
    (p for p in de_candidates if p.exists()),
    None
)

if DE_FILE is None:
    raise FileNotFoundError(
        "Cannot find combined 10-cell-type DE file."
    )


# ------------------------------------------------------------------
# Locate L6 pseudobulk file automatically
# ------------------------------------------------------------------

pb_candidates = [
    p for p in PB_DIR.glob("*.csv*")
    if "l6" in p.name.lower()
    and "ex" in p.name.lower()
    and "gene_info" not in p.name.lower()
]

if not pb_candidates:
    raise FileNotFoundError(
        f"No L6 Ex pseudobulk CSV found in {PB_DIR}"
    )

print("=" * 110)
print("L6 EX ROBUSTNESS DIAGNOSTICS")
print("=" * 110)

print("\nL6 pseudobulk candidates:")
for p in pb_candidates:
    print(" ", p.name)


# Pick candidate that looks like donor x gene matrix
PB_FILE = None
pb = None

for p in pb_candidates:

    try:
        x = read_csv_safe(p)

        if 20 <= len(x) <= 30 and x.shape[1] > 1000:
            PB_FILE = p
            pb = x
            break

    except Exception as e:
        print(
            "Skipping:",
            p.name,
            "|",
            str(e)
        )


if PB_FILE is None:
    raise RuntimeError(
        "Could not identify the 23-donor L6 pseudobulk matrix."
    )


print("\nSelected pseudobulk:")
print(PB_FILE)

print("Shape:", pb.shape)


# ==================================================================
# Load DE results
# ==================================================================

de = read_csv_safe(DE_FILE)

l6_de = de[
    de["cell_type"] == "L6 Ex"
].copy()

if "gene_source" in l6_de.columns:
    l6_de = l6_de[
        l6_de["gene_source"] == "envi_imputed"
    ].copy()

print("\nL6 DE genes:", len(l6_de))

if len(l6_de) != 34687:
    print(
        "WARNING: expected 34,687 imputed L6 genes, found",
        len(l6_de)
    )


# ==================================================================
# Load 141 significant genes
# ==================================================================

sig = pd.read_csv(SIG_FILE)

if len(sig) != 141:
    raise RuntimeError(
        f"Expected 141 significant L6 genes, found {len(sig)}"
    )

sig = sig.sort_values(
    [
        "FDR_within_celltype",
        "P.value",
    ]
).reset_index(drop=True)

print("L6 within-FDR genes:", len(sig))


# ==================================================================
# P-VALUE HISTOGRAM
# ==================================================================

pvals = (
    l6_de["P.value"]
    .dropna()
    .astype(float)
)

plt.figure(figsize=(8, 6))

plt.hist(
    pvals,
    bins=np.linspace(0, 1, 21),
    edgecolor="black"
)

plt.xlabel("Nominal P-value")
plt.ylabel("Number of genes")
plt.title(
    "L6 Ex imputed-gene DE: P-value distribution"
)

plt.tight_layout()

plt.savefig(
    OUT / "L6_Ex_pvalue_histogram.png",
    dpi=300
)

plt.close()


# ==================================================================
# QQ PLOT
# ==================================================================

observed = np.sort(
    np.clip(
        pvals.to_numpy(),
        1e-300,
        1
    )
)

n = len(observed)

expected = (
    np.arange(1, n + 1) - 0.5
) / n

x = -np.log10(expected)
y = -np.log10(observed)

limit = max(
    x.max(),
    y.max()
)

plt.figure(figsize=(7, 7))

plt.scatter(
    x,
    y,
    s=9,
    alpha=0.65
)

plt.plot(
    [0, limit],
    [0, limit],
    linewidth=1.5
)

plt.xlabel("Expected -log10(P)")
plt.ylabel("Observed -log10(P)")
plt.title(
    "L6 Ex imputed-gene DE: QQ plot"
)

plt.tight_layout()

plt.savefig(
    OUT / "L6_Ex_QQ_plot.png",
    dpi=300
)

plt.close()


# ==================================================================
# IDENTIFY DONOR COLUMN
# ==================================================================

if pb.columns[0].startswith("Unnamed"):
    pb = pb.rename(
        columns={
            pb.columns[0]: "donor"
        }
    )


donor_candidates = [
    "donor",
    "BrNum",
    "sample",
    "sample_id",
]

DONOR_COL = next(
    (
        c for c in donor_candidates
        if c in pb.columns
    ),
    None
)

if DONOR_COL is None:
    raise RuntimeError(
        "Could not identify donor column in pseudobulk."
    )

pb[DONOR_COL] = (
    pb[DONOR_COL]
    .astype(str)
    .str.strip()
)


# ==================================================================
# METADATA
# ==================================================================

meta = pd.read_csv(META_FILE)

print("\nMetadata columns:")
print(meta.columns.tolist())

# Metadata uses:
# Brnumbr = donor ID
# AGE     = age
# SEX     = sex
# Dx      = diagnosis

required_meta = [
    "Brnumbr",
    "Dx",
    "AGE",
    "SEX",
]

missing_meta = [
    c for c in required_meta
    if c not in meta.columns
]

if missing_meta:
    raise RuntimeError(
        f"Missing metadata columns: {missing_meta}. "
        f"Available columns: {meta.columns.tolist()}"
    )

meta = meta[
    [
        "Brnumbr",
        "Dx",
        "AGE",
        "SEX",
    ]
].copy()

meta["Brnumbr"] = (
    meta["Brnumbr"]
    .astype(str)
    .str.strip()
)

# Standardize column names for downstream code
meta = meta.rename(
    columns={
        "Brnumbr": "BrNum",
        "AGE": "Age",
        "SEX": "Sex",
    }
)

pb = pb.merge(
    meta,
    left_on=DONOR_COL,
    right_on="BrNum",
    how="left",
    validate="one_to_one"
)

if pb["Dx"].isna().any():
    print(
        "\nWARNING: missing metadata for donors:"
    )
    print(
        pb.loc[
            pb["Dx"].isna(),
            DONOR_COL
        ].tolist()
    )
    raise RuntimeError(
        "Some pseudobulk donors did not match metadata."
    )

print("\nDiagnosis:")
print(
    pb["Dx"]
    .value_counts()
    .to_string()
)


# ==================================================================
# GENE COLUMNS
# ==================================================================

gene_info = pd.read_csv(
    GENE_INFO
)

all_genes = (
    gene_info["gene"]
    .astype(str)
    .tolist()
)

gene_cols = [
    g for g in all_genes
    if g in pb.columns
]

print("\nGene columns in pseudobulk:", len(gene_cols))

if len(gene_cols) != 34987:
    print(
        "WARNING: expected 34,987 total genes."
    )


# ==================================================================
# Reproduce log2(CPM + 1)
# ==================================================================

counts = (
    pb[gene_cols]
    .to_numpy(dtype=float)
)

library_size = counts.sum(axis=1)

if np.any(library_size <= 0):
    raise RuntimeError(
        "At least one donor has zero pseudobulk library size."
    )

cpm = (
    counts
    / library_size[:, None]
    * 1e6
)

logcpm = np.log2(
    cpm + 1
)

logcpm_df = pd.DataFrame(
    logcpm,
    index=pb[DONOR_COL].values,
    columns=gene_cols
)


# ==================================================================
# PCA
# ==================================================================

variance = logcpm_df.var(axis=0)

top_var_genes = (
    variance
    .sort_values(ascending=False)
    .head(min(2000, len(variance)))
    .index
)

X = (
    logcpm_df[top_var_genes]
    .to_numpy()
)

# Center genes
X = X - X.mean(axis=0)

pca = PCA(
    n_components=2
)

pcs = pca.fit_transform(X)

pca_df = pd.DataFrame(
    {
        "donor": pb[DONOR_COL].values,
        "Dx": pb["Dx"].values,
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1],
    }
)

pca_df.to_csv(
    OUT / "L6_Ex_donor_PCA_coordinates.csv",
    index=False
)


plt.figure(figsize=(8, 7))

for dx, group in pca_df.groupby("Dx"):

    plt.scatter(
        group["PC1"],
        group["PC2"],
        s=60,
        label=dx
    )

    for _, r in group.iterrows():

        plt.annotate(
            r["donor"],
            (
                r["PC1"],
                r["PC2"]
            ),
            fontsize=7,
            xytext=(3, 3),
            textcoords="offset points"
        )


plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}%)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}%)"
)

plt.title(
    "L6 Ex donor pseudobulk PCA"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    OUT / "L6_Ex_donor_PCA.png",
    dpi=300
)

plt.close()


print("\nPCA:")
print(
    "PC1 variance:",
    f"{pca.explained_variance_ratio_[0] * 100:.2f}%"
)

print(
    "PC2 variance:",
    f"{pca.explained_variance_ratio_[1] * 100:.2f}%"
)


# ==================================================================
# DESIGN MATRIX FOR LEAVE-ONE-DONOR-OUT CHECK
# ==================================================================

analysis_meta = pb[
    [
        DONOR_COL,
        "Dx",
        "Age",
        "Sex",
    ]
].copy()


analysis_meta["Dx_SCZ"] = (
    analysis_meta["Dx"]
    .astype(str)
    .str.upper()
    .eq("SCZ")
    .astype(float)
)


analysis_meta["Age_numeric"] = pd.to_numeric(
    analysis_meta["Age"],
    errors="raise"
)

analysis_meta["Age_centered"] = (
    analysis_meta["Age_numeric"]
    -
    analysis_meta["Age_numeric"].mean()
)


sex_values = sorted(
    analysis_meta["Sex"]
    .astype(str)
    .unique()
)

if len(sex_values) != 2:
    raise RuntimeError(
        f"Expected 2 Sex categories, found {sex_values}"
    )

analysis_meta["Sex_dummy"] = (
    analysis_meta["Sex"]
    .astype(str)
    .eq(sex_values[1])
    .astype(float)
)


X_full = np.column_stack(
    [
        np.ones(len(analysis_meta)),
        analysis_meta["Dx_SCZ"],
        analysis_meta["Age_centered"],
        analysis_meta["Sex_dummy"],
    ]
)


def fit_dx_beta(X, y):

    beta = np.linalg.lstsq(
        X,
        y,
        rcond=None
    )[0]

    return float(beta[1])


# ==================================================================
# LEAVE-ONE-DONOR-OUT FOR ALL 141 GENES
# ==================================================================

loo_rows = []

for _, r in sig.iterrows():

    gene = r["gene"]

    if gene not in logcpm_df.columns:
        continue

    y = logcpm_df[
        gene
    ].to_numpy()

    full_beta = fit_dx_beta(
        X_full,
        y
    )

    loo_betas = []
    donor_names = []


    for i in range(len(y)):

        keep = np.ones(
            len(y),
            dtype=bool
        )

        keep[i] = False

        beta = fit_dx_beta(
            X_full[keep],
            y[keep]
        )

        loo_betas.append(beta)

        donor_names.append(
            analysis_meta.iloc[i][DONOR_COL]
        )


    loo_betas = np.asarray(
        loo_betas
    )

    original_sign = np.sign(
        full_beta
    )

    flip_mask = (
        np.sign(loo_betas)
        != original_sign
    )

    delta = np.abs(
        loo_betas - full_beta
    )

    influential_i = int(
        np.argmax(delta)
    )


    loo_rows.append(
        {
            "gene": gene,

            "reported_log2FC":
                r["log2FC"],

            "recomputed_log2FC":
                full_beta,

            "FDR_within_celltype":
                r["FDR_within_celltype"],

            "loo_min_log2FC":
                loo_betas.min(),

            "loo_max_log2FC":
                loo_betas.max(),

            "loo_mean_log2FC":
                loo_betas.mean(),

            "loo_sd_log2FC":
                loo_betas.std(ddof=1),

            "n_direction_flips":
                int(flip_mask.sum()),

            "direction_stable_all_23":
                bool(
                    flip_mask.sum()
                    == 0
                ),

            "max_abs_beta_change":
                delta.max(),

            "most_influential_donor":
                donor_names[
                    influential_i
                ],
        }
    )


loo = pd.DataFrame(
    loo_rows
)

loo = loo.sort_values(
    [
        "n_direction_flips",
        "max_abs_beta_change",
    ],
    ascending=[
        True,
        True,
    ]
)

loo.to_csv(
    OUT / "L6_Ex_141_leave_one_donor_out.csv",
    index=False
)


# ==================================================================
# TOP 10 DONOR-LEVEL PLOTS
# ==================================================================

top10 = (
    sig.head(10)
    .copy()
)


plot_summary = []


for rank, (_, r) in enumerate(
    top10.iterrows(),
    start=1
):

    gene = r["gene"]

    if gene not in logcpm_df.columns:
        continue


    vals = pd.DataFrame(
        {
            "donor": pb[DONOR_COL].values,
            "Dx": pb["Dx"].values,
            "expression": logcpm_df[
                gene
            ].values,
        }
    )


    ntc = vals[
        vals["Dx"] == "NTC"
    ]["expression"].to_numpy()

    scz = vals[
        vals["Dx"] == "SCZ"
    ]["expression"].to_numpy()


    plt.figure(
        figsize=(7, 6)
    )

    plt.boxplot(
        [
            ntc,
            scz,
        ],
        labels=[
            "NTC",
            "SCZ",
        ],
        showfliers=False
    )


    rng = np.random.default_rng(
        12345
    )


    for pos, dx in [
        (1, "NTC"),
        (2, "SCZ"),
    ]:

        z = vals[
            vals["Dx"] == dx
        ]

        jitter = rng.normal(
            0,
            0.035,
            len(z)
        )

        plt.scatter(
            np.full(len(z), pos)
            + jitter,
            z["expression"],
            s=35
        )

        for xx, (_, rr) in zip(
            np.full(len(z), pos) + jitter,
            z.iterrows()
        ):

            plt.annotate(
                rr["donor"],
                (
                    xx,
                    rr["expression"]
                ),
                fontsize=6,
                xytext=(2, 2),
                textcoords="offset points"
            )


    plt.ylabel(
        "Donor pseudobulk log2(CPM + 1)"
    )

    plt.title(
        f"{gene} | L6 Ex\n"
        f"log2FC={r['log2FC']:.3f}, "
        f"FDR={r['FDR_within_celltype']:.4g}"
    )

    plt.tight_layout()

    plt.savefig(
        OUT
        / f"L6_top{rank:02d}_{gene}_donor_expression.png",
        dpi=300
    )

    plt.close()


    plot_summary.append(
        {
            "rank": rank,
            "gene": gene,
            "log2FC": r["log2FC"],
            "FDR": r["FDR_within_celltype"],
        }
    )


pd.DataFrame(
    plot_summary
).to_csv(
    OUT / "L6_Ex_top10_plotted_genes.csv",
    index=False
)


# ==================================================================
# FINAL SUMMARY
# ==================================================================

n_stable = int(
    loo["direction_stable_all_23"].sum()
)

n_unstable = len(loo) - n_stable

print("\n" + "=" * 110)
print("LEAVE-ONE-DONOR-OUT SUMMARY")
print("=" * 110)

print(
    "Genes tested                     :",
    len(loo)
)

print(
    "Direction stable after every LOO :",
    n_stable
)

print(
    "At least one direction flip      :",
    n_unstable
)


print("\nMost donor-sensitive genes:")

print(
    loo[
        [
            "gene",
            "reported_log2FC",
            "loo_min_log2FC",
            "loo_max_log2FC",
            "n_direction_flips",
            "max_abs_beta_change",
            "most_influential_donor",
        ]
    ]
    .sort_values(
        "max_abs_beta_change",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)


# Recomputed vs reported coefficient check

merged_beta = loo.dropna(
    subset=[
        "reported_log2FC",
        "recomputed_log2FC",
    ]
)

max_beta_diff = (
    merged_beta[
        "reported_log2FC"
    ]
    -
    merged_beta[
        "recomputed_log2FC"
    ]
).abs().max()


print(
    "\nMax |reported log2FC - recomputed log2FC|:",
    max_beta_diff
)


print("\nOutputs:")
print(OUT)

print("\nGenerated:")
print("  L6_Ex_pvalue_histogram.png")
print("  L6_Ex_QQ_plot.png")
print("  L6_Ex_donor_PCA.png")
print("  L6_Ex_141_leave_one_donor_out.csv")
print("  10 donor-level top-gene plots")

print("\nSUCCESS")
