#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

INFILE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/"
    / "ENVI_L6_Ex_all141_withinFDR05.csv"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/L6_validation"
)

OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INFILE)

required = {
    "gene",
    "P.value",
    "FDR_within_celltype",
    "FDR_global_10celltypes",
    "log2FC",
    "direction",
}

missing = required - set(df.columns)

if missing:
    raise RuntimeError(
        f"Missing columns: {sorted(missing)}"
    )

if len(df) != 141:
    raise RuntimeError(
        f"Expected 141 L6 genes, found {len(df)}"
    )

# ------------------------------------------------------------------
# Explicit direction based on SCZ - NTC coefficient
# ------------------------------------------------------------------

df["DE_direction"] = np.where(
    df["log2FC"] > 0,
    "Higher_in_SCZ",
    "Higher_in_NTC"
)

scz = (
    df[df["log2FC"] > 0]
    .copy()
    .sort_values(
        ["FDR_within_celltype", "P.value"]
    )
)

ntc = (
    df[df["log2FC"] < 0]
    .copy()
    .sort_values(
        ["FDR_within_celltype", "P.value"]
    )
)

assert len(scz) == 77
assert len(ntc) == 64


# ------------------------------------------------------------------
# Add effect-size ranking
# ------------------------------------------------------------------

for x in [df, scz, ntc]:
    x["abs_log2FC"] = x["log2FC"].abs()


# ------------------------------------------------------------------
# Save complete tables
# ------------------------------------------------------------------

scz.to_csv(
    OUT / "L6_Ex_77_higher_in_SCZ.csv",
    index=False
)

ntc.to_csv(
    OUT / "L6_Ex_64_higher_in_NTC.csv",
    index=False
)


# ------------------------------------------------------------------
# Plain gene lists for enrichment
# ------------------------------------------------------------------

scz["gene"].to_csv(
    OUT / "L6_Ex_77_higher_in_SCZ_genes.txt",
    index=False,
    header=False
)

ntc["gene"].to_csv(
    OUT / "L6_Ex_64_higher_in_NTC_genes.txt",
    index=False,
    header=False
)


# ------------------------------------------------------------------
# Rank all 141 by statistical evidence
# ------------------------------------------------------------------

ranked_fdr = df.sort_values(
    ["FDR_within_celltype", "P.value"]
).copy()

ranked_fdr.insert(
    0,
    "validation_rank",
    range(1, len(ranked_fdr) + 1)
)

ranked_fdr.to_csv(
    OUT / "L6_Ex_all141_ranked_by_FDR.csv",
    index=False
)


# ------------------------------------------------------------------
# Rank by effect size
# ------------------------------------------------------------------

ranked_effect = df.sort_values(
    "abs_log2FC",
    ascending=False
).copy()

ranked_effect.insert(
    0,
    "effect_rank",
    range(1, len(ranked_effect) + 1)
)

ranked_effect.to_csv(
    OUT / "L6_Ex_all141_ranked_by_effect.csv",
    index=False
)


# ------------------------------------------------------------------
# Advisor-facing table
# ------------------------------------------------------------------

advisor = ranked_fdr[
    [
        "gene",
        "P.value",
        "FDR_within_celltype",
        "log2FC",
        "DE_direction",
    ]
].copy()

advisor.columns = [
    "gene",
    "P-value",
    "FDR",
    "log2FC",
    "direction",
]

advisor.to_csv(
    OUT / "L6_Ex_141_advisor_table.csv",
    index=False
)


# ------------------------------------------------------------------
# Print summary
# ------------------------------------------------------------------

print("=" * 100)
print("L6 EX BIOLOGICAL VALIDATION INPUT")
print("=" * 100)

print("\nTotal significant genes :", len(df))
print("Higher in SCZ          :", len(scz))
print("Higher in NTC          :", len(ntc))

print("\nNOTE:")
print(
    "All 141 satisfy within-L6 FDR < 0.05, "
    "but none satisfy global 10-cell-type FDR < 0.05."
)

print("\n" + "=" * 100)
print("TOP 20 SCZ-HIGHER GENES")
print("=" * 100)

print(
    scz[
        [
            "gene",
            "P.value",
            "FDR_within_celltype",
            "log2FC",
        ]
    ]
    .head(20)
    .to_string(index=False)
)

print("\n" + "=" * 100)
print("TOP 20 NTC-HIGHER GENES")
print("=" * 100)

print(
    ntc[
        [
            "gene",
            "P.value",
            "FDR_within_celltype",
            "log2FC",
        ]
    ]
    .head(20)
    .to_string(index=False)
)

print("\nSaved:")
print(OUT)
