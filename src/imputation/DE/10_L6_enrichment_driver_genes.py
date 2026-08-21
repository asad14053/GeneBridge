#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BASE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/L6_validation"
)

ENRICH = (
    BASE
    / "pathway_enrichment"
    / "L6_Ex_pathway_enrichment_FDR05.csv"
)

# Allow either filename, since both versions were used during extraction.
candidate_files = [
    BASE / "L6_Ex_all141_ranked_by_FDR.csv",
    BASE.parent / "ENVI_L6_Ex_all141_withinFDR05.csv",
    BASE.parent / "ENVI_L6_Ex_withinFDR05_all.csv",
]

GENES = next(
    (p for p in candidate_files if p.exists()),
    None
)

if GENES is None:
    raise FileNotFoundError(
        "Could not find the L6 141-gene table.\nTried:\n"
        + "\n".join(str(p) for p in candidate_files)
    )

if not ENRICH.exists():
    raise FileNotFoundError(f"Missing enrichment file: {ENRICH}")

OUT = BASE / "driver_genes"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 110)
print("L6 EX ENRICHMENT DRIVER-GENE ANALYSIS")
print("=" * 110)

print("\nGene table:")
print(GENES)

print("\nSignificant pathways:")
print(ENRICH)


# ---------------------------------------------------------------------
# Read data
# ---------------------------------------------------------------------

genes = pd.read_csv(GENES)
enrich = pd.read_csv(ENRICH)

print("\nL6 candidate genes:", len(genes))
print("Significant pathway rows:", len(enrich))


# ---------------------------------------------------------------------
# Normalize direction labels in gene-level table
# ---------------------------------------------------------------------

if "DE_direction" in genes.columns:
    genes["direction_clean"] = genes["DE_direction"]

elif "direction" in genes.columns:
    genes["direction_clean"] = genes["direction"].replace(
        {
            "Higher SCZ": "Higher_in_SCZ",
            "Higher NTC": "Higher_in_NTC",
            "Higher_in_SCZ": "Higher_in_SCZ",
            "Higher_in_NTC": "Higher_in_NTC",
        }
    )

else:
    genes["direction_clean"] = genes["log2FC"].apply(
        lambda x: "Higher_in_SCZ" if x > 0 else "Higher_in_NTC"
    )


# ---------------------------------------------------------------------
# Expand pathway -> gene memberships
# ---------------------------------------------------------------------

rows = []

for _, r in enrich.iterrows():

    gene_list = [
        g.strip()
        for g in str(r["genes"]).split("/")
        if g.strip()
    ]

    for gene in gene_list:

        rows.append(
            {
                "gene": gene,
                "direction": r["direction"],
                "database": r["database"],
                "pathway": r["pathway"],
                "pathway_P": r["P-value"],
                "pathway_FDR": r["FDR"],
            }
        )

memberships = pd.DataFrame(rows)

if memberships.empty:
    raise RuntimeError("No pathway-gene memberships were recovered.")


# ---------------------------------------------------------------------
# Verify direction consistency
# ---------------------------------------------------------------------

direction_map = (
    genes[["gene", "direction_clean"]]
    .drop_duplicates("gene")
    .set_index("gene")["direction_clean"]
    .to_dict()
)

memberships["gene_direction"] = memberships["gene"].map(direction_map)

bad = memberships[
    memberships["gene_direction"].notna()
    &
    (memberships["direction"] != memberships["gene_direction"])
]

if len(bad):
    print("\nWARNING: direction mismatches found:")
    print(
        bad[
            ["gene", "direction", "gene_direction", "pathway"]
        ].to_string(index=False)
    )
else:
    print("\nDirection consistency: PASS")


# ---------------------------------------------------------------------
# Count pathway recurrence for each gene
# ---------------------------------------------------------------------

driver = (
    memberships
    .groupby(["gene", "direction"], as_index=False)
    .agg(
        n_significant_pathways=("pathway", "nunique"),
        n_databases=("database", "nunique"),
        min_pathway_P=("pathway_P", "min"),
        min_pathway_FDR=("pathway_FDR", "min"),
    )
)


# Database-specific counts
db_counts = (
    memberships
    .groupby(["gene", "direction", "database"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

for col in ["GO_BP", "Reactome", "Hallmark"]:
    if col not in db_counts.columns:
        db_counts[col] = 0

db_counts = db_counts.rename(
    columns={
        "GO_BP": "n_GO_BP",
        "Reactome": "n_Reactome",
        "Hallmark": "n_Hallmark",
    }
)

driver = driver.merge(
    db_counts[
        [
            "gene",
            "direction",
            "n_GO_BP",
            "n_Reactome",
            "n_Hallmark",
        ]
    ],
    on=["gene", "direction"],
    how="left",
)


# ---------------------------------------------------------------------
# Add pathway names
# ---------------------------------------------------------------------

pathway_names = (
    memberships
    .groupby(["gene", "direction"])["pathway"]
    .apply(lambda x: " | ".join(sorted(set(x))))
    .reset_index(name="significant_pathways")
)

driver = driver.merge(
    pathway_names,
    on=["gene", "direction"],
    how="left",
)


# ---------------------------------------------------------------------
# Add gene-level DE statistics
# ---------------------------------------------------------------------

gene_cols = ["gene", "log2FC"]

for col in [
    "P.value",
    "FDR_within_celltype",
    "FDR_global_10celltypes",
]:
    if col in genes.columns:
        gene_cols.append(col)

gene_stats = (
    genes[gene_cols]
    .drop_duplicates("gene")
)

driver = driver.merge(
    gene_stats,
    on="gene",
    how="left",
)

driver["abs_log2FC"] = driver["log2FC"].abs()


# ---------------------------------------------------------------------
# Rank genes
#
# Primary ranking:
# 1. Number of significant pathways
# 2. Best pathway FDR
# 3. Gene-level within-L6 FDR
# 4. Absolute log2FC
# ---------------------------------------------------------------------

sort_cols = [
    "n_significant_pathways",
    "min_pathway_FDR",
]

ascending = [
    False,
    True,
]

if "FDR_within_celltype" in driver.columns:
    sort_cols.append("FDR_within_celltype")
    ascending.append(True)

sort_cols.append("abs_log2FC")
ascending.append(False)

driver = (
    driver
    .sort_values(sort_cols, ascending=ascending)
    .reset_index(drop=True)
)

driver.insert(
    0,
    "driver_rank",
    range(1, len(driver) + 1)
)


# ---------------------------------------------------------------------
# Save all driver genes
# ---------------------------------------------------------------------

driver.to_csv(
    OUT / "L6_Ex_significant_pathway_driver_genes.csv",
    index=False
)

scz = driver[
    driver["direction"] == "Higher_in_SCZ"
].copy()

ntc = driver[
    driver["direction"] == "Higher_in_NTC"
].copy()

scz.to_csv(
    OUT / "L6_Ex_SCZ_higher_driver_genes.csv",
    index=False
)

ntc.to_csv(
    OUT / "L6_Ex_NTC_higher_driver_genes.csv",
    index=False
)


# ---------------------------------------------------------------------
# Save pathway-gene membership table
# ---------------------------------------------------------------------

memberships.to_csv(
    OUT / "L6_Ex_significant_pathway_gene_memberships.csv",
    index=False
)


# ---------------------------------------------------------------------
# Print compact advisor-facing tables
# ---------------------------------------------------------------------

show_cols = [
    "gene",
    "n_significant_pathways",
    "n_GO_BP",
    "n_Reactome",
    "min_pathway_FDR",
    "log2FC",
]

if "FDR_within_celltype" in driver.columns:
    show_cols.append("FDR_within_celltype")


print("\n" + "=" * 110)
print("SCZ-HIGHER — TOP DRIVER GENES")
print("=" * 110)

print(
    scz[show_cols]
    .head(30)
    .to_string(index=False)
)


print("\n" + "=" * 110)
print("NTC-HIGHER / SCZ-LOWER — TOP DRIVER GENES")
print("=" * 110)

print(
    ntc[show_cols]
    .head(30)
    .to_string(index=False)
)


# ---------------------------------------------------------------------
# Basic summary
# ---------------------------------------------------------------------

print("\n" + "=" * 110)
print("SUMMARY")
print("=" * 110)

print(
    "Unique genes driving significant SCZ-higher pathways:",
    scz["gene"].nunique()
)

print(
    "Unique genes driving significant NTC-higher pathways:",
    ntc["gene"].nunique()
)

print(
    "Total unique pathway-driver genes:",
    driver["gene"].nunique()
)

print("\nSCZ-higher most recurrent:")
print(
    scz[
        ["gene", "n_significant_pathways"]
    ]
    .head(15)
    .to_string(index=False)
)

print("\nNTC-higher / SCZ-lower most recurrent:")
print(
    ntc[
        ["gene", "n_significant_pathways"]
    ]
    .head(15)
    .to_string(index=False)
)

print("\nSaved to:")
print(OUT)

print("\nSUCCESS")
