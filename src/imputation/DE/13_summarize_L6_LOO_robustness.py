#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


# =============================================================================
# Paths
# =============================================================================

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

LOO_FILE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/L6_validation/robustness_qc/"
    / "L6_Ex_141_leave_one_donor_out.csv"
)

DRIVER_FILE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/L6_validation/driver_genes/"
    / "L6_Ex_significant_pathway_driver_genes.csv"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/L6_validation/robustness_qc/"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# =============================================================================
# Check inputs
# =============================================================================

for p in [LOO_FILE, DRIVER_FILE]:
    if not p.exists():
        raise FileNotFoundError(
            f"Missing input file: {p}"
        )


# =============================================================================
# Load
# =============================================================================

loo = pd.read_csv(LOO_FILE)
drivers = pd.read_csv(DRIVER_FILE)

driver_set = set(
    drivers["gene"]
    .dropna()
    .astype(str)
)

loo["gene"] = (
    loo["gene"]
    .astype(str)
    .str.strip()
)

loo["pathway_driver"] = (
    loo["gene"].isin(driver_set)
)


# =============================================================================
# Robust boolean parsing
# =============================================================================

def parse_bool(x):

    if isinstance(x, (bool, np.bool_)):
        return bool(x)

    if pd.isna(x):
        return False

    return str(x).strip().lower() in {
        "true",
        "1",
        "yes",
        "y",
        "t",
    }


loo["direction_stable"] = (
    loo["direction_stable_all_23"]
    .apply(parse_bool)
)


# =============================================================================
# Sanity checks
# =============================================================================

required = [
    "gene",
    "reported_log2FC",
    "recomputed_log2FC",
    "loo_min_log2FC",
    "loo_max_log2FC",
    "n_direction_flips",
    "max_abs_beta_change",
    "most_influential_donor",
]

missing = [
    c for c in required
    if c not in loo.columns
]

if missing:
    raise RuntimeError(
        f"Missing LOO columns: {missing}"
    )


print("=" * 115)
print("L6 EX — LEAVE-ONE-DONOR-OUT ROBUSTNESS SUMMARY")
print("=" * 115)

print("\nInput:")
print(LOO_FILE)

print("\nGenes tested:", len(loo))

if len(loo) != 141:
    print(
        "WARNING: expected 141 genes, found",
        len(loo)
    )


# =============================================================================
# Main robustness summary
# =============================================================================

n_total = len(loo)

n_stable = int(
    loo["direction_stable"].sum()
)

n_unstable = (
    n_total - n_stable
)

pct_stable = (
    100 * n_stable / n_total
    if n_total > 0
    else np.nan
)


print("\n" + "=" * 115)
print("ALL 141 L6 GENES")
print("=" * 115)

print(
    "Genes tested                              :",
    n_total
)

print(
    "Direction stable after removing every donor:",
    n_stable
)

print(
    "At least one direction flip               :",
    n_unstable
)

print(
    "Percent direction-stable                  :",
    f"{pct_stable:.1f}%"
)


# =============================================================================
# Pathway-driver subset
# =============================================================================

drv = (
    loo[
        loo["pathway_driver"]
    ]
    .copy()
)

n_drv = len(drv)

n_drv_stable = int(
    drv["direction_stable"].sum()
)

n_drv_unstable = (
    n_drv - n_drv_stable
)

pct_drv_stable = (
    100 * n_drv_stable / n_drv
    if n_drv > 0
    else np.nan
)


print("\n" + "=" * 115)
print("37 PATHWAY-DRIVER GENES")
print("=" * 115)

print(
    "Driver genes tested                       :",
    n_drv
)

print(
    "Stable after every donor removal           :",
    n_drv_stable
)

print(
    "At least one direction flip               :",
    n_drv_unstable
)

print(
    "Percent direction-stable                  :",
    f"{pct_drv_stable:.1f}%"
)


# =============================================================================
# Reproduction of original log2FC
# =============================================================================

loo["beta_difference"] = (
    loo["reported_log2FC"]
    -
    loo["recomputed_log2FC"]
).abs()


print("\n" + "=" * 115)
print("COEFFICIENT REPRODUCTION CHECK")
print("=" * 115)

print(
    "Max |reported - recomputed log2FC| :",
    loo["beta_difference"].max()
)

print(
    "Median |reported - recomputed|     :",
    loo["beta_difference"].median()
)


# =============================================================================
# Genes that flip direction
# =============================================================================

flipped = (
    loo[
        loo["n_direction_flips"] > 0
    ]
    .copy()
)

flipped = flipped.sort_values(
    [
        "n_direction_flips",
        "max_abs_beta_change",
    ],
    ascending=[
        False,
        False,
    ]
)


print("\n" + "=" * 115)
print("GENES WITH AT LEAST ONE DIRECTION FLIP")
print("=" * 115)

flip_cols = [
    "gene",
    "reported_log2FC",
    "recomputed_log2FC",
    "loo_min_log2FC",
    "loo_max_log2FC",
    "n_direction_flips",
    "max_abs_beta_change",
    "most_influential_donor",
    "pathway_driver",
]


if len(flipped) == 0:

    print(
        "NONE — all tested genes preserve "
        "their SCZ-vs-NTC direction after "
        "removing every donor individually."
    )

else:

    print(
        flipped[
            flip_cols
        ].to_string(
            index=False
        )
    )


flipped.to_csv(
    OUT
    / "L6_Ex_genes_with_LOO_direction_flip.csv",
    index=False
)


# =============================================================================
# Most donor-sensitive genes
# =============================================================================

sensitive = (
    loo.sort_values(
        "max_abs_beta_change",
        ascending=False
    )
    .copy()
)


print("\n" + "=" * 115)
print("TOP 25 MOST DONOR-SENSITIVE GENES")
print("=" * 115)

sensitive_cols = [
    "gene",
    "reported_log2FC",
    "loo_min_log2FC",
    "loo_max_log2FC",
    "n_direction_flips",
    "max_abs_beta_change",
    "most_influential_donor",
    "pathway_driver",
]

print(
    sensitive[
        sensitive_cols
    ]
    .head(25)
    .to_string(
        index=False
    )
)


sensitive.to_csv(
    OUT
    / "L6_Ex_141_ranked_by_donor_sensitivity.csv",
    index=False
)


# =============================================================================
# Which donors are most influential?
# =============================================================================

donor_counts = (
    loo[
        "most_influential_donor"
    ]
    .value_counts()
    .rename_axis(
        "donor"
    )
    .reset_index(
        name="n_genes_most_influential"
    )
)


print("\n" + "=" * 115)
print("MOST FREQUENTLY INFLUENTIAL DONORS")
print("=" * 115)

print(
    donor_counts.to_string(
        index=False
    )
)


donor_counts.to_csv(
    OUT
    / "L6_Ex_most_influential_donor_counts.csv",
    index=False
)


# =============================================================================
# Specifically inspect Br2039
# =============================================================================

br2039_n = int(
    (
        loo["most_influential_donor"]
        .astype(str)
        == "Br2039"
    ).sum()
)


print("\nBr2039 was the most influential donor for:")
print(
    f"{br2039_n} / {n_total} genes "
    f"({100 * br2039_n / n_total:.1f}%)"
)


# =============================================================================
# Pathway-driver details
# =============================================================================

drv = drv.sort_values(
    "max_abs_beta_change",
    ascending=False
)


print("\n" + "=" * 115)
print("37 PATHWAY-DRIVER GENES — LOO DETAILS")
print("=" * 115)

print(
    drv[
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
    .to_string(
        index=False
    )
)


drv.to_csv(
    OUT
    / "L6_Ex_37_pathway_drivers_LOO_robustness.csv",
    index=False
)


# =============================================================================
# Summary table
# =============================================================================

summary = pd.DataFrame(
    [
        {
            "gene_set": "All_141_L6_genes",
            "genes_tested": n_total,
            "direction_stable": n_stable,
            "direction_unstable": n_unstable,
            "percent_stable": pct_stable,
        },
        {
            "gene_set": "Pathway_driver_genes",
            "genes_tested": n_drv,
            "direction_stable": n_drv_stable,
            "direction_unstable": n_drv_unstable,
            "percent_stable": pct_drv_stable,
        },
    ]
)


summary.to_csv(
    OUT
    / "L6_Ex_LOO_robustness_summary.csv",
    index=False
)


# =============================================================================
# Final
# =============================================================================

print("\n" + "=" * 115)
print("FINAL SUMMARY")
print("=" * 115)

print(
    f"All L6 genes: "
    f"{n_stable}/{n_total} "
    f"({pct_stable:.1f}%) direction-stable"
)

print(
    f"Pathway drivers: "
    f"{n_drv_stable}/{n_drv} "
    f"({pct_drv_stable:.1f}%) direction-stable"
)

print(
    f"Genes with at least one direction flip: "
    f"{n_unstable}"
)

print(
    f"Br2039 most influential for: "
    f"{br2039_n}/{n_total} genes"
)

print("\nSaved outputs:")
print(OUT)

print("\nSUCCESS")
