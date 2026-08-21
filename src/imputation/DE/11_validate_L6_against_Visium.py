#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------

VISIUM_FILE = (
    ROOT
    / "data/external/scz_validation/"
    / "Kwon2026_SuppTable5_layer_restricted_DE.xlsx"
)

L6_FILE = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/results/"
    / "significant_candidates/L6_validation/"
    / "L6_Ex_all141_ranked_by_FDR.csv"
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
    / "significant_candidates/L6_validation/"
    / "visium_validation"
)

OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Check files
# ---------------------------------------------------------------------

for p in [VISIUM_FILE, L6_FILE, DRIVER_FILE]:
    if not p.exists():
        raise FileNotFoundError(f"Missing file: {p}")


print("=" * 110)
print("L6 EX — ENVI vs VISIUM SpD03-L6 VALIDATION")
print("=" * 110)


# =====================================================================
# Read ENVI L6 genes
# =====================================================================

envi = pd.read_csv(L6_FILE)

if len(envi) != 141:
    raise RuntimeError(
        f"Expected 141 ENVI L6 genes, found {len(envi)}"
    )

required = {
    "gene",
    "log2FC",
    "P.value",
    "FDR_within_celltype",
}

missing = required - set(envi.columns)

if missing:
    raise RuntimeError(
        f"Missing ENVI columns: {sorted(missing)}"
    )


envi = envi.copy()

envi["ENVI_direction"] = np.where(
    envi["log2FC"] > 0,
    "Higher_in_SCZ",
    "Higher_in_NTC"
)


# =====================================================================
# Read pathway-driver genes
# =====================================================================

drivers = pd.read_csv(DRIVER_FILE)

driver_genes = set(
    drivers["gene"].dropna().astype(str)
)

print("\nENVI L6 significant genes :", len(envi))
print("Pathway-driver genes      :", len(driver_genes))


# =====================================================================
# Read Visium layer-restricted DE
# =====================================================================

visium = pd.read_excel(
    VISIUM_FILE,
    sheet_name="Layer-res.DEGs"
)

print("\nFull Visium table:", visium.shape)

print("\nSpatial domains:")
print(
    visium["PRECAST_spd"]
    .value_counts()
    .sort_index()
    .to_string()
)


# ---------------------------------------------------------------------
# Restrict to cortical L6
# ---------------------------------------------------------------------

l6_visium = (
    visium[
        visium["PRECAST_spd"] == "SpD03-L6"
    ]
    .copy()
)

print("\nVisium SpD03-L6 genes:", len(l6_visium))

if len(l6_visium) == 0:
    raise RuntimeError(
        "No SpD03-L6 rows found."
    )


# ---------------------------------------------------------------------
# Clean symbols
# ---------------------------------------------------------------------

l6_visium = l6_visium[
    l6_visium["gene"].notna()
].copy()

l6_visium["gene"] = (
    l6_visium["gene"]
    .astype(str)
    .str.strip()
)

envi["gene"] = (
    envi["gene"]
    .astype(str)
    .str.strip()
)


# Check duplicate symbols
dup = l6_visium["gene"].duplicated(keep=False)

if dup.any():

    print(
        "\nWARNING: duplicate Visium gene symbols:",
        l6_visium.loc[dup, "gene"].nunique()
    )

    # If duplicate symbols occur, retain the row with
    # smallest nominal p-value.
    l6_visium = (
        l6_visium
        .sort_values("P.Value")
        .drop_duplicates("gene", keep="first")
    )


# =====================================================================
# Rename Visium columns
# =====================================================================

visium_keep = l6_visium[
    [
        "ensembl_id",
        "gene",
        "logFC",
        "AverageExpression",
        "t-statistics",
        "P.Value",
        "adj.P.Val",
        "sig_10",
        "layer_specificity",
    ]
].copy()

visium_keep = visium_keep.rename(
    columns={
        "ensembl_id": "Visium_ensembl_id",
        "logFC": "Visium_log2FC",
        "AverageExpression": "Visium_AverageExpression",
        "t-statistics": "Visium_t",
        "P.Value": "Visium_P",
        "adj.P.Val": "Visium_FDR",
        "sig_10": "Visium_FDR10",
        "layer_specificity": "Visium_layer_specific",
    }
)


# =====================================================================
# Merge ENVI and measured Visium
# =====================================================================

merged = envi.merge(
    visium_keep,
    on="gene",
    how="left",
    validate="one_to_one"
)


merged["present_in_Visium_L6"] = (
    merged["Visium_log2FC"].notna()
)


merged["Visium_direction"] = np.where(
    merged["Visium_log2FC"].isna(),
    "Not_available",
    np.where(
        merged["Visium_log2FC"] > 0,
        "Higher_in_SCZ",
        np.where(
            merged["Visium_log2FC"] < 0,
            "Higher_in_NTC",
            "No_change"
        )
    )
)


merged["same_direction"] = (
    merged["present_in_Visium_L6"]
    &
    (
        merged["ENVI_direction"]
        ==
        merged["Visium_direction"]
    )
)


merged["Visium_nominal_P05"] = (
    merged["Visium_P"] < 0.05
)

merged["Visium_FDR05"] = (
    merged["Visium_FDR"] < 0.05
)

merged["Visium_FDR10_calc"] = (
    merged["Visium_FDR"] < 0.10
)


# ---------------------------------------------------------------------
# Concordance categories
# ---------------------------------------------------------------------

def classify(row):

    if not row["present_in_Visium_L6"]:
        return "Not_available_in_Visium"

    same = row["same_direction"]

    if row["Visium_FDR05"]:
        return (
            "FDR05_same_direction"
            if same
            else "FDR05_opposite_direction"
        )

    if row["Visium_FDR10_calc"]:
        return (
            "FDR10_same_direction"
            if same
            else "FDR10_opposite_direction"
        )

    if row["Visium_nominal_P05"]:
        return (
            "Nominal_P05_same_direction"
            if same
            else "Nominal_P05_opposite_direction"
        )

    return (
        "Same_direction_not_significant"
        if same
        else "Opposite_direction_not_significant"
    )


merged["validation_category"] = merged.apply(
    classify,
    axis=1
)


merged["pathway_driver"] = (
    merged["gene"].isin(driver_genes)
)


# =====================================================================
# Save full 141-gene comparison
# =====================================================================

merged.to_csv(
    OUT / "L6_Ex_141_ENVI_vs_Visium_L6.csv",
    index=False
)


# =====================================================================
# Summary helper
# =====================================================================

def summarize(df, label):

    matched = df[
        df["present_in_Visium_L6"]
    ].copy()

    n_total = len(df)
    n_matched = len(matched)

    n_same = int(
        matched["same_direction"].sum()
    )

    n_opposite = n_matched - n_same


    nominal = matched[
        matched["Visium_nominal_P05"]
    ].copy()

    fdr10 = matched[
        matched["Visium_FDR10_calc"]
    ].copy()

    fdr05 = matched[
        matched["Visium_FDR05"]
    ].copy()


    def same_n(x):
        return int(
            x["same_direction"].sum()
        )


    print("\n" + "=" * 110)
    print(label)
    print("=" * 110)

    print("ENVI genes                  :", n_total)
    print("Present in Visium L6        :", n_matched)

    print(
        "Overall same direction     :",
        n_same,
        f"({100*n_same/n_matched:.1f}%)"
        if n_matched else ""
    )

    print(
        "Overall opposite direction :",
        n_opposite
    )

    print("\nVisium nominal P < 0.05     :", len(nominal))
    print("  Same direction            :", same_n(nominal))
    print(
        "  Opposite direction        :",
        len(nominal) - same_n(nominal)
    )

    print("\nVisium FDR < 0.10           :", len(fdr10))
    print("  Same direction            :", same_n(fdr10))
    print(
        "  Opposite direction        :",
        len(fdr10) - same_n(fdr10)
    )

    print("\nVisium FDR < 0.05           :", len(fdr05))
    print("  Same direction            :", same_n(fdr05))
    print(
        "  Opposite direction        :",
        len(fdr05) - same_n(fdr05)
    )


    # ---------------------------------------------------------------
    # Effect-size concordance
    # ---------------------------------------------------------------

    if n_matched >= 3:

        rho, p_rho = spearmanr(
            matched["log2FC"],
            matched["Visium_log2FC"],
            nan_policy="omit"
        )

        print("\nENVI vs Visium logFC:")
        print(f"  Spearman rho              : {rho:.4f}")
        print(f"  P-value                   : {p_rho:.4g}")


    # ---------------------------------------------------------------
    # Simple sign-concordance test against 50%
    # Descriptive because cohorts are related/not independent.
    # ---------------------------------------------------------------

    if n_matched > 0:

        b = binomtest(
            n_same,
            n_matched,
            p=0.5,
            alternative="greater"
        )

        print("\nSign concordance vs 50%:")
        print(f"  Binomial P                : {b.pvalue:.4g}")


    return {
        "set": label,
        "n_ENVI": n_total,
        "n_present_Visium_L6": n_matched,
        "n_same_direction": n_same,
        "pct_same_direction": (
            100 * n_same / n_matched
            if n_matched else np.nan
        ),
        "n_opposite_direction": n_opposite,

        "Visium_nominal_P05": len(nominal),
        "Visium_nominal_P05_same_direction": same_n(nominal),

        "Visium_FDR10": len(fdr10),
        "Visium_FDR10_same_direction": same_n(fdr10),

        "Visium_FDR05": len(fdr05),
        "Visium_FDR05_same_direction": same_n(fdr05),
    }


# =====================================================================
# Summaries
# =====================================================================

summary_rows = []

summary_rows.append(
    summarize(
        merged,
        "ALL 141 L6 GENES"
    )
)


driver_merged = merged[
    merged["pathway_driver"]
].copy()

summary_rows.append(
    summarize(
        driver_merged,
        "37 PATHWAY-DRIVER GENES"
    )
)


# ---------------------------------------------------------------------
# SCZ-higher and NTC-higher separately
# ---------------------------------------------------------------------

summary_rows.append(
    summarize(
        merged[
            merged["ENVI_direction"]
            == "Higher_in_SCZ"
        ],
        "77 ENVI SCZ-HIGHER GENES"
    )
)

summary_rows.append(
    summarize(
        merged[
            merged["ENVI_direction"]
            == "Higher_in_NTC"
        ],
        "64 ENVI NTC-HIGHER / SCZ-LOWER GENES"
    )
)


summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT / "L6_Ex_ENVI_vs_Visium_validation_summary.csv",
    index=False
)


# =====================================================================
# High-confidence concordant candidates
# =====================================================================

concordant_nominal = merged[
    merged["present_in_Visium_L6"]
    &
    merged["same_direction"]
    &
    merged["Visium_nominal_P05"]
].copy()


concordant_fdr10 = merged[
    merged["present_in_Visium_L6"]
    &
    merged["same_direction"]
    &
    merged["Visium_FDR10_calc"]
].copy()


concordant_fdr05 = merged[
    merged["present_in_Visium_L6"]
    &
    merged["same_direction"]
    &
    merged["Visium_FDR05"]
].copy()


concordant_nominal.to_csv(
    OUT / "L6_Ex_concordant_Visium_nominalP05.csv",
    index=False
)

concordant_fdr10.to_csv(
    OUT / "L6_Ex_concordant_Visium_FDR10.csv",
    index=False
)

concordant_fdr05.to_csv(
    OUT / "L6_Ex_concordant_Visium_FDR05.csv",
    index=False
)


# =====================================================================
# Print strongest concordant genes
# =====================================================================

print("\n" + "=" * 110)
print("CONCORDANT ENVI + VISIUM L6 GENES")
print("=" * 110)

show = [
    "gene",
    "ENVI_direction",
    "log2FC",
    "FDR_within_celltype",
    "Visium_log2FC",
    "Visium_P",
    "Visium_FDR",
    "Visium_layer_specific",
    "pathway_driver",
]


if len(concordant_nominal):

    display = (
        concordant_nominal
        .sort_values(
            [
                "Visium_FDR",
                "FDR_within_celltype",
            ]
        )
    )

    print(
        display[show]
        .to_string(index=False)
    )

else:
    print(
        "No same-direction genes with "
        "Visium nominal P < 0.05."
    )


print("\nSaved to:")
print(OUT)

print("\nSUCCESS")
