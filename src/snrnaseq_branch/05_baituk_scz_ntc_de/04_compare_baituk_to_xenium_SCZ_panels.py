#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

# ======================================================================================
# INPUTS
# ======================================================================================

BAITUK_DE = (
    ROOT
    / "outputs/snrnaseq_branch/05_baituk_scz_ntc_de"
    / "baituk_donor_level_SCZ_vs_NTC_DE_all.csv"
)

# Original measured Xenium: validated paper-compatible WITH-slide analysis
XENIUM_DE = (
    ROOT
    / "outputs/imputation_full/DE/task2_3_preimputation_baseline"
    / "domain_de/with_slide"
    / "original_xenium_domain_SCZ_DE_dx_results.csv"
)

PAPER_FILE = (
    ROOT
    / "data/reference/paper_xenium/paper_concordance_de.xlsx"
)

OUT_DIR = (
    ROOT
    / "outputs/snrnaseq_branch/05_baituk_scz_ntc_de"
    / "targeted_xenium_panel_comparison"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

OUT_300 = (
    OUT_DIR
    / "baituk_vs_xenium_300gene_targeted_DE.csv"
)

OUT_19 = (
    OUT_DIR
    / "baituk_vs_paper19_targeted_DE.csv"
)

OUT_SUMMARY = (
    OUT_DIR
    / "baituk_vs_xenium_targeted_panel_summary.csv"
)


# ======================================================================================
# HELPERS
# ======================================================================================

def section(title):
    print()
    print("=" * 115)
    print(title)
    print("=" * 115)


def norm_gene(x):
    return (
        x.astype(str)
        .str.strip()
        .str.upper()
    )


def bh_adjust(pvalues):
    """
    Benjamini-Hochberg FDR.
    NaN p-values remain NaN.
    """

    p = np.asarray(pvalues, dtype=float)

    out = np.full(
        len(p),
        np.nan,
        dtype=float
    )

    valid = np.isfinite(p)

    if valid.sum() == 0:
        return out

    pv = p[valid]

    order = np.argsort(pv)
    ranked = pv[order]

    n = len(ranked)

    q = ranked * n / np.arange(
        1,
        n + 1
    )

    # Enforce monotonicity from largest rank backward
    q = np.minimum.accumulate(
        q[::-1]
    )[::-1]

    q = np.minimum(
        q,
        1.0
    )

    restored = np.empty(
        n,
        dtype=float
    )

    restored[order] = q

    out[valid] = restored

    return out


def safe_corr(df, x, y, method):
    z = df[
        [x, y]
    ].dropna()

    if len(z) < 3:
        return np.nan

    return z[x].corr(
        z[y],
        method=method
    )


# ======================================================================================
# LOAD BAITUK DONOR-LEVEL DE
# ======================================================================================

section("LOAD BAITUK DONOR-LEVEL SCZ vs NTC DE")

baituk = pd.read_csv(
    BAITUK_DE
)

required = {
    "gene",
    "logFC",
    "PValue",
    "FDR"
}

missing = required - set(
    baituk.columns
)

if missing:
    raise RuntimeError(
        f"Baituk DE missing columns: {sorted(missing)}"
    )


baituk["gene"] = (
    baituk["gene"]
    .astype(str)
    .str.strip()
)

baituk["gene_std"] = norm_gene(
    baituk["gene"]
)

baituk["logFC"] = pd.to_numeric(
    baituk["logFC"],
    errors="coerce"
)

baituk["PValue"] = pd.to_numeric(
    baituk["PValue"],
    errors="coerce"
)

baituk["FDR"] = pd.to_numeric(
    baituk["FDR"],
    errors="coerce"
)


if baituk["gene_std"].duplicated().any():
    dup = baituk.loc[
        baituk["gene_std"].duplicated(False),
        "gene"
    ].tolist()

    raise RuntimeError(
        f"Duplicate Baituk genes after normalization: {dup[:20]}"
    )


print("Baituk genes tested:", len(baituk))
print(
    "Genome-wide FDR < 0.10:",
    int((baituk["FDR"] < 0.10).sum())
)


# ======================================================================================
# LOAD ORIGINAL MEASURED XENIUM DE
# ======================================================================================

section("LOAD ORIGINAL MEASURED XENIUM 300-GENE DE")

xenium = pd.read_csv(
    XENIUM_DE
)

required_x = {
    "gene",
    "logFC_SCZ",
    "p_value_SCZ",
    "fdr_SCZ"
}

missing = required_x - set(
    xenium.columns
)

if missing:
    raise RuntimeError(
        f"Xenium DE missing columns: {sorted(missing)}"
    )


xenium["gene"] = (
    xenium["gene"]
    .astype(str)
    .str.strip()
)

xenium["gene_std"] = norm_gene(
    xenium["gene"]
)

xenium["logFC_SCZ"] = pd.to_numeric(
    xenium["logFC_SCZ"],
    errors="coerce"
)

xenium["p_value_SCZ"] = pd.to_numeric(
    xenium["p_value_SCZ"],
    errors="coerce"
)

xenium["fdr_SCZ"] = pd.to_numeric(
    xenium["fdr_SCZ"],
    errors="coerce"
)


if len(xenium) != 300:
    raise RuntimeError(
        f"Expected 300 Xenium genes; found {len(xenium)}"
    )


if xenium["gene_std"].duplicated().any():
    raise RuntimeError(
        "Duplicate Xenium genes after normalization."
    )


print("Xenium genes:", len(xenium))
print(
    "Xenium FDR < 0.10:",
    int((xenium["fdr_SCZ"] < 0.10).sum())
)


# ======================================================================================
# LOAD PAPER 19-GENE REFERENCE
# ======================================================================================

section("LOAD PAPER DONOR-DOMAIN 19-GENE REFERENCE")

paper = pd.read_excel(
    PAPER_FILE,
    sheet_name="donor_domain_level",
    engine="openpyxl"
)

for col in [
    "gene",
    "logFC_SCZ",
    "fdr_SCZ"
]:
    if col not in paper.columns:
        raise RuntimeError(
            f"Paper sheet missing column: {col}"
        )


paper["gene"] = (
    paper["gene"]
    .astype(str)
    .str.strip()
)

paper["gene_std"] = norm_gene(
    paper["gene"]
)

paper["logFC_SCZ"] = pd.to_numeric(
    paper["logFC_SCZ"],
    errors="coerce"
)

paper["fdr_SCZ"] = pd.to_numeric(
    paper["fdr_SCZ"],
    errors="coerce"
)


paper19 = paper.loc[
    paper["fdr_SCZ"] < 0.10
].copy()


if len(paper19) != 19:
    raise RuntimeError(
        f"Expected exactly 19 paper genes; "
        f"found {len(paper19)}"
    )


print("Paper FDR < 0.10 genes:", len(paper19))


# ======================================================================================
# PANEL A — XENIUM 300
# ======================================================================================

section("TARGETED PANEL A — XENIUM 300 GENES")

panel300 = xenium[
    [
        "gene",
        "gene_std",
        "logFC_SCZ",
        "p_value_SCZ",
        "fdr_SCZ"
    ]
].copy()

panel300 = panel300.rename(
    columns={
        "gene": "xenium_gene",
        "logFC_SCZ": "xenium_logFC_SCZ",
        "p_value_SCZ": "xenium_PValue",
        "fdr_SCZ": "xenium_FDR"
    }
)


baituk_small = baituk[
    [
        "gene_std",
        "gene",
        "logFC",
        "PValue",
        "FDR"
    ]
].copy()

baituk_small = baituk_small.rename(
    columns={
        "gene": "baituk_gene",
        "logFC": "baituk_logFC_SCZ",
        "PValue": "baituk_PValue",
        "FDR": "baituk_genomewide_FDR"
    }
)


panel300 = panel300.merge(
    baituk_small,
    on="gene_std",
    how="left",
    validate="one_to_one"
)


panel300["baituk_tested"] = (
    panel300["baituk_PValue"]
    .notna()
)


# Panel-specific correction among the externally fixed
# Xenium genes that survived Baituk expression filtering.
panel300["baituk_panel300_FDR"] = bh_adjust(
    panel300["baituk_PValue"]
)


panel300["same_direction"] = np.where(
    panel300["baituk_logFC_SCZ"].notna(),
    np.sign(
        panel300["baituk_logFC_SCZ"]
    )
    ==
    np.sign(
        panel300["xenium_logFC_SCZ"]
    ),
    np.nan
)


panel300["baituk_nominal_P005"] = (
    panel300["baituk_PValue"] < 0.05
)

panel300["baituk_panel_FDR010"] = (
    panel300["baituk_panel300_FDR"] < 0.10
)


panel300.to_csv(
    OUT_300,
    index=False
)


tested300 = panel300.loc[
    panel300["baituk_tested"]
].copy()

n_tested300 = len(
    tested300
)

nominal300 = int(
    (
        tested300["baituk_PValue"] < 0.05
    ).sum()
)

fdr300 = int(
    (
        tested300["baituk_panel300_FDR"] < 0.10
    ).sum()
)

same300 = int(
    (
        tested300["same_direction"] == True
    ).sum()
)

spearman300 = safe_corr(
    tested300,
    "baituk_logFC_SCZ",
    "xenium_logFC_SCZ",
    "spearman"
)

pearson300 = safe_corr(
    tested300,
    "baituk_logFC_SCZ",
    "xenium_logFC_SCZ",
    "pearson"
)


print("Xenium panel size              :", 300)
print("Tested in Baituk              :", n_tested300)
print("Not tested after Baituk filter:", 300 - n_tested300)
print("Baituk nominal P < 0.05       :", nominal300)
print("Baituk panel FDR < 0.10       :", fdr300)

print(
    "Direction concordance         :",
    f"{same300}/{n_tested300}"
)

print(
    "Direction concordance %       :",
    f"{100 * same300 / n_tested300:.2f}%"
    if n_tested300 else "NA"
)

print(
    "Spearman logFC correlation    :",
    f"{spearman300:.4f}"
    if pd.notna(spearman300)
    else "NA"
)

print(
    "Pearson logFC correlation     :",
    f"{pearson300:.4f}"
    if pd.notna(pearson300)
    else "NA"
)


# ======================================================================================
# PANEL B — EXACT PAPER 19
# ======================================================================================

section("TARGETED PANEL B — PAPER DOMAIN 19 GENES")

panel19 = paper19[
    [
        "gene",
        "gene_std",
        "logFC_SCZ",
        "fdr_SCZ"
    ]
].copy()

panel19 = panel19.rename(
    columns={
        "gene": "paper_gene",
        "logFC_SCZ": "paper_logFC_SCZ",
        "fdr_SCZ": "paper_FDR"
    }
)


# Add original measured Xenium values
x19 = xenium[
    [
        "gene_std",
        "logFC_SCZ",
        "p_value_SCZ",
        "fdr_SCZ"
    ]
].copy()

x19 = x19.rename(
    columns={
        "logFC_SCZ": "xenium_logFC_SCZ",
        "p_value_SCZ": "xenium_PValue",
        "fdr_SCZ": "xenium_FDR"
    }
)


panel19 = panel19.merge(
    x19,
    on="gene_std",
    how="left",
    validate="one_to_one"
)

panel19 = panel19.merge(
    baituk_small,
    on="gene_std",
    how="left",
    validate="one_to_one"
)


panel19["baituk_tested"] = (
    panel19["baituk_PValue"]
    .notna()
)


panel19["baituk_panel19_FDR"] = bh_adjust(
    panel19["baituk_PValue"]
)


panel19["same_direction_vs_paper"] = np.where(
    panel19["baituk_logFC_SCZ"].notna(),
    np.sign(
        panel19["baituk_logFC_SCZ"]
    )
    ==
    np.sign(
        panel19["paper_logFC_SCZ"]
    ),
    np.nan
)


panel19["same_direction_vs_xenium"] = np.where(
    panel19["baituk_logFC_SCZ"].notna(),
    np.sign(
        panel19["baituk_logFC_SCZ"]
    )
    ==
    np.sign(
        panel19["xenium_logFC_SCZ"]
    ),
    np.nan
)


panel19["baituk_nominal_P005"] = (
    panel19["baituk_PValue"] < 0.05
)

panel19["baituk_panel_FDR010"] = (
    panel19["baituk_panel19_FDR"] < 0.10
)


panel19 = panel19.sort_values(
    [
        "baituk_PValue",
        "paper_FDR"
    ],
    na_position="last"
)


panel19.to_csv(
    OUT_19,
    index=False
)


tested19 = panel19.loc[
    panel19["baituk_tested"]
].copy()

n_tested19 = len(
    tested19
)

nominal19 = int(
    (
        tested19["baituk_PValue"] < 0.05
    ).sum()
)

fdr19 = int(
    (
        tested19["baituk_panel19_FDR"] < 0.10
    ).sum()
)

same19 = int(
    (
        tested19["same_direction_vs_paper"] == True
    ).sum()
)

spearman19 = safe_corr(
    tested19,
    "baituk_logFC_SCZ",
    "paper_logFC_SCZ",
    "spearman"
)

pearson19 = safe_corr(
    tested19,
    "baituk_logFC_SCZ",
    "paper_logFC_SCZ",
    "pearson"
)


print("Paper panel size               :", 19)
print("Tested in Baituk              :", n_tested19)
print("Not tested after Baituk filter:", 19 - n_tested19)
print("Baituk nominal P < 0.05       :", nominal19)
print("Baituk panel FDR < 0.10       :", fdr19)

print(
    "Direction concordance         :",
    f"{same19}/{n_tested19}"
)

print(
    "Direction concordance %       :",
    f"{100 * same19 / n_tested19:.2f}%"
    if n_tested19 else "NA"
)

print(
    "Spearman logFC correlation    :",
    f"{spearman19:.4f}"
    if pd.notna(spearman19)
    else "NA"
)

print(
    "Pearson logFC correlation     :",
    f"{pearson19:.4f}"
    if pd.notna(pearson19)
    else "NA"
)


# ======================================================================================
# PRINT ALL 19 GENES
# ======================================================================================

section("PAPER-19 GENE-BY-GENE RESULT")

print_cols = [
    "paper_gene",
    "paper_logFC_SCZ",
    "paper_FDR",
    "baituk_logFC_SCZ",
    "baituk_PValue",
    "baituk_genomewide_FDR",
    "baituk_panel19_FDR",
    "same_direction_vs_paper"
]

print(
    panel19[
        print_cols
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.6g}"
    )
)


# ======================================================================================
# SUMMARY
# ======================================================================================

summary = pd.DataFrame(
    [
        {
            "panel": "Xenium_300",
            "panel_size": 300,
            "baituk_tested": n_tested300,
            "baituk_not_tested": 300 - n_tested300,
            "baituk_nominal_P_lt_005": nominal300,
            "baituk_panel_FDR_lt_010": fdr300,
            "same_direction": same300,
            "direction_concordance_percent":
                100 * same300 / n_tested300
                if n_tested300 else np.nan,
            "logFC_spearman": spearman300,
            "logFC_pearson": pearson300,
        },
        {
            "panel": "Paper_domain_19",
            "panel_size": 19,
            "baituk_tested": n_tested19,
            "baituk_not_tested": 19 - n_tested19,
            "baituk_nominal_P_lt_005": nominal19,
            "baituk_panel_FDR_lt_010": fdr19,
            "same_direction": same19,
            "direction_concordance_percent":
                100 * same19 / n_tested19
                if n_tested19 else np.nan,
            "logFC_spearman": spearman19,
            "logFC_pearson": pearson19,
        }
    ]
)


summary.to_csv(
    OUT_SUMMARY,
    index=False
)


section("FINAL SUMMARY")

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}"
    )
)


print("\nSaved Xenium-300 table:")
print(OUT_300)

print("\nSaved paper-19 table:")
print(OUT_19)

print("\nSaved summary:")
print(OUT_SUMMARY)

print(
    "\nFINAL STATUS: PASS — BAITUK TARGETED XENIUM SCZ PANEL COMPARISON COMPLETE"
)
