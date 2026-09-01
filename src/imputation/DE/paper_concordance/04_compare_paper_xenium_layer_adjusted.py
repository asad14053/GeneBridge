#!/usr/bin/env python3

"""
Compare GeneBridge 23-donor Xenium layer-adjusted SCZ DE
against the published 24-donor Xenium layer-adjusted DE results.

Inputs
------
Paper:
    data/reference/paper_xenium/
        paper_xenium_layer_adjusted_24donor.csv

GeneBridge:
    outputs/imputation_full/DE/paper_concordance/
        layer_adjusted_de/original/
        original_layer_adjusted_SCZ_DE_all.csv

Outputs
-------
    outputs/imputation_full/DE/paper_concordance/
        layer_adjusted_de/paper_comparison/

        paper24_vs_genebridge23_layer_adjusted_all300.csv
        paper24_vs_genebridge23_reproduced_FDR10.csv
        paper24_only_FDR10.csv
        genebridge23_only_FDR10.csv
        paper24_vs_genebridge23_summary.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


# =============================================================================
# Paths
# =============================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

PAPER_FILE = (
    ROOT
    / "data/reference/paper_xenium/"
      "paper_xenium_layer_adjusted_24donor.csv"
)

OURS_FILE = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "layer_adjusted_de/original/"
      "original_layer_adjusted_SCZ_DE_all.csv"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "layer_adjusted_de/paper_comparison"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# =============================================================================
# Helpers
# =============================================================================

def section(title):
    print(
        "\n"
        + "=" * 100
        + f"\n{title}\n"
        + "=" * 100
    )


def safe_pearson(x, y):
    if len(x) < 3:
        return np.nan

    return float(
        pearsonr(x, y).statistic
    )


def safe_spearman(x, y):
    if len(x) < 3:
        return np.nan

    return float(
        spearmanr(x, y).statistic
    )


# =============================================================================
# Load
# =============================================================================

section("PAPER vs GENEBRIDGE XENIUM LAYER-ADJUSTED SCZ DE")

print("Paper file:")
print(PAPER_FILE)

print("\nGeneBridge file:")
print(OURS_FILE)


if not PAPER_FILE.exists():
    raise FileNotFoundError(
        f"Paper result file not found:\n{PAPER_FILE}"
    )

if not OURS_FILE.exists():
    raise FileNotFoundError(
        f"GeneBridge result file not found:\n{OURS_FILE}"
    )


paper = pd.read_csv(
    PAPER_FILE
)

ours = pd.read_csv(
    OURS_FILE
)


print("\nRaw paper shape:", paper.shape)
print("Raw GeneBridge shape:", ours.shape)


# =============================================================================
# Validate required columns
# =============================================================================

required = [
    "gene",
    "logFC_SCZ",
    "t_stat_SCZ",
    "p_value_SCZ",
    "fdr_SCZ",
]


for name, df in [
    ("paper", paper),
    ("GeneBridge", ours),
]:

    missing = [
        x
        for x in required
        if x not in df.columns
    ]

    if missing:
        raise RuntimeError(
            f"{name}: missing required columns: {missing}"
        )


# =============================================================================
# Reduce to relevant columns
# =============================================================================

paper = (
    paper[
        required
    ]
    .copy()
    .rename(
        columns={
            "logFC_SCZ": "paper_logFC",
            "t_stat_SCZ": "paper_t",
            "p_value_SCZ": "paper_p",
            "fdr_SCZ": "paper_FDR",
        }
    )
)


ours = (
    ours[
        required
    ]
    .copy()
    .rename(
        columns={
            "logFC_SCZ": "our_logFC",
            "t_stat_SCZ": "our_t",
            "p_value_SCZ": "our_p",
            "fdr_SCZ": "our_FDR",
        }
    )
)


paper["gene"] = (
    paper["gene"]
    .astype(str)
    .str.strip()
)

ours["gene"] = (
    ours["gene"]
    .astype(str)
    .str.strip()
)


if paper["gene"].duplicated().any():
    dup = (
        paper.loc[
            paper["gene"].duplicated(),
            "gene"
        ]
        .tolist()
    )

    raise RuntimeError(
        f"Duplicate genes in paper output: {dup[:20]}"
    )


if ours["gene"].duplicated().any():
    dup = (
        ours.loc[
            ours["gene"].duplicated(),
            "gene"
        ]
        .tolist()
    )

    raise RuntimeError(
        f"Duplicate genes in GeneBridge output: {dup[:20]}"
    )


# =============================================================================
# Merge
# =============================================================================

x = paper.merge(
    ours,
    on="gene",
    how="inner",
    validate="one_to_one",
)


section("GENE OVERLAP")

print("Paper genes       :", len(paper))
print("GeneBridge genes  :", len(ours))
print("Shared genes      :", len(x))


paper_only_genes = sorted(
    set(paper["gene"])
    - set(ours["gene"])
)

ours_only_genes = sorted(
    set(ours["gene"])
    - set(paper["gene"])
)


print(
    "Paper-only genes  :",
    len(paper_only_genes),
)

print(
    "GeneBridge-only   :",
    len(ours_only_genes),
)


if paper_only_genes:
    print(
        "\nPaper-only gene examples:",
        paper_only_genes[:20],
    )


if ours_only_genes:
    print(
        "\nGeneBridge-only gene examples:",
        ours_only_genes[:20],
    )


if len(x) == 0:
    raise RuntimeError(
        "No shared genes between paper and GeneBridge."
    )


# =============================================================================
# Significance and direction
# =============================================================================

x["paper_sig10"] = (
    x["paper_FDR"] < 0.10
)

x["our_sig10"] = (
    x["our_FDR"] < 0.10
)


x["paper_sig05"] = (
    x["paper_FDR"] < 0.05
)

x["our_sig05"] = (
    x["our_FDR"] < 0.05
)


x["same_direction"] = (
    np.sign(x["paper_logFC"])
    == np.sign(x["our_logFC"])
)


x["abs_logFC_difference"] = np.abs(
    x["paper_logFC"]
    - x["our_logFC"]
)


x["status"] = np.select(
    [
        x["paper_sig10"]
        & x["our_sig10"],

        x["paper_sig10"]
        & ~x["our_sig10"],

        ~x["paper_sig10"]
        & x["our_sig10"],
    ],
    [
        "reproduced",
        "paper_only",
        "genebridge_only",
    ],
    default="not_sig",
)


paper_sig = x.loc[
    x["paper_sig10"]
].copy()


ours_sig = x.loc[
    x["our_sig10"]
].copy()


overlap = x.loc[
    x["paper_sig10"]
    & x["our_sig10"]
].copy()


paper_only = x.loc[
    x["status"] == "paper_only"
].copy()


ours_only = x.loc[
    x["status"] == "genebridge_only"
].copy()


# =============================================================================
# Correlations
# =============================================================================

section("EFFECT-SIZE CONCORDANCE")


pearson_all = safe_pearson(
    x["paper_logFC"],
    x["our_logFC"],
)


spearman_all = safe_spearman(
    x["paper_logFC"],
    x["our_logFC"],
)


pearson_t_all = safe_pearson(
    x["paper_t"],
    x["our_t"],
)


spearman_t_all = safe_spearman(
    x["paper_t"],
    x["our_t"],
)


direction_all = float(
    x["same_direction"].mean()
)


print(
    "Pearson logFC, all shared genes :",
    f"{pearson_all:.4f}",
)

print(
    "Spearman logFC, all shared genes:",
    f"{spearman_all:.4f}",
)

print(
    "Pearson t-stat, all shared genes :",
    f"{pearson_t_all:.4f}",
)

print(
    "Spearman t-stat, all shared genes:",
    f"{spearman_t_all:.4f}",
)

print(
    "Direction concordance, all genes :",
    f"{direction_all:.2%}",
)


# =============================================================================
# FDR concordance
# =============================================================================

section("FDR < 0.10 CONCORDANCE")


n_paper_sig = len(
    paper_sig
)

n_ours_sig = len(
    ours_sig
)

n_overlap = len(
    overlap
)


paper_recovery = (
    n_overlap / n_paper_sig
    if n_paper_sig > 0
    else np.nan
)


ours_confirmation = (
    n_overlap / n_ours_sig
    if n_ours_sig > 0
    else np.nan
)


print(
    "Paper FDR < 0.10:",
    n_paper_sig,
)

print(
    "GeneBridge FDR < 0.10:",
    n_ours_sig,
)

print(
    "Overlap:",
    n_overlap,
)

print(
    "Paper DEG recovery:",
    (
        f"{n_overlap}/{n_paper_sig} "
        f"({paper_recovery:.1%})"
        if n_paper_sig > 0
        else "NA"
    ),
)

print(
    "GeneBridge DEG confirmation:",
    (
        f"{n_overlap}/{n_ours_sig} "
        f"({ours_confirmation:.1%})"
        if n_ours_sig > 0
        else "NA"
    ),
)


if n_overlap > 0:

    overlap_direction = float(
        overlap[
            "same_direction"
        ].mean()
    )

else:

    overlap_direction = np.nan


print(
    "Direction concordance among reproduced DEGs:",
    (
        f"{overlap_direction:.1%}"
        if np.isfinite(overlap_direction)
        else "NA"
    ),
)


pearson_overlap = safe_pearson(
    overlap["paper_logFC"],
    overlap["our_logFC"],
)


spearman_overlap = safe_spearman(
    overlap["paper_logFC"],
    overlap["our_logFC"],
)


print(
    "Pearson logFC among reproduced DEGs:",
    (
        f"{pearson_overlap:.4f}"
        if np.isfinite(pearson_overlap)
        else "NA"
    ),
)

print(
    "Spearman logFC among reproduced DEGs:",
    (
        f"{spearman_overlap:.4f}"
        if np.isfinite(spearman_overlap)
        else "NA"
    ),
)


# =============================================================================
# DEG tables
# =============================================================================

show_cols = [
    "gene",
    "paper_logFC",
    "our_logFC",
    "paper_FDR",
    "our_FDR",
    "same_direction",
]


section("REPRODUCED PAPER DEGs")

if overlap.empty:

    print("NONE")

else:

    print(
        overlap[
            show_cols
        ]
        .sort_values(
            "paper_FDR"
        )
        .to_string(
            index=False
        )
    )


section("PAPER FDR < 0.10 BUT GENEBRIDGE FDR >= 0.10")

if paper_only.empty:

    print("NONE")

else:

    print(
        paper_only[
            show_cols
        ]
        .sort_values(
            "paper_FDR"
        )
        .to_string(
            index=False
        )
    )


section("GENEBRIDGE FDR < 0.10 BUT PAPER FDR >= 0.10")

if ours_only.empty:

    print("NONE")

else:

    print(
        ours_only[
            show_cols
        ]
        .sort_values(
            "our_FDR"
        )
        .to_string(
            index=False
        )
    )


# =============================================================================
# Summary
# =============================================================================

summary = pd.DataFrame(
    [
        {
            "paper_donors": 24,
            "genebridge_donors": 23,

            "paper_genes":
                len(paper),

            "genebridge_genes":
                len(ours),

            "shared_genes":
                len(x),

            "pearson_logFC_all":
                pearson_all,

            "spearman_logFC_all":
                spearman_all,

            "pearson_t_all":
                pearson_t_all,

            "spearman_t_all":
                spearman_t_all,

            "direction_concordance_all":
                direction_all,

            "paper_fdr10":
                n_paper_sig,

            "genebridge_fdr10":
                n_ours_sig,

            "fdr10_overlap":
                n_overlap,

            "paper_deg_recovery":
                paper_recovery,

            "genebridge_deg_confirmation":
                ours_confirmation,

            "direction_concordance_reproduced":
                overlap_direction,

            "pearson_logFC_reproduced":
                pearson_overlap,

            "spearman_logFC_reproduced":
                spearman_overlap,
        }
    ]
)


# =============================================================================
# Save
# =============================================================================

section("SAVE OUTPUTS")


all_file = (
    OUT
    / "paper24_vs_genebridge23_layer_adjusted_all300.csv"
)

overlap_file = (
    OUT
    / "paper24_vs_genebridge23_reproduced_FDR10.csv"
)

paper_only_file = (
    OUT
    / "paper24_only_FDR10.csv"
)

ours_only_file = (
    OUT
    / "genebridge23_only_FDR10.csv"
)

summary_file = (
    OUT
    / "paper24_vs_genebridge23_summary.csv"
)


x.to_csv(
    all_file,
    index=False,
)

overlap.to_csv(
    overlap_file,
    index=False,
)

paper_only.to_csv(
    paper_only_file,
    index=False,
)

ours_only.to_csv(
    ours_only_file,
    index=False,
)

summary.to_csv(
    summary_file,
    index=False,
)


print(all_file)
print(overlap_file)
print(paper_only_file)
print(ours_only_file)
print(summary_file)


# =============================================================================
# Final
# =============================================================================

section("FINAL SUMMARY")

print(
    summary.to_string(
        index=False
    )
)


print(
    "\nFINAL STATUS: PASS — "
    "published 24-donor Xenium layer-adjusted DE "
    "successfully compared against GeneBridge 23-donor results."
)
