#!/usr/bin/env python3

"""
TASK 2.3 — PAPER CONCORDANCE
============================

Canonical paper reference:
    data/reference/paper_xenium/paper_concordance_de.xlsx

Sheets:
    donor_domain_level
        -> expected 19 genes with fdr_SCZ < 0.10
        -> compare against donor/domain-adjusted DE

    donor_cell_type_level
        -> expected 15 genes with fdr_SCZ < 0.10
        -> compare against donor/cell-type-adjusted DE

For each of the 6 imputation models:
    - genes tested
    - model genes FDR < 0.10
    - model UP / DOWN
    - paper significant genes
    - paper UP / DOWN
    - paper genes actually tested by model
    - significant overlap
    - same-direction overlap
    - recovery %
    - matched gene names
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd


# ======================================================================================
# CONFIG
# ======================================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

TASK_ROOT = (
    ROOT
    / "outputs/imputation_full/DE/task2_3_ntc_dysregulation"
)

WORKBOOK = (
    ROOT
    / "data/reference/paper_xenium/paper_concordance_de.xlsx"
)

CELLTYPE_WITH_SLIDE_ROOT = (
    TASK_ROOT
    / "celltype_adjusted_de_nofilter"
)

CELLTYPE_NOSLIDE_ROOT = (
    TASK_ROOT
    / "celltype_adjusted_de_nofilter_noslide"
)

MODELS = [
    "VISTA",
    "gimVI",
    "ENVI",
    "TransImp",
    "Tangram",
    "SpaGE",
]

FDR = 0.10


# ======================================================================================
# OUTPUTS
# ======================================================================================

REF_OUT = (
    ROOT
    / "data/reference/paper_xenium"
)

DOMAIN_REF_ALL = (
    REF_OUT
    / "paper_concordance_donor_domain_level.csv"
)

DOMAIN_REF_SIG = (
    REF_OUT
    / "paper_concordance_donor_domain_FDR010_19genes.csv"
)

CELL_REF_ALL = (
    REF_OUT
    / "paper_concordance_donor_cell_type_level.csv"
)

CELL_REF_SIG = (
    REF_OUT
    / "paper_concordance_donor_cell_type_FDR010_15genes.csv"
)

DOMAIN_SLIDE_SUMMARY = (
    TASK_ROOT
    / "paper_concordance_domain_DE_with_slide_summary.csv"
)

DOMAIN_SLIDE_MATCHED = (
    TASK_ROOT
    / "paper_concordance_domain_DE_with_slide_matched_genes.csv"
)

DOMAIN_NOSLIDE_SUMMARY = (
    TASK_ROOT
    / "paper_concordance_domain_DE_noslide_summary.csv"
)

DOMAIN_NOSLIDE_MATCHED = (
    TASK_ROOT
    / "paper_concordance_domain_DE_noslide_matched_genes.csv"
)

CELL_SLIDE_SUMMARY = (
    CELLTYPE_WITH_SLIDE_ROOT
    / "paper_concordance_celltype_DE_with_slide_summary.csv"
)

CELL_SLIDE_MATCHED = (
    CELLTYPE_WITH_SLIDE_ROOT
    / "paper_concordance_celltype_DE_with_slide_matched_genes.csv"
)

CELL_NOSLIDE_SUMMARY = (
    CELLTYPE_NOSLIDE_ROOT
    / "paper_concordance_celltype_DE_noslide_summary.csv"
)

CELL_NOSLIDE_MATCHED = (
    CELLTYPE_NOSLIDE_ROOT
    / "paper_concordance_celltype_DE_noslide_matched_genes.csv"
)


# ======================================================================================
# HELPERS
# ======================================================================================

def section(title):
    print()
    print("=" * 120)
    print(title)
    print("=" * 120)


def norm(s):
    return re.sub(
        r"[^a-z0-9]",
        "",
        str(s).lower(),
    )


def find_col(df, candidates, label):

    normalized = {
        norm(c): c
        for c in df.columns
    }

    for candidate in candidates:

        key = norm(candidate)

        if key in normalized:
            return normalized[key]

    raise RuntimeError(
        f"{label}: could not identify column from "
        f"{candidates}\n"
        f"Available columns: {list(df.columns)}"
    )


def standardize(df, label):

    gene_col = find_col(
        df,
        [
            "gene",
            "Symbol",
            "gene_name",
            "Gene",
        ],
        label,
    )

    fdr_col = find_col(
        df,
        [
            "fdr_SCZ",
            "FDR_SCZ",
            "fdr",
            "FDR",
            "adj.P.Val",
        ],
        label,
    )

    logfc_col = find_col(
        df,
        [
            "logFC_SCZ",
            "logFC",
            "log2FC",
            "log2FoldChange",
        ],
        label,
    )

    x = df.copy()

    x["gene_std"] = (
        x[gene_col]
        .astype(str)
        .str.strip()
    )

    x["fdr_std"] = pd.to_numeric(
        x[fdr_col],
        errors="coerce",
    )

    x["logFC_std"] = pd.to_numeric(
        x[logfc_col],
        errors="coerce",
    )

    if x["gene_std"].duplicated().any():

        dup = (
            x.loc[
                x["gene_std"].duplicated(keep=False),
                "gene_std",
            ]
            .unique()
            .tolist()
        )

        raise RuntimeError(
            f"{label}: duplicated genes detected: {dup}"
        )

    return x


# ======================================================================================
# LOAD PAPER WORKBOOK DIRECTLY
# ======================================================================================

section("READ CANONICAL PAPER WORKBOOK")

if not WORKBOOK.exists():
    raise FileNotFoundError(WORKBOOK)

excel = pd.ExcelFile(
    WORKBOOK,
    engine="openpyxl",
)

print("Workbook:")
print(WORKBOOK)

print("\nSheets:")
for s in excel.sheet_names:
    print("  ", s)


required_sheets = {
    "donor_domain_level",
    "donor_cell_type_level",
}

missing_sheets = (
    required_sheets
    - set(excel.sheet_names)
)

if missing_sheets:
    raise RuntimeError(
        f"Missing required sheets: {missing_sheets}"
    )


domain_raw = pd.read_excel(
    WORKBOOK,
    sheet_name="donor_domain_level",
    engine="openpyxl",
)

cell_raw = pd.read_excel(
    WORKBOOK,
    sheet_name="donor_cell_type_level",
    engine="openpyxl",
)


domain = standardize(
    domain_raw,
    "donor_domain_level",
)

cell = standardize(
    cell_raw,
    "donor_cell_type_level",
)


# ======================================================================================
# EXTRACT PAPER SIGNIFICANT SETS
# ======================================================================================

domain_sig = (
    domain.loc[
        domain["fdr_std"] < FDR
    ]
    .copy()
    .sort_values(
        ["fdr_std", "gene_std"]
    )
)

cell_sig = (
    cell.loc[
        cell["fdr_std"] < FDR
    ]
    .copy()
    .sort_values(
        ["fdr_std", "gene_std"]
    )
)


section("VALIDATE PAPER REFERENCE SETS")

print(
    "donor_domain_level FDR < 0.10:",
    len(domain_sig),
)

print(
    "donor_cell_type_level FDR < 0.10:",
    len(cell_sig),
)


if len(domain_sig) != 19:
    raise RuntimeError(
        "STOP: donor_domain_level should contain exactly "
        f"19 genes with FDR < 0.10; found {len(domain_sig)}."
    )


if len(cell_sig) != 15:
    raise RuntimeError(
        "STOP: donor_cell_type_level should contain exactly "
        f"15 genes with FDR < 0.10; found {len(cell_sig)}."
    )


def paper_direction_summary(sig):

    up = int(
        (sig["logFC_std"] > 0).sum()
    )

    down = int(
        (sig["logFC_std"] < 0).sum()
    )

    return up, down


domain_up, domain_down = paper_direction_summary(
    domain_sig
)

cell_up, cell_down = paper_direction_summary(
    cell_sig
)


print()
print(
    f"DOMAIN PAPER:    19 total | "
    f"UP={domain_up} | DOWN={domain_down}"
)

print(
    f"CELL-TYPE PAPER: 15 total | "
    f"UP={cell_up} | DOWN={cell_down}"
)


print("\nDOMAIN 19 genes:")

print(
    domain_sig[
        [
            "gene_std",
            "logFC_std",
            "fdr_std",
        ]
    ].to_string(
        index=False
    )
)


print("\nCELL-TYPE 15 genes:")

print(
    cell_sig[
        [
            "gene_std",
            "logFC_std",
            "fdr_std",
        ]
    ].to_string(
        index=False
    )
)


# ======================================================================================
# SAVE CANONICAL EXTRACTS
# ======================================================================================

domain_raw.to_csv(
    DOMAIN_REF_ALL,
    index=False,
)

cell_raw.to_csv(
    CELL_REF_ALL,
    index=False,
)


domain_sig[
    [
        "gene_std",
        "logFC_std",
        "fdr_std",
    ]
].rename(
    columns={
        "gene_std": "gene",
        "logFC_std": "logFC_SCZ",
        "fdr_std": "fdr_SCZ",
    }
).to_csv(
    DOMAIN_REF_SIG,
    index=False,
)


cell_sig[
    [
        "gene_std",
        "logFC_std",
        "fdr_std",
    ]
].rename(
    columns={
        "gene_std": "gene",
        "logFC_std": "logFC_SCZ",
        "fdr_std": "fdr_SCZ",
    }
).to_csv(
    CELL_REF_SIG,
    index=False,
)


# ======================================================================================
# FIND DOMAIN DE RESULT
#
# Existing Task 2.3 domain workflow has both filtered/no-filter history.
# Our primary fixed-300-gene comparison used the NO-FILTER branch.
# Prefer a path containing "nofilter" / "no_filter" when multiple files exist.
# ======================================================================================


def find_domain_de_with_slide(model):

    """
    Domain DE with slide_id adjustment.
    """

    f = (
        TASK_ROOT
        / "layer_adjusted_de_nofilter"
        / model
        / f"{model}_layer_adjusted_SCZ_DE_dx_results.csv"
    )

    if not f.exists():
        raise FileNotFoundError(f)

    return f


def find_domain_de_noslide(model):

    """
    Domain DE without slide_id adjustment.
    """

    f = (
        TASK_ROOT
        / "layer_adjusted_de_nofilter_noslide"
        / model
        / f"{model}_layer_adjusted_SCZ_DE_dx_results.csv"
    )

    if not f.exists():
        raise FileNotFoundError(f)

    return f


# ======================================================================================
# EXACT CELL-TYPE DE PATH
# ======================================================================================

def find_cell_de_with_slide(model):

    f = (
        CELLTYPE_WITH_SLIDE_ROOT
        / model
        / f"{model}_celltype_adjusted_SCZ_DE_dx_results.csv"
    )

    if not f.exists():
        raise FileNotFoundError(f)

    return f


def find_cell_de_noslide(model):

    f = (
        CELLTYPE_NOSLIDE_ROOT
        / model
        / f"{model}_celltype_adjusted_SCZ_DE_dx_results.csv"
    )

    if not f.exists():
        raise FileNotFoundError(f)

    return f


# ======================================================================================
# COMPARISON FUNCTION
# ======================================================================================

def compare_level(
    level_name,
    paper_sig,
    expected_n,
    de_file_func,
):

    paper_genes = set(
        paper_sig["gene_std"]
    )

    paper_lookup = (
        paper_sig[
            [
                "gene_std",
                "logFC_std",
                "fdr_std",
            ]
        ]
        .rename(
            columns={
                "gene_std": "gene",
                "logFC_std": "paper_logFC_SCZ",
                "fdr_std": "paper_fdr_SCZ",
            }
        )
    )


    paper_up = int(
        (paper_sig["logFC_std"] > 0).sum()
    )

    paper_down = int(
        (paper_sig["logFC_std"] < 0).sum()
    )


    summaries = []
    matches = []


    for model in MODELS:

        section(
            f"{level_name.upper()} — {model}"
        )

        de_file = de_file_func(
            model
        )

        print("DE file:")
        print(de_file)


        de_raw = pd.read_csv(
            de_file
        )

        de = standardize(
            de_raw,
            f"{level_name}/{model}",
        )


        sig = de.loc[
            de["fdr_std"] < FDR
        ].copy()


        n_up = int(
            (sig["logFC_std"] > 0).sum()
        )

        n_down = int(
            (sig["logFC_std"] < 0).sum()
        )


        tested = set(
            de["gene_std"]
        )

        paper_tested = (
            paper_genes
            & tested
        )

        paper_not_tested = (
            paper_genes
            - tested
        )


        matched_gene_set = (
            set(sig["gene_std"])
            & paper_genes
        )


        matched = (
            sig.loc[
                sig["gene_std"].isin(
                    matched_gene_set
                ),
                [
                    "gene_std",
                    "logFC_std",
                    "fdr_std",
                ],
            ]
            .rename(
                columns={
                    "gene_std": "gene",
                    "logFC_std": "model_logFC_SCZ",
                    "fdr_std": "model_fdr_SCZ",
                }
            )
            .merge(
                paper_lookup,
                on="gene",
                how="left",
                validate="one_to_one",
            )
        )


        matched["model"] = model

        matched["analysis"] = (
            level_name
        )


        matched["model_direction"] = np.where(
            matched["model_logFC_SCZ"] > 0,
            "UP",
            np.where(
                matched["model_logFC_SCZ"] < 0,
                "DOWN",
                "ZERO",
            ),
        )


        matched["paper_direction"] = np.where(
            matched["paper_logFC_SCZ"] > 0,
            "UP",
            np.where(
                matched["paper_logFC_SCZ"] < 0,
                "DOWN",
                "ZERO",
            ),
        )


        matched["same_direction"] = (
            np.sign(
                matched["model_logFC_SCZ"]
            )
            ==
            np.sign(
                matched["paper_logFC_SCZ"]
            )
        )


        n_match = len(
            matched
        )

        n_same = int(
            matched["same_direction"].sum()
        )


        recovery = (
            100.0 * n_match / expected_n
        )


        tested_recovery = (
            100.0 * n_match / len(paper_tested)
            if len(paper_tested) > 0
            else np.nan
        )


        matched_genes = sorted(
            matched["gene"].tolist()
        )


        same_direction_genes = sorted(
            matched.loc[
                matched["same_direction"],
                "gene",
            ].tolist()
        )


        print()
        print(
            f"Genes tested                       : {len(de)}"
        )

        print(
            f"Model genes FDR < 0.10            : {len(sig)}"
        )

        print(
            f"  Model UP                        : {n_up}"
        )

        print(
            f"  Model DOWN                      : {n_down}"
        )

        print(
            f"Paper genes FDR < 0.10            : {expected_n}"
        )

        print(
            f"  Paper UP                        : {paper_up}"
        )

        print(
            f"  Paper DOWN                      : {paper_down}"
        )

        print(
            f"Paper significant genes tested    : {len(paper_tested)}"
        )

        print(
            f"Paper significant genes not tested: {len(paper_not_tested)}"
        )

        print(
            f"Significant overlap               : {n_match}/{expected_n}"
        )

        print(
            f"Same-direction overlap            : {n_same}/{expected_n}"
        )

        print(
            f"Recovery                          : {recovery:.2f}%"
        )

        print(
            "Matched genes                     : "
            + (
                ", ".join(matched_genes)
                if matched_genes
                else "NONE"
            )
        )

        print(
            "Same-direction matched genes      : "
            + (
                ", ".join(same_direction_genes)
                if same_direction_genes
                else "NONE"
            )
        )


        if paper_not_tested:

            print(
                "Paper genes not tested            : "
                + ", ".join(
                    sorted(paper_not_tested)
                )
            )


        summaries.append(
            {
                "analysis":
                    level_name,

                "model":
                    model,

                "genes_tested":
                    len(de),

                "model_FDR_lt_010":
                    len(sig),

                "model_up_FDR010":
                    n_up,

                "model_down_FDR010":
                    n_down,

                "paper_FDR_lt_010":
                    expected_n,

                "paper_up_FDR010":
                    paper_up,

                "paper_down_FDR010":
                    paper_down,

                "paper_genes_tested":
                    len(paper_tested),

                "paper_genes_not_tested":
                    len(paper_not_tested),

                "significant_overlap":
                    n_match,

                "same_direction_overlap":
                    n_same,

                "recovery_percent":
                    recovery,

                "recovery_among_tested_percent":
                    tested_recovery,

                "matched_genes":
                    ";".join(
                        matched_genes
                    ),

                "same_direction_genes":
                    ";".join(
                        same_direction_genes
                    ),

                "paper_not_tested_genes":
                    ";".join(
                        sorted(
                            paper_not_tested
                        )
                    ),

                "de_file":
                    str(de_file),
            }
        )


        if len(matched) > 0:
            matches.append(
                matched
            )


    summary = pd.DataFrame(
        summaries
    )


    if matches:

        matched_all = pd.concat(
            matches,
            ignore_index=True,
        )

    else:

        matched_all = pd.DataFrame()


    return summary, matched_all


# ======================================================================================
# DOMAIN COMPARISON — PAPER 19
# ======================================================================================

section(
    "DOMAIN-LEVEL PAPER CONCORDANCE — EXPECT 19"
)

domain_slide_summary, domain_slide_matched = compare_level(
    level_name="donor_domain_level_WITH_SLIDE",
    paper_sig=domain_sig,
    expected_n=19,
    de_file_func=find_domain_de_with_slide,
)


domain_noslide_summary, domain_noslide_matched = compare_level(
    level_name="donor_domain_level_NO_SLIDE",
    paper_sig=domain_sig,
    expected_n=19,
    de_file_func=find_domain_de_noslide,
)


# ======================================================================================
# CELL-TYPE COMPARISON — PAPER 15
# ======================================================================================

section(
    "CELL-TYPE PAPER CONCORDANCE — EXPECT 15"
)

cell_slide_summary, cell_slide_matched = compare_level(
    level_name="donor_cell_type_level_WITH_SLIDE",
    paper_sig=cell_sig,
    expected_n=15,
    de_file_func=find_cell_de_with_slide,
)


cell_noslide_summary, cell_noslide_matched = compare_level(
    level_name="donor_cell_type_level_NO_SLIDE",
    paper_sig=cell_sig,
    expected_n=15,
    de_file_func=find_cell_de_noslide,
)


# ======================================================================================
# SAVE
# ======================================================================================

domain_slide_summary.to_csv(
    DOMAIN_SLIDE_SUMMARY,
    index=False,
)

domain_slide_matched.to_csv(
    DOMAIN_SLIDE_MATCHED,
    index=False,
)

domain_noslide_summary.to_csv(
    DOMAIN_NOSLIDE_SUMMARY,
    index=False,
)

domain_noslide_matched.to_csv(
    DOMAIN_NOSLIDE_MATCHED,
    index=False,
)

cell_slide_summary.to_csv(
    CELL_SLIDE_SUMMARY,
    index=False,
)

cell_slide_matched.to_csv(
    CELL_SLIDE_MATCHED,
    index=False,
)

cell_noslide_summary.to_csv(
    CELL_NOSLIDE_SUMMARY,
    index=False,
)

cell_noslide_matched.to_csv(
    CELL_NOSLIDE_MATCHED,
    index=False,
)


# ======================================================================================
# FINAL TABLES
# ======================================================================================

cols = [
    "model",
    "genes_tested",
    "model_FDR_lt_010",
    "model_up_FDR010",
    "model_down_FDR010",
    "paper_FDR_lt_010",
    "significant_overlap",
    "same_direction_overlap",
    "recovery_percent",
]


section(
    "FINAL DOMAIN-LEVEL SUMMARY — WITH SLIDE — PAPER 19"
)

print(
    domain_slide_summary[
        cols
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}",
    )
)


section(
    "FINAL DOMAIN-LEVEL SUMMARY — NO SLIDE — PAPER 19"
)

print(
    domain_noslide_summary[
        cols
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}",
    )
)


section(
    "FINAL CELL-TYPE SUMMARY — WITH SLIDE — PAPER 15"
)

print(
    cell_slide_summary[
        cols
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}",
    )
)


section(
    "FINAL CELL-TYPE SUMMARY — NO SLIDE — PAPER 15"
)

print(
    cell_noslide_summary[
        cols
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}",
    )
)


section(
    "OUTPUT FILES"
)

print("Paper domain all:")
print(DOMAIN_REF_ALL)

print("\nPaper domain FDR<0.10:")
print(DOMAIN_REF_SIG)

print("\nPaper cell-type all:")
print(CELL_REF_ALL)

print("\nPaper cell-type FDR<0.10:")
print(CELL_REF_SIG)

print("\nDomain WITH-SLIDE comparison:")
print(DOMAIN_SLIDE_SUMMARY)

print("\nDomain WITH-SLIDE matched genes:")
print(DOMAIN_SLIDE_MATCHED)

print("\nDomain NO-SLIDE comparison:")
print(DOMAIN_NOSLIDE_SUMMARY)

print("\nDomain NO-SLIDE matched genes:")
print(DOMAIN_NOSLIDE_MATCHED)

print("\nCell-type WITH-SLIDE comparison:")
print(CELL_SLIDE_SUMMARY)

print("\nCell-type WITH-SLIDE matched genes:")
print(CELL_SLIDE_MATCHED)

print("\nCell-type NO-SLIDE comparison:")
print(CELL_NOSLIDE_SUMMARY)

print("\nCell-type NO-SLIDE matched genes:")
print(CELL_NOSLIDE_MATCHED)


section(
    "FINAL STATUS"
)

print(
    "PASS — PAPER CONCORDANCE READ DIRECTLY FROM CORRECT XLSX SHEETS"
)
