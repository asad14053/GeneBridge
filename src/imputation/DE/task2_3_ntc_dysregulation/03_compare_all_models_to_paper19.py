#!/usr/bin/env python3

"""
03_compare_all_models_to_paper19.py

TASK 2.3
========
Compare SCZ differential-expression results from all six NTC-trained
OOF imputation models against the 19 SCZ-associated genes reported
in the reference Xenium analysis.

MODELS
------
1. VISTA
2. gimVI
3. ENVI
4. TransImp
5. Tangram
6. SpaGE

INPUT DE DESIGN
---------------
- 24 donors
    - 12 NTC
    - 12 SCZ
- 7 spatial domains
- 168 donor × domain pseudobulks
- 300 OOF-imputed genes
- no-expression-filter DE branch

For every model, report:
    1. Number of genes tested
    2. Number of genes with SCZ FDR < 0.10
    3. Number upregulated at FDR < 0.10
    4. Number downregulated at FDR < 0.10
    5. Number of exact overlaps with the paper's FDR < 0.10 genes
    6. Number of overlapping genes with matching direction
    7. Names of recovered paper genes

Paper reference:
    paper_xenium_layer_adjusted_24donor.csv

The paper significant set is defined independently as:
    paper fdr_SCZ < 0.10

A paper gene is considered recovered only when:
    paper fdr_SCZ < 0.10
AND
    model fdr_SCZ < 0.10
"""

from pathlib import Path

import numpy as np
import pandas as pd


# ======================================================================================
# CONFIGURATION
# ======================================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

PAPER_FILE = (
    ROOT
    / "data"
    / "reference"
    / "paper_xenium"
    / "paper_xenium_layer_adjusted_24donor.csv"
)

DE_ROOT = (
    ROOT
    / "outputs"
    / "imputation_full"
    / "DE"
    / "task2_3_ntc_dysregulation"
    / "layer_adjusted_de_nofilter"
)

MODELS = [
    "VISTA",
    "gimVI",
    "ENVI",
    "TransImp",
    "Tangram",
    "SpaGE",
]

FDR_THRESHOLD = 0.10


# ======================================================================================
# OUTPUT FILES
# ======================================================================================

SUMMARY_FILE = (
    DE_ROOT
    / "all_models_vs_paper19_DE_summary.csv"
)

MATCHED_GENES_FILE = (
    DE_ROOT
    / "all_models_vs_paper19_matched_genes.csv"
)

PAPER19_FILE = (
    DE_ROOT
    / "paper_reference_19genes_FDR010.csv"
)


# ======================================================================================
# HELPERS
# ======================================================================================

def section(title):
    print()
    print("=" * 110)
    print(title)
    print("=" * 110)


def check_file(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Missing required file:\n{path}"
        )


def prepare_de_table(df, label):
    """
    Validate and standardize a DE result table.
    """

    required = {
        "gene",
        "fdr_SCZ",
        "logFC_SCZ",
    }

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            f"{label}: missing required columns: "
            f"{sorted(missing)}\n"
            f"Available columns:\n{list(df.columns)}"
        )

    df = df.copy()

    df["gene"] = (
        df["gene"]
        .astype(str)
        .str.strip()
    )

    df["fdr_SCZ"] = pd.to_numeric(
        df["fdr_SCZ"],
        errors="coerce",
    )

    df["logFC_SCZ"] = pd.to_numeric(
        df["logFC_SCZ"],
        errors="coerce",
    )

    if df["gene"].duplicated().any():

        duplicated = (
            df.loc[
                df["gene"].duplicated(keep=False),
                "gene",
            ]
            .unique()
            .tolist()
        )

        raise RuntimeError(
            f"{label}: duplicated genes detected:\n"
            + ", ".join(duplicated)
        )

    return df


# ======================================================================================
# LOAD PAPER REFERENCE
# ======================================================================================

section(
    "LOAD PAPER MEDIA-10 / REFERENCE XENIUM DE"
)

check_file(
    PAPER_FILE
)

paper = pd.read_csv(
    PAPER_FILE
)

paper = prepare_de_table(
    paper,
    "Paper reference",
)


# ======================================================================================
# DEFINE PAPER SIGNIFICANT GENE SET
# ======================================================================================

paper_sig = (
    paper.loc[
        paper["fdr_SCZ"] < FDR_THRESHOLD
    ]
    .copy()
    .sort_values(
        "fdr_SCZ"
    )
)

paper_sig_genes = set(
    paper_sig["gene"]
)

print(
    "Paper file:",
    PAPER_FILE,
)

print(
    "Paper genes tested:",
    len(paper),
)

print(
    "Paper genes with FDR < 0.10:",
    len(paper_sig),
)


if len(paper_sig) != 19:

    print(
        "\nWARNING:"
        f" Expected 19 paper genes,"
        f" but found {len(paper_sig)}."
    )


# Save the exact reference gene set used.

paper_sig[
    [
        "gene",
        "logFC_SCZ",
        "fdr_SCZ",
    ]
].to_csv(
    PAPER19_FILE,
    index=False,
)


print(
    "\nPaper FDR < 0.10 genes:"
)

print(
    paper_sig[
        [
            "gene",
            "logFC_SCZ",
            "fdr_SCZ",
        ]
    ].to_string(
        index=False
    )
)


# ======================================================================================
# PROCESS ALL SIX MODELS
# ======================================================================================

summary_rows = []
matched_rows = []


for model in MODELS:

    section(
        f"MODEL: {model}"
    )


    de_file = (
        DE_ROOT
        / model
        / f"{model}_layer_adjusted_SCZ_DE_dx_results.csv"
    )


    print(
        "DE file:",
        de_file,
    )


    check_file(
        de_file
    )


    de = pd.read_csv(
        de_file
    )

    de = prepare_de_table(
        de,
        model,
    )


    # ==================================================================================
    # MODEL SIGNIFICANT GENES
    # ==================================================================================

    sig = (
        de.loc[
            de["fdr_SCZ"] < FDR_THRESHOLD
        ]
        .copy()
    )


    up = sig.loc[
        sig["logFC_SCZ"] > 0
    ].copy()


    down = sig.loc[
        sig["logFC_SCZ"] < 0
    ].copy()


    zero_logfc = sig.loc[
        sig["logFC_SCZ"] == 0
    ].copy()


    # ==================================================================================
    # PAPER-19 OVERLAP
    #
    # Both conditions must be true:
    #
    # paper fdr_SCZ < 0.10
    # model fdr_SCZ < 0.10
    # ==================================================================================

    matched = (
        sig.loc[
            sig["gene"].isin(
                paper_sig_genes
            )
        ]
        .copy()
    )


    # Add paper statistics.

    paper_lookup = (
        paper_sig[
            [
                "gene",
                "logFC_SCZ",
                "fdr_SCZ",
            ]
        ]
        .rename(
            columns={
                "logFC_SCZ":
                    "paper_logFC_SCZ",

                "fdr_SCZ":
                    "paper_fdr_SCZ",
            }
        )
    )


    matched = matched.merge(
        paper_lookup,
        on="gene",
        how="left",
        validate="one_to_one",
    )


    matched = matched.rename(
        columns={
            "logFC_SCZ":
                "model_logFC_SCZ",

            "fdr_SCZ":
                "model_fdr_SCZ",
        }
    )


    # ==================================================================================
    # DIRECTION CONCORDANCE
    # ==================================================================================

    matched["direction_match"] = (
        np.sign(
            matched["model_logFC_SCZ"]
        )
        ==
        np.sign(
            matched["paper_logFC_SCZ"]
        )
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


    matched["model_direction"] = np.where(
        matched["model_logFC_SCZ"] > 0,
        "UP",
        np.where(
            matched["model_logFC_SCZ"] < 0,
            "DOWN",
            "ZERO",
        ),
    )


    matched["model"] = model


    matched = matched[
        [
            "model",
            "gene",
            "paper_logFC_SCZ",
            "paper_fdr_SCZ",
            "paper_direction",
            "model_logFC_SCZ",
            "model_fdr_SCZ",
            "model_direction",
            "direction_match",
        ]
    ]


    matched = matched.sort_values(
        [
            "model_fdr_SCZ",
            "gene",
        ]
    )


    matched_genes = sorted(
        matched["gene"]
        .tolist()
    )


    matched_same_direction_genes = sorted(
        matched.loc[
            matched["direction_match"],
            "gene",
        ].tolist()
    )


    n_matched = len(
        matched
    )


    n_same_direction = int(
        matched["direction_match"]
        .sum()
    )


    recovery_rate = (
        n_matched
        /
        len(paper_sig)
        if len(paper_sig) > 0
        else np.nan
    )


    # ==================================================================================
    # PRINT MODEL SUMMARY
    # ==================================================================================

    print()

    print(
        f"Genes tested                     : {len(de)}"
    )

    print(
        f"FDR < 0.10                      : {len(sig)}"
    )

    print(
        f"  Upregulated                   : {len(up)}"
    )

    print(
        f"  Downregulated                 : {len(down)}"
    )

    print(
        f"  Zero logFC                    : {len(zero_logfc)}"
    )

    print(
        f"Paper genes FDR < 0.10          : {len(paper_sig)}"
    )

    print(
        f"Exact paper-19 overlap          : {n_matched}"
    )

    print(
        f"Overlap + same direction        : {n_same_direction}"
    )

    print(
        f"Recovery rate                   : {100 * recovery_rate:.2f}%"
    )

    print(
        "Matched genes                   :",
        (
            ", ".join(
                matched_genes
            )
            if matched_genes
            else "NONE"
        ),
    )

    print(
        "Matched + same-direction genes  :",
        (
            ", ".join(
                matched_same_direction_genes
            )
            if matched_same_direction_genes
            else "NONE"
        ),
    )


    # ==================================================================================
    # STORE MODEL SUMMARY
    # ==================================================================================

    summary_rows.append(
        {
            "model":
                model,

            "genes_tested":
                len(de),

            "FDR_lt_010":
                len(sig),

            "up_FDR010":
                len(up),

            "down_FDR010":
                len(down),

            "zero_logFC_FDR010":
                len(zero_logfc),

            "paper_FDR_lt_010":
                len(paper_sig),

            "matched_paper19_FDR010":
                n_matched,

            "matched_same_direction":
                n_same_direction,

            "recovery_rate":
                recovery_rate,

            "recovery_percent":
                100 * recovery_rate,

            "matched_genes":
                ";".join(
                    matched_genes
                ),

            "matched_same_direction_genes":
                ";".join(
                    matched_same_direction_genes
                ),
        }
    )


    if len(matched) > 0:

        matched_rows.append(
            matched
        )


# ======================================================================================
# COMBINED SUMMARY
# ======================================================================================

section(
    "FINAL SIX-MODEL SUMMARY"
)


summary = pd.DataFrame(
    summary_rows
)


summary_columns = [
    "model",
    "genes_tested",
    "FDR_lt_010",
    "up_FDR010",
    "down_FDR010",
    "paper_FDR_lt_010",
    "matched_paper19_FDR010",
    "matched_same_direction",
    "recovery_percent",
]


print(
    summary[
        summary_columns
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}",
    )
)


# ======================================================================================
# PRINT MATCHED GENE NAMES
# ======================================================================================

section(
    "PAPER-19 GENES RECOVERED BY EACH MODEL"
)


for _, row in summary.iterrows():

    genes = (
        row["matched_genes"]
        if row["matched_genes"]
        else "NONE"
    )

    print(
        f"{row['model']:<10}: {genes}"
    )


# ======================================================================================
# SAVE RESULTS
# ======================================================================================

summary.to_csv(
    SUMMARY_FILE,
    index=False,
)


if matched_rows:

    matched_all = pd.concat(
        matched_rows,
        ignore_index=True,
    )

else:

    matched_all = pd.DataFrame(
        columns=[
            "model",
            "gene",
            "paper_logFC_SCZ",
            "paper_fdr_SCZ",
            "paper_direction",
            "model_logFC_SCZ",
            "model_fdr_SCZ",
            "model_direction",
            "direction_match",
        ]
    )


matched_all.to_csv(
    MATCHED_GENES_FILE,
    index=False,
)


# ======================================================================================
# FINAL OUTPUT
# ======================================================================================

section(
    "OUTPUT FILES"
)


print(
    "Paper 19-gene reference:"
)

print(
    PAPER19_FILE
)


print(
    "\nSix-model summary:"
)

print(
    SUMMARY_FILE
)


print(
    "\nMatched-gene details:"
)

print(
    MATCHED_GENES_FILE
)


section(
    "FINAL STATUS"
)

print(
    "PASS — six-model DE comparison to paper FDR < 0.10 genes complete."
)
