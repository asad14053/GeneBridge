#!/usr/bin/env python3

"""
06_compare_celltype_DE_to_media10.py

TASK 2.3 — CELL-TYPE-ADJUSTED DE vs MEDIA-10
=============================================

Purpose
-------
Compare SCZ-vs-NTC cell-type-adjusted DE results from all six
NTC-trained OOF imputation models against the significant genes
in the Media-10 reference file.

Reference significant gene definition:
    fdr_SCZ < 0.10

For each model report:
    - genes tested
    - number of DE genes with FDR < 0.10
    - upregulated genes with FDR < 0.10
    - downregulated genes with FDR < 0.10
    - number of Media-10 significant genes present in tested DE set
    - exact overlap:
          model FDR < 0.10 AND Media-10 FDR < 0.10
    - overlap with same logFC direction
    - matched gene names
    - recovery percentage

Expected Media-10 significant set:
    15 genes

IMPORTANT
---------
The reference set is derived directly from the file. It is NOT
hard-coded. If the reference file does not contain 15 genes at
FDR < 0.10, the script prints a warning.
"""

from pathlib import Path
import numpy as np
import pandas as pd


# ======================================================================================
# CONFIG
# ======================================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

REFERENCE_FILE = (
    ROOT
    / "data/reference/paper_xenium"
    / "paper_xenium_layer_adjusted_24donor.csv"
)

DE_ROOT = (
    ROOT
    / "outputs/imputation_full/DE/task2_3_ntc_dysregulation"
    / "celltype_adjusted_de"
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
EXPECTED_REFERENCE_GENES = 15


SUMMARY_FILE = (
    DE_ROOT
    / "all_models_celltype_DE_vs_media10_summary.csv"
)

MATCHED_FILE = (
    DE_ROOT
    / "all_models_celltype_DE_vs_media10_matched_genes.csv"
)

REFERENCE_SIG_FILE = (
    DE_ROOT
    / "media10_reference_FDR010_genes.csv"
)


# ======================================================================================
# HELPERS
# ======================================================================================

def section(title):
    print()
    print("=" * 110)
    print(title)
    print("=" * 110)


def standardize(df, label):

    required = {
        "gene",
        "fdr_SCZ",
        "logFC_SCZ",
    }

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            f"{label}: missing columns {sorted(missing)}\n"
            f"Available columns: {list(df.columns)}"
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
            f"{label}: duplicated genes: {duplicated}"
        )

    return df


# ======================================================================================
# LOAD MEDIA-10 REFERENCE
# ======================================================================================

section("LOAD MEDIA-10 REFERENCE")

if not REFERENCE_FILE.exists():
    raise FileNotFoundError(
        REFERENCE_FILE
    )

ref = pd.read_csv(
    REFERENCE_FILE
)

ref = standardize(
    ref,
    "Media-10 reference",
)


# Strict criterion requested:
# fdr_SCZ < 0.10

ref_sig = (
    ref.loc[
        ref["fdr_SCZ"] < FDR_THRESHOLD
    ]
    .copy()
    .sort_values(
        ["fdr_SCZ", "gene"]
    )
)


ref_sig_genes = set(
    ref_sig["gene"]
)


print("Reference file:")
print(REFERENCE_FILE)

print(
    "\nTotal reference genes:",
    len(ref)
)

print(
    "Reference genes FDR < 0.10:",
    len(ref_sig)
)


# ======================================================================================
# CHECK EXPECTED 15 GENES
# ======================================================================================

if len(ref_sig) != EXPECTED_REFERENCE_GENES:

    print()
    print("!" * 110)

    print(
        f"WARNING: expected {EXPECTED_REFERENCE_GENES} "
        f"Media-10 genes at FDR < 0.10, "
        f"but this file contains {len(ref_sig)}."
    )

    print(
        "Comparison will use the actual significant genes "
        "found in the supplied reference file."
    )

    print("!" * 110)


print("\nMedia-10 FDR < 0.10 genes:")

print(
    ref_sig[
        [
            "gene",
            "logFC_SCZ",
            "fdr_SCZ",
        ]
    ].to_string(
        index=False
    )
)


ref_sig[
    [
        "gene",
        "logFC_SCZ",
        "fdr_SCZ",
    ]
].to_csv(
    REFERENCE_SIG_FILE,
    index=False,
)


# ======================================================================================
# SIX-MODEL COMPARISON
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
        / f"{model}_celltype_adjusted_SCZ_DE_dx_results.csv"
    )


    if not de_file.exists():
        raise FileNotFoundError(
            f"{model}: missing DE file:\n{de_file}"
        )


    de = pd.read_csv(
        de_file
    )

    de = standardize(
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


    # ==================================================================================
    # HOW MANY OF THE REFERENCE GENES WERE TESTED?
    #
    # Important because filterByExpr may remove genes from some models.
    # ==================================================================================

    tested_gene_set = set(
        de["gene"]
    )


    reference_tested = (
        ref_sig_genes
        & tested_gene_set
    )


    reference_not_tested = (
        ref_sig_genes
        - tested_gene_set
    )


    # ==================================================================================
    # EXACT SIGNIFICANT OVERLAP
    #
    # Both:
    #
    #   Media-10 fdr_SCZ < 0.10
    #
    # AND
    #
    #   model fdr_SCZ < 0.10
    # ==================================================================================

    matched = (
        sig.loc[
            sig["gene"].isin(
                ref_sig_genes
            )
        ]
        .copy()
    )


    # ==================================================================================
    # ADD MEDIA-10 EFFECT / FDR
    # ==================================================================================

    ref_lookup = (
        ref_sig[
            [
                "gene",
                "logFC_SCZ",
                "fdr_SCZ",
            ]
        ]
        .rename(
            columns={
                "logFC_SCZ":
                    "media10_logFC_SCZ",

                "fdr_SCZ":
                    "media10_fdr_SCZ",
            }
        )
    )


    matched = matched.merge(
        ref_lookup,
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
    # DIRECTION
    # ==================================================================================

    matched["direction_match"] = (
        np.sign(
            matched["model_logFC_SCZ"]
        )
        ==
        np.sign(
            matched["media10_logFC_SCZ"]
        )
    )


    matched["media10_direction"] = np.where(
        matched["media10_logFC_SCZ"] > 0,
        "UP",
        np.where(
            matched["media10_logFC_SCZ"] < 0,
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


    # ==================================================================================
    # COUNTS
    # ==================================================================================

    n_sig = len(
        sig
    )

    n_up = len(
        up
    )

    n_down = len(
        down
    )

    n_match = len(
        matched
    )

    n_same_direction = int(
        matched["direction_match"].sum()
    )


    recovery_rate = (
        n_match / len(ref_sig)
        if len(ref_sig) > 0
        else np.nan
    )


    matched_genes = sorted(
        matched["gene"].tolist()
    )


    same_direction_genes = sorted(
        matched.loc[
            matched["direction_match"],
            "gene",
        ].tolist()
    )


    # ==================================================================================
    # PRINT
    # ==================================================================================

    print("DE file:")
    print(de_file)

    print()
    print(
        f"Genes tested                       : {len(de)}"
    )

    print(
        f"Genes FDR < 0.10                  : {n_sig}"
    )

    print(
        f"  Upregulated                     : {n_up}"
    )

    print(
        f"  Downregulated                   : {n_down}"
    )

    print(
        f"Media-10 reference FDR < 0.10     : {len(ref_sig)}"
    )

    print(
        f"Reference genes tested in model   : {len(reference_tested)}"
    )

    print(
        f"Reference genes filtered/not tested: {len(reference_not_tested)}"
    )

    print(
        f"Exact significant overlap         : {n_match}"
    )

    print(
        f"Overlap + same direction          : {n_same_direction}"
    )

    print(
        f"Recovery                          : "
        f"{n_match}/{len(ref_sig)} "
        f"({100 * recovery_rate:.2f}%)"
    )

    print(
        "Matched genes                     :",
        ", ".join(matched_genes)
        if matched_genes
        else "NONE"
    )

    print(
        "Matched + same-direction genes    :",
        ", ".join(same_direction_genes)
        if same_direction_genes
        else "NONE"
    )


    if reference_not_tested:

        print(
            "Reference genes not tested       :",
            ", ".join(
                sorted(
                    reference_not_tested
                )
            )
        )


    # ==================================================================================
    # SUMMARY ROW
    # ==================================================================================

    summary_rows.append(
        {
            "model":
                model,

            "genes_tested":
                len(de),

            "FDR_lt_010":
                n_sig,

            "up_FDR010":
                n_up,

            "down_FDR010":
                n_down,

            "media10_FDR_lt_010":
                len(ref_sig),

            "media10_genes_tested":
                len(reference_tested),

            "media10_genes_not_tested":
                len(reference_not_tested),

            "matched_media10_FDR010":
                n_match,

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
                    same_direction_genes
                ),

            "reference_not_tested_genes":
                ";".join(
                    sorted(
                        reference_not_tested
                    )
                ),
        }
    )


    # ==================================================================================
    # MATCHED-GENE TABLE
    # ==================================================================================

    if len(matched) > 0:

        matched_rows.append(
            matched[
                [
                    "model",
                    "gene",
                    "media10_logFC_SCZ",
                    "media10_fdr_SCZ",
                    "media10_direction",
                    "model_logFC_SCZ",
                    "model_fdr_SCZ",
                    "model_direction",
                    "direction_match",
                ]
            ]
        )


# ======================================================================================
# FINAL SUMMARY
# ======================================================================================

section(
    "FINAL SIX-MODEL CELL-TYPE DE vs MEDIA-10 SUMMARY"
)


summary = pd.DataFrame(
    summary_rows
)


display_cols = [
    "model",
    "genes_tested",
    "FDR_lt_010",
    "up_FDR010",
    "down_FDR010",
    "media10_FDR_lt_010",
    "matched_media10_FDR010",
    "matched_same_direction",
    "recovery_percent",
]


print(
    summary[
        display_cols
    ].to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}",
    )
)


# ======================================================================================
# MATCHED GENE NAMES
# ======================================================================================

section(
    "MATCHED MEDIA-10 GENES BY MODEL"
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
# SAVE
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
            "media10_logFC_SCZ",
            "media10_fdr_SCZ",
            "media10_direction",
            "model_logFC_SCZ",
            "model_fdr_SCZ",
            "model_direction",
            "direction_match",
        ]
    )


matched_all.to_csv(
    MATCHED_FILE,
    index=False,
)


# ======================================================================================
# FINAL STATUS
# ======================================================================================

section(
    "OUTPUT FILES"
)


print(
    "Reference significant genes:"
)

print(
    REFERENCE_SIG_FILE
)


print(
    "\nSix-model comparison:"
)

print(
    SUMMARY_FILE
)


print(
    "\nMatched genes:"
)

print(
    MATCHED_FILE
)


section(
    "FINAL STATUS"
)

print(
    "PASS — CELL-TYPE-ADJUSTED DE vs MEDIA-10 COMPARISON COMPLETE"
)
