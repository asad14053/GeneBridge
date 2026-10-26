#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

PAPER_FILE = (
    ROOT
    / "data/reference/paper_xenium/paper_concordance_de.xlsx"
)

DE_ROOT = (
    ROOT
    / "outputs/imputation_full/DE/task2_3_preimputation_baseline/domain_de"
)

OUT_FILE = (
    DE_ROOT
    / "original_xenium_domain_vs_paper19_summary.csv"
)

MATCH_FILE = (
    DE_ROOT
    / "original_xenium_domain_vs_paper19_matched_genes.csv"
)

FDR = 0.10


# ======================================================================================
# PAPER REFERENCE
# ======================================================================================

print("=" * 110)
print("PAPER DONOR-DOMAIN REFERENCE")
print("=" * 110)

paper = pd.read_excel(
    PAPER_FILE,
    sheet_name="donor_domain_level",
    engine="openpyxl"
)

required = ["gene", "logFC_SCZ", "fdr_SCZ"]

missing = [
    c for c in required
    if c not in paper.columns
]

if missing:
    raise RuntimeError(
        f"Paper reference missing columns: {missing}\n"
        f"Available columns: {list(paper.columns)}"
    )

paper["gene"] = (
    paper["gene"]
    .astype(str)
    .str.strip()
)

paper["logFC_SCZ"] = pd.to_numeric(
    paper["logFC_SCZ"],
    errors="coerce"
)

paper["fdr_SCZ"] = pd.to_numeric(
    paper["fdr_SCZ"],
    errors="coerce"
)

paper_sig = (
    paper.loc[
        paper["fdr_SCZ"] < FDR
    ]
    .copy()
    .sort_values("fdr_SCZ")
)

if len(paper_sig) != 19:
    raise RuntimeError(
        f"Expected 19 donor-domain paper genes at FDR<0.10; "
        f"found {len(paper_sig)}"
    )

paper_genes = set(
    paper_sig["gene"]
)

paper_lookup = (
    paper_sig
    .set_index("gene")
)

print("Paper genes FDR < 0.10:", len(paper_sig))
print()
print(
    paper_sig[
        ["gene", "logFC_SCZ", "fdr_SCZ"]
    ].to_string(index=False)
)


# ======================================================================================
# COMPARE BOTH ORIGINAL-XENIUM DE MODELS
# ======================================================================================

rows = []
match_rows = []

branches = [
    ("WITH_SLIDE", "with_slide"),
    ("NO_SLIDE", "no_slide"),
]

for label, folder in branches:

    print()
    print("=" * 110)
    print(f"ORIGINAL XENIUM — {label}")
    print("=" * 110)

    de_file = (
        DE_ROOT
        / folder
        / "original_xenium_domain_SCZ_DE_dx_results.csv"
    )

    if not de_file.exists():
        raise FileNotFoundError(de_file)

    de = pd.read_csv(de_file)

    for c in ["gene", "logFC_SCZ", "fdr_SCZ"]:
        if c not in de.columns:
            raise RuntimeError(
                f"{label}: missing {c}; columns={list(de.columns)}"
            )

    de["gene"] = (
        de["gene"]
        .astype(str)
        .str.strip()
    )

    de["logFC_SCZ"] = pd.to_numeric(
        de["logFC_SCZ"],
        errors="coerce"
    )

    de["fdr_SCZ"] = pd.to_numeric(
        de["fdr_SCZ"],
        errors="coerce"
    )

    sig = de.loc[
        de["fdr_SCZ"] < FDR
    ].copy()

    model_sig_genes = set(
        sig["gene"]
    )

    # EXACT definition requested:
    #
    # paper gene name with paper FDR < 0.10
    #        INTERSECT
    # model gene name with model FDR < 0.10
    #
    overlap = sorted(
        paper_genes
        & model_sig_genes
    )

    same_direction = []

    for gene in overlap:

        model_fc = float(
            sig.loc[
                sig["gene"] == gene,
                "logFC_SCZ"
            ].iloc[0]
        )

        paper_fc = float(
            paper_lookup.loc[
                gene,
                "logFC_SCZ"
            ]
        )

        if np.sign(model_fc) == np.sign(paper_fc):
            same_direction.append(gene)

        match_rows.append({
            "analysis": label,
            "gene": gene,
            "paper_logFC_SCZ": paper_fc,
            "paper_fdr_SCZ": float(
                paper_lookup.loc[
                    gene,
                    "fdr_SCZ"
                ]
            ),
            "model_logFC_SCZ": model_fc,
            "model_fdr_SCZ": float(
                sig.loc[
                    sig["gene"] == gene,
                    "fdr_SCZ"
                ].iloc[0]
            ),
            "same_direction":
                np.sign(model_fc) == np.sign(paper_fc),
        })

    n_up = int(
        (sig["logFC_SCZ"] > 0).sum()
    )

    n_down = int(
        (sig["logFC_SCZ"] < 0).sum()
    )

    print("Genes tested                    :", len(de))
    print("Model genes FDR < 0.10         :", len(sig))
    print("  Model UP                     :", n_up)
    print("  Model DOWN                   :", n_down)

    print("Paper genes FDR < 0.10         : 19")

    print(
        "Significant overlap            :",
        f"{len(overlap)}/19"
    )

    print(
        "Same-direction overlap         :",
        f"{len(same_direction)}/19"
    )

    print(
        "Recovery                       :",
        f"{100 * len(overlap) / 19:.2f}%"
    )

    print(
        "Matched genes                  :",
        ", ".join(overlap)
        if overlap else "NONE"
    )

    print(
        "Same-direction matched genes   :",
        ", ".join(same_direction)
        if same_direction else "NONE"
    )

    rows.append({
        "analysis": label,
        "genes_tested": len(de),
        "model_FDR_lt_010": len(sig),
        "model_up_FDR010": n_up,
        "model_down_FDR010": n_down,
        "paper_FDR_lt_010": 19,
        "significant_overlap": len(overlap),
        "same_direction_overlap": len(same_direction),
        "recovery_percent":
            100 * len(overlap) / 19,
        "matched_genes":
            ";".join(overlap),
        "same_direction_genes":
            ";".join(same_direction),
    })


# ======================================================================================
# FINAL SUMMARY
# ======================================================================================

summary = pd.DataFrame(rows)

print()
print("=" * 110)
print("PRE-IMPUTATION ORIGINAL XENIUM vs PAPER 19 — FINAL SUMMARY")
print("=" * 110)

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.2f}"
    )
)

summary.to_csv(
    OUT_FILE,
    index=False
)

pd.DataFrame(
    match_rows
).to_csv(
    MATCH_FILE,
    index=False
)

print("\nSaved summary:")
print(OUT_FILE)

print("\nSaved matched genes:")
print(MATCH_FILE)

print("\nFINAL STATUS: PASS — ORIGINAL XENIUM vs PAPER DOMAIN-19 COMPARISON COMPLETE")
