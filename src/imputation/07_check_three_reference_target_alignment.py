#!/usr/bin/env python3

"""
Validate three Huuki snRNA references against the same Br8667 Xenium target.

This script does not modify any H5AD file.

Checks
------
1. The raw and QC Xenium files have the same gene panel and gene order.
2. Each Xenium target gene maps uniquely to each reference.
3. Target-gene coverage for each reference.
4. Duplicate gene names.
5. Reference cohort composition:
   - Experiment 5: Br8667 only
   - Experiment 5.1: all 10 Huuki brains
   - Experiment 5.2/5.3: 9 brains, excluding Br8667
6. All Br8667 variants are absent from the nine-brain reference.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

TARGET_RAW = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "spatial_data_xenium_Br8667_vista.h5ad"
)

TARGET_QC = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "spatial_data_xenium_Br8667_vista_qc.h5ad"
)

REFERENCE_EX5 = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "seq_data_huuki_snrna_Br8667_vista.h5ad"
)

REFERENCE_EX51 = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

REFERENCE_NINE_BRAIN = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "seq_data_huuki_snrna_9samples_excluding_xenium_overlap.h5ad"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667/reference_alignment_check"
)

SUMMARY_CSV = OUTPUT_DIR / "reference_target_alignment_summary.csv"
TARGET_COMPARISON_CSV = OUTPUT_DIR / "xenium_raw_vs_qc_comparison.csv"
COHORT_SUMMARY_CSV = OUTPUT_DIR / "reference_cohort_summary.csv"

MINIMUM_TARGET_COVERAGE = 0.98


def log(message: object = "") -> None:
    print(message, flush=True)


def section(title: str) -> None:
    log()
    log("=" * 100)
    log(title)
    log("=" * 100)



def gene_checksum(genes: list[str]) -> str:
    joined = "\n".join(genes).encode("utf-8")

    return hashlib.sha256(joined).hexdigest()


def extract_brain_id(value: object) -> str | None:
    if pd.isna(value):
        return None

    match = re.search(
        r"Br\d+",
        str(value),
        flags=re.IGNORECASE,
    )

    if match is None:
        return None

    digits = re.sub(
        r"\D",
        "",
        match.group(0),
    )

    return f"Br{digits}"


def unique_nonempty(values) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()

    for value in values:
        if pd.isna(value):
            continue

        text = str(value).strip()

        if not text or text.lower() in {
            "nan",
            "none",
            "na",
        }:
            continue

        if text not in seen:
            result.append(text)
            seen.add(text)

    return result


def get_brain_series(adata: ad.AnnData) -> pd.Series:
    if "BrNum" in adata.obs.columns:
        brains = adata.obs["BrNum"].map(
            extract_brain_id
        )
    else:
        brains = pd.Series(
            [None] * adata.n_obs,
            index=adata.obs_names,
            dtype="object",
        )

    if "Sample" in adata.obs.columns:
        missing = brains.isna()

        if missing.any():
            brains.loc[missing] = (
                adata.obs.loc[
                    missing,
                    "Sample",
                ]
                .map(extract_brain_id)
            )

    return brains


def get_variants(adata: ad.AnnData) -> list[str]:
    if "Sample" not in adata.obs.columns:
        return []

    return sorted(
        unique_nonempty(
            adata.obs["Sample"]
        )
    )


def target_raw_qc_check() -> dict[str, object]:
    section("CHECKING RAW AND QC XENIUM TARGETS")

    if not TARGET_RAW.exists():
        raise FileNotFoundError(TARGET_RAW)

    if not TARGET_QC.exists():
        raise FileNotFoundError(TARGET_QC)

    raw = ad.read_h5ad(
        TARGET_RAW,
        backed="r",
    )

    qc = ad.read_h5ad(
        TARGET_QC,
        backed="r",
    )

    raw_genes = raw.var_names.astype(str).tolist()
    qc_genes = qc.var_names.astype(str).tolist()

    same_gene_set = set(raw_genes) == set(qc_genes)
    same_gene_order = raw_genes == qc_genes

    raw_duplicates = int(
        pd.Index(raw_genes).duplicated().sum()
    )

    qc_duplicates = int(
        pd.Index(qc_genes).duplicated().sum()
    )

    result = {
        "raw_file": str(TARGET_RAW),
        "qc_file": str(TARGET_QC),
        "raw_cells": raw.n_obs,
        "qc_cells": qc.n_obs,
        "raw_genes": raw.n_vars,
        "qc_genes": qc.n_vars,
        "same_gene_set": same_gene_set,
        "same_gene_order": same_gene_order,
        "raw_duplicate_genes": raw_duplicates,
        "qc_duplicate_genes": qc_duplicates,
        "raw_gene_checksum": gene_checksum(raw_genes),
        "qc_gene_checksum": gene_checksum(qc_genes),
    }

    raw.file.close()
    qc.file.close()

    pd.DataFrame(
        [result]
    ).to_csv(
        TARGET_COMPARISON_CSV,
        index=False,
    )

    log(pd.DataFrame([result]).to_string(index=False))

    if not same_gene_set:
        raise RuntimeError(
            "Raw and QC Xenium files do not contain the same gene set. "
            "The three experiments must use one canonical target."
        )

    if not same_gene_order:
        raise RuntimeError(
            "Raw and QC Xenium files contain the same genes but in "
            "different orders. Use only the QC target for all experiments."
        )

    if raw_duplicates > 0 or qc_duplicates > 0:
        raise RuntimeError(
            "The Xenium target contains duplicated gene names."
        )

    log()
    log("Canonical target for all experiments:")
    log(TARGET_QC)

    return result


def inspect_reference(
    experiment: str,
    reference_path: Path,
    target_genes: list[str],
    expected_brains: int,
    should_contain_br8667: bool,
    br8667_only: bool = False,
) -> tuple[dict[str, object], dict[str, object]]:

    section(f"CHECKING {experiment}")

    if not reference_path.exists():
        raise FileNotFoundError(
            f"Reference not found:\n{reference_path}"
        )

    reference = ad.read_h5ad(
        reference_path,
        backed="r",
    )

    reference_genes = (
        reference.var_names.astype(str).tolist()
    )

    reference_index = pd.Index(
        reference_genes
    )

    target_index = pd.Index(
        target_genes
    )

    reference_duplicate_count = int(
        reference_index.duplicated().sum()
    )

    target_duplicate_count = int(
        target_index.duplicated().sum()
    )

    mapping_indices = reference_index.get_indexer(
        target_index
    )

    found_mask = mapping_indices >= 0

    shared_genes = target_index[
        found_mask
    ].astype(str).tolist()

    missing_target_genes = target_index[
        ~found_mask
    ].astype(str).tolist()

    target_coverage = (
        len(shared_genes) / len(target_genes)
        if target_genes
        else 0.0
    )

    reindexable_to_target_order = (
        reference_duplicate_count == 0
        and target_duplicate_count == 0
        and len(missing_target_genes) == 0
    )

    mapping_table = pd.DataFrame(
        {
            "target_gene_order": np.arange(
                1,
                len(target_genes) + 1,
            ),
            "target_gene": target_genes,
            "present_in_reference": found_mask,
            "reference_var_index": mapping_indices,
        }
    )

    safe_experiment = (
        experiment.lower()
        .replace(" ", "_")
        .replace("/", "_")
        .replace(".", "_")
    )

    mapping_table.to_csv(
        OUTPUT_DIR
        / f"{safe_experiment}_target_gene_mapping.csv",
        index=False,
    )

    pd.DataFrame(
        {
            "missing_target_gene": missing_target_genes,
        }
    ).to_csv(
        OUTPUT_DIR
        / f"{safe_experiment}_missing_target_genes.csv",
        index=False,
    )

    pd.DataFrame(
        {
            "shared_target_gene": shared_genes,
        }
    ).to_csv(
        OUTPUT_DIR
        / f"{safe_experiment}_shared_target_genes.csv",
        index=False,
    )

    brain_series = get_brain_series(
        reference
    )

    brains = sorted(
        unique_nonempty(brain_series)
    )

    variants = get_variants(
        reference
    )

    contains_br8667 = "Br8667" in brains

    br8667_variants = [
        variant
        for variant in variants
        if extract_brain_id(variant) == "Br8667"
    ]

    brain_count_pass = (
        len(brains) == expected_brains
    )

    br8667_presence_pass = (
        contains_br8667 == should_contain_br8667
    )

    br8667_only_pass = True

    if br8667_only:
        br8667_only_pass = (
            brains == ["Br8667"]
        )

    cohort_pass = (
        brain_count_pass
        and br8667_presence_pass
        and br8667_only_pass
    )

    coverage_pass = (
        target_coverage >= MINIMUM_TARGET_COVERAGE
    )

    duplicate_pass = (
        reference_duplicate_count == 0
        and target_duplicate_count == 0
    )

    overall_pass = (
        coverage_pass
        and duplicate_pass
        and cohort_pass
    )

    alignment_status = (
        "COMPLETE_100_PERCENT"
        if target_coverage == 1.0
        and duplicate_pass
        else (
            "PASS_AT_LEAST_98_PERCENT"
            if coverage_pass and duplicate_pass
            else "REVIEW_REQUIRED"
        )
    )

    result = {
        "experiment": experiment,
        "reference_file": str(reference_path),
        "target_file": str(TARGET_QC),
        "reference_cells": reference.n_obs,
        "reference_genes": reference.n_vars,
        "target_genes": len(target_genes),
        "shared_target_genes": len(shared_genes),
        "missing_target_genes": len(missing_target_genes),
        "target_coverage_percent": round(
            target_coverage * 100,
            4,
        ),
        "reference_duplicate_genes": reference_duplicate_count,
        "target_duplicate_genes": target_duplicate_count,
        "reindexable_to_exact_target_order": (
            reindexable_to_target_order
        ),
        "reference_gene_checksum": gene_checksum(
            reference_genes
        ),
        "target_gene_checksum": gene_checksum(
            target_genes
        ),
        "alignment_status": alignment_status,
        "overall_pass": overall_pass,
        "layers": ";".join(
            reference.layers.keys()
        ),
    }

    cohort_result = {
        "experiment": experiment,
        "reference_file": str(reference_path),
        "expected_brains": expected_brains,
        "observed_brains": len(brains),
        "brain_ids": ";".join(brains),
        "observed_sample_variants": len(variants),
        "sample_variants": ";".join(variants),
        "contains_Br8667": contains_br8667,
        "expected_Br8667_presence": should_contain_br8667,
        "Br8667_variants": ";".join(br8667_variants),
        "brain_count_pass": brain_count_pass,
        "Br8667_presence_pass": br8667_presence_pass,
        "Br8667_only_pass": br8667_only_pass,
        "cohort_pass": cohort_pass,
    }

    log()
    log("Alignment:")
    log(pd.DataFrame([result]).to_string(index=False))

    log()
    log("Cohort:")
    log(pd.DataFrame([cohort_result]).to_string(index=False))

    if missing_target_genes:
        log()
        log("Missing Xenium target genes:")
        log(missing_target_genes[:30])

    reference.file.close()

    return result, cohort_result


def main() -> None:
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    target_raw_qc_check()

    target = ad.read_h5ad(
        TARGET_QC,
        backed="r",
    )

    target_genes = (
        target.var_names.astype(str).tolist()
    )

    target.file.close()

    nine_brain_reference = REFERENCE_NINE_BRAIN

    reference_definitions = [
        {
            "experiment": "Experiment_5_Br8667_only",
            "path": REFERENCE_EX5,
            "expected_brains": 1,
            "should_contain_br8667": True,
            "br8667_only": True,
        },
        {
            "experiment": "Experiment_5_1_All_10_Huuki",
            "path": REFERENCE_EX51,
            "expected_brains": 10,
            "should_contain_br8667": True,
            "br8667_only": False,
        },
        {
            "experiment": "Experiment_5_2_or_5_3_Nine_Huuki",
            "path": nine_brain_reference,
            "expected_brains": 9,
            "should_contain_br8667": False,
            "br8667_only": False,
        },
    ]

    alignment_results = []
    cohort_results = []

    for definition in reference_definitions:
        alignment, cohort = inspect_reference(
            experiment=definition["experiment"],
            reference_path=definition["path"],
            target_genes=target_genes,
            expected_brains=definition["expected_brains"],
            should_contain_br8667=(
                definition["should_contain_br8667"]
            ),
            br8667_only=definition["br8667_only"],
        )

        alignment_results.append(alignment)
        cohort_results.append(cohort)

    alignment_summary = pd.DataFrame(
        alignment_results
    )

    cohort_summary = pd.DataFrame(
        cohort_results
    )

    alignment_summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    cohort_summary.to_csv(
        COHORT_SUMMARY_CSV,
        index=False,
    )

    section("FINAL ALIGNMENT SUMMARY")

    log(
        alignment_summary[
            [
                "experiment",
                "reference_cells",
                "reference_genes",
                "target_genes",
                "shared_target_genes",
                "missing_target_genes",
                "target_coverage_percent",
                "reindexable_to_exact_target_order",
                "alignment_status",
                "overall_pass",
            ]
        ].to_string(index=False)
    )

    section("FINAL COHORT SUMMARY")

    log(
        cohort_summary[
            [
                "experiment",
                "expected_brains",
                "observed_brains",
                "observed_sample_variants",
                "contains_Br8667",
                "cohort_pass",
            ]
        ].to_string(index=False)
    )

    failed = alignment_summary.loc[
        ~alignment_summary["overall_pass"]
    ]

    if not failed.empty:
        raise RuntimeError(
            "One or more experiments failed alignment or cohort "
            "validation. Review the generated CSV reports."
        )

    section("ALL THREE EXPERIMENTS PASSED")

    log(f"Alignment summary:\n{SUMMARY_CSV}")
    log()
    log(f"Cohort summary:\n{COHORT_SUMMARY_CSV}")
    log()
    log(
        "The same QC Br8667 Xenium target is suitable for all "
        "three reference experiments."
    )


if __name__ == "__main__":
    main()
