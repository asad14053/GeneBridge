#!/usr/bin/env python3

"""
Experiment 5.3

Create a Huuki snRNA-seq reference containing all brains except the brain
represented by the Br8667 Xenium dataset.

Huuki structure
---------------
obs["BrNum"]  = brain/donor identifier, for example Br8667
obs["Sample"] = tissue variant, for example Br8667_ant or Br8667_mid

Therefore, all Sample variants associated with the Xenium brain are removed.

Expected:
    10 brains -> 9 brains
    Br8667_ant and Br8667_mid are both removed
    all genes and layers are retained
"""

from __future__ import annotations

import gc
import os
import re
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

HUUKI_INPUT = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

XENIUM_INPUT = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "spatial_data_xenium_Br8667_vista.h5ad"
)

OUTPUT_H5AD = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "seq_data_huuki_snrna_9brains_excluding_Br8667.h5ad"
)

# Compatibility link for the earlier filename.
COMPATIBILITY_LINK = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "seq_data_huuki_snrna_9samples_excluding_xenium_overlap.h5ad"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667/experiment_5_3"
)

SUMMARY_CSV = OUTPUT_DIR / "experiment_5_3_summary.csv"
BEFORE_COUNTS_CSV = OUTPUT_DIR / "brain_variant_counts_before.csv"
AFTER_COUNTS_CSV = OUTPUT_DIR / "brain_variant_counts_after.csv"
REMOVED_VARIANTS_CSV = OUTPUT_DIR / "removed_brain_variants.csv"

EXPECTED_INPUT_BRAINS = 10
EXPECTED_OUTPUT_BRAINS = 9

BRAIN_COLUMN = "BrNum"
VARIANT_COLUMN = "Sample"


def log(message: object = "") -> None:
    print(message, flush=True)


def section(title: str) -> None:
    log()
    log("=" * 100)
    log(title)
    log("=" * 100)


def extract_brain_id(value: object) -> str | None:
    """Extract identifiers such as Br8667 from arbitrary text."""

    if pd.isna(value):
        return None

    match = re.search(
        r"\bBr\d+\b",
        str(value),
        flags=re.IGNORECASE,
    )

    if match is None:
        return None

    value = match.group(0)

    return "Br" + re.sub(
        r"\D",
        "",
        value,
    )


def unique_nonempty(values) -> list[str]:
    output: list[str] = []
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
            output.append(text)
            seen.add(text)

    return output


def identify_xenium_brain(
    xenium: ad.AnnData,
) -> tuple[str, str]:

    candidate_columns = [
        "BrNum",
        "brain_id",
        "BrainID",
        "donor",
        "donor_id",
        "subject_id",
        "Sample",
        "sample",
        "sample_id",
        "SAMPLE_ID",
    ]

    identifiers: set[str] = set()

    for column in candidate_columns:
        if column not in xenium.obs.columns:
            continue

        for value in xenium.obs[column]:
            brain_id = extract_brain_id(value)

            if brain_id is not None:
                identifiers.add(brain_id)

    # Also infer from the input path. This identifies Br8667 even when
    # the Xenium obs metadata does not contain a donor column.
    path_brain = extract_brain_id(str(XENIUM_INPUT))

    if path_brain is not None:
        identifiers.add(path_brain)

    if len(identifiers) != 1:
        raise RuntimeError(
            "Expected exactly one Xenium brain identifier, but found: "
            f"{sorted(identifiers)}"
        )

    brain_id = next(iter(identifiers))

    if path_brain == brain_id:
        method = "Xenium metadata and/or input path"
    else:
        method = "Xenium obs metadata"

    return brain_id, method


def make_count_table(
    brain_values: pd.Series,
    variant_values: pd.Series,
) -> pd.DataFrame:

    table = pd.DataFrame(
        {
            "BrNum": brain_values.astype(str),
            "Sample": variant_values.astype(str),
        }
    )

    return (
        table.groupby(
            ["BrNum", "Sample"],
            observed=True,
        )
        .size()
        .reset_index(name="n_cells")
        .sort_values(
            ["BrNum", "Sample"]
        )
        .reset_index(drop=True)
    )


def remove_unused_categories(
    adata: ad.AnnData,
) -> None:

    for column in adata.obs.columns:
        if isinstance(
            adata.obs[column].dtype,
            pd.CategoricalDtype,
        ):
            adata.obs[column] = (
                adata.obs[column]
                .cat.remove_unused_categories()
            )


def main() -> None:

    start_time = time.time()

    section(
        "EXPERIMENT 5.3: REMOVE ALL HUUKI VARIANTS "
        "ASSOCIATED WITH THE XENIUM BRAIN"
    )

    for path in [
        HUUKI_INPUT,
        XENIUM_INPUT,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required input not found:\n{path}"
            )

    OUTPUT_H5AD.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    log(f"Huuki input:\n{HUUKI_INPUT}")
    log()
    log(f"Xenium input:\n{XENIUM_INPUT}")
    log()
    log(f"Output:\n{OUTPUT_H5AD}")

    # ------------------------------------------------------------------
    # Open both files in backed mode to inspect metadata.
    # ------------------------------------------------------------------

    section("READING METADATA")

    huuki = ad.read_h5ad(
        HUUKI_INPUT,
        backed="r",
    )

    xenium = ad.read_h5ad(
        XENIUM_INPUT,
        backed="r",
    )

    log(f"Huuki shape: {huuki.shape}")
    log(f"Xenium shape: {xenium.shape}")
    log(f"Huuki layers: {list(huuki.layers.keys())}")

    for required_column in [
        BRAIN_COLUMN,
        VARIANT_COLUMN,
    ]:
        if required_column not in huuki.obs.columns:
            raise KeyError(
                f"Huuki obs column '{required_column}' is missing.\n"
                f"Available columns:\n{list(huuki.obs.columns)}"
            )

    xenium_brain, identification_method = (
        identify_xenium_brain(xenium)
    )

    log()
    log(f"Xenium brain identified as: {xenium_brain}")
    log(f"Identification method: {identification_method}")

    # Normalize BrNum using either BrNum itself or the Sample label.
    huuki_brain_ids = huuki.obs[BRAIN_COLUMN].map(
        extract_brain_id
    )

    missing_brain_mask = huuki_brain_ids.isna()

    if missing_brain_mask.any():
        huuki_brain_ids.loc[missing_brain_mask] = (
            huuki.obs.loc[
                missing_brain_mask,
                VARIANT_COLUMN,
            ]
            .map(extract_brain_id)
        )

    if huuki_brain_ids.isna().any():
        raise RuntimeError(
            f"Could not determine BrNum for "
            f"{int(huuki_brain_ids.isna().sum())} Huuki cells."
        )

    huuki_variants = (
        huuki.obs[VARIANT_COLUMN]
        .astype(str)
    )

    input_brains = unique_nonempty(
        huuki_brain_ids
    )

    input_variants = unique_nonempty(
        huuki_variants
    )

    before_counts = make_count_table(
        huuki_brain_ids,
        huuki_variants,
    )

    before_counts.to_csv(
        BEFORE_COUNTS_CSV,
        index=False,
    )

    section("HUUKI BRAIN AND VARIANT SUMMARY")

    log(before_counts.to_string(index=False))
    log()
    log(f"Unique brains: {len(input_brains)}")
    log(f"Unique Sample variants: {len(input_variants)}")

    if len(input_brains) != EXPECTED_INPUT_BRAINS:
        raise RuntimeError(
            f"Expected {EXPECTED_INPUT_BRAINS} Huuki brains, "
            f"but found {len(input_brains)}:\n{input_brains}"
        )

    if xenium_brain not in input_brains:
        raise RuntimeError(
            f"Xenium brain {xenium_brain} was not found "
            "in the Huuki BrNum values."
        )

    # ------------------------------------------------------------------
    # Identify every Sample variant associated with the Xenium brain.
    # ------------------------------------------------------------------

    removal_mask = (
        huuki_brain_ids.astype(str)
        == xenium_brain
    ).to_numpy()

    keep_mask = ~removal_mask

    removed_variants = sorted(
        unique_nonempty(
            huuki_variants.loc[removal_mask]
        )
    )

    removed_cells = int(
        removal_mask.sum()
    )

    retained_cells = int(
        keep_mask.sum()
    )

    section("VARIANTS TO REMOVE")

    log(f"Brain to remove: {xenium_brain}")
    log(f"Cells to remove: {removed_cells:,}")
    log(f"Associated Sample variants: {removed_variants}")

    if removed_cells == 0:
        raise RuntimeError(
            f"No Huuki cells were found for {xenium_brain}."
        )

    if len(removed_variants) == 0:
        raise RuntimeError(
            f"No Sample variants were found for {xenium_brain}."
        )

    pd.DataFrame(
        {
            "excluded_brain": xenium_brain,
            "excluded_variant": removed_variants,
        }
    ).to_csv(
        REMOVED_VARIANTS_CSV,
        index=False,
    )

    input_shape = huuki.shape
    input_gene_names = (
        huuki.var_names.astype(str).to_numpy()
    )
    input_layers = list(
        huuki.layers.keys()
    )

    # ------------------------------------------------------------------
    # Select retained cells while the Huuki file is open in backed mode.
    #
    # Row indices are sorted, making this substantially safer and faster
    # than arbitrary backed gene-column filtering.
    # ------------------------------------------------------------------

    section("LOADING THE NINE-BRAIN HUUKI SUBSET")

    retained_indices = np.flatnonzero(
        keep_mask
    )

    subset_start = time.time()

    huuki_filtered = (
        huuki[
            retained_indices,
            :,
        ]
        .to_memory()
    )

    huuki.file.close()
    xenium.file.close()

    del huuki
    del xenium
    gc.collect()

    remove_unused_categories(
        huuki_filtered
    )

    log(
        f"Subset loaded in "
        f"{(time.time() - subset_start) / 60:.2f} minutes."
    )

    log(f"Subset shape: {huuki_filtered.shape}")
    log(f"Layers retained: {list(huuki_filtered.layers.keys())}")

    # ------------------------------------------------------------------
    # Validate the resulting brain and variant structure.
    # ------------------------------------------------------------------

    section("PRE-WRITE VALIDATION")

    output_brain_ids = huuki_filtered.obs[
        BRAIN_COLUMN
    ].map(extract_brain_id)

    missing_output_brains = output_brain_ids.isna()

    if missing_output_brains.any():
        output_brain_ids.loc[missing_output_brains] = (
            huuki_filtered.obs.loc[
                missing_output_brains,
                VARIANT_COLUMN,
            ]
            .map(extract_brain_id)
        )

    output_variants = (
        huuki_filtered.obs[VARIANT_COLUMN]
        .astype(str)
    )

    output_brains = unique_nonempty(
        output_brain_ids
    )

    output_variant_labels = unique_nonempty(
        output_variants
    )

    after_counts = make_count_table(
        output_brain_ids,
        output_variants,
    )

    after_counts.to_csv(
        AFTER_COUNTS_CSV,
        index=False,
    )

    log(after_counts.to_string(index=False))
    log()
    log(f"Input brains: {len(input_brains)}")
    log(f"Output brains: {len(output_brains)}")
    log(f"Input Sample variants: {len(input_variants)}")
    log(f"Output Sample variants: {len(output_variant_labels)}")
    log(f"Removed variants: {removed_variants}")
    log(f"Removed cells: {removed_cells:,}")
    log(f"Retained cells: {retained_cells:,}")

    if len(output_brains) != EXPECTED_OUTPUT_BRAINS:
        raise RuntimeError(
            f"Expected {EXPECTED_OUTPUT_BRAINS} output brains, "
            f"but found {len(output_brains)}:\n{output_brains}"
        )

    if xenium_brain in output_brains:
        raise RuntimeError(
            f"{xenium_brain} is still present after filtering."
        )

    remaining_removed_variants = sorted(
        set(removed_variants).intersection(
            output_variant_labels
        )
    )

    if remaining_removed_variants:
        raise RuntimeError(
            "Some excluded variants remain in the output: "
            f"{remaining_removed_variants}"
        )

    if huuki_filtered.n_obs != retained_cells:
        raise RuntimeError(
            "Retained cell count does not match the expected count."
        )

    if huuki_filtered.n_vars != input_shape[1]:
        raise RuntimeError(
            "Gene count changed during brain-level filtering."
        )

    if not np.array_equal(
        huuki_filtered.var_names.astype(str).to_numpy(),
        input_gene_names,
    ):
        raise RuntimeError(
            "Gene names or gene order changed."
        )

    if list(huuki_filtered.layers.keys()) != input_layers:
        raise RuntimeError(
            "The Huuki layer list changed during filtering."
        )

    huuki_filtered.uns[
        "experiment_5_3"
    ] = {
        "description": (
            "Huuki reference after removing every tissue variant "
            "associated with the Xenium brain"
        ),
        "xenium_brain": xenium_brain,
        "removed_variants": removed_variants,
        "brain_column": BRAIN_COLUMN,
        "variant_column": VARIANT_COLUMN,
        "input_brains": len(input_brains),
        "output_brains": len(output_brains),
        "input_variants": len(input_variants),
        "output_variants": len(output_variant_labels),
        "removed_cells": removed_cells,
        "retained_cells": retained_cells,
        "huuki_input": str(HUUKI_INPUT),
        "xenium_input": str(XENIUM_INPUT),
    }

    # ------------------------------------------------------------------
    # Write output.
    # ------------------------------------------------------------------

    section("WRITING EXPERIMENT 5.3 OUTPUT")

    if OUTPUT_H5AD.exists() or OUTPUT_H5AD.is_symlink():
        OUTPUT_H5AD.unlink()

    if (
        COMPATIBILITY_LINK.exists()
        or COMPATIBILITY_LINK.is_symlink()
    ):
        COMPATIBILITY_LINK.unlink()

    write_start = time.time()

    huuki_filtered.write_h5ad(
        OUTPUT_H5AD,
        compression="lzf",
    )

    log(
        f"Write completed in "
        f"{(time.time() - write_start) / 60:.2f} minutes."
    )

    expected_shape = huuki_filtered.shape

    del huuki_filtered
    gc.collect()

    # Create a relative symbolic link using the earlier filename.
    os.symlink(
        OUTPUT_H5AD.name,
        COMPATIBILITY_LINK,
    )

    # ------------------------------------------------------------------
    # Read-back validation.
    # ------------------------------------------------------------------

    section("READ-BACK VALIDATION")

    check = ad.read_h5ad(
        OUTPUT_H5AD,
        backed="r",
    )

    check_brains = check.obs[
        BRAIN_COLUMN
    ].map(extract_brain_id)

    missing_check_brains = check_brains.isna()

    if missing_check_brains.any():
        check_brains.loc[missing_check_brains] = (
            check.obs.loc[
                missing_check_brains,
                VARIANT_COLUMN,
            ]
            .map(extract_brain_id)
        )

    check_variants = unique_nonempty(
        check.obs[VARIANT_COLUMN]
    )

    check_brain_list = unique_nonempty(
        check_brains
    )

    if check.shape != expected_shape:
        raise RuntimeError(
            f"Read-back shape mismatch: {check.shape} "
            f"versus {expected_shape}"
        )

    if len(check_brain_list) != EXPECTED_OUTPUT_BRAINS:
        raise RuntimeError(
            f"Read-back contains {len(check_brain_list)} brains, "
            f"not {EXPECTED_OUTPUT_BRAINS}."
        )

    if xenium_brain in check_brain_list:
        raise RuntimeError(
            f"Read-back still contains {xenium_brain}."
        )

    if set(removed_variants).intersection(check_variants):
        raise RuntimeError(
            "One or more removed variants are present after read-back."
        )

    if not np.array_equal(
        check.var_names.astype(str).to_numpy(),
        input_gene_names,
    ):
        raise RuntimeError(
            "Read-back gene names or order changed."
        )

    summary = pd.DataFrame(
        [
            {
                "experiment": "5.3",
                "huuki_input": str(HUUKI_INPUT),
                "xenium_input": str(XENIUM_INPUT),
                "output_h5ad": str(OUTPUT_H5AD),
                "excluded_brain": xenium_brain,
                "excluded_variants": ";".join(removed_variants),
                "input_brains": len(input_brains),
                "output_brains": len(check_brain_list),
                "input_sample_variants": len(input_variants),
                "output_sample_variants": len(check_variants),
                "input_cells": input_shape[0],
                "removed_cells": removed_cells,
                "output_cells": check.n_obs,
                "input_genes": input_shape[1],
                "output_genes": check.n_vars,
                "layers": ";".join(check.layers.keys()),
                "output_size_GB": round(
                    OUTPUT_H5AD.stat().st_size / 1024**3,
                    4,
                ),
                "runtime_minutes": round(
                    (time.time() - start_time) / 60,
                    2,
                ),
            }
        ]
    )

    check.file.close()

    summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    section("EXPERIMENT 5.3 COMPLETED")

    log(summary.to_string(index=False))
    log()
    log("Primary output:")
    log(OUTPUT_H5AD)
    log()
    log("Compatibility link:")
    log(COMPATIBILITY_LINK)
    log()
    log("Removed variants:")
    log(removed_variants)


if __name__ == "__main__":
    main()
