#!/usr/bin/env python3

"""
Create matched Baituk and Huuki H5AD files containing the exact same genes
in the exact same order.

Inputs
------
Baituk full:
    data/processed/snrnaseq/Baituk/
    baituk_snrna_reference_full_allgenes.h5ad

Huuki full:
    data/processed/snrnaseq/sce_DLPFC_annotated/
    huuki_snrna_reference_full_allgenes.h5ad

Outputs
-------
Baituk filtered:
    data/processed/snrnaseq/Baituk/
    baituk_snrna_reference_huuki_overlap.h5ad

Huuki filtered:
    data/processed/snrnaseq/sce_DLPFC_annotated/
    huuki_snrna_reference_baituk_overlap.h5ad
"""

from __future__ import annotations

import gc
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BAITUK_FULL = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "Baituk"
    / "baituk_snrna_reference_full_allgenes.h5ad"
)

HUUKI_FULL = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

BAITUK_FILTERED = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "Baituk"
    / "baituk_snrna_reference_huuki_overlap.h5ad"
)

HUUKI_FILTERED = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "sce_DLPFC_annotated"
    / "huuki_snrna_reference_baituk_overlap.h5ad"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "snrnaseq_branch"
    / "03_batiuk_huuki_gene_overlap"
)

SHARED_GENES_CSV = (
    OUTPUT_DIR
    / "baituk_huuki_shared_genes_final.csv"
)

SOURCE_COMPARISON_CSV = (
    OUTPUT_DIR
    / "huuki_gene_source_for_filtering.csv"
)

SUMMARY_CSV = (
    OUTPUT_DIR
    / "baituk_huuki_filtered_h5ad_summary.csv"
)


def clean_gene(value: object) -> str | None:
    if pd.isna(value):
        return None

    value = str(value).strip()

    if value == "" or value.lower() in {
        "nan",
        "none",
        "na",
    }:
        return None

    return value


def first_index_map(
    values,
) -> tuple[dict[str, int], int]:

    mapping: dict[str, int] = {}
    duplicates = 0

    for index, value in enumerate(values):
        gene = clean_gene(value)

        if gene is None:
            continue

        if gene in mapping:
            duplicates += 1
            continue

        mapping[gene] = index

    return mapping, duplicates


def get_huuki_candidates(
    huuki: ad.AnnData,
) -> dict[str, pd.Series]:

    candidates: dict[str, pd.Series] = {
        "var_names": pd.Series(
            huuki.var_names.astype(str),
            dtype="string",
        )
    }

    for column in [
        "gene_name",
        "gene_symbol",
        "symbol",
        "Symbol",
        "gene_id",
    ]:
        if column in huuki.var.columns:
            candidates[column] = pd.Series(
                huuki.var[column].to_numpy(),
                dtype="string",
            )

    return candidates


def validate_input(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Required input does not exist:\n{path}"
        )


def main() -> None:

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    BAITUK_FILTERED.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    HUUKI_FILTERED.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    validate_input(BAITUK_FULL)
    validate_input(HUUKI_FULL)

    print("=" * 100)
    print("READING FULL BAITUK AND HUUKI REFERENCES")
    print("=" * 100)

    baituk = ad.read_h5ad(
        BAITUK_FULL,
        backed="r",
    )

    huuki = ad.read_h5ad(
        HUUKI_FULL,
        backed="r",
    )

    print("\nBaituk:")
    print(baituk)

    print("\nHuuki:")
    print(huuki)

    baituk_map, baituk_duplicates = first_index_map(
        baituk.var_names
    )

    baituk_gene_order = list(
        baituk_map.keys()
    )

    print("\nUnique Baituk genes:")
    print(len(baituk_map))

    print("\nDuplicated Baituk gene names skipped:")
    print(baituk_duplicates)

    ###########################################################################
    # Find the Huuki gene field with maximum exact overlap
    ###########################################################################

    candidate_records = []
    candidate_maps: dict[str, dict[str, int]] = {}

    for source_name, source_values in get_huuki_candidates(
        huuki
    ).items():

        source_map, duplicate_count = first_index_map(
            source_values
        )

        candidate_maps[source_name] = source_map

        shared_count = sum(
            gene in source_map
            for gene in baituk_gene_order
        )

        candidate_records.append(
            {
                "huuki_gene_source": source_name,
                "n_unique_huuki_values": len(source_map),
                "duplicate_values_skipped": duplicate_count,
                "shared_exact_genes": shared_count,
            }
        )

    source_comparison = (
        pd.DataFrame(candidate_records)
        .sort_values(
            [
                "shared_exact_genes",
                "n_unique_huuki_values",
            ],
            ascending=False,
        )
        .reset_index(drop=True)
    )

    source_comparison.to_csv(
        SOURCE_COMPARISON_CSV,
        index=False,
    )

    selected_source = str(
        source_comparison.loc[
            0,
            "huuki_gene_source",
        ]
    )

    huuki_map = candidate_maps[
        selected_source
    ]

    ###########################################################################
    # Shared genes in Baituk order
    ###########################################################################

    shared_genes = [
        gene
        for gene in baituk_gene_order
        if gene in huuki_map
    ]

    if len(shared_genes) == 0:
        raise RuntimeError(
            "No exact shared genes were found."
        )

    baituk_indices = np.asarray(
        [
            baituk_map[gene]
            for gene in shared_genes
        ],
        dtype=np.int64,
    )

    huuki_indices = np.asarray(
        [
            huuki_map[gene]
            for gene in shared_genes
        ],
        dtype=np.int64,
    )

    print("\n" + "=" * 100)
    print("GENE OVERLAP")
    print("=" * 100)

    print(source_comparison.to_string(index=False))

    print("\nSelected Huuki gene source:")
    print(selected_source)

    print("\nShared exact genes:")
    print(len(shared_genes))

    pd.DataFrame(
        {
            "gene_order": np.arange(
                1,
                len(shared_genes) + 1,
            ),
            "gene": shared_genes,
            "baituk_var_index": baituk_indices,
            "huuki_var_index": huuki_indices,
        }
    ).to_csv(
        SHARED_GENES_CSV,
        index=False,
    )

    ###########################################################################
    # Filter Baituk
    #
    # Baituk indices are naturally increasing because shared genes retain
    # Baituk's original order.
    ###########################################################################

    print("\n" + "=" * 100)
    print("FILTERING BAITUK")
    print("=" * 100)

    # Backed column slicing is extremely slow for sparse H5AD matrices.
    # Close the backed object, load the full file sequentially, and then
    # perform the gene filtering in memory.
    baituk.file.close()
    del baituk
    gc.collect()

    print("Loading complete Baituk H5AD sequentially into RAM...")
    baituk_memory = ad.read_h5ad(BAITUK_FULL)

    print("Full Baituk shape:", baituk_memory.shape)
    print("Subsetting Baituk genes in memory...")

    baituk_subset = baituk_memory[
        :,
        baituk_indices,
    ].copy()

    del baituk_memory
    gc.collect()

    baituk_subset.var[
        "original_var_name"
    ] = baituk_subset.var_names.astype(str)

    baituk_subset.var[
        "alignment_gene"
    ] = shared_genes

    baituk_subset.var_names = pd.Index(
        shared_genes
    )

    baituk_subset.uns[
        "gene_alignment"
    ] = {
        "paired_dataset": "Huuki",
        "gene_source": "Baituk var_names",
        "n_shared_genes": len(shared_genes),
        "gene_order": "Baituk original order",
    }

    if BAITUK_FILTERED.exists():
        BAITUK_FILTERED.unlink()

    print("Writing:")
    print(BAITUK_FILTERED)

    baituk_subset.write_h5ad(
        BAITUK_FILTERED,
        compression="lzf",
    )

    baituk_cells = baituk_subset.n_obs
    baituk_genes = baituk_subset.n_vars

    del baituk_subset
    gc.collect()

    ###########################################################################
    # Filter Huuki
    #
    # h5py requires backed fancy indices to be sorted.
    # We read sorted Huuki columns first, then restore Baituk gene order
    # after loading the subset into memory.
    ###########################################################################

    print("\n" + "=" * 100)
    print("FILTERING HUUKI")
    print("=" * 100)

    # Load the complete Huuki file sequentially and filter in memory.
    # In-memory indexing does not require sorted HDF5 indices.
    huuki.file.close()
    del huuki
    gc.collect()

    print("Loading complete Huuki H5AD sequentially into RAM...")
    huuki_memory = ad.read_h5ad(HUUKI_FULL)

    print("Full Huuki shape:", huuki_memory.shape)
    print("Subsetting Huuki genes in memory...")

    huuki_subset = huuki_memory[
        :,
        huuki_indices,
    ].copy()

    del huuki_memory
    gc.collect()

    huuki_subset.var[
        "original_var_name"
    ] = huuki_subset.var_names.astype(str)

    huuki_subset.var[
        "alignment_gene"
    ] = shared_genes

    huuki_subset.var_names = pd.Index(
        shared_genes
    )

    huuki_subset.uns[
        "gene_alignment"
    ] = {
        "paired_dataset": "Baituk",
        "huuki_gene_source": selected_source,
        "n_shared_genes": len(shared_genes),
        "gene_order": "Baituk original order",
    }

    if HUUKI_FILTERED.exists():
        HUUKI_FILTERED.unlink()

    print("Writing:")
    print(HUUKI_FILTERED)

    huuki_subset.write_h5ad(
        HUUKI_FILTERED,
        compression="lzf",
    )

    huuki_cells = huuki_subset.n_obs
    huuki_genes = huuki_subset.n_vars

    del huuki_subset
    gc.collect()

    ###########################################################################
    # Validate both outputs
    ###########################################################################

    print("\n" + "=" * 100)
    print("VALIDATING FILTERED FILES")
    print("=" * 100)

    baituk_check = ad.read_h5ad(
        BAITUK_FILTERED,
        backed="r",
    )

    huuki_check = ad.read_h5ad(
        HUUKI_FILTERED,
        backed="r",
    )

    if baituk_check.n_vars != huuki_check.n_vars:
        raise ValueError(
            "Filtered Baituk and Huuki gene counts differ."
        )

    if not np.array_equal(
        baituk_check.var_names.to_numpy(),
        huuki_check.var_names.to_numpy(),
    ):
        raise ValueError(
            "Filtered Baituk and Huuki gene order differs."
        )

    print("Baituk shape:", baituk_check.shape)
    print("Huuki shape:", huuki_check.shape)
    print("Shared genes:", baituk_check.n_vars)
    print("Gene order identical: TRUE")

    baituk_check.file.close()
    huuki_check.file.close()

    ###########################################################################
    # Summary
    ###########################################################################

    summary = pd.DataFrame(
        [
            {
                "dataset": "Baituk",
                "input_h5ad": str(BAITUK_FULL),
                "output_h5ad": str(BAITUK_FILTERED),
                "n_cells": baituk_cells,
                "n_genes": baituk_genes,
                "gene_source": "var_names",
            },
            {
                "dataset": "Huuki",
                "input_h5ad": str(HUUKI_FULL),
                "output_h5ad": str(HUUKI_FILTERED),
                "n_cells": huuki_cells,
                "n_genes": huuki_genes,
                "gene_source": selected_source,
            },
        ]
    )

    summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    print("\n" + "=" * 100)
    print("DONE")
    print("=" * 100)

    print(summary.to_string(index=False))

    print("\nShared gene list:")
    print(SHARED_GENES_CSV)

    print("\nFiltered Baituk:")
    print(BAITUK_FILTERED)

    print("\nFiltered Huuki:")
    print(HUUKI_FILTERED)

    print("\nSummary:")
    print(SUMMARY_CSV)


if __name__ == "__main__":
    main()
