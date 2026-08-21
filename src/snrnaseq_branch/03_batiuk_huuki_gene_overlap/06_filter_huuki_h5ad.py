#!/usr/bin/env python3

"""
Filter the full Huuki H5AD using the exact genes and gene order in the
filtered Baituk H5AD.

Input
-----
Huuki full:
data/processed/snrnaseq/sce_DLPFC_annotated/
huuki_snrna_reference_full_allgenes.h5ad

Filtered Baituk:
data/processed/snrnaseq/Baituk/
baituk_snrna_reference_huuki_overlap.h5ad

Output
------
data/processed/snrnaseq/sce_DLPFC_annotated/
huuki_snrna_reference_baituk_overlap.h5ad
"""

from __future__ import annotations

import gc
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

BAITUK_FILTERED = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/Baituk"
    / "baituk_snrna_reference_huuki_overlap.h5ad"
)

HUUKI_OUTPUT = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
    / "huuki_snrna_reference_baituk_overlap.h5ad"
)

SUMMARY_CSV = (
    PROJECT_ROOT
    / "outputs/snrnaseq_branch/03_batiuk_huuki_gene_overlap"
    / "huuki_filtered_h5ad_summary.csv"
)

EXPECTED_SHARED_GENES = 36119


def log(message: object = "") -> None:
    print(message, flush=True)


def section(title: str) -> None:
    log()
    log("=" * 100)
    log(title)
    log("=" * 100)


def main() -> None:
    section("FILTER FULL HUUKI H5AD")

    for path in [
        HUUKI_INPUT,
        BAITUK_FILTERED,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required input not found:\n{path}"
            )

    section("READING FILTERED BAITUK GENE ORDER")

    baituk = ad.read_h5ad(
        BAITUK_FILTERED,
        backed="r",
    )

    shared_genes = baituk.var_names.astype(str).tolist()

    log(f"Filtered Baituk shape: {baituk.shape}")
    log(f"Shared genes: {len(shared_genes):,}")

    baituk.file.close()
    del baituk

    if len(shared_genes) != EXPECTED_SHARED_GENES:
        raise RuntimeError(
            f"Expected {EXPECTED_SHARED_GENES:,} shared genes, "
            f"but filtered Baituk has {len(shared_genes):,}."
        )

    section("LOADING FULL HUUKI INTO MEMORY")

    start_time = time.time()

    huuki = ad.read_h5ad(
        HUUKI_INPUT,
    )

    log(f"Full Huuki: {huuki}")
    log(
        f"Load time: {(time.time() - start_time) / 60:.1f} minutes"
    )

    if not huuki.var_names.is_unique:
        raise RuntimeError(
            "Huuki var_names are not unique."
        )

    section("MATCHING SHARED GENES")

    huuki_indices = huuki.var_names.get_indexer(
        shared_genes
    )

    if np.any(huuki_indices < 0):
        missing = np.asarray(shared_genes)[
            huuki_indices < 0
        ]

        raise RuntimeError(
            f"Huuki is missing {len(missing)} genes. "
            f"Examples: {missing[:20].tolist()}"
        )

    log(
        f"Matched {len(huuki_indices):,} genes in the exact "
        "filtered Baituk order."
    )

    section("SUBSETTING HUUKI")

    subset_start = time.time()

    huuki_filtered = huuki[
        :,
        huuki_indices,
    ].copy()

    del huuki
    gc.collect()

    huuki_filtered.var[
        "alignment_gene"
    ] = shared_genes

    huuki_filtered.var_names = pd.Index(
        shared_genes,
        name=huuki_filtered.var_names.name,
    )

    huuki_filtered.uns[
        "gene_alignment"
    ] = {
        "paired_dataset": "Baituk",
        "paired_file": str(BAITUK_FILTERED),
        "n_shared_genes": len(shared_genes),
        "gene_order": "filtered Baituk var_names order",
    }

    log(f"Filtered Huuki: {huuki_filtered}")
    log(f"Layers retained: {list(huuki_filtered.layers.keys())}")
    log(
        f"Subset time: {(time.time() - subset_start) / 60:.1f} minutes"
    )

    section("WRITING FILTERED HUUKI")

    HUUKI_OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    SUMMARY_CSV.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if HUUKI_OUTPUT.exists():
        HUUKI_OUTPUT.unlink()

    write_start = time.time()

    huuki_filtered.write_h5ad(
        HUUKI_OUTPUT,
        compression="lzf",
    )

    expected_shape = huuki_filtered.shape
    expected_layers = list(
        huuki_filtered.layers.keys()
    )

    del huuki_filtered
    gc.collect()

    log(
        f"Write time: {(time.time() - write_start) / 60:.1f} minutes"
    )

    section("VALIDATING BOTH FILTERED FILES")

    baituk_check = ad.read_h5ad(
        BAITUK_FILTERED,
        backed="r",
    )

    huuki_check = ad.read_h5ad(
        HUUKI_OUTPUT,
        backed="r",
    )

    if huuki_check.shape != expected_shape:
        raise RuntimeError(
            f"Unexpected Huuki output shape: {huuki_check.shape}"
        )

    if baituk_check.n_vars != huuki_check.n_vars:
        raise RuntimeError(
            "Filtered Baituk and Huuki gene counts differ."
        )

    identical_order = np.array_equal(
        baituk_check.var_names.to_numpy(),
        huuki_check.var_names.to_numpy(),
    )

    if not identical_order:
        raise RuntimeError(
            "Filtered Baituk and Huuki gene order differs."
        )

    log(f"Baituk filtered shape: {baituk_check.shape}")
    log(f"Huuki filtered shape: {huuki_check.shape}")
    log(f"Huuki layers: {list(huuki_check.layers.keys())}")
    log("Gene names and order identical: TRUE")

    summary = pd.DataFrame(
        [
            {
                "dataset": "Baituk",
                "file": str(BAITUK_FILTERED),
                "n_cells": baituk_check.n_obs,
                "n_genes": baituk_check.n_vars,
                "layers": ";".join(
                    baituk_check.layers.keys()
                ),
            },
            {
                "dataset": "Huuki",
                "file": str(HUUKI_OUTPUT),
                "n_cells": huuki_check.n_obs,
                "n_genes": huuki_check.n_vars,
                "layers": ";".join(
                    huuki_check.layers.keys()
                ),
            },
        ]
    )

    baituk_check.file.close()
    huuki_check.file.close()

    summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    section("DONE")

    log(summary.to_string(index=False))
    log()
    log("Created filtered Huuki H5AD:")
    log(HUUKI_OUTPUT)
    log()
    log("Retained Huuki layers:")
    log(expected_layers)


if __name__ == "__main__":
    main()
