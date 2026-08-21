#!/usr/bin/env python3

"""
Filter the already-created full Baituk H5AD to Baituk–Huuki shared genes.

Input
-----
data/processed/snrnaseq/Baituk/
baituk_snrna_reference_full_allgenes.h5ad

Output
------
data/processed/snrnaseq/Baituk/
baituk_snrna_reference_huuki_overlap.h5ad

The sparse X matrix is filtered directly on disk. The full matrix is not
loaded into memory.
"""

from __future__ import annotations

import copy
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

INPUT_H5AD = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/Baituk"
    / "baituk_snrna_reference_full_allgenes.h5ad"
)

SHARED_GENE_FILES = [
    (
        PROJECT_ROOT
        / "outputs/snrnaseq_branch/03_batiuk_huuki_gene_overlap"
        / "baituk_huuki_shared_genes_final.csv"
    ),
    (
        PROJECT_ROOT
        / "outputs/snrnaseq_branch/03_batiuk_huuki_gene_overlap"
        / "batiuk_huuki_shared_genes.csv"
    ),
]

OUTPUT_H5AD = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/Baituk"
    / "baituk_snrna_reference_huuki_overlap.h5ad"
)

SUMMARY_CSV = (
    PROJECT_ROOT
    / "outputs/snrnaseq_branch/03_batiuk_huuki_gene_overlap"
    / "baituk_filtered_h5ad_summary.csv"
)

EXPECTED_SHARED_GENES = 36119
CSR_ROW_CHUNK = 2000


def log(message: object = "") -> None:
    print(message, flush=True)


def section(title: str) -> None:
    log()
    log("=" * 100)
    log(title)
    log("=" * 100)


def decode_attribute(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def read_shared_genes() -> tuple[list[str], Path]:
    shared_file = next(
        (path for path in SHARED_GENE_FILES if path.exists()),
        None,
    )

    if shared_file is None:
        checked = "\n".join(str(path) for path in SHARED_GENE_FILES)
        raise FileNotFoundError(
            f"No shared-gene CSV was found. Checked:\n{checked}"
        )

    table = pd.read_csv(shared_file)

    candidate_columns = [
        "gene",
        "shared_gene",
        "gene_name",
        "genes",
    ]

    gene_column = next(
        (
            column
            for column in candidate_columns
            if column in table.columns
        ),
        table.columns[0],
    )

    genes: list[str] = []
    seen: set[str] = set()

    for value in table[gene_column]:
        if pd.isna(value):
            continue

        gene = str(value).strip()

        if not gene or gene.lower() in {"nan", "none", "na"}:
            continue

        if gene not in seen:
            genes.append(gene)
            seen.add(gene)

    if not genes:
        raise RuntimeError(
            f"The shared-gene CSV contains no usable genes: {shared_file}"
        )

    return genes, shared_file


def contiguous_runs(indices: np.ndarray):
    """Return position ranges for consecutive source-column indices."""

    run_start = 0

    for position in range(1, len(indices)):
        if indices[position] != indices[position - 1] + 1:
            yield run_start, position
            run_start = position

    yield run_start, len(indices)


def initialize_sparse_group(
    destination: h5py.File,
    encoding_type: str,
    shape: tuple[int, int],
) -> h5py.Group:
    if "X" in destination:
        del destination["X"]

    group = destination.create_group("X")
    group.attrs["encoding-type"] = encoding_type
    group.attrs["encoding-version"] = "0.1.0"
    group.attrs["shape"] = np.asarray(shape, dtype=np.int64)

    return group


def write_csc_subset(
    source_group: h5py.Group,
    destination: h5py.File,
    selected_columns: np.ndarray,
    n_obs: int,
) -> None:
    section("Streaming CSC matrix")

    source_data = source_group["data"]
    source_indices = source_group["indices"]
    source_indptr = np.asarray(
        source_group["indptr"][:],
        dtype=np.int64,
    )

    column_nnz = (
        source_indptr[selected_columns + 1]
        - source_indptr[selected_columns]
    )

    output_indptr = np.zeros(
        len(selected_columns) + 1,
        dtype=np.int64,
    )

    output_indptr[1:] = np.cumsum(
        column_nnz,
        dtype=np.int64,
    )

    total_nnz = int(output_indptr[-1])

    log(f"Selected columns: {len(selected_columns):,}")
    log(f"Output nonzero values: {total_nnz:,}")

    output_group = initialize_sparse_group(
        destination,
        "csc_matrix",
        (n_obs, len(selected_columns)),
    )

    chunk_length = min(
        max(100_000, 1),
        max(total_nnz, 1),
    )

    output_data = output_group.create_dataset(
        "data",
        shape=(total_nnz,),
        dtype=source_data.dtype,
        chunks=(min(1_000_000, max(total_nnz, 1)),),
        compression="lzf",
    )

    output_indices = output_group.create_dataset(
        "indices",
        shape=(total_nnz,),
        dtype=source_indices.dtype,
        chunks=(min(1_000_000, max(total_nnz, 1)),),
        compression="lzf",
    )

    output_group.create_dataset(
        "indptr",
        data=output_indptr,
        compression="lzf",
    )

    processed_columns = 0
    start_time = time.time()

    for run_number, (position_start, position_end) in enumerate(
        contiguous_runs(selected_columns),
        start=1,
    ):
        source_first_column = int(selected_columns[position_start])
        source_last_column = int(selected_columns[position_end - 1])

        source_start = int(
            source_indptr[source_first_column]
        )

        source_end = int(
            source_indptr[source_last_column + 1]
        )

        output_start = int(
            output_indptr[position_start]
        )

        output_end = int(
            output_indptr[position_end]
        )

        output_data[output_start:output_end] = (
            source_data[source_start:source_end]
        )

        output_indices[output_start:output_end] = (
            source_indices[source_start:source_end]
        )

        processed_columns = position_end

        if (
            run_number == 1
            or processed_columns % 1000 == 0
            or processed_columns == len(selected_columns)
        ):
            elapsed = time.time() - start_time
            log(
                f"Copied {processed_columns:,}/"
                f"{len(selected_columns):,} genes "
                f"({elapsed / 60:.1f} minutes)"
            )


def write_csr_subset(
    source_group: h5py.Group,
    destination: h5py.File,
    selected_columns: np.ndarray,
    n_obs: int,
    n_vars: int,
) -> None:
    section("Streaming CSR matrix")

    source_data = source_group["data"]
    source_indices = source_group["indices"]
    source_indptr = np.asarray(
        source_group["indptr"][:],
        dtype=np.int64,
    )

    column_mapping = np.full(
        n_vars,
        -1,
        dtype=np.int64,
    )

    column_mapping[selected_columns] = np.arange(
        len(selected_columns),
        dtype=np.int64,
    )

    output_group = initialize_sparse_group(
        destination,
        "csr_matrix",
        (n_obs, len(selected_columns)),
    )

    output_data = output_group.create_dataset(
        "data",
        shape=(0,),
        maxshape=(None,),
        dtype=source_data.dtype,
        chunks=(1_000_000,),
        compression="lzf",
    )

    output_indices = output_group.create_dataset(
        "indices",
        shape=(0,),
        maxshape=(None,),
        dtype=source_indices.dtype,
        chunks=(1_000_000,),
        compression="lzf",
    )

    output_indptr = np.zeros(
        n_obs + 1,
        dtype=np.int64,
    )

    written_nnz = 0
    start_time = time.time()

    for row_start in range(0, n_obs, CSR_ROW_CHUNK):
        row_end = min(
            row_start + CSR_ROW_CHUNK,
            n_obs,
        )

        source_start = int(source_indptr[row_start])
        source_end = int(source_indptr[row_end])

        chunk_data = np.asarray(
            source_data[source_start:source_end]
        )

        chunk_indices = np.asarray(
            source_indices[source_start:source_end]
        )

        mapped_indices = column_mapping[chunk_indices]
        keep = mapped_indices >= 0

        kept_data = chunk_data[keep]
        kept_indices = mapped_indices[keep].astype(
            source_indices.dtype,
            copy=False,
        )

        previous_nnz = written_nnz
        written_nnz += len(kept_data)

        output_data.resize((written_nnz,))
        output_indices.resize((written_nnz,))

        output_data[previous_nnz:written_nnz] = kept_data
        output_indices[previous_nnz:written_nnz] = kept_indices

        local_indptr = (
            source_indptr[row_start : row_end + 1]
            - source_start
        )

        keep_prefix = np.empty(
            len(keep) + 1,
            dtype=np.int64,
        )

        keep_prefix[0] = 0
        np.cumsum(
            keep.astype(np.int64),
            out=keep_prefix[1:],
        )

        row_counts = (
            keep_prefix[local_indptr[1:]]
            - keep_prefix[local_indptr[:-1]]
        )

        output_indptr[row_start + 1 : row_end + 1] = (
            previous_nnz
            + np.cumsum(row_counts, dtype=np.int64)
        )

        elapsed = time.time() - start_time

        log(
            f"Processed cells {row_start:,}:{row_end:,} / "
            f"{n_obs:,}; output nnz={written_nnz:,}; "
            f"elapsed={elapsed / 60:.1f} minutes"
        )

    output_group.create_dataset(
        "indptr",
        data=output_indptr,
        compression="lzf",
    )


def write_dense_subset(
    source_dataset: h5py.Dataset,
    destination: h5py.File,
    selected_columns: np.ndarray,
    n_obs: int,
) -> None:
    section("Streaming dense matrix")

    if "X" in destination:
        del destination["X"]

    output = destination.create_dataset(
        "X",
        shape=(n_obs, len(selected_columns)),
        dtype=source_dataset.dtype,
        chunks=(min(500, n_obs), min(1000, len(selected_columns))),
        compression="lzf",
    )

    output.attrs["encoding-type"] = "array"
    output.attrs["encoding-version"] = "0.2.0"

    row_chunk = 500

    for row_start in range(0, n_obs, row_chunk):
        row_end = min(row_start + row_chunk, n_obs)

        complete_block = np.asarray(
            source_dataset[row_start:row_end, :]
        )

        output[row_start:row_end, :] = complete_block[
            :,
            selected_columns,
        ]

        log(
            f"Processed cells {row_start:,}:{row_end:,} / {n_obs:,}"
        )


def main() -> None:
    section("FILTER EXISTING FULL BAITUK H5AD")

    if not INPUT_H5AD.exists():
        raise FileNotFoundError(
            f"Full Baituk H5AD not found:\n{INPUT_H5AD}"
        )

    shared_genes, shared_gene_file = read_shared_genes()

    log(f"Input: {INPUT_H5AD}")
    log(f"Shared-gene CSV: {shared_gene_file}")
    log(f"Output: {OUTPUT_H5AD}")
    log(f"Shared genes from CSV: {len(shared_genes):,}")

    if len(shared_genes) != EXPECTED_SHARED_GENES:
        raise RuntimeError(
            f"Expected {EXPECTED_SHARED_GENES:,} shared genes, "
            f"but found {len(shared_genes):,}."
        )

    section("READING BAITUK METADATA")

    backed = ad.read_h5ad(
        INPUT_H5AD,
        backed="r",
    )

    log(backed)

    if not backed.var_names.is_unique:
        raise RuntimeError(
            "Baituk var_names are not unique."
        )

    selected_columns = backed.var_names.get_indexer(
        shared_genes
    )

    if np.any(selected_columns < 0):
        missing = np.asarray(shared_genes)[
            selected_columns < 0
        ]

        raise RuntimeError(
            f"Baituk is missing {len(missing)} shared genes. "
            f"Examples: {missing[:20].tolist()}"
        )

    n_obs = backed.n_obs
    n_vars = backed.n_vars

    obs = backed.obs.copy()
    var = backed.var.iloc[selected_columns].copy()

    var["alignment_gene"] = shared_genes
    var.index = pd.Index(
        shared_genes,
        name=backed.var_names.name,
    )

    original_uns = copy.deepcopy(
        dict(backed.uns)
    )

    layer_names = list(backed.layers.keys())

    if layer_names:
        raise RuntimeError(
            "The full Baituk H5AD unexpectedly contains layers: "
            f"{layer_names}. This script intentionally filters X only."
        )

    backed.file.close()

    section("DETECTING MATRIX STORAGE")

    with h5py.File(INPUT_H5AD, "r") as source:
        source_x = source["X"]

        if isinstance(source_x, h5py.Group):
            encoding_type = decode_attribute(
                source_x.attrs.get(
                    "encoding-type",
                    "",
                )
            )

            matrix_dtype = source_x["data"].dtype

        elif isinstance(source_x, h5py.Dataset):
            encoding_type = decode_attribute(
                source_x.attrs.get(
                    "encoding-type",
                    "array",
                )
            )

            matrix_dtype = source_x.dtype

        else:
            raise RuntimeError(
                "Unsupported H5AD X representation."
            )

    log(f"Input X encoding: {encoding_type}")
    log(f"Input X dtype: {matrix_dtype}")
    log(f"Input shape: {n_obs:,} cells × {n_vars:,} genes")
    log(
        f"Output shape: {n_obs:,} cells × "
        f"{len(shared_genes):,} genes"
    )

    section("CREATING OUTPUT METADATA")

    OUTPUT_H5AD.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    SUMMARY_CSV.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if OUTPUT_H5AD.exists():
        OUTPUT_H5AD.unlink()

    placeholder = sparse.csr_matrix(
        (n_obs, len(shared_genes)),
        dtype=matrix_dtype,
    )

    output_adata = ad.AnnData(
        X=placeholder,
        obs=obs,
        var=var,
        uns=original_uns,
    )

    output_adata.uns["gene_alignment"] = {
        "paired_dataset": "Huuki",
        "shared_gene_file": str(shared_gene_file),
        "n_shared_genes": len(shared_genes),
        "gene_order": "Baituk var_names order",
    }

    output_adata.write_h5ad(
        OUTPUT_H5AD,
        compression="lzf",
    )

    del output_adata
    del placeholder

    section("STREAMING FILTERED X MATRIX")

    with h5py.File(INPUT_H5AD, "r") as source:
        with h5py.File(OUTPUT_H5AD, "r+") as destination:
            source_x = source["X"]

            if isinstance(source_x, h5py.Group):
                encoding_type = decode_attribute(
                    source_x.attrs.get(
                        "encoding-type",
                        "",
                    )
                )

                if encoding_type == "csc_matrix":
                    write_csc_subset(
                        source_x,
                        destination,
                        selected_columns,
                        n_obs,
                    )

                elif encoding_type == "csr_matrix":
                    write_csr_subset(
                        source_x,
                        destination,
                        selected_columns,
                        n_obs,
                        n_vars,
                    )

                else:
                    raise RuntimeError(
                        f"Unsupported sparse encoding: {encoding_type}"
                    )

            elif isinstance(source_x, h5py.Dataset):
                write_dense_subset(
                    source_x,
                    destination,
                    selected_columns,
                    n_obs,
                )

    section("VALIDATION")

    check = ad.read_h5ad(
        OUTPUT_H5AD,
        backed="r",
    )

    if check.shape != (
        n_obs,
        len(shared_genes),
    ):
        raise RuntimeError(
            f"Unexpected output shape: {check.shape}"
        )

    if not np.array_equal(
        check.var_names.to_numpy(),
        np.asarray(shared_genes),
    ):
        raise RuntimeError(
            "Output Baituk gene names or order are incorrect."
        )

    log(check)
    log("Baituk gene names and order validated: TRUE")

    check.file.close()

    summary = pd.DataFrame(
        [
            {
                "dataset": "Baituk",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad": str(OUTPUT_H5AD),
                "shared_gene_file": str(shared_gene_file),
                "n_cells": n_obs,
                "n_input_genes": n_vars,
                "n_output_genes": len(shared_genes),
                "input_X_encoding": encoding_type,
                "output_size_GB": round(
                    OUTPUT_H5AD.stat().st_size / 1024**3,
                    4,
                ),
            }
        ]
    )

    summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    section("DONE")

    log(summary.to_string(index=False))
    log()
    log("Created filtered Baituk H5AD:")
    log(OUTPUT_H5AD)


if __name__ == "__main__":
    main()
