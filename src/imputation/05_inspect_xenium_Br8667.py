#!/usr/bin/env python

"""
05_inspect_xenium_Br8667.py

Inspect:

1. Xenium .X
2. Xenium layers['counts']
3. Xenium layers['log1p_norm']
4. Whether .X is raw or normalized
5. Spatial coordinates
6. Cell-type and layer annotations
7. Gene overlap with existing Br8667 snRNA-seq
8. Preliminary VISTA readiness
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

############## If Xenium is layer-cell annotated Xenium  #########################

XENIUM_PATH = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/"
      "xenium_Br8667_annotated.h5ad"
)

# ############## If Xenium is raw Xenium  #########################

# XENIUM_PATH = (
#     PROJECT_ROOT
#     / "data/processed/imputation_beta/Br8667/"
#       "xenium_Br8667_annotated.h5ad"
# )

############## If ScRNA is Hukki-scRNA  #########################

SNRNA_PATH = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/"
      "seq_data_huuki_snrna_Br8667_shared.h5ad"
)
############## If ScRNA is Hukki-visium  #########################
# SNRNA_PATH = (
#     PROJECT_ROOT
#     / "data/processed/imputation_beta/Br8667/"
#       "spatial_data_huuki_visium_Br8667.h5ad"
# )

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667"
)

MATRIX_SUMMARY_CSV = (
    OUTPUT_DIR
    / "05_xenium_Br8667_matrix_summary.csv"
)

GENE_OVERLAP_CSV = (
    OUTPUT_DIR
    / "05_xenium_snrna_gene_overlap.csv"
)

MISSING_GENES_CSV = (
    OUTPUT_DIR
    / "05_xenium_genes_missing_from_snrna.csv"
)

READINESS_CSV = (
    OUTPUT_DIR
    / "05_vista_Br8667_readiness_summary.csv"
)

CELLTYPE_COUNTS_CSV = (
    OUTPUT_DIR
    / "05_xenium_Br8667_celltype_counts.csv"
)

LAYER_COUNTS_CSV = (
    OUTPUT_DIR
    / "05_xenium_Br8667_layer_counts.csv"
)


def get_values(
    matrix,
    maximum_values: int = 1_000_000,
) -> np.ndarray:
    """
    Obtain stored matrix values for diagnostics.
    """

    if sparse.issparse(matrix):
        values = np.asarray(matrix.data)
    else:
        values = np.asarray(matrix).reshape(-1)

    if values.size == 0:
        return np.array([], dtype=np.float64)

    if values.size > maximum_values:
        positions = np.linspace(
            0,
            values.size - 1,
            maximum_values,
            dtype=int,
        )
        values = values[positions]

    return values.astype(
        np.float64,
        copy=False,
    )


def summarize_matrix(name: str, matrix) -> dict:
    values = get_values(matrix)

    total_elements = int(
        matrix.shape[0] * matrix.shape[1]
    )

    if sparse.issparse(matrix):
        nnz = int(matrix.nnz)
        matrix_type = type(matrix).__name__
        dtype = str(matrix.dtype)
    else:
        array = np.asarray(matrix)
        nnz = int(np.count_nonzero(array))
        matrix_type = type(matrix).__name__
        dtype = str(array.dtype)

    density = (
        nnz / total_elements
        if total_elements
        else 0.0
    )

    if values.size:
        integer_like_fraction = float(
            np.mean(
                np.isclose(
                    values,
                    np.round(values),
                    atol=1e-6,
                )
            )
        )

        minimum = float(np.min(values))
        maximum = float(np.max(values))
        mean = float(np.mean(values))
        has_negative = bool(np.any(values < 0))
    else:
        integer_like_fraction = 1.0
        minimum = 0.0
        maximum = 0.0
        mean = 0.0
        has_negative = False

    return {
        "matrix": name,
        "matrix_type": matrix_type,
        "shape": f"{matrix.shape[0]} x {matrix.shape[1]}",
        "dtype": dtype,
        "nonzero_values": nnz,
        "density": density,
        "sparsity": 1.0 - density,
        "sampled_minimum": minimum,
        "sampled_maximum": maximum,
        "sampled_mean": mean,
        "integer_like_fraction": integer_like_fraction,
        "has_negative_values": has_negative,
    }


def get_matrix_block(
    matrix,
    rows: np.ndarray,
    columns: np.ndarray,
) -> np.ndarray:
    block = matrix[rows, :][:, columns]

    if sparse.issparse(block):
        return block.toarray()

    return np.asarray(block)


def sampled_matrices_equal(
    first,
    second,
) -> bool:
    if first.shape != second.shape:
        return False

    n_rows, n_columns = first.shape

    rows = np.unique(
        np.linspace(
            0,
            n_rows - 1,
            min(200, n_rows),
            dtype=int,
        )
    )

    columns = np.unique(
        np.linspace(
            0,
            n_columns - 1,
            min(300, n_columns),
            dtype=int,
        )
    )

    first_block = get_matrix_block(
        first,
        rows,
        columns,
    )

    second_block = get_matrix_block(
        second,
        rows,
        columns,
    )

    return bool(
        np.allclose(
            first_block,
            second_block,
            atol=1e-6,
            equal_nan=True,
        )
    )


def get_gene_symbols(
    adata: ad.AnnData,
    dataset_name: str,
) -> tuple[pd.Index, str]:
    """
    Choose the most appropriate gene-symbol column.
    """

    candidates = [
        "gene_symbol",
        "Symbol",
        "gene_name",
        "symbol",
    ]

    for column in candidates:
        if column in adata.var.columns:
            symbols = (
                adata.var[column]
                .astype(str)
                .str.strip()
            )

            print(
                f"{dataset_name}: using var['{column}'] "
                "for gene symbols."
            )

            return pd.Index(symbols), column

    print(
        f"{dataset_name}: using var_names "
        "for gene symbols."
    )

    return (
        pd.Index(
            adata.var_names.astype(str)
        ),
        "var_names",
    )


def classify_x(
    x_matches_counts: bool,
    x_matches_log1p: bool,
) -> str:
    if x_matches_counts and not x_matches_log1p:
        return "raw_counts"

    if x_matches_log1p and not x_matches_counts:
        return "log1p_normalized"

    if x_matches_counts and x_matches_log1p:
        return "matches_both_layers"

    return "does_not_match_tested_layers"


def main() -> None:
    print("=" * 100)
    print("Task-02: Inspect Xenium Br8667 and snRNA overlap")
    print("=" * 100)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    for path in [XENIUM_PATH, SNRNA_PATH]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing required input:\n{path}"
            )

    print("\nXenium input:")
    print(XENIUM_PATH)

    print("\nsnRNA input:")
    print(SNRNA_PATH)

    print("\nLoading Xenium Br8667...")
    xenium = ad.read_h5ad(XENIUM_PATH)

    print("Loading snRNA Br8667...")
    snrna = ad.read_h5ad(SNRNA_PATH)

    print("\nXenium Br8667:")
    print(xenium)

    print("\nsnRNA Br8667:")
    print(snrna)

    print("\nXenium layers:")
    print(list(xenium.layers.keys()))

    print("\nsnRNA layers:")
    print(list(snrna.layers.keys()))

    matrix_records = [
        summarize_matrix(
            "xenium:X",
            xenium.X,
        )
    ]

    for layer_name in xenium.layers.keys():
        matrix_records.append(
            summarize_matrix(
                f"xenium:layer:{layer_name}",
                xenium.layers[layer_name],
            )
        )

    matrix_records.append(
        summarize_matrix(
            "snrna:X",
            snrna.X,
        )
    )

    if "counts" in snrna.layers:
        matrix_records.append(
            summarize_matrix(
                "snrna:layer:counts",
                snrna.layers["counts"],
            )
        )

    matrix_summary = pd.DataFrame(
        matrix_records
    )

    matrix_summary.to_csv(
        MATRIX_SUMMARY_CSV,
        index=False,
    )

    print("\nMatrix summary:")
    print(
        matrix_summary.to_string(
            index=False
        )
    )

    x_matches_counts = False
    x_matches_log1p = False

    if "counts" in xenium.layers:
        x_matches_counts = sampled_matrices_equal(
            xenium.X,
            xenium.layers["counts"],
        )

    if "log1p_norm" in xenium.layers:
        x_matches_log1p = sampled_matrices_equal(
            xenium.X,
            xenium.layers["log1p_norm"],
        )

    x_classification = classify_x(
        x_matches_counts,
        x_matches_log1p,
    )

    print("\nXenium X comparison:")
    print(
        "X matches layers['counts']:",
        x_matches_counts,
    )
    print(
        "X matches layers['log1p_norm']:",
        x_matches_log1p,
    )
    print(
        "X classification:",
        x_classification,
    )

    xenium_symbols, xenium_gene_source = (
        get_gene_symbols(
            xenium,
            "Xenium",
        )
    )

    snrna_symbols, snrna_gene_source = (
        get_gene_symbols(
            snrna,
            "snRNA",
        )
    )

    snrna_symbol_set = set(
        snrna_symbols
    )

    overlap_table = pd.DataFrame(
        {
            "xenium_var_name": (
                xenium.var_names.astype(str)
            ),
            "xenium_gene_symbol": (
                xenium_symbols.astype(str)
            ),
        }
    )

    overlap_table["present_in_snrna"] = (
        overlap_table["xenium_gene_symbol"]
        .isin(snrna_symbol_set)
    )

    overlap_table.to_csv(
        GENE_OVERLAP_CSV,
        index=False,
    )

    missing_table = overlap_table.loc[
        ~overlap_table["present_in_snrna"]
    ].copy()

    missing_table.to_csv(
        MISSING_GENES_CSV,
        index=False,
    )

    shared_gene_count = int(
        overlap_table["present_in_snrna"].sum()
    )

    missing_gene_count = int(
        (~overlap_table["present_in_snrna"]).sum()
    )

    all_xenium_symbols_in_snrna = (
        missing_gene_count == 0
    )

    raw_var_names_subset = (
        set(xenium.var_names.astype(str))
        <= set(snrna.var_names.astype(str))
    )

    print("\nGene overlap:")
    print(
        f"Xenium genes: {xenium.n_vars:,}"
    )
    print(
        f"snRNA genes: {snrna.n_vars:,}"
    )
    print(
        f"Shared by gene symbol: {shared_gene_count:,}"
    )
    print(
        f"Xenium genes missing from snRNA: "
        f"{missing_gene_count:,}"
    )
    print(
        "All Xenium symbols present in snRNA:",
        all_xenium_symbols_in_snrna,
    )
    print(
        "Raw Xenium var_names subset of snRNA var_names:",
        raw_var_names_subset,
    )

    if "cell_type_annotation" in xenium.obs.columns:
        celltype_counts = (
            xenium.obs["cell_type_annotation"]
            .astype(str)
            .value_counts(dropna=False)
            .rename_axis("cell_type_annotation")
            .reset_index(name="n_cells")
        )

        celltype_counts.to_csv(
            CELLTYPE_COUNTS_CSV,
            index=False,
        )

        print("\nXenium cell-type counts:")
        print(
            celltype_counts.to_string(
                index=False
            )
        )

    if "layer_annotation" in xenium.obs.columns:
        layer_counts = (
            xenium.obs["layer_annotation"]
            .astype(str)
            .value_counts(dropna=False)
            .rename_axis("layer_annotation")
            .reset_index(name="n_cells")
        )

        layer_counts.to_csv(
            LAYER_COUNTS_CSV,
            index=False,
        )

        print("\nXenium cortical-layer counts:")
        print(
            layer_counts.to_string(
                index=False
            )
        )

    has_spatial = (
        "spatial" in xenium.obsm
    )

    if has_spatial:
        spatial = np.asarray(
            xenium.obsm["spatial"]
        )

        print("\nXenium spatial coordinates:")
        print("Shape:", spatial.shape)
        print(
            "All finite:",
            bool(np.isfinite(spatial).all()),
        )
        print(
            "X coordinate range:",
            float(spatial[:, 0].min()),
            "to",
            float(spatial[:, 0].max()),
        )
        print(
            "Y coordinate range:",
            float(spatial[:, 1].min()),
            "to",
            float(spatial[:, 1].max()),
        )

    xenium_counts_valid = False
    snrna_counts_valid = False

    if "counts" in xenium.layers:
        xenium_count_record = summarize_matrix(
            "xenium_counts_check",
            xenium.layers["counts"],
        )

        xenium_counts_valid = (
            not xenium_count_record["has_negative_values"]
            and xenium_count_record[
                "integer_like_fraction"
            ] > 0.999
        )

    if "counts" in snrna.layers:
        snrna_count_record = summarize_matrix(
            "snrna_counts_check",
            snrna.layers["counts"],
        )

        snrna_counts_valid = (
            not snrna_count_record["has_negative_values"]
            and snrna_count_record[
                "integer_like_fraction"
            ] > 0.999
        )

    readiness = pd.DataFrame(
        [
            {
                "check": "xenium_has_counts_layer",
                "value": "counts" in xenium.layers,
            },
            {
                "check": "xenium_counts_integer_nonnegative",
                "value": xenium_counts_valid,
            },
            {
                "check": "snrna_has_counts_layer",
                "value": "counts" in snrna.layers,
            },
            {
                "check": "snrna_counts_integer_nonnegative",
                "value": snrna_counts_valid,
            },
            {
                "check": "xenium_has_spatial",
                "value": has_spatial,
            },
            {
                "check": "xenium_has_names",
                "value": "names" in xenium.obs.columns,
            },
            {
                "check": "xenium_has_batch",
                "value": "batch" in xenium.obs.columns,
            },
            {
                "check": "xenium_has_scClassify",
                "value": "scClassify" in xenium.obs.columns,
            },
            {
                "check": "snrna_has_scClassify",
                "value": "scClassify" in snrna.obs.columns,
            },
            {
                "check": "xenium_symbols_subset_of_snrna",
                "value": all_xenium_symbols_in_snrna,
            },
            {
                "check": "raw_var_names_subset",
                "value": raw_var_names_subset,
            },
            {
                "check": "xenium_X_classification",
                "value": x_classification,
            },
            {
                "check": "xenium_gene_source",
                "value": xenium_gene_source,
            },
            {
                "check": "snrna_gene_source",
                "value": snrna_gene_source,
            },
            {
                "check": "shared_gene_count",
                "value": shared_gene_count,
            },
            {
                "check": "missing_gene_count",
                "value": missing_gene_count,
            },
        ]
    )

    readiness.to_csv(
        READINESS_CSV,
        index=False,
    )

    print("\n" + "=" * 100)
    print("VISTA PRELIMINARY DECISION")
    print("=" * 100)

    print(
        "Use Xenium layers['counts'] as the "
        "spatial count matrix."
    )

    print(
        "Use snRNA layers['counts'] as the "
        "reference count matrix."
    )

    print(
        "Do not use normalized/log-transformed "
        "values as the VISTA training counts."
    )

    print("\nPreliminary readiness:")
    print(
        readiness.to_string(
            index=False
        )
    )

    print("\nCreated reports:")
    for path in [
        MATRIX_SUMMARY_CSV,
        GENE_OVERLAP_CSV,
        MISSING_GENES_CSV,
        READINESS_CSV,
        CELLTYPE_COUNTS_CSV,
        LAYER_COUNTS_CSV,
    ]:
        if path.exists():
            print(path)

    print("\nDONE: Task-02 completed.")


if __name__ == "__main__":
    main()
