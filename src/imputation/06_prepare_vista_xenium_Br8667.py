#!/usr/bin/env python

"""
06_prepare_vista_xenium_Br8667.py

Prepare final VISTA inputs for Br8667.

Reference:
    Huuki snRNA-seq Br8667

Spatial target:
    Xenium Br8667

Inputs:
    data/processed/imputation_beta/Br8667/
        seq_data_huuki_snrna_Br8667_shared.h5ad
        xenium_Br8667_annotated.h5ad

Outputs:
    data/processed/imputation_beta/Br8667/
        seq_data_huuki_snrna_Br8667_vista.h5ad
        spatial_data_xenium_Br8667_vista.h5ad

    outputs/imputation_beta/Br8667/
        06_vista_preparation_summary.csv
        06_vista_celltype_counts.csv
        06_vista_gene_overlap.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

INPUT_DIR = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
)

SNRNA_INPUT = (
    INPUT_DIR
    / "seq_data_huuki_snrna_Br8667_shared.h5ad"
)

XENIUM_INPUT = (
    INPUT_DIR
    / "xenium_Br8667_annotated.h5ad"
)

SNRNA_OUTPUT = (
    INPUT_DIR
    / "seq_data_huuki_snrna_Br8667_vista.h5ad"
)

XENIUM_OUTPUT = (
    INPUT_DIR
    / "spatial_data_xenium_Br8667_vista.h5ad"
)

REPORT_DIR = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667"
)

SUMMARY_OUTPUT = (
    REPORT_DIR
    / "06_vista_preparation_summary.csv"
)

CELLTYPE_OUTPUT = (
    REPORT_DIR
    / "06_vista_celltype_counts.csv"
)

GENE_OVERLAP_OUTPUT = (
    REPORT_DIR
    / "06_vista_gene_overlap.csv"
)


SNRNA_LABEL_CANDIDATES = [
    "cellType_broad_k",
    "cellType_broad_hc",
    "cellType_k",
    "cellType_hc",
    "cellType_layer",
    "layer_annotation",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare Br8667 snRNA and Xenium data for VISTA."
    )

    parser.add_argument(
        "--snrna-label-column",
        default="auto",
        help=(
            "snRNA obs column used to create scClassify. "
            "Default: auto."
        ),
    )

    parser.add_argument(
        "--compression",
        choices=["gzip", "none"],
        default="gzip",
        help="Output h5ad compression. Default: gzip.",
    )

    return parser.parse_args()


def print_created_file(label: str, path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Expected output was not created: {path}"
        )

    size_mb = path.stat().st_size / (1024**2)

    print(f"\nCreated {label}:")
    print(path)
    print(f"Size: {size_mb:.2f} MB")


def choose_snrna_label_column(
    adata: ad.AnnData,
    requested_column: str,
) -> str:
    if requested_column != "auto":
        if requested_column not in adata.obs.columns:
            raise KeyError(
                f"Requested snRNA label column "
                f"'{requested_column}' was not found.\n"
                f"Available obs columns:\n"
                f"{adata.obs.columns.tolist()}"
            )

        return requested_column

    for column in SNRNA_LABEL_CANDIDATES:
        if column in adata.obs.columns:
            return column

    raise KeyError(
        "Could not automatically identify an snRNA cell-type column.\n"
        f"Checked: {SNRNA_LABEL_CANDIDATES}\n"
        f"Available obs columns:\n"
        f"{adata.obs.columns.tolist()}"
    )


def harmonize_cell_type(value) -> str:
    """
    Convert detailed snRNA/Xenium annotations into comparable
    broad cortical cell-type classes.
    """

    if pd.isna(value):
        return "Other"

    original = str(value).strip()
    label = original.lower()

    if label in {"", "nan", "none", "na", "unknown"}:
        return "Other"

    # Ambiguous or low-confidence annotations
    if label in {
        "microoligo",
        "micro oligo",
        "drop",
        "ambig/oligo",
        "ambig/in/endo",
    }:
        return "Other"

    # Excitatory neurons
    excitatory_terms = [
        "excit",
        "glut",
        "pyr",
        "l2/3 ex",
        "l3 ex",
        "l4 ex",
        "l4/5 ex",
        "l5 ex",
        "l6 ex",
    ]

    if any(term in label for term in excitatory_terms):
        return "Excit"

    # Inhibitory neurons
    inhibitory_terms = [
        "inhib",
        "gaba",
        "interneuron",
        "mge",
        "cge",
        "sst",
        "pvalb",
        "pv ",
        "vip",
        "lamp5",
    ]

    if any(term in label for term in inhibitory_terms):
        return "Inhib"

    if "ambig/in/endo" in label:
        return "Other"

    # OPC must be checked before Oligo
    if (
        "opc" in label
        or "oligodendrocyte precursor" in label
        or "precursor" in label
    ):
        return "OPC"

    # Mature oligodendrocytes
    if "oligo" in label:
        return "Oligo"

    # Astrocytes
    if (
        "astro" in label
        or label == "ast"
        or label.startswith("ast ")
    ):
        return "Astro"

    # Microglia / immune
    if (
        "micro" in label
        or label == "mic"
        or label.startswith("mic ")
        or "immune" in label
        or "macrophage" in label
    ):
        return "Micro"

    # Endothelial / vascular / mural cells
    vascular_terms = [
        "endo",
        "vascular",
        "pericyte",
        "mural",
        "smooth muscle",
        "vlmc",
    ]

    if any(term in label for term in vascular_terms):
        return "EndoMural"

    return "Other"


def get_gene_symbols(
    adata: ad.AnnData,
    dataset_name: str,
) -> tuple[pd.Series, str]:
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
                "as gene symbols."
            )

            return symbols, column

    print(
        f"{dataset_name}: using var_names as gene symbols."
    )

    return (
        pd.Series(
            adata.var_names.astype(str),
            index=adata.var_names,
        ),
        "var_names",
    )


def standardize_genes(
    adata: ad.AnnData,
    dataset_name: str,
) -> tuple[ad.AnnData, str, int]:
    symbols, source_column = get_gene_symbols(
        adata,
        dataset_name,
    )

    symbols = (
        symbols
        .astype(str)
        .str.strip()
        .reset_index(drop=True)
    )

    gene_table = pd.DataFrame(
        {
            "gene_symbol": symbols,
            "position": np.arange(adata.n_vars),
        }
    )

    invalid_values = {
        "",
        "nan",
        "NaN",
        "None",
        "NA",
    }

    valid = ~gene_table["gene_symbol"].isin(
        invalid_values
    )

    gene_table = gene_table.loc[valid].copy()

    duplicate_count = int(
        gene_table["gene_symbol"].duplicated().sum()
    )

    gene_table = (
        gene_table
        .drop_duplicates(
            subset="gene_symbol",
            keep="first",
        )
        .copy()
    )

    positions = gene_table["position"].to_numpy()
    final_symbols = gene_table["gene_symbol"].tolist()

    standardized = adata[:, positions].copy()

    standardized.var["original_var_name"] = (
        standardized.var_names.astype(str)
    )

    standardized.var["vista_gene_symbol"] = (
        final_symbols
    )

    standardized.var_names = pd.Index(
        final_symbols,
        dtype="str",
    )

    standardized.var.index.name = None
    standardized.obs.index.name = None

    if not standardized.var_names.is_unique:
        raise ValueError(
            f"{dataset_name} gene names are not unique "
            "after standardization."
        )

    print(f"\n{dataset_name} gene standardization:")
    print(f"  Input genes: {adata.n_vars:,}")
    print(f"  Final genes: {standardized.n_vars:,}")
    print(f"  Duplicates removed: {duplicate_count:,}")
    print(f"  Gene source: {source_column}")

    return (
        standardized,
        source_column,
        duplicate_count,
    )


def use_raw_counts_as_x(
    adata: ad.AnnData,
    dataset_name: str,
) -> ad.AnnData:
    if "counts" not in adata.layers:
        raise KeyError(
            f"{dataset_name} does not contain layers['counts']."
        )

    counts = adata.layers["counts"]

    if sparse.issparse(counts):
        counts_x = counts.copy()
    else:
        counts_x = np.asarray(counts).copy()

    adata.X = counts_x

    # Final VISTA files only need raw counts in X.
    # Removing layers prevents duplication of the large count matrix.
    for layer_name in list(adata.layers.keys()):
        del adata.layers[layer_name]

    adata.uns["vista_expression_matrix"] = (
        "X contains raw counts copied from input layers['counts']"
    )

    print(
        f"{dataset_name}: copied layers['counts'] into X "
        "and removed duplicated layers."
    )

    return adata


def matrix_is_integer_nonnegative(matrix) -> bool:
    if sparse.issparse(matrix):
        values = np.asarray(matrix.data)
    else:
        values = np.asarray(matrix).reshape(-1)

    if values.size == 0:
        return True

    if values.size > 1_000_000:
        indices = np.linspace(
            0,
            values.size - 1,
            1_000_000,
            dtype=int,
        )
        values = values[indices]

    nonnegative = bool(np.all(values >= 0))

    integer_like = bool(
        np.allclose(
            values,
            np.round(values),
            atol=1e-6,
        )
    )

    return nonnegative and integer_like


def write_h5ad(
    adata: ad.AnnData,
    path: Path,
    compression: str,
) -> None:
    if path.exists():
        print("\nRemoving old output:")
        print(path)
        path.unlink()

    if compression == "gzip":
        adata.write_h5ad(
            path,
            compression="gzip",
        )
    else:
        adata.write_h5ad(path)


def main() -> None:
    args = parse_args()

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    INPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 100)
    print("Task-03: Prepare final Br8667 VISTA inputs")
    print("=" * 100)

    print("\nsnRNA input:")
    print(SNRNA_INPUT)

    print("\nXenium input:")
    print(XENIUM_INPUT)

    for path in [SNRNA_INPUT, XENIUM_INPUT]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing required input:\n{path}"
            )

    print("\nLoading snRNA Br8667...")
    snrna = ad.read_h5ad(SNRNA_INPUT)

    print("Loading Xenium Br8667...")
    xenium = ad.read_h5ad(XENIUM_INPUT)

    print("\nOriginal snRNA:")
    print(snrna)

    print("\nOriginal Xenium:")
    print(xenium)

    # ------------------------------------------------------------------
    # Cell-type labels
    # ------------------------------------------------------------------

    snrna_label_column = choose_snrna_label_column(
        snrna,
        args.snrna_label_column,
    )

    print("\nsnRNA label source:")
    print(snrna_label_column)

    snrna.obs["scClassify_original"] = (
        snrna.obs[snrna_label_column]
        .astype(str)
    )

    snrna.obs["scClassify"] = (
        snrna.obs[snrna_label_column]
        .map(harmonize_cell_type)
        .astype("category")
    )

    xenium_label_column = "cell_type_annotation"

    if xenium_label_column not in xenium.obs.columns:
        if "scClassify" in xenium.obs.columns:
            xenium_label_column = "scClassify"
        else:
            raise KeyError(
                "Xenium has neither "
                "obs['cell_type_annotation'] nor "
                "obs['scClassify']."
            )

    xenium.obs["scClassify_original"] = (
        xenium.obs[xenium_label_column]
        .astype(str)
    )

    xenium.obs["scClassify"] = (
        xenium.obs[xenium_label_column]
        .map(harmonize_cell_type)
        .astype("category")
    )

    print("\nHarmonized snRNA cell types:")
    print(
        snrna.obs["scClassify"]
        .value_counts(dropna=False)
        .to_string()
    )

    print("\nHarmonized Xenium cell types:")
    print(
        xenium.obs["scClassify"]
        .value_counts(dropna=False)
        .to_string()
    )

    snrna_types = set(
        snrna.obs["scClassify"].astype(str)
    )

    xenium_types = set(
        xenium.obs["scClassify"].astype(str)
    )

    common_cell_types = sorted(
        snrna_types & xenium_types
    )

    print("\nCommon harmonized cell types:")
    print(common_cell_types)

    if len(common_cell_types) < 2:
        raise ValueError(
            "Fewer than two common scClassify labels were "
            "found between snRNA and Xenium."
        )

    # ------------------------------------------------------------------
    # Required metadata
    # ------------------------------------------------------------------

    snrna.obs["names"] = (
        snrna.obs_names.astype(str)
    )

    xenium.obs["names"] = (
        xenium.obs_names.astype(str)
    )

    snrna.obs["batch"] = "Br8667_snrna"
    xenium.obs["batch"] = "Br8667_xenium"

    snrna.obs_names_make_unique()
    xenium.obs_names_make_unique()

    if "spatial" not in xenium.obsm:
        raise KeyError(
            "Xenium does not contain obsm['spatial']."
        )

    xenium_spatial = np.asarray(
        xenium.obsm["spatial"]
    )

    if (
        xenium_spatial.ndim != 2
        or xenium_spatial.shape[1] < 2
    ):
        raise ValueError(
            "Xenium obsm['spatial'] must be an "
            "n_cells x 2 matrix."
        )

    xenium.obsm["spatial"] = (
        xenium_spatial[:, :2]
        .astype(np.float32, copy=False)
    )

    # ------------------------------------------------------------------
    # Standardize genes
    # ------------------------------------------------------------------

    snrna, snrna_gene_source, snrna_duplicates = (
        standardize_genes(
            snrna,
            "snRNA",
        )
    )

    xenium, xenium_gene_source, xenium_duplicates = (
        standardize_genes(
            xenium,
            "Xenium",
        )
    )

    snrna_gene_set = set(
        snrna.var_names.astype(str)
    )

    gene_overlap = pd.DataFrame(
        {
            "xenium_gene": (
                xenium.var_names.astype(str)
            )
        }
    )

    gene_overlap["present_in_snrna"] = (
        gene_overlap["xenium_gene"]
        .isin(snrna_gene_set)
    )

    gene_overlap.to_csv(
        GENE_OVERLAP_OUTPUT,
        index=False,
    )

    missing_genes = gene_overlap.loc[
        ~gene_overlap["present_in_snrna"],
        "xenium_gene",
    ].tolist()

    print("\nGene overlap:")
    print(f"snRNA genes: {snrna.n_vars:,}")
    print(f"Xenium genes: {xenium.n_vars:,}")
    print(
        "Xenium genes present in snRNA:",
        int(gene_overlap["present_in_snrna"].sum()),
    )
    print(
        "Xenium genes missing from snRNA:",
        len(missing_genes),
    )

    if missing_genes:
        raise ValueError(
            "The following Xenium genes are missing "
            f"from snRNA:\n{missing_genes}"
        )

    # ------------------------------------------------------------------
    # Put raw counts into X
    # ------------------------------------------------------------------

    snrna = use_raw_counts_as_x(
        snrna,
        "snRNA",
    )

    xenium = use_raw_counts_as_x(
        xenium,
        "Xenium",
    )

    snrna_counts_valid = matrix_is_integer_nonnegative(
        snrna.X
    )

    xenium_counts_valid = matrix_is_integer_nonnegative(
        xenium.X
    )

    print("\nCount validation:")
    print(
        "snRNA X integer and nonnegative:",
        snrna_counts_valid,
    )
    print(
        "Xenium X integer and nonnegative:",
        xenium_counts_valid,
    )

    if not snrna_counts_valid:
        raise ValueError(
            "snRNA X does not look like raw count data."
        )

    if not xenium_counts_valid:
        raise ValueError(
            "Xenium X does not look like raw count data."
        )

    # ------------------------------------------------------------------
    # Final metadata cleanup
    # ------------------------------------------------------------------

    for dataset in [snrna, xenium]:
        dataset.obs["names"] = (
            dataset.obs_names.astype(str)
        )

        dataset.obs["names"] = np.asarray(
            dataset.obs_names.astype(str),
            dtype=object,
        )

        dataset.obs["batch"] = (
            dataset.obs["batch"].astype("category")
        )

        dataset.obs["scClassify"] = (
            dataset.obs["scClassify"].astype("category")
        )

        dataset.obs.index.name = None
        dataset.var.index.name = None

    if not snrna.obs["names"].is_unique:
        raise ValueError(
            "snRNA obs['names'] is not unique."
        )

    if not xenium.obs["names"].is_unique:
        raise ValueError(
            "Xenium obs['names'] is not unique."
        )

    if not set(xenium.var_names) <= set(snrna.var_names):
        raise ValueError(
            "Final Xenium genes are not a subset "
            "of final snRNA genes."
        )

    # ------------------------------------------------------------------
    # Reports
    # ------------------------------------------------------------------

    snrna_celltypes = (
        snrna.obs["scClassify"]
        .value_counts(dropna=False)
        .rename_axis("scClassify")
        .reset_index(name="n_observations")
    )

    snrna_celltypes.insert(
        0,
        "dataset",
        "snRNA",
    )

    xenium_celltypes = (
        xenium.obs["scClassify"]
        .value_counts(dropna=False)
        .rename_axis("scClassify")
        .reset_index(name="n_observations")
    )

    xenium_celltypes.insert(
        0,
        "dataset",
        "Xenium",
    )

    celltype_report = pd.concat(
        [
            snrna_celltypes,
            xenium_celltypes,
        ],
        ignore_index=True,
    )

    celltype_report.to_csv(
        CELLTYPE_OUTPUT,
        index=False,
    )

    summary = pd.DataFrame(
        [
            {
                "metric": "snrna_input",
                "value": str(SNRNA_INPUT),
            },
            {
                "metric": "xenium_input",
                "value": str(XENIUM_INPUT),
            },
            {
                "metric": "snrna_output",
                "value": str(SNRNA_OUTPUT),
            },
            {
                "metric": "xenium_output",
                "value": str(XENIUM_OUTPUT),
            },
            {
                "metric": "snrna_cells",
                "value": snrna.n_obs,
            },
            {
                "metric": "xenium_cells",
                "value": xenium.n_obs,
            },
            {
                "metric": "snrna_genes",
                "value": snrna.n_vars,
            },
            {
                "metric": "xenium_genes",
                "value": xenium.n_vars,
            },
            {
                "metric": "xenium_genes_in_snrna",
                "value": int(
                    gene_overlap[
                        "present_in_snrna"
                    ].sum()
                ),
            },
            {
                "metric": "missing_xenium_genes",
                "value": len(missing_genes),
            },
            {
                "metric": "snrna_label_source",
                "value": snrna_label_column,
            },
            {
                "metric": "xenium_label_source",
                "value": xenium_label_column,
            },
            {
                "metric": "common_cell_types",
                "value": ",".join(common_cell_types),
            },
            {
                "metric": "snrna_gene_source",
                "value": snrna_gene_source,
            },
            {
                "metric": "xenium_gene_source",
                "value": xenium_gene_source,
            },
            {
                "metric": "snrna_duplicate_genes_removed",
                "value": snrna_duplicates,
            },
            {
                "metric": "xenium_duplicate_genes_removed",
                "value": xenium_duplicates,
            },
            {
                "metric": "snrna_X_is_raw_counts",
                "value": snrna_counts_valid,
            },
            {
                "metric": "xenium_X_is_raw_counts",
                "value": xenium_counts_valid,
            },
            {
                "metric": "xenium_has_spatial",
                "value": "spatial" in xenium.obsm,
            },
            {
                "metric": "xenium_spatial_shape",
                "value": str(
                    xenium.obsm["spatial"].shape
                ),
            },
        ]
    )

    summary.to_csv(
        SUMMARY_OUTPUT,
        index=False,
    )

    # ------------------------------------------------------------------
    # Save final VISTA inputs
    # ------------------------------------------------------------------

    print("\nFinal snRNA VISTA input:")
    print(snrna)

    print("\nFinal Xenium VISTA input:")
    print(xenium)

    print("\nWriting final snRNA VISTA input...")
    write_h5ad(
        snrna,
        SNRNA_OUTPUT,
        args.compression,
    )

    print_created_file(
        "snRNA VISTA input",
        SNRNA_OUTPUT,
    )

    print("\nWriting final Xenium VISTA input...")
    write_h5ad(
        xenium,
        XENIUM_OUTPUT,
        args.compression,
    )

    print_created_file(
        "Xenium VISTA input",
        XENIUM_OUTPUT,
    )

    print_created_file(
        "preparation summary",
        SUMMARY_OUTPUT,
    )

    print_created_file(
        "cell-type report",
        CELLTYPE_OUTPUT,
    )

    print_created_file(
        "gene-overlap report",
        GENE_OVERLAP_OUTPUT,
    )

    # ------------------------------------------------------------------
    # Reload validation
    # ------------------------------------------------------------------

    print("\nReloading final files for validation...")

    snrna_check = ad.read_h5ad(
        SNRNA_OUTPUT
    )

    xenium_check = ad.read_h5ad(
        XENIUM_OUTPUT
    )

    final_checks = {
        "snrna_has_scClassify":
            "scClassify" in snrna_check.obs,
        "xenium_has_scClassify":
            "scClassify" in xenium_check.obs,
        "snrna_has_names":
            "names" in snrna_check.obs,
        "xenium_has_names":
            "names" in xenium_check.obs,
        "xenium_has_spatial":
            "spatial" in xenium_check.obsm,
        "xenium_genes_subset_snrna":
            set(xenium_check.var_names)
            <= set(snrna_check.var_names),
        "snrna_X_raw_counts":
            matrix_is_integer_nonnegative(
                snrna_check.X
            ),
        "xenium_X_raw_counts":
            matrix_is_integer_nonnegative(
                xenium_check.X
            ),
    }

    print("\nFinal validation:")
    for check, value in final_checks.items():
        print(f"{check}: {value}")

    if not all(final_checks.values()):
        failed = [
            check
            for check, value
            in final_checks.items()
            if not value
        ]

        raise ValueError(
            f"Final VISTA validation failed: {failed}"
        )

    print("\n" + "=" * 100)
    print("PASS: Br8667 VISTA inputs are ready")
    print("=" * 100)

    print("\nReference input:")
    print(SNRNA_OUTPUT)

    print("\nSpatial target input:")
    print(XENIUM_OUTPUT)


if __name__ == "__main__":
    main()
