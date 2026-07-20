#!/usr/bin/env python
"""
11_prepare_holdout_data_Br8667.py

Create benchmark datasets from the frozen Br8667 210/45/45 gene split.

Input Xenium
-------------
data/processed/imputation_beta/Br8667/
    spatial_data_xenium_Br8667_vista.h5ad

Frozen split files
------------------
data/processed/imputation_beta/Br8667/benchmark/splits/
    train_genes_210.csv
    validation_genes_45.csv
    test_genes_45.csv

Outputs
-------
data/processed/imputation_beta/Br8667/benchmark/inputs/
    spatial_train_210genes.h5ad

data/processed/imputation_beta/Br8667/benchmark/truth/
    validation_truth_45genes.h5ad
    test_truth_45genes.h5ad

Reports
-------
outputs/imputation_beta/Br8667/benchmark/
    11_holdout_dataset_summary.csv
    11_holdout_leakage_checks.csv

Scientific use
--------------
- spatial_train_210genes.h5ad:
    This is the only Xenium expression input used during model development.

- validation_truth_45genes.h5ad:
    Never pass this file to model training. Use it only for validation and
    model/hyperparameter selection.

- test_truth_45genes.h5ad:
    Keep untouched until the final model configuration is frozen. Use it once
    for the final blind test.

The Huuki snRNA reference remains unchanged and retains all 28,907 genes.
"""

from __future__ import annotations

import argparse
import gc
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

DEFAULT_XENIUM = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "spatial_data_xenium_Br8667_vista.h5ad"
)

DEFAULT_SPLIT_DIR = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "benchmark/splits"
)

DEFAULT_INPUT_DIR = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "benchmark/inputs"
)

DEFAULT_TRUTH_DIR = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "benchmark/truth"
)

DEFAULT_REPORT_DIR = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667"
    / "benchmark"
)

EXPECTED_TOTAL_GENES = 300
EXPECTED_TRAIN_GENES = 210
EXPECTED_VALIDATION_GENES = 45
EXPECTED_TEST_GENES = 45
EXPECTED_CELLS = 66164
SPLIT_SEED = 8667


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create Br8667 training, validation-truth, "
            "and test-truth h5ad files from the frozen gene split."
        )
    )

    parser.add_argument(
        "--xenium",
        type=Path,
        default=DEFAULT_XENIUM,
        help="Original 300-gene Br8667 Xenium h5ad.",
    )

    parser.add_argument(
        "--split-dir",
        type=Path,
        default=DEFAULT_SPLIT_DIR,
        help="Directory containing the frozen split CSV files.",
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Output directory for model-input h5ad files.",
    )

    parser.add_argument(
        "--truth-dir",
        type=Path,
        default=DEFAULT_TRUTH_DIR,
        help="Output directory for validation/test truth h5ad files.",
    )

    parser.add_argument(
        "--report-dir",
        type=Path,
        default=DEFAULT_REPORT_DIR,
        help="Output directory for validation reports.",
    )

    parser.add_argument(
        "--compression",
        choices=["gzip", "none"],
        default="gzip",
        help="h5ad compression. Default: gzip.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Replace existing derived h5ad files. "
            "The frozen split CSV files are never modified."
        ),
    )

    return parser.parse_args()


def read_gene_list(
    path: Path,
    expected_count: int,
    split_name: str,
) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {split_name} split file:\n{path}"
        )

    table = pd.read_csv(path)

    if "gene" not in table.columns:
        raise KeyError(
            f"{path} does not contain a 'gene' column."
        )

    genes = (
        table["gene"]
        .astype(str)
        .str.strip()
        .tolist()
    )

    if len(genes) != expected_count:
        raise ValueError(
            f"{split_name} contains {len(genes)} genes; "
            f"expected {expected_count}."
        )

    if len(set(genes)) != len(genes):
        duplicated = (
            pd.Series(genes)[
                pd.Series(genes).duplicated(keep=False)
            ]
            .tolist()
        )

        raise ValueError(
            f"{split_name} contains duplicated genes: "
            f"{duplicated[:20]}"
        )

    return genes


def validate_raw_counts(
    matrix,
    dataset_name: str,
) -> None:
    if sparse.issparse(matrix):
        values = np.asarray(matrix.data)
    else:
        values = np.asarray(matrix).reshape(-1)

    if values.size == 0:
        raise ValueError(
            f"{dataset_name} expression matrix is empty."
        )

    if not np.isfinite(values).all():
        raise ValueError(
            f"{dataset_name} contains NaN or infinite values."
        )

    if np.any(values < 0):
        raise ValueError(
            f"{dataset_name} contains negative values."
        )

    if not np.allclose(
        values,
        np.round(values),
        atol=1e-6,
    ):
        raise ValueError(
            f"{dataset_name} does not look like raw integer counts."
        )


def write_h5ad(
    adata: ad.AnnData,
    path: Path,
    compression: str,
) -> None:
    if compression == "gzip":
        adata.write_h5ad(
            path,
            compression="gzip",
        )
    else:
        adata.write_h5ad(path)


def add_benchmark_metadata(
    adata: ad.AnnData,
    *,
    role: str,
    split_name: str,
    source_path: Path,
    split_file: Path,
) -> None:
    adata.uns["benchmark_role"] = role
    adata.uns["benchmark_split_name"] = split_name
    adata.uns["benchmark_split_seed"] = int(SPLIT_SEED)
    adata.uns["benchmark_source_xenium"] = str(source_path)
    adata.uns["benchmark_gene_split_file"] = str(split_file)
    adata.uns["benchmark_gene_count"] = int(adata.n_vars)

    if split_name == "train":
        adata.uns["benchmark_usage"] = (
            "Allowed as Xenium expression input during model development."
        )
    elif split_name == "validation":
        adata.uns["benchmark_usage"] = (
            "Ground truth only. Do not pass to model training. "
            "Use for model and hyperparameter selection."
        )
    elif split_name == "test":
        adata.uns["benchmark_usage"] = (
            "Blind ground truth only. Keep untouched until the final "
            "configuration is frozen."
        )


def create_subset(
    xenium: ad.AnnData,
    genes: list[str],
    *,
    role: str,
    split_name: str,
    source_path: Path,
    split_file: Path,
) -> ad.AnnData:
    """
    Create a gene subset while preserving original Xenium gene order.
    """
    gene_set = set(genes)

    ordered_genes = [
        gene
        for gene in xenium.var_names.astype(str)
        if gene in gene_set
    ]

    if len(ordered_genes) != len(genes):
        missing = sorted(
            gene_set
            - set(xenium.var_names.astype(str))
        )

        raise ValueError(
            f"{split_name} genes missing from Xenium: "
            f"{missing[:20]}"
        )

    subset = xenium[:, ordered_genes].copy()

    subset.obs.index.name = None
    subset.var.index.name = None

    add_benchmark_metadata(
        subset,
        role=role,
        split_name=split_name,
        source_path=source_path,
        split_file=split_file,
    )

    validate_raw_counts(
        subset.X,
        f"{split_name} subset",
    )

    return subset


def validate_split_sets(
    original_genes: set[str],
    train_genes: list[str],
    validation_genes: list[str],
    test_genes: list[str],
) -> dict[str, bool]:
    train_set = set(train_genes)
    validation_set = set(validation_genes)
    test_set = set(test_genes)

    checks = {
        "original_gene_count_is_300":
            len(original_genes) == EXPECTED_TOTAL_GENES,
        "train_gene_count_is_210":
            len(train_set) == EXPECTED_TRAIN_GENES,
        "validation_gene_count_is_45":
            len(validation_set) == EXPECTED_VALIDATION_GENES,
        "test_gene_count_is_45":
            len(test_set) == EXPECTED_TEST_GENES,
        "train_validation_disjoint":
            train_set.isdisjoint(validation_set),
        "train_test_disjoint":
            train_set.isdisjoint(test_set),
        "validation_test_disjoint":
            validation_set.isdisjoint(test_set),
        "all_300_original_genes_assigned":
            (
                train_set
                | validation_set
                | test_set
            )
            == original_genes,
    }

    return checks


def validate_saved_file(
    path: Path,
    *,
    expected_cells: int,
    expected_genes: int,
    expected_gene_set: set[str],
    expected_obs_names: pd.Index,
) -> dict[str, object]:
    saved = ad.read_h5ad(
        path,
        backed="r",
    )

    checks = {
        "file_exists": path.exists(),
        "cell_count_correct":
            saved.n_obs == expected_cells,
        "gene_count_correct":
            saved.n_vars == expected_genes,
        "gene_set_correct":
            set(saved.var_names.astype(str))
            == expected_gene_set,
        "cell_order_matches_original":
            saved.obs_names.equals(expected_obs_names),
        "cell_ids_unique":
            saved.obs_names.is_unique,
        "gene_names_unique":
            saved.var_names.is_unique,
        "has_spatial":
            "spatial" in saved.obsm,
        "spatial_shape_correct":
            (
                "spatial" in saved.obsm
                and saved.obsm["spatial"].shape
                == (expected_cells, 2)
            ),
    }

    saved.file.close()

    return checks


def refuse_unintended_overwrite(
    paths: list[Path],
    overwrite: bool,
) -> None:
    existing = [
        path
        for path in paths
        if path.exists()
    ]

    if existing and not overwrite:
        raise FileExistsError(
            "Derived benchmark files already exist and were not overwritten:\n"
            + "\n".join(str(path) for path in existing)
            + "\nUse --overwrite only when you intentionally "
            "want to rebuild them from the same frozen split."
        )


def main() -> None:
    args = parse_args()

    train_split_path = (
        args.split_dir
        / "train_genes_210.csv"
    )

    validation_split_path = (
        args.split_dir
        / "validation_genes_45.csv"
    )

    test_split_path = (
        args.split_dir
        / "test_genes_45.csv"
    )

    train_output = (
        args.input_dir
        / "spatial_train_210genes.h5ad"
    )

    validation_output = (
        args.truth_dir
        / "validation_truth_45genes.h5ad"
    )

    test_output = (
        args.truth_dir
        / "test_truth_45genes.h5ad"
    )

    output_paths = [
        train_output,
        validation_output,
        test_output,
    ]

    refuse_unintended_overwrite(
        output_paths,
        args.overwrite,
    )

    args.input_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.truth_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.report_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 100)
    print(
        "Day 2: Prepare Br8667 train, validation-truth, "
        "and test-truth h5ad files"
    )
    print("=" * 100)

    print("\nOriginal Xenium:")
    print(args.xenium)

    if not args.xenium.exists():
        raise FileNotFoundError(args.xenium)

    train_genes = read_gene_list(
        train_split_path,
        EXPECTED_TRAIN_GENES,
        "train",
    )

    validation_genes = read_gene_list(
        validation_split_path,
        EXPECTED_VALIDATION_GENES,
        "validation",
    )

    test_genes = read_gene_list(
        test_split_path,
        EXPECTED_TEST_GENES,
        "test",
    )

    xenium = ad.read_h5ad(
        args.xenium
    )

    print("\nLoaded Xenium:")
    print(xenium)

    if xenium.n_obs != EXPECTED_CELLS:
        raise ValueError(
            f"Expected {EXPECTED_CELLS} cells, "
            f"found {xenium.n_obs}."
        )

    if xenium.n_vars != EXPECTED_TOTAL_GENES:
        raise ValueError(
            f"Expected {EXPECTED_TOTAL_GENES} genes, "
            f"found {xenium.n_vars}."
        )

    if "spatial" not in xenium.obsm:
        raise KeyError(
            "Original Xenium does not contain obsm['spatial']."
        )

    original_gene_set = set(
        xenium.var_names.astype(str)
    )

    split_checks = validate_split_sets(
        original_gene_set,
        train_genes,
        validation_genes,
        test_genes,
    )

    print("\nFrozen split validation:")
    for name, value in split_checks.items():
        print(f"{name}: {value}")

    if not all(split_checks.values()):
        failed = [
            name
            for name, value in split_checks.items()
            if not value
        ]

        raise ValueError(
            f"Frozen split validation failed: {failed}"
        )

    dataset_specs = [
        {
            "split_name": "train",
            "role": "spatial_training_input",
            "genes": train_genes,
            "split_file": train_split_path,
            "output": train_output,
            "expected_genes": EXPECTED_TRAIN_GENES,
        },
        {
            "split_name": "validation",
            "role": "validation_ground_truth",
            "genes": validation_genes,
            "split_file": validation_split_path,
            "output": validation_output,
            "expected_genes": EXPECTED_VALIDATION_GENES,
        },
        {
            "split_name": "test",
            "role": "blind_test_ground_truth",
            "genes": test_genes,
            "split_file": test_split_path,
            "output": test_output,
            "expected_genes": EXPECTED_TEST_GENES,
        },
    ]

    summary_rows: list[dict] = []
    saved_validation_rows: list[dict] = []

    for spec in dataset_specs:
        print(
            f"\nCreating {spec['split_name']} dataset..."
        )

        subset = create_subset(
            xenium=xenium,
            genes=spec["genes"],
            role=spec["role"],
            split_name=spec["split_name"],
            source_path=args.xenium,
            split_file=spec["split_file"],
        )

        print(subset)

        write_h5ad(
            subset,
            spec["output"],
            args.compression,
        )

        size_mb = (
            spec["output"].stat().st_size
            / (1024**2)
        )

        summary_rows.append(
            {
                "dataset": spec["split_name"],
                "role": spec["role"],
                "path": str(spec["output"]),
                "n_cells": subset.n_obs,
                "n_genes": subset.n_vars,
                "file_size_mb": size_mb,
                "has_spatial":
                    "spatial" in subset.obsm,
                "raw_counts_integer_nonnegative":
                    True,
            }
        )

        validation = validate_saved_file(
            spec["output"],
            expected_cells=EXPECTED_CELLS,
            expected_genes=spec["expected_genes"],
            expected_gene_set=set(spec["genes"]),
            expected_obs_names=xenium.obs_names,
        )

        for check_name, check_value in validation.items():
            saved_validation_rows.append(
                {
                    "dataset": spec["split_name"],
                    "check": check_name,
                    "value": check_value,
                }
            )

        if not all(validation.values()):
            failed = [
                name
                for name, value in validation.items()
                if not value
            ]

            raise ValueError(
                f"Saved {spec['split_name']} file failed "
                f"validation: {failed}"
            )

        print(
            f"Created: {spec['output']}\n"
            f"Size: {size_mb:.2f} MB"
        )

        del subset
        gc.collect()

    leakage_checks = {
        "validation_absent_from_train":
            set(validation_genes).isdisjoint(
                train_genes
            ),
        "test_absent_from_train":
            set(test_genes).isdisjoint(
                train_genes
            ),
        "test_absent_from_validation":
            set(test_genes).isdisjoint(
                validation_genes
            ),
        "train_plus_validation_plus_test_is_300":
            len(
                set(train_genes)
                | set(validation_genes)
                | set(test_genes)
            )
            == EXPECTED_TOTAL_GENES,
    }

    for name, value in leakage_checks.items():
        saved_validation_rows.append(
            {
                "dataset": "cross_dataset",
                "check": name,
                "value": value,
            }
        )

    if not all(leakage_checks.values()):
        failed = [
            name
            for name, value in leakage_checks.items()
            if not value
        ]

        raise ValueError(
            f"Gene leakage detected: {failed}"
        )

    summary_path = (
        args.report_dir
        / "11_holdout_dataset_summary.csv"
    )

    validation_path = (
        args.report_dir
        / "11_holdout_leakage_checks.csv"
    )

    pd.DataFrame(
        summary_rows
    ).to_csv(
        summary_path,
        index=False,
    )

    pd.DataFrame(
        saved_validation_rows
    ).to_csv(
        validation_path,
        index=False,
    )

    print("\n" + "=" * 100)
    print(
        "PASS: Br8667 holdout datasets created "
        "without gene leakage"
    )
    print("=" * 100)

    print("\nCreated model input:")
    print(train_output)

    print("\nCreated ground-truth files:")
    print(validation_output)
    print(test_output)

    print("\nCreated reports:")
    print(summary_path)
    print(validation_path)

    print(
        "\nUsage rule:\n"
        "  Train with spatial_train_210genes.h5ad only.\n"
        "  Use validation_truth_45genes.h5ad only for development evaluation.\n"
        "  Do not inspect test_truth_45genes.h5ad performance until the "
        "final configuration is frozen."
    )


if __name__ == "__main__":
    main()
