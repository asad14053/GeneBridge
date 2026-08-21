#!/usr/bin/env python3

"""
Experiment 5.1 preflight validation.

Experiment 5.1:
    Reference:
        All 10 Huuki snRNA-seq brains, including Br8667.

    Target:
        Br8667 Xenium QC dataset containing the full 300-gene panel.

    Cross-validation:
        Reuse the existing three gene-wise folds:
            200 observed genes
            100 held-out genes

This script is read-only with respect to the source H5AD files.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

REFERENCE_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

TARGET_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "imputation_beta"
    / "Br8667"
    / "spatial_data_xenium_Br8667_vista_qc.h5ad"
)

FOLD_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "imputation_beta"
    / "Br8667"
    / "gene_folds_200_100"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
    / "experiment_5_1"
    / "input_validation"
)

MANIFEST_FILE = (
    OUTPUT_DIR
    / "experiment_5_1_manifest.json"
)

FOLD_VALIDATION_FILE = (
    OUTPUT_DIR
    / "fold_validation.csv"
)

COHORT_SUMMARY_FILE = (
    OUTPUT_DIR
    / "reference_cohort_summary.csv"
)

PREFLIGHT_REPORT_FILE = (
    OUTPUT_DIR
    / "experiment_5_1_preflight_report.txt"
)


def section(title: str) -> None:
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(
            f"Required file was not found:\n{path}"
        )


def require_directory(path: Path) -> None:
    if not path.is_dir():
        raise FileNotFoundError(
            f"Required directory was not found:\n{path}"
        )


def clean_gene_list(values: Any) -> list[str]:
    genes = [
        str(value).strip()
        for value in values
        if pd.notna(value) and str(value).strip()
    ]

    return genes


def read_gene_file(path: Path) -> list[str]:
    suffix = path.suffix.lower()

    if suffix == ".h5ad":
        adata = ad.read_h5ad(
            path,
            backed="r",
        )

        try:
            return clean_gene_list(
                adata.var_names.tolist()
            )
        finally:
            adata.file.close()

    if suffix in {".txt", ".tsv"}:
        separator = "\t" if suffix == ".tsv" else None

        frame = pd.read_csv(
            path,
            sep=separator,
            header=None,
            comment="#",
            engine="python",
        )

        return clean_gene_list(
            frame.iloc[:, 0].tolist()
        )

    if suffix == ".csv":
        frame = pd.read_csv(path)

        preferred_columns = [
            "gene",
            "gene_name",
            "gene_symbol",
            "genes",
            "var_names",
        ]

        selected_column = None

        for column in preferred_columns:
            if column in frame.columns:
                selected_column = column
                break

        if selected_column is None:
            selected_column = frame.columns[0]

        return clean_gene_list(
            frame[selected_column].tolist()
        )

    if suffix == ".json":
        with path.open("r", encoding="utf-8") as handle:
            content = json.load(handle)

        if isinstance(content, list):
            return clean_gene_list(content)

        if isinstance(content, dict):
            for key in (
                "genes",
                "gene_list",
                "observed_genes",
                "heldout_genes",
                "held_out_genes",
            ):
                if key in content and isinstance(
                    content[key],
                    list,
                ):
                    return clean_gene_list(content[key])

        raise ValueError(
            f"Could not identify a gene list in JSON file:\n{path}"
        )

    raise ValueError(
        f"Unsupported gene file format:\n{path}"
    )


def classify_fold_file(path: Path) -> tuple[int, str] | None:
    name = path.name.lower()

    fold_match = re.search(
        r"fold[_-]?([123])",
        name,
    )

    if fold_match is None:
        return None

    fold = int(
        fold_match.group(1)
    )

    if any(
        token in name
        for token in (
            "heldout",
            "held_out",
            "holdout",
            "test_genes",
            "evaluation_genes",
        )
    ):
        return fold, "heldout"

    if any(
        token in name
        for token in (
            "observed",
            "training",
            "train_genes",
            "input_genes",
        )
    ):
        return fold, "observed"

    return None


def discover_fold_files(
    fold_dir: Path,
) -> dict[int, dict[str, Path]]:
    supported_suffixes = {
        ".h5ad",
        ".csv",
        ".txt",
        ".tsv",
        ".json",
    }

    fold_files: dict[int, dict[str, Path]] = {
        1: {},
        2: {},
        3: {},
    }

    for path in sorted(
        fold_dir.rglob("*")
    ):
        if (
            not path.is_file()
            or path.suffix.lower()
            not in supported_suffixes
        ):
            continue

        classification = classify_fold_file(path)

        if classification is None:
            continue

        fold, role = classification

        if role in fold_files[fold]:
            existing = fold_files[fold][role]

            raise RuntimeError(
                f"Multiple candidate files were found for "
                f"fold {fold} {role}:\n"
                f"1. {existing}\n"
                f"2. {path}"
            )

        fold_files[fold][role] = path

    missing = []

    for fold in (1, 2, 3):
        for role in ("observed", "heldout"):
            if role not in fold_files[fold]:
                missing.append(
                    f"fold {fold}: {role}"
                )

    if missing:
        available = "\n".join(
            str(path)
            for path in sorted(
                fold_dir.rglob("*")
            )
            if path.is_file()
        )

        raise RuntimeError(
            "Could not automatically locate every fold gene file.\n"
            f"Missing: {', '.join(missing)}\n\n"
            "Files currently present in the fold directory:\n"
            f"{available}"
        )

    return fold_files


def sample_x_statistics(
    adata: ad.AnnData,
) -> dict[str, Any]:
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(
            f"Invalid AnnData shape: {adata.shape}"
        )

    row_indices = sorted(
        {
            0,
            min(99, adata.n_obs - 1),
            adata.n_obs // 2,
            max(0, adata.n_obs - 100),
            adata.n_obs - 1,
        }
    )

    sampled_blocks: list[np.ndarray] = []

    for row in row_indices:
        block = adata.X[
            row : row + 1,
            :,
        ]

        if sparse.issparse(block):
            block = block.toarray()

        sampled_blocks.append(
            np.asarray(
                block,
                dtype=np.float64,
            ).ravel()
        )

    values = np.concatenate(
        sampled_blocks
    )

    finite = np.isfinite(values)

    if not finite.all():
        raise ValueError(
            "The sampled matrix contains NaN or infinite values."
        )

    nonzero = values[
        values != 0
    ]

    minimum = float(
        values.min()
    )

    maximum = float(
        values.max()
    )

    integer_like_fraction = float(
        np.isclose(
            values,
            np.round(values),
            atol=1e-6,
        ).mean()
    )

    nonnegative = bool(
        np.all(values >= 0)
    )

    if (
        nonnegative
        and integer_like_fraction >= 0.999
    ):
        classification = "likely_raw_counts"
    elif (
        nonnegative
        and integer_like_fraction < 0.99
        and maximum < 30
    ):
        classification = "likely_log_normalized"
    else:
        classification = "uncertain"

    return {
        "sampled_rows": row_indices,
        "sampled_values": int(values.size),
        "minimum": minimum,
        "maximum": maximum,
        "mean_sampled_nonzero": (
            float(nonzero.mean())
            if nonzero.size
            else 0.0
        ),
        "integer_like_fraction": integer_like_fraction,
        "nonnegative": nonnegative,
        "classification": classification,
    }


def summarize_reference_cohort(
    reference: ad.AnnData,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    obs = reference.obs.copy()

    summary_rows: list[dict[str, Any]] = []

    if "BrNum" in obs.columns:
        brain_counts = (
            obs["BrNum"]
            .astype(str)
            .value_counts(dropna=False)
            .sort_index()
        )

        for brain, count in brain_counts.items():
            summary_rows.append(
                {
                    "summary_level": "BrNum",
                    "label": brain,
                    "n_cells": int(count),
                }
            )

        unique_brains = int(
            obs["BrNum"]
            .astype(str)
            .nunique()
        )

        includes_br8667 = bool(
            (
                obs["BrNum"]
                .astype(str)
                == "Br8667"
            ).any()
        )
    else:
        unique_brains = None
        includes_br8667 = None

    if "Sample" in obs.columns:
        sample_counts = (
            obs["Sample"]
            .astype(str)
            .value_counts(dropna=False)
            .sort_index()
        )

        for sample, count in sample_counts.items():
            summary_rows.append(
                {
                    "summary_level": "Sample",
                    "label": sample,
                    "n_cells": int(count),
                }
            )

        unique_samples = int(
            obs["Sample"]
            .astype(str)
            .nunique()
        )

        sample_includes_br8667 = bool(
            obs["Sample"]
            .astype(str)
            .str.startswith("Br8667")
            .any()
        )

        if includes_br8667 is None:
            includes_br8667 = sample_includes_br8667
    else:
        unique_samples = None

    cohort_frame = pd.DataFrame(
        summary_rows,
        columns=[
            "summary_level",
            "label",
            "n_cells",
        ],
    )

    cohort_metadata = {
        "BrNum_column_present": "BrNum" in obs.columns,
        "Sample_column_present": "Sample" in obs.columns,
        "n_unique_brains": unique_brains,
        "n_unique_samples": unique_samples,
        "includes_Br8667": includes_br8667,
    }

    return cohort_frame, cohort_metadata


def main() -> None:
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    require_file(REFERENCE_FILE)
    require_file(TARGET_FILE)
    require_directory(FOLD_DIR)

    section("EXPERIMENT 5.1 PREFLIGHT")

    print("Reference:")
    print(REFERENCE_FILE)

    print("\nTarget:")
    print(TARGET_FILE)

    print("\nFold directory:")
    print(FOLD_DIR)

    reference = ad.read_h5ad(
        REFERENCE_FILE,
        backed="r",
    )

    target = ad.read_h5ad(
        TARGET_FILE,
        backed="r",
    )

    try:
        section("REFERENCE VALIDATION")

        print(
            "Reference shape [cells x genes]:",
            reference.shape,
        )
        print(
            "Reference obs columns:",
            reference.obs.columns.tolist(),
        )
        print(
            "Reference layers:",
            list(reference.layers.keys()),
        )
        print(
            "Unique cell IDs:",
            reference.obs_names.is_unique,
        )
        print(
            "Unique gene names:",
            reference.var_names.is_unique,
        )

        reference_x_audit = sample_x_statistics(
            reference
        )

        print(
            "Reference X classification:",
            reference_x_audit["classification"],
        )

        cohort_frame, cohort_metadata = (
            summarize_reference_cohort(
                reference
            )
        )

        print(
            "Unique brains:",
            cohort_metadata["n_unique_brains"],
        )
        print(
            "Unique samples:",
            cohort_metadata["n_unique_samples"],
        )
        print(
            "Includes Br8667:",
            cohort_metadata["includes_Br8667"],
        )

        section("TARGET VALIDATION")

        print(
            "Target shape [cells x genes]:",
            target.shape,
        )
        print(
            "Target layers:",
            list(target.layers.keys()),
        )
        print(
            "Target obsm keys:",
            list(target.obsm.keys()),
        )
        print(
            "Unique cell IDs:",
            target.obs_names.is_unique,
        )
        print(
            "Unique gene names:",
            target.var_names.is_unique,
        )

        target_x_audit = sample_x_statistics(
            target
        )

        print(
            "Target X classification:",
            target_x_audit["classification"],
        )

        target_genes = [
            str(gene)
            for gene in target.var_names
        ]

        reference_gene_set = set(
            str(gene)
            for gene in reference.var_names
        )

        target_gene_set = set(
            target_genes
        )

        target_genes_missing_reference = sorted(
            target_gene_set
            - reference_gene_set
        )

        has_spatial = (
            "spatial" in target.obsm
        )

        spatial_shape = (
            list(target.obsm["spatial"].shape)
            if has_spatial
            else None
        )

        section("FOLD VALIDATION")

        fold_files = discover_fold_files(
            FOLD_DIR
        )

        fold_rows: list[dict[str, Any]] = []
        fold_manifest: list[dict[str, Any]] = []

        all_fold_checks_pass = True

        for fold in (1, 2, 3):
            observed_file = (
                fold_files[fold]["observed"]
            )

            heldout_file = (
                fold_files[fold]["heldout"]
            )

            observed_genes = read_gene_file(
                observed_file
            )

            heldout_genes = read_gene_file(
                heldout_file
            )

            observed_set = set(
                observed_genes
            )

            heldout_set = set(
                heldout_genes
            )

            overlap = sorted(
                observed_set
                & heldout_set
            )

            union = (
                observed_set
                | heldout_set
            )

            observed_missing_reference = sorted(
                observed_set
                - reference_gene_set
            )

            heldout_missing_reference = sorted(
                heldout_set
                - reference_gene_set
            )

            observed_missing_target = sorted(
                observed_set
                - target_gene_set
            )

            heldout_missing_target = sorted(
                heldout_set
                - target_gene_set
            )

            missing_from_fold_union = sorted(
                target_gene_set
                - union
            )

            extra_in_fold_union = sorted(
                union
                - target_gene_set
            )

            fold_pass = all(
                [
                    len(observed_genes) == 200,
                    len(heldout_genes) == 100,
                    len(observed_set) == 200,
                    len(heldout_set) == 100,
                    len(overlap) == 0,
                    len(union) == 300,
                    len(observed_missing_reference) == 0,
                    len(heldout_missing_reference) == 0,
                    len(observed_missing_target) == 0,
                    len(heldout_missing_target) == 0,
                    len(missing_from_fold_union) == 0,
                    len(extra_in_fold_union) == 0,
                ]
            )

            all_fold_checks_pass = (
                all_fold_checks_pass
                and fold_pass
            )

            fold_rows.append(
                {
                    "fold": fold,
                    "observed_file": str(
                        observed_file
                    ),
                    "heldout_file": str(
                        heldout_file
                    ),
                    "n_observed": len(
                        observed_genes
                    ),
                    "n_observed_unique": len(
                        observed_set
                    ),
                    "n_heldout": len(
                        heldout_genes
                    ),
                    "n_heldout_unique": len(
                        heldout_set
                    ),
                    "n_overlap": len(overlap),
                    "n_union": len(union),
                    "observed_missing_reference": len(
                        observed_missing_reference
                    ),
                    "heldout_missing_reference": len(
                        heldout_missing_reference
                    ),
                    "observed_missing_target": len(
                        observed_missing_target
                    ),
                    "heldout_missing_target": len(
                        heldout_missing_target
                    ),
                    "target_genes_missing_union": len(
                        missing_from_fold_union
                    ),
                    "extra_genes_in_union": len(
                        extra_in_fold_union
                    ),
                    "fold_pass": fold_pass,
                }
            )

            fold_manifest.append(
                {
                    "fold": fold,
                    "observed_file": str(
                        observed_file
                    ),
                    "heldout_file": str(
                        heldout_file
                    ),
                    "observed_genes": observed_genes,
                    "heldout_genes": heldout_genes,
                    "validation_pass": fold_pass,
                }
            )

            print(
                f"Fold {fold}: "
                f"observed={len(observed_genes)}, "
                f"heldout={len(heldout_genes)}, "
                f"overlap={len(overlap)}, "
                f"pass={fold_pass}"
            )

        fold_validation = pd.DataFrame(
            fold_rows
        )

        overall_pass = all(
            [
                reference.obs_names.is_unique,
                reference.var_names.is_unique,
                target.obs_names.is_unique,
                target.var_names.is_unique,
                reference_x_audit[
                    "classification"
                ]
                == "likely_raw_counts",
                target_x_audit[
                    "classification"
                ]
                == "likely_raw_counts",
                target.n_vars == 300,
                len(
                    target_genes_missing_reference
                )
                == 0,
                has_spatial,
                spatial_shape is not None,
                spatial_shape[0] == target.n_obs,
                spatial_shape[1] >= 2,
                cohort_metadata[
                    "n_unique_brains"
                ]
                == 10,
                cohort_metadata[
                    "includes_Br8667"
                ]
                is True,
                all_fold_checks_pass,
            ]
        )

        manifest = {
            "experiment": "Experiment_5_1",
            "description": (
                "All 10 Huuki snRNA-seq brains, including "
                "Br8667, used as the reference for Br8667 "
                "Xenium imputation."
            ),
            "reference": {
                "path": str(
                    REFERENCE_FILE
                ),
                "shape": [
                    int(reference.n_obs),
                    int(reference.n_vars),
                ],
                "matrix_source": "X",
                "matrix_audit": reference_x_audit,
                "layers": list(
                    reference.layers.keys()
                ),
                "cohort": cohort_metadata,
            },
            "target": {
                "path": str(
                    TARGET_FILE
                ),
                "shape": [
                    int(target.n_obs),
                    int(target.n_vars),
                ],
                "matrix_source": "X",
                "matrix_audit": target_x_audit,
                "layers": list(
                    target.layers.keys()
                ),
                "has_spatial": has_spatial,
                "spatial_shape": spatial_shape,
                "n_target_genes": len(
                    target_genes
                ),
                "target_genes_missing_reference": (
                    target_genes_missing_reference
                ),
            },
            "cross_validation": {
                "scheme": (
                    "3-fold gene-wise cross-validation"
                ),
                "observed_genes_per_fold": 200,
                "heldout_genes_per_fold": 100,
                "fold_directory": str(
                    FOLD_DIR
                ),
                "reuse_existing_folds": True,
                "folds": fold_manifest,
            },
            "methods": [
                "gimVI",
                "VISTA",
                "Tangram",
                "ENVI",
                "TransImp",
                "SpaGE",
            ],
            "overall_pass": overall_pass,
        }

        fold_validation.to_csv(
            FOLD_VALIDATION_FILE,
            index=False,
        )

        cohort_frame.to_csv(
            COHORT_SUMMARY_FILE,
            index=False,
        )

        with MANIFEST_FILE.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                manifest,
                handle,
                indent=2,
            )

        report_lines = [
            "EXPERIMENT 5.1 PREFLIGHT REPORT",
            "=" * 100,
            "",
            f"Reference: {REFERENCE_FILE}",
            (
                "Reference shape: "
                f"{reference.n_obs} cells x "
                f"{reference.n_vars} genes"
            ),
            (
                "Reference X classification: "
                f"{reference_x_audit['classification']}"
            ),
            (
                "Reference unique brains: "
                f"{cohort_metadata['n_unique_brains']}"
            ),
            (
                "Reference unique samples: "
                f"{cohort_metadata['n_unique_samples']}"
            ),
            (
                "Reference includes Br8667: "
                f"{cohort_metadata['includes_Br8667']}"
            ),
            "",
            f"Target: {TARGET_FILE}",
            (
                "Target shape: "
                f"{target.n_obs} cells x "
                f"{target.n_vars} genes"
            ),
            (
                "Target X classification: "
                f"{target_x_audit['classification']}"
            ),
            f"Target has spatial coordinates: {has_spatial}",
            f"Target spatial shape: {spatial_shape}",
            (
                "Target genes missing from reference: "
                f"{len(target_genes_missing_reference)}"
            ),
            "",
            "Fold results:",
        ]

        for row in fold_rows:
            report_lines.append(
                f"  Fold {row['fold']}: "
                f"observed={row['n_observed']}, "
                f"heldout={row['n_heldout']}, "
                f"overlap={row['n_overlap']}, "
                f"pass={row['fold_pass']}"
            )

        report_lines.extend(
            [
                "",
                f"OVERALL PASS: {overall_pass}",
                "",
                f"Manifest: {MANIFEST_FILE}",
                (
                    "Fold validation: "
                    f"{FOLD_VALIDATION_FILE}"
                ),
                (
                    "Cohort summary: "
                    f"{COHORT_SUMMARY_FILE}"
                ),
            ]
        )

        PREFLIGHT_REPORT_FILE.write_text(
            "\n".join(report_lines) + "\n",
            encoding="utf-8",
        )

        section("FINAL RESULT")

        print(
            "Reference X:",
            reference_x_audit[
                "classification"
            ],
        )
        print(
            "Target X:",
            target_x_audit[
                "classification"
            ],
        )
        print(
            "Reference brains:",
            cohort_metadata[
                "n_unique_brains"
            ],
        )
        print(
            "Reference samples:",
            cohort_metadata[
                "n_unique_samples"
            ],
        )
        print(
            "Includes Br8667:",
            cohort_metadata[
                "includes_Br8667"
            ],
        )
        print(
            "Target genes:",
            target.n_vars,
        )
        print(
            "Target genes missing from reference:",
            len(
                target_genes_missing_reference
            ),
        )
        print(
            "All fold checks pass:",
            all_fold_checks_pass,
        )
        print(
            "OVERALL PASS:",
            overall_pass,
        )

        print("\nSaved:")
        print(MANIFEST_FILE)
        print(FOLD_VALIDATION_FILE)
        print(COHORT_SUMMARY_FILE)
        print(PREFLIGHT_REPORT_FILE)

        if not overall_pass:
            raise RuntimeError(
                "Experiment 5.1 preflight validation failed. "
                "Review the generated report before running models."
            )

    finally:
        reference.file.close()
        target.file.close()


if __name__ == "__main__":
    main()
