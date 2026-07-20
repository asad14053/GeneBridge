#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
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
    / "data/processed/imputation_beta/Br8667/benchmark/splits"
)
DEFAULT_REPORT_DIR = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667/gene_split"
)

EXPECTED_GENES = 300
SEED = 8667
N_STRATA = 5
PER_STRATUM = {"train": 42, "validation": 9, "test": 9}
EXPECTED_SPLIT = {"train": 210, "validation": 45, "test": 45}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the fixed Br8667 210/45/45 gene split."
    )
    parser.add_argument("--xenium", type=Path, default=DEFAULT_XENIUM)
    parser.add_argument("--split-dir", type=Path, default=DEFAULT_SPLIT_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing frozen split files.",
    )
    return parser.parse_args()


def calculate_gene_statistics(adata: ad.AnnData) -> pd.DataFrame:
    matrix = adata.X

    if sparse.issparse(matrix):
        values = np.asarray(matrix.data)
        means = np.asarray(matrix.mean(axis=0)).ravel()
        mean_squared = np.asarray(
            matrix.multiply(matrix).mean(axis=0)
        ).ravel()
        detected = np.asarray(matrix.getnnz(axis=0)).ravel()
    else:
        values = np.asarray(matrix).reshape(-1)
        dense = np.asarray(matrix)
        means = dense.mean(axis=0)
        mean_squared = np.square(dense).mean(axis=0)
        detected = np.count_nonzero(dense > 0, axis=0)

    if values.size == 0:
        raise ValueError("Xenium expression matrix is empty.")
    if not np.isfinite(values).all():
        raise ValueError("Xenium expression contains NaN or infinity.")
    if np.any(values < 0):
        raise ValueError("Xenium expression contains negative values.")

    variances = np.maximum(mean_squared - means**2, 0)
    detection_fraction = detected / float(adata.n_obs)

    stats = pd.DataFrame(
        {
            "gene": adata.var_names.astype(str),
            "mean_expression": means,
            "variance_expression": variances,
            "detected_cell_count": detected,
            "detection_fraction": detection_fraction,
        }
    )

    if not stats["gene"].is_unique:
        raise ValueError("Xenium gene names are not unique.")

    stats["log1p_mean_expression"] = np.log1p(
        stats["mean_expression"]
    )
    stats["log1p_variance_expression"] = np.log1p(
        stats["variance_expression"]
    )
    stats["mean_percentile"] = stats[
        "log1p_mean_expression"
    ].rank(method="average", pct=True)
    stats["variance_percentile"] = stats[
        "log1p_variance_expression"
    ].rank(method="average", pct=True)
    stats["detection_percentile"] = stats[
        "detection_fraction"
    ].rank(method="average", pct=True)

    stats["composite_expression_score"] = stats[
        [
            "mean_percentile",
            "variance_percentile",
            "detection_percentile",
        ]
    ].mean(axis=1)

    return stats


def assign_strata(stats: pd.DataFrame) -> pd.DataFrame:
    result = stats.copy()

    # method="first" is used only to break ties deterministically.
    score_rank = result[
        "composite_expression_score"
    ].rank(method="first")

    result["expression_stratum_number"] = (
        pd.qcut(score_rank, q=N_STRATA, labels=False).astype(int) + 1
    )

    labels = {
        1: "S1_lowest",
        2: "S2_low",
        3: "S3_medium",
        4: "S4_high",
        5: "S5_highest",
    }
    result["expression_stratum"] = result[
        "expression_stratum_number"
    ].map(labels)

    counts = result["expression_stratum"].value_counts()

    if len(counts) != N_STRATA or not counts.eq(60).all():
        raise ValueError(
            "Expected five strata with 60 genes each:\n"
            + counts.to_string()
        )

    return result


def assign_splits(
    stratified: pd.DataFrame,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    parts = []

    strata = [
        "S1_lowest",
        "S2_low",
        "S3_medium",
        "S4_high",
        "S5_highest",
    ]

    for stratum in strata:
        group = stratified.loc[
            stratified["expression_stratum"] == stratum
        ].copy()

        shuffled = rng.permutation(group.index.to_numpy())
        train_end = PER_STRATUM["train"]
        validation_end = train_end + PER_STRATUM["validation"]

        group["split"] = ""
        group.loc[shuffled[:train_end], "split"] = "train"
        group.loc[
            shuffled[train_end:validation_end],
            "split",
        ] = "validation"
        group.loc[shuffled[validation_end:], "split"] = "test"
        parts.append(group)

    assignments = pd.concat(parts).reset_index(drop=True)
    return assignments


def validate_assignments(
    assignments: pd.DataFrame,
    original_genes: set[str],
) -> None:
    sets = {
        split: set(
            assignments.loc[
                assignments["split"] == split,
                "gene",
            ]
        )
        for split in EXPECTED_SPLIT
    }

    checks = {
        "total_genes_is_300": len(assignments) == EXPECTED_GENES,
        "train_genes_is_210": len(sets["train"]) == 210,
        "validation_genes_is_45": len(sets["validation"]) == 45,
        "test_genes_is_45": len(sets["test"]) == 45,
        "train_validation_disjoint": sets["train"].isdisjoint(
            sets["validation"]
        ),
        "train_test_disjoint": sets["train"].isdisjoint(
            sets["test"]
        ),
        "validation_test_disjoint": sets["validation"].isdisjoint(
            sets["test"]
        ),
        "all_genes_assigned_once": (
            sets["train"] | sets["validation"] | sets["test"]
        )
        == original_genes,
        "gene_names_unique": assignments["gene"].is_unique,
    }

    counts = assignments.groupby(
        ["expression_stratum", "split"]
    ).size().unstack(fill_value=0)

    checks["every_stratum_has_42_9_9"] = (
        counts["train"].eq(42).all()
        and counts["validation"].eq(9).all()
        and counts["test"].eq(9).all()
    )

    print("\nSplit validation:")
    for name, passed in checks.items():
        print(f"{name}: {passed}")

    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError(f"Split validation failed: {failed}")


def expression_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for split in ["train", "validation", "test"]:
        subset = assignments.loc[assignments["split"] == split]
        rows.append(
            {
                "split": split,
                "n_genes": len(subset),
                "mean_expression_mean": subset[
                    "mean_expression"
                ].mean(),
                "mean_expression_median": subset[
                    "mean_expression"
                ].median(),
                "detection_fraction_mean": subset[
                    "detection_fraction"
                ].mean(),
                "detection_fraction_median": subset[
                    "detection_fraction"
                ].median(),
                "variance_mean": subset[
                    "variance_expression"
                ].mean(),
                "variance_median": subset[
                    "variance_expression"
                ].median(),
                "composite_score_mean": subset[
                    "composite_expression_score"
                ].mean(),
            }
        )

    return pd.DataFrame(rows)


def save_plot(counts: pd.DataFrame, path: Path) -> None:
    pivot = counts.pivot(
        index="expression_stratum",
        columns="split",
        values="n_genes",
    ).fillna(0)

    pivot = pivot.reindex(
        [
            "S1_lowest",
            "S2_low",
            "S3_medium",
            "S4_high",
            "S5_highest",
        ]
    )[
        ["train", "validation", "test"]
    ]

    axis = pivot.plot(kind="bar", figsize=(9, 5))
    axis.set_xlabel("Expression stratum")
    axis.set_ylabel("Number of genes")
    axis.set_title("Br8667 gene split within expression strata")
    axis.tick_params(axis="x", rotation=0)

    figure = axis.get_figure()
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)


def main() -> None:
    args = parse_args()

    if not args.xenium.exists():
        raise FileNotFoundError(args.xenium)

    args.split_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {
        "train": args.split_dir / "train_genes_210.csv",
        "validation": args.split_dir / "validation_genes_45.csv",
        "test": args.split_dir / "test_genes_45.csv",
        "all": args.split_dir / "all_gene_split_assignments.csv",
        "summary": args.split_dir / "gene_split_summary.csv",
        "manifest": args.split_dir / "split_manifest.csv",
    }

    existing = [path for path in output_paths.values() if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Frozen split files already exist and were not overwritten:\n"
            + "\n".join(str(path) for path in existing)
            + "\nUse --overwrite only intentionally."
        )

    print("=" * 100)
    print("Day 2: Create fixed Br8667 train/validation/test gene split")
    print("=" * 100)

    xenium = ad.read_h5ad(args.xenium)
    print("\nInput Xenium:")
    print(xenium)

    if xenium.n_vars != EXPECTED_GENES:
        raise ValueError(
            f"Expected 300 genes, found {xenium.n_vars}."
        )

    stats = calculate_gene_statistics(xenium)
    stratified = assign_strata(stats)
    assignments = assign_splits(stratified, args.seed)

    validate_assignments(
        assignments,
        set(xenium.var_names.astype(str)),
    )

    columns = [
        "gene",
        "split",
        "expression_stratum",
        "expression_stratum_number",
        "mean_expression",
        "variance_expression",
        "detected_cell_count",
        "detection_fraction",
        "log1p_mean_expression",
        "log1p_variance_expression",
        "mean_percentile",
        "variance_percentile",
        "detection_percentile",
        "composite_expression_score",
    ]

    assignments[columns].to_csv(output_paths["all"], index=False)

    for split in ["train", "validation", "test"]:
        assignments.loc[
            assignments["split"] == split,
            columns,
        ].to_csv(output_paths[split], index=False)

    summary = expression_summary(assignments)
    summary.to_csv(output_paths["summary"], index=False)

    manifest = pd.DataFrame(
        [
            {"field": "input_xenium", "value": str(args.xenium)},
            {"field": "random_seed", "value": args.seed},
            {"field": "total_genes", "value": 300},
            {"field": "train_genes", "value": 210},
            {"field": "validation_genes", "value": 45},
            {"field": "test_genes", "value": 45},
            {"field": "number_of_strata", "value": 5},
            {"field": "genes_per_stratum", "value": 60},
            {"field": "train_per_stratum", "value": 42},
            {"field": "validation_per_stratum", "value": 9},
            {"field": "test_per_stratum", "value": 9},
            {
                "field": "stratification",
                "value": (
                    "percentile ranks of log1p mean, log1p variance, "
                    "and detection fraction"
                ),
            },
            {"field": "split_status", "value": "FROZEN"},
        ]
    )
    manifest.to_csv(output_paths["manifest"], index=False)

    stratum_counts = assignments.groupby(
        ["expression_stratum", "split"]
    ).size().rename("n_genes").reset_index()

    stratum_counts_path = (
        args.report_dir / "10_split_stratum_counts.csv"
    )
    expression_summary_path = (
        args.report_dir / "10_split_expression_summary.csv"
    )
    figure_path = (
        args.report_dir / "10_split_stratum_counts.png"
    )

    stratum_counts.to_csv(stratum_counts_path, index=False)
    summary.to_csv(expression_summary_path, index=False)
    save_plot(stratum_counts, figure_path)

    print("\nSplit sizes:")
    print(
        assignments["split"]
        .value_counts()
        .reindex(["train", "validation", "test"])
        .to_string()
    )

    print("\nStratum allocation:")
    print(
        stratum_counts.pivot(
            index="expression_stratum",
            columns="split",
            values="n_genes",
        )
        .fillna(0)
        .astype(int)
        .to_string()
    )

    print("\nExpression summary:")
    print(summary.to_string(index=False))

    print("\n" + "=" * 100)
    print("PASS: fixed Br8667 210/45/45 gene split created")
    print("=" * 100)

    print("\nPermanent split files:")
    for path in output_paths.values():
        print(path)

    print("\nQC reports:")
    print(stratum_counts_path)
    print(expression_summary_path)
    print(figure_path)

    print(
        "\nDo not regenerate this split separately for "
        "VISTA, gimVI, Tangram, or other methods."
    )


if __name__ == "__main__":
    main()
