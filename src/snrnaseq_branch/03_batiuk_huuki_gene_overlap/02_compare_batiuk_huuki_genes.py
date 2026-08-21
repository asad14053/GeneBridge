#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

# All Sample - uncomment
HUUKI_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

# Single Sample - uncomment
#HUUKI_H5AD = (
#    PROJECT_ROOT
#    / "data"
#    / "processed"
#    / "imputation_beta"
#    / "Br8667"
#    / "seq_data_huuki_snrna_Br8667_vista.h5ad"
#)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "snrnaseq_branch"
    / "03_batiuk_huuki_gene_overlap"
)

BATIUK_UNION_CSV = (
    OUTPUT_DIR
    / "batiuk_gene_union.csv"
)

BATIUK_INTERSECTION_CSV = (
    OUTPUT_DIR
    / "batiuk_gene_intersection.csv"
)

SOURCE_COMPARISON_CSV = (
    OUTPUT_DIR
    / "huuki_gene_source_comparison.csv"
)

SUMMARY_CSV = (
    OUTPUT_DIR
    / "batiuk_huuki_gene_overlap_summary.csv"
)

SHARED_GENES_CSV = (
    OUTPUT_DIR
    / "batiuk_huuki_shared_genes.csv"
)

BATIUK_ONLY_CSV = (
    OUTPUT_DIR
    / "genes_only_in_batiuk.csv"
)

HUUKI_ONLY_CSV = (
    OUTPUT_DIR
    / "genes_only_in_huuki.csv"
)

SHARED_CASE_INSENSITIVE_CSV = (
    OUTPUT_DIR
    / "batiuk_huuki_shared_genes_case_insensitive.csv"
)


def clean_gene_values(values) -> pd.Index:
    series = pd.Series(
        values,
        dtype="string",
    )

    series = series.str.strip()

    series = series[
        series.notna()
        & series.ne("")
        & series.str.lower().ne("nan")
    ]

    return pd.Index(
        pd.unique(
            series.astype(str)
        )
    )


def add_summary(
    records: list[dict],
    metric: str,
    value,
) -> None:
    records.append(
        {
            "metric": metric,
            "value": value,
        }
    )


def main() -> None:
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    for path in [
        HUUKI_H5AD,
        BATIUK_UNION_CSV,
        BATIUK_INTERSECTION_CSV,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Required input not found:\n{path}"
            )

    batiuk_union = clean_gene_values(
        pd.read_csv(
            BATIUK_UNION_CSV
        )["gene"]
    )

    batiuk_intersection = clean_gene_values(
        pd.read_csv(
            BATIUK_INTERSECTION_CSV
        )["gene"]
    )

    batiuk_union_set = set(
        batiuk_union
    )

    batiuk_intersection_set = set(
        batiuk_intersection
    )

    print("=" * 80)
    print("Reading Huuki H5AD in backed mode")
    print("=" * 80)

    huuki = ad.read_h5ad(
        HUUKI_H5AD,
        backed="r",
    )

    print(huuki)
    print(
        "Huuki shape [cells x genes]:",
        huuki.shape,
    )

    candidate_sources: dict[str, pd.Index] = {
        "var_names": clean_gene_values(
            huuki.var_names
        )
    }

    candidate_columns = [
        "gene_name",
        "gene_symbol",
        "symbol",
        "Symbol",
        "gene_id",
    ]

    for column in candidate_columns:
        if column in huuki.var.columns:
            candidate_sources[column] = (
                clean_gene_values(
                    huuki.var[column]
                )
            )

    source_records = []

    for source_name, huuki_genes in candidate_sources.items():

        huuki_set = set(
            huuki_genes
        )

        shared_exact = (
            batiuk_union_set
            & huuki_set
        )

        batiuk_upper = {
            gene.upper()
            for gene in batiuk_union_set
        }

        huuki_upper = {
            gene.upper()
            for gene in huuki_set
        }

        shared_case_insensitive = (
            batiuk_upper
            & huuki_upper
        )

        source_records.append(
            {
                "huuki_gene_source": source_name,
                "n_unique_huuki_values": len(
                    huuki_set
                ),
                "shared_exact_with_batiuk_union": len(
                    shared_exact
                ),
                "shared_case_insensitive_with_batiuk_union": len(
                    shared_case_insensitive
                ),
            }
        )

    source_comparison = pd.DataFrame(
        source_records
    ).sort_values(
        by=[
            "shared_exact_with_batiuk_union",
            "shared_case_insensitive_with_batiuk_union",
        ],
        ascending=False,
    )

    source_comparison.to_csv(
        SOURCE_COMPARISON_CSV,
        index=False,
    )

    selected_source = str(
        source_comparison.iloc[0][
            "huuki_gene_source"
        ]
    )

    huuki_genes = candidate_sources[
        selected_source
    ]

    huuki.file.close()

    huuki_set = set(
        huuki_genes
    )

    shared_exact = sorted(
        batiuk_union_set
        & huuki_set
    )

    shared_intersection_exact = sorted(
        batiuk_intersection_set
        & huuki_set
    )

    batiuk_only = sorted(
        batiuk_union_set
        - huuki_set
    )

    huuki_only = sorted(
        huuki_set
        - batiuk_union_set
    )

    batiuk_upper = {
        gene.upper()
        for gene in batiuk_union_set
    }

    huuki_upper = {
        gene.upper()
        for gene in huuki_set
    }

    shared_case_insensitive = sorted(
        batiuk_upper
        & huuki_upper
    )

    pd.DataFrame(
        {"gene": shared_exact}
    ).to_csv(
        SHARED_GENES_CSV,
        index=False,
    )

    pd.DataFrame(
        {"gene": shared_case_insensitive}
    ).to_csv(
        SHARED_CASE_INSENSITIVE_CSV,
        index=False,
    )

    pd.DataFrame(
        {"gene": batiuk_only}
    ).to_csv(
        BATIUK_ONLY_CSV,
        index=False,
    )

    pd.DataFrame(
        {"gene": huuki_only}
    ).to_csv(
        HUUKI_ONLY_CSV,
        index=False,
    )

    summary_records: list[dict] = []

    add_summary(
        summary_records,
        "selected_huuki_gene_source",
        selected_source,
    )

    add_summary(
        summary_records,
        "batiuk_union_genes",
        len(batiuk_union_set),
    )

    add_summary(
        summary_records,
        "batiuk_intersection_genes",
        len(batiuk_intersection_set),
    )

    add_summary(
        summary_records,
        "huuki_unique_genes",
        len(huuki_set),
    )

    add_summary(
        summary_records,
        "shared_exact_batiuk_union_vs_huuki",
        len(shared_exact),
    )

    add_summary(
        summary_records,
        "shared_exact_batiuk_intersection_vs_huuki",
        len(shared_intersection_exact),
    )

    add_summary(
        summary_records,
        "shared_case_insensitive",
        len(shared_case_insensitive),
    )

    add_summary(
        summary_records,
        "genes_only_in_batiuk",
        len(batiuk_only),
    )

    add_summary(
        summary_records,
        "genes_only_in_huuki",
        len(huuki_only),
    )

    add_summary(
        summary_records,
        "percent_batiuk_union_shared_with_huuki",
        round(
            100
            * len(shared_exact)
            / len(batiuk_union_set),
            4,
        ),
    )

    add_summary(
        summary_records,
        "percent_huuki_shared_with_batiuk",
        round(
            100
            * len(shared_exact)
            / len(huuki_set),
            4,
        ),
    )

    summary = pd.DataFrame(
        summary_records
    )

    summary.to_csv(
        SUMMARY_CSV,
        index=False,
    )

    print("\n" + "=" * 80)
    print("HUUKI–BATIUK GENE OVERLAP")
    print("=" * 80)

    print(
        source_comparison.to_string(
            index=False
        )
    )

    print("\nSelected Huuki gene source:")
    print(selected_source)

    print("\nFinal summary:")
    print(
        summary.to_string(
            index=False
        )
    )

    print("\nCreated:")
    print(SUMMARY_CSV)
    print(SHARED_GENES_CSV)
    print(BATIUK_ONLY_CSV)
    print(HUUKI_ONLY_CSV)


if __name__ == "__main__":
    main()
