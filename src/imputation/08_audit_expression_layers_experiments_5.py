#!/usr/bin/env python3

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

FILES = {
    "Experiment_5_Br8667_only": (
        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667"
        / "seq_data_huuki_snrna_Br8667_vista.h5ad"
    ),
    "Experiment_5_1_All_10_Huuki": (
        PROJECT_ROOT
        / "data/processed/snrnaseq/sce_DLPFC_annotated"
        / "huuki_snrna_reference_full_allgenes.h5ad"
    ),
    "Experiment_5_3_Nine_Huuki": (
        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667"
        / "seq_data_huuki_snrna_9samples_excluding_xenium_overlap.h5ad"
    ),
    "Target_Xenium_Br8667_QC": (
        PROJECT_ROOT
        / "data/processed/imputation_beta/Br8667"
        / "spatial_data_xenium_Br8667_vista_qc.h5ad"
    ),
}

OUTPUT = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667/reference_alignment_check"
    / "expression_matrix_audit.csv"
)


def summarize_matrix(matrix) -> dict:
    if sparse.issparse(matrix):
        values = matrix.data
        nnz = matrix.nnz
        total = matrix.shape[0] * matrix.shape[1]
    else:
        array = np.asarray(matrix)
        values = array.ravel()
        nnz = np.count_nonzero(array)
        total = array.size

    values = np.asarray(values)

    if values.size == 0:
        return {
            "minimum": 0.0,
            "maximum": 0.0,
            "mean_nonzero": 0.0,
            "integer_like_fraction": 1.0,
            "nonnegative": True,
            "nonzero_fraction": 0.0,
        }

    integer_like = np.isclose(
        values,
        np.round(values),
        atol=1e-6,
    )

    return {
        "minimum": float(values.min()),
        "maximum": float(values.max()),
        "mean_nonzero": float(values.mean()),
        "integer_like_fraction": float(integer_like.mean()),
        "nonnegative": bool(np.all(values >= 0)),
        "nonzero_fraction": float(nnz / total),
    }


def read_sample(matrix, n_obs: int):
    blocks = []

    ranges = [
        (0, min(100, n_obs)),
        (
            max(0, n_obs // 2 - 50),
            min(n_obs, n_obs // 2 + 50),
        ),
        (
            max(0, n_obs - 100),
            n_obs,
        ),
    ]

    for start, end in ranges:
        if end <= start:
            continue

        block = matrix[start:end, :]

        if sparse.issparse(block):
            block = block.tocsr()

        blocks.append(block)

    if not blocks:
        raise RuntimeError("No matrix rows were sampled.")

    if sparse.issparse(blocks[0]):
        return sparse.vstack(blocks, format="csr")

    return np.vstack(
        [np.asarray(block) for block in blocks]
    )


def classify_matrix(stats: dict) -> str:
    if (
        stats["nonnegative"]
        and stats["integer_like_fraction"] >= 0.999
    ):
        return "likely_raw_counts"

    if (
        stats["nonnegative"]
        and stats["integer_like_fraction"] < 0.95
        and stats["maximum"] < 30
    ):
        return "likely_log_normalized"

    return "uncertain_review_required"


def main() -> None:
    records = []

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    for experiment, path in FILES.items():
        print("\n" + "=" * 100)
        print(experiment)
        print("=" * 100)
        print(path)

        if not path.exists():
            raise FileNotFoundError(path)

        adata = ad.read_h5ad(
            path,
            backed="r",
        )

        print("Shape:", adata.shape)
        print("Layers:", list(adata.layers.keys()))

        matrix_sources = {
            "X": adata.X,
        }

        for layer in adata.layers.keys():
            matrix_sources[f"layer:{layer}"] = (
                adata.layers[layer]
            )

        for source_name, matrix in matrix_sources.items():
            print("\nSampling:", source_name)

            sampled = read_sample(
                matrix,
                adata.n_obs,
            )

            stats = summarize_matrix(sampled)
            classification = classify_matrix(stats)

            record = {
                "experiment": experiment,
                "file": str(path),
                "n_obs": adata.n_obs,
                "n_vars": adata.n_vars,
                "matrix_source": source_name,
                **stats,
                "classification": classification,
            }

            records.append(record)

            print(
                pd.DataFrame([record]).to_string(
                    index=False
                )
            )

        adata.file.close()

    results = pd.DataFrame(records)

    results.to_csv(
        OUTPUT,
        index=False,
    )

    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        results[
            [
                "experiment",
                "matrix_source",
                "minimum",
                "maximum",
                "mean_nonzero",
                "integer_like_fraction",
                "classification",
            ]
        ].to_string(index=False)
    )

    print("\nSaved:")
    print(OUTPUT)


if __name__ == "__main__":
    main()
