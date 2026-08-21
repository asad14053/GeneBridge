#!/usr/bin/env python3

from pathlib import Path

import h5py
import numpy as np
import pandas as pd


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
}

OUTPUT = (
    PROJECT_ROOT
    / "outputs/imputation_beta/Br8667/reference_alignment_check"
    / "expression_matrix_audit_lowmem.csv"
)

MAX_SAMPLE_VALUES = 300_000


def decode(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")

    return str(value)


def sample_dataset(
    dataset: h5py.Dataset,
    max_values: int = MAX_SAMPLE_VALUES,
) -> np.ndarray:
    """
    Read small contiguous sections from the beginning, middle, and end.
    This avoids loading the complete compressed HDF5 dataset.
    """

    n_values = int(dataset.shape[0])

    if n_values == 0:
        return np.asarray([], dtype=float)

    section_size = min(
        max_values // 3,
        n_values,
    )

    starts = [
        0,
        max(0, n_values // 2 - section_size // 2),
        max(0, n_values - section_size),
    ]

    arrays = []

    for start in starts:
        end = min(
            start + section_size,
            n_values,
        )

        if end > start:
            arrays.append(
                np.asarray(dataset[start:end])
            )

    if not arrays:
        return np.asarray([], dtype=float)

    return np.concatenate(arrays)


def sample_dense_matrix(
    dataset: h5py.Dataset,
) -> np.ndarray:
    n_obs = dataset.shape[0]

    row_indices = sorted(
        set(
            [
                0,
                n_obs // 2,
                n_obs - 1,
            ]
        )
    )

    arrays = []

    for row in row_indices:
        arrays.append(
            np.asarray(dataset[row : row + 1, :]).ravel()
        )

    return np.concatenate(arrays)


def summarize_values(
    values: np.ndarray,
    total_values: int,
    nnz: int,
) -> dict:
    values = np.asarray(values)

    if values.size == 0:
        return {
            "minimum": np.nan,
            "maximum": np.nan,
            "mean_sampled_nonzero": np.nan,
            "integer_like_fraction": np.nan,
            "nonnegative": True,
            "nonzero_fraction": 0.0,
        }

    nonzero_values = values[
        values != 0
    ]

    if nonzero_values.size == 0:
        mean_nonzero = 0.0
    else:
        mean_nonzero = float(
            nonzero_values.mean()
        )

    integer_like = np.isclose(
        values,
        np.round(values),
        atol=1e-6,
    )

    return {
        "minimum": float(
            min(0.0, float(values.min()))
        ),
        "maximum": float(values.max()),
        "mean_sampled_nonzero": mean_nonzero,
        "integer_like_fraction": float(
            integer_like.mean()
        ),
        "nonnegative": bool(
            np.all(values >= 0)
        ),
        "nonzero_fraction": float(
            nnz / total_values
        ),
    }


def classify(stats: dict) -> str:
    if (
        stats["nonnegative"]
        and stats["integer_like_fraction"] >= 0.999
    ):
        return "likely_raw_counts"

    if (
        stats["nonnegative"]
        and stats["integer_like_fraction"] < 0.99
        and stats["maximum"] < 30
    ):
        return "likely_log_normalized"

    return "uncertain_review_required"


def inspect_matrix(
    h5file: h5py.File,
    matrix_path: str,
) -> dict:
    obj = h5file[matrix_path]

    if isinstance(obj, h5py.Group):
        encoding = decode(
            obj.attrs.get(
                "encoding-type",
                "unknown_sparse",
            )
        )

        shape = tuple(
            int(value)
            for value in obj.attrs["shape"]
        )

        data = obj["data"]

        values = sample_dataset(data)
        nnz = int(data.shape[0])
        total_values = int(
            shape[0] * shape[1]
        )

    elif isinstance(obj, h5py.Dataset):
        encoding = decode(
            obj.attrs.get(
                "encoding-type",
                "array",
            )
        )

        shape = tuple(
            int(value)
            for value in obj.shape
        )

        values = sample_dense_matrix(obj)
        nnz = int(
            np.count_nonzero(values)
        )

        sampled_fraction = (
            values.size
            / int(shape[0] * shape[1])
        )

        estimated_nnz = (
            nnz / sampled_fraction
            if sampled_fraction > 0
            else 0
        )

        nnz = int(estimated_nnz)
        total_values = int(
            shape[0] * shape[1]
        )

    else:
        raise TypeError(
            f"Unsupported HDF5 object: {matrix_path}"
        )

    stats = summarize_values(
        values=values,
        total_values=total_values,
        nnz=nnz,
    )

    return {
        "matrix_source": matrix_path,
        "encoding": encoding,
        "n_obs": shape[0],
        "n_vars": shape[1],
        "sampled_values": values.size,
        **stats,
        "classification": classify(stats),
    }


def main() -> None:
    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    records = []

    for experiment, path in FILES.items():
        print("\n" + "=" * 100)
        print(experiment)
        print("=" * 100)
        print(path)

        if not path.exists():
            raise FileNotFoundError(path)

        with h5py.File(path, "r") as h5file:
            matrix_paths = ["X"]

            if "layers" in h5file:
                for layer_name in h5file["layers"].keys():
                    matrix_paths.append(
                        f"layers/{layer_name}"
                    )

            print("Matrix sources:", matrix_paths)

            for matrix_path in matrix_paths:
                result = inspect_matrix(
                    h5file,
                    matrix_path,
                )

                result = {
                    "experiment": experiment,
                    "file": str(path),
                    **result,
                }

                records.append(result)

                print()
                print(
                    pd.DataFrame(
                        [result]
                    ).to_string(index=False)
                )

    results = pd.DataFrame(records)

    results.to_csv(
        OUTPUT,
        index=False,
    )

    print("\n" + "=" * 100)
    print("FINAL SUMMARY")
    print("=" * 100)

    print(
        results[
            [
                "experiment",
                "matrix_source",
                "encoding",
                "minimum",
                "maximum",
                "mean_sampled_nonzero",
                "integer_like_fraction",
                "classification",
            ]
        ].to_string(index=False)
    )

    print("\nSaved:")
    print(OUTPUT)


if __name__ == "__main__":
    main()
