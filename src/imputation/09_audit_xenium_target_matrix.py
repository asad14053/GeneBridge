#!/usr/bin/env python3

"""
Audit the Br8667 Xenium target H5AD and determine whether adata.X
contains raw counts, log-normalized values, or an uncertain matrix.

Default input:
    data/processed/imputation_beta/Br8667/
    spatial_data_xenium_Br8667_vista_qc.h5ad
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

DEFAULT_INPUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "imputation_beta"
    / "Br8667"
    / "spatial_data_xenium_Br8667_vista_qc.h5ad"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit whether Xenium target adata.X contains raw counts."
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input Xenium H5AD file.",
    )

    return parser.parse_args()


def classify_matrix(
    nonnegative: bool,
    integer_like_fraction: float,
    maximum: float,
) -> str:
    if nonnegative and integer_like_fraction >= 0.999:
        return "RAW COUNTS"

    if (
        nonnegative
        and integer_like_fraction < 0.99
        and maximum < 30
    ):
        return "LIKELY LOG-NORMALIZED"

    return "UNCERTAIN — REVIEW REQUIRED"


def main() -> None:
    args = parse_arguments()
    input_path = args.input.resolve()

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input H5AD was not found:\n{input_path}"
        )

    adata = ad.read_h5ad(
        input_path,
        backed="r",
    )

    try:
        print("=" * 90)
        print("XENIUM TARGET MATRIX AUDIT")
        print("=" * 90)
        print("File:", input_path)
        print("Shape [cells x genes]:", adata.shape)
        print("Layers:", list(adata.layers.keys()))
        print("Has adata.raw:", adata.raw is not None)
        print("obsm:", list(adata.obsm.keys()))
        print("uns keys:", list(adata.uns.keys()))
        print("Unique obs names:", adata.obs_names.is_unique)
        print("Unique var names:", adata.var_names.is_unique)

        if adata.n_obs == 0 or adata.n_vars == 0:
            raise RuntimeError(
                f"Target has an invalid shape: {adata.shape}"
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

        blocks: list[np.ndarray] = []

        for row in row_indices:
            block = adata.X[
                row : row + 1,
                :,
            ]

            if sparse.issparse(block):
                block = block.toarray()

            blocks.append(
                np.asarray(block).ravel()
            )

        values = np.concatenate(blocks)

        if values.size == 0:
            raise RuntimeError(
                "No values were read from adata.X."
            )

        nonzero = values[
            values != 0
        ]

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

        minimum = float(
            values.min()
        )

        maximum = float(
            values.max()
        )

        mean_nonzero = (
            float(nonzero.mean())
            if nonzero.size
            else 0.0
        )

        zero_fraction = float(
            np.mean(values == 0)
        )

        classification = classify_matrix(
            nonnegative=nonnegative,
            integer_like_fraction=integer_like_fraction,
            maximum=maximum,
        )

        print()
        print("Sampled X statistics")
        print("-" * 90)
        print("Sampled rows:", row_indices)
        print("Sampled values:", values.size)
        print("Minimum:", minimum)
        print("Maximum:", maximum)
        print("Nonnegative:", nonnegative)
        print(
            "Integer-like fraction:",
            integer_like_fraction,
        )
        print("Mean nonzero:", mean_nonzero)
        print("Zero fraction:", zero_fraction)

        if "log1p" in adata.uns:
            print()
            print("Warning: adata.uns contains log1p metadata:")
            print(adata.uns["log1p"])

        print()
        print("=" * 90)
        print("CLASSIFICATION:", classification)
        print("=" * 90)

        if classification == "RAW COUNTS":
            print(
                "Interpretation: adata.X is suitable as raw-count input "
                "for VISTA, gimVI, ENVI, and TransImp."
            )

        elif classification == "LIKELY LOG-NORMALIZED":
            print(
                "Interpretation: adata.X should not be used as raw-count "
                "input for count-based models."
            )

        else:
            print(
                "Interpretation: inspect the data-generation workflow "
                "before using this matrix."
            )

    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
