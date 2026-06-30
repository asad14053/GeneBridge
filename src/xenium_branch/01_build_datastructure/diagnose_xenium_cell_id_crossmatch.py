#!/usr/bin/env python3

from pathlib import Path
import gzip
import io
import h5py
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "xenium"

OUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "xenium_branch"
    / "01_build_datastructure"
    / "cell_id_crossmatch"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)


def read_h5_barcodes(matrix_path: Path) -> pd.Index:
    with h5py.File(matrix_path, "r") as f:
        barcodes = f["matrix"]["barcodes"][:]
    barcodes = [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in barcodes
    ]
    return pd.Index(barcodes.astype(str) if hasattr(barcodes, "astype") else barcodes)


def read_cells_metadata(cells_path: Path) -> pd.DataFrame:
    name = cells_path.name.lower()

    if name.endswith(".parquet.gz"):
        try:
            return pd.read_parquet(cells_path)
        except Exception:
            with gzip.open(cells_path, "rb") as f:
                raw = f.read()
            return pd.read_parquet(io.BytesIO(raw))

    if name.endswith(".parquet"):
        return pd.read_parquet(cells_path)

    if name.endswith(".csv.gz") or name.endswith(".csv"):
        return pd.read_csv(cells_path)

    raise ValueError(f"Unsupported metadata file: {cells_path}")


def find_sample_files():
    rows = []

    for sample_dir in sorted(RAW_DIR.glob("Br*")):
        if not sample_dir.is_dir():
            continue

        brnum = sample_dir.name

        matrix_files = sorted(sample_dir.glob("*cell_feature_matrix.h5"))
        cells_files = sorted(
            list(sample_dir.glob("*cells*.parquet*"))
            + list(sample_dir.glob("*cells*.csv*"))
        )

        for matrix_path in matrix_files:
            for cells_path in cells_files:
                rows.append(
                    {
                        "BrNum": brnum,
                        "matrix_path": matrix_path,
                        "cells_path": cells_path,
                    }
                )

    return pd.DataFrame(rows)


def main():
    print("=" * 100)
    print("Xenium matrix-vs-cells cell ID crossmatch")
    print("=" * 100)

    files = find_sample_files()

    if files.empty:
        raise RuntimeError(f"No Xenium files found under {RAW_DIR}")

    print("\nFound sample file pairs:")
    print(files[["BrNum", "matrix_path", "cells_path"]].to_string(index=False))

    # Load all matrix barcodes.
    matrix_ids = {}
    matrix_sizes = {}

    for _, row in files.drop_duplicates("matrix_path").iterrows():
        matrix_path = Path(row["matrix_path"])
        matrix_brnum = Path(matrix_path).parent.name

        ids = read_h5_barcodes(matrix_path)
        matrix_ids[matrix_brnum] = set(ids.astype(str))
        matrix_sizes[matrix_brnum] = len(ids)

        print(f"Loaded matrix {matrix_brnum}: {len(ids)} barcodes")

    # Load all metadata cell IDs.
    metadata_ids = {}
    metadata_sizes = {}
    metadata_paths = {}

    for _, row in files.drop_duplicates("cells_path").iterrows():
        cells_path = Path(row["cells_path"])
        metadata_brnum = Path(cells_path).parent.name

        try:
            cells = read_cells_metadata(cells_path)
        except Exception as e:
            print(f"FAILED reading metadata {metadata_brnum}: {cells_path}")
            print(e)
            continue

        if "cell_id" not in cells.columns:
            print(f"WARNING: no cell_id column in {cells_path}")
            print("Columns:", cells.columns.tolist())
            continue

        ids = set(cells["cell_id"].astype(str))
        metadata_ids[metadata_brnum] = ids
        metadata_sizes[metadata_brnum] = len(ids)
        metadata_paths[metadata_brnum] = str(cells_path)

        print(f"Loaded metadata {metadata_brnum}: {len(ids)} unique cell_ids")

    # Crossmatch every matrix against every metadata file.
    rows = []

    for matrix_brnum, h5_set in matrix_ids.items():
        for metadata_brnum, meta_set in metadata_ids.items():
            exact = len(h5_set & meta_set)
            frac_matrix = exact / len(h5_set) if len(h5_set) else 0

            rows.append(
                {
                    "matrix_BrNum": matrix_brnum,
                    "metadata_BrNum": metadata_brnum,
                    "matrix_cells": matrix_sizes[matrix_brnum],
                    "metadata_unique_cell_ids": metadata_sizes[metadata_brnum],
                    "exact_matches": exact,
                    "fraction_of_matrix_matched": frac_matrix,
                    "metadata_path": metadata_paths[metadata_brnum],
                }
            )

    report = pd.DataFrame(rows)

    report = report.sort_values(
        ["matrix_BrNum", "exact_matches"],
        ascending=[True, False],
    )

    long_csv = OUT_DIR / "xenium_matrix_metadata_crossmatch_long.csv"
    report.to_csv(long_csv, index=False)

    pivot = report.pivot(
        index="matrix_BrNum",
        columns="metadata_BrNum",
        values="exact_matches",
    ).fillna(0).astype(int)

    pivot_csv = OUT_DIR / "xenium_matrix_metadata_crossmatch_pivot.csv"
    pivot.to_csv(pivot_csv)

    print("\nBest metadata match for each matrix:")
    best = report.sort_values(
        ["matrix_BrNum", "exact_matches"],
        ascending=[True, False],
    ).groupby("matrix_BrNum").head(1)

    print(
        best[
            [
                "matrix_BrNum",
                "metadata_BrNum",
                "matrix_cells",
                "metadata_unique_cell_ids",
                "exact_matches",
                "fraction_of_matrix_matched",
                "metadata_path",
            ]
        ].to_string(index=False)
    )

    print("\nSaved:")
    print(long_csv)
    print(pivot_csv)

    print("\nBr6432-specific results:")
    print(
        report[report["matrix_BrNum"] == "Br6432"][
            [
                "matrix_BrNum",
                "metadata_BrNum",
                "matrix_cells",
                "metadata_unique_cell_ids",
                "exact_matches",
                "fraction_of_matrix_matched",
                "metadata_path",
            ]
        ].sort_values("exact_matches", ascending=False).head(10).to_string(index=False)
    )


if __name__ == "__main__":
    main()
