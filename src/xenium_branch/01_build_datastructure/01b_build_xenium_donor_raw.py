#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]

ORIGINAL_SCRIPT = (
    PROJECT_ROOT
    / "src/xenium_branch/01_build_datastructure/01_build_xenium_N24_raw.py"
)

RAW_ROOT = PROJECT_ROOT / "data/raw/xenium"

DONOR_OUT_DIR = (
    PROJECT_ROOT
    / "data/processed/xenium/raw_by_donor"
)

SUMMARY_DIR = (
    PROJECT_ROOT
    / "outputs/xenium_branch/01_build_datastructure/by_donor"
)

DONOR_OUT_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_DIR.mkdir(parents=True, exist_ok=True)


def section(title):
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def load_original_module():
    spec = importlib.util.spec_from_file_location(
        "xenium_n24_original",
        ORIGINAL_SCRIPT,
    )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    return module


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample",
        required=True,
        help="Single donor, e.g. Br1113",
    )

    args = parser.parse_args()

    donor = args.sample

    section(f"BUILD SINGLE XENIUM DONOR: {donor}")

    mod = load_original_module()

    donor_dir = RAW_ROOT / donor

    if not donor_dir.exists():
        raise FileNotFoundError(
            f"Missing donor directory: {donor_dir}"
        )

    matrix_path = mod.find_first(
        donor_dir,
        [
            "cell_feature_matrix.h5",
            "*cell_feature_matrix.h5",
        ],
    )

    cells_path = mod.find_first(
        donor_dir,
        [
            "cells.csv.gz",
            "*cells.csv.gz",
            "cells.parquet.gz",
            "*cells.parquet.gz",
            "cells.parquet",
            "*cells.parquet",
            "cells.csv",
            "*cells.csv",
        ],
    )

    cell_boundaries_path = mod.find_first(
        donor_dir,
        [
            "cell_boundaries.csv.gz",
            "*cell_boundaries.csv.gz",
            "cell_boundaries.csv",
            "*cell_boundaries.csv",
        ],
    )

    nucleus_boundaries_path = mod.find_first(
        donor_dir,
        [
            "nucleus_boundaries.csv.gz",
            "*nucleus_boundaries.csv.gz",
            "nucleus_boundaries.csv",
            "*nucleus_boundaries.csv",
        ],
    )

    if not matrix_path:
        raise FileNotFoundError(
            f"No cell_feature_matrix.h5 for {donor}"
        )

    if not cells_path:
        raise FileNotFoundError(
            f"No cells metadata for {donor}"
        )

    row = pd.Series(
        {
            "BrNum": donor,
            "path": str(donor_dir),
            "matrix_path": matrix_path,
            "cells_path": cells_path,
            "cell_boundaries_path": cell_boundaries_path,
            "nucleus_boundaries_path": nucleus_boundaries_path,
            "has_matrix": True,
            "has_cells": True,
            "has_cell_boundaries": bool(cell_boundaries_path),
            "has_nucleus_boundaries": bool(nucleus_boundaries_path),
        }
    )

    adata, summary = mod.load_one_xenium_sample(row)

    # ---------------------------------------------------------
    # Critical validation
    # ---------------------------------------------------------

    observed_donors = set(
        adata.obs["BrNum"].astype(str).unique()
    )

    if observed_donors != {donor}:
        raise RuntimeError(
            f"Expected only {donor}; found {observed_donors}"
        )

    if not adata.obs_names.is_unique:
        raise RuntimeError(
            f"Non-unique cell IDs in {donor}"
        )

    if "counts" not in adata.layers:
        raise RuntimeError(
            f"Missing counts layer in {donor}"
        )

    if "spatial" not in adata.obsm:
        raise RuntimeError(
            f"Missing spatial coordinates in {donor}"
        )

    # ---------------------------------------------------------
    # Save donor object
    # ---------------------------------------------------------

    out_h5ad = DONOR_OUT_DIR / f"{donor}.h5ad"

    section("SAVE DONOR")

    print(f"Output:   {out_h5ad}")
    print(f"Cells:    {adata.n_obs}")
    print(f"Features: {adata.n_vars}")

    adata.write_h5ad(
        out_h5ad,
        compression="gzip",
    )

    summary_df = pd.DataFrame([summary])

    summary_csv = (
        SUMMARY_DIR
        / f"{donor}_summary.csv"
    )

    summary_df.to_csv(
        summary_csv,
        index=False,
    )

    section("PASS")

    print(f"Donor:    {donor}")
    print(f"Cells:    {adata.n_obs}")
    print(f"Features: {adata.n_vars}")
    print(f"H5AD:     {out_h5ad}")


if __name__ == "__main__":
    main()
