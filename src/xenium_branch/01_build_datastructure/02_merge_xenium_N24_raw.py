#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]

DONOR_DIR = (
    PROJECT_ROOT
    / "data/processed/xenium/raw_by_donor"
)

OUT_H5AD = (
    PROJECT_ROOT
    / "data/processed/xenium/xenium_N24_raw.h5ad"
)

OUT_DIR = (
    PROJECT_ROOT
    / "outputs/xenium_branch/01_build_datastructure"
)

TABLE_DIR = OUT_DIR / "tables"

TABLE_DIR.mkdir(parents=True, exist_ok=True)


ORIGINAL_SCRIPT = (
    PROJECT_ROOT
    / "src/xenium_branch/01_build_datastructure/01_build_xenium_N24_raw.py"
)


DONORS = [
    "Br1113",
    "Br1139",
    "Br2039",
    "Br2421",
    "Br2719",
    "Br5314",
    "Br5373",
    "Br5400",
    "Br5436",
    "Br5588",
    "Br5590",
    "Br5622",
    "Br5639",
    "Br5746",
    "Br5931",
    "Br5973",
    "Br6032",
    "Br6389",
    "Br6432",
    "Br6437",
    "Br6496",
    "Br8433",
    "Br8667",
    "Br8772",
]


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

    section("VERIFY 24 DONOR FILES")

    missing = []

    for donor in DONORS:

        path = DONOR_DIR / f"{donor}.h5ad"

        if not path.exists():
            missing.append(str(path))
        else:
            print(f"FOUND: {path}")

    if missing:

        print()
        print("Missing donor files:")

        for x in missing:
            print(x)

        raise RuntimeError(
            f"{len(missing)} donor h5ads are missing."
        )

    # ---------------------------------------------------------
    # Load donors
    # ---------------------------------------------------------

    section("LOAD DONORS")

    adata_list = []
    summaries = []

    reference_var_names = None

    total_cells = 0

    for donor in DONORS:

        path = DONOR_DIR / f"{donor}.h5ad"

        print()
        print(f"Loading {donor}")
        print(path)

        adata = ad.read_h5ad(path)

        actual = set(
            adata.obs["BrNum"]
            .astype(str)
            .unique()
        )

        if actual != {donor}:
            raise RuntimeError(
                f"{donor}: BrNum mismatch: {actual}"
            )

        if not adata.obs_names.is_unique:
            raise RuntimeError(
                f"{donor}: duplicate obs_names"
            )

        if "counts" not in adata.layers:
            raise RuntimeError(
                f"{donor}: missing counts layer"
            )

        if "spatial" not in adata.obsm:
            raise RuntimeError(
                f"{donor}: missing spatial"
            )

        if reference_var_names is None:

            reference_var_names = (
                adata.var_names.astype(str).tolist()
            )

        else:

            current = (
                adata.var_names.astype(str).tolist()
            )

            if current != reference_var_names:

                raise RuntimeError(
                    f"{donor}: feature names/order "
                    "do not exactly match reference donor."
                )

        summaries.append(
            {
                "BrNum": donor,
                "n_cells": int(adata.n_obs),
                "n_features": int(adata.n_vars),
                "has_counts": "counts" in adata.layers,
                "has_spatial": "spatial" in adata.obsm,
            }
        )

        total_cells += adata.n_obs

        print(
            f"  cells={adata.n_obs:,}; "
            f"features={adata.n_vars}"
        )

        adata_list.append(adata)

    # ---------------------------------------------------------
    # Merge
    # ---------------------------------------------------------

    section("CONCATENATE 24 DONORS")

    print(
        f"Expected cells from donor sum: "
        f"{total_cells:,}"
    )

    adata_all = ad.concat(
        adata_list,
        axis=0,
        join="outer",
        merge="first",
        uns_merge="first",
        index_unique=None,
    )

    # ---------------------------------------------------------
    # Global metadata
    # ---------------------------------------------------------

    adata_all.uns["stage"] = (
        "01_build_datastructure"
    )

    adata_all.uns["object_name"] = (
        "xenium_N24_raw"
    )

    adata_all.uns["description"] = (
        "Raw Xenium AnnData built from "
        "24 donor-level h5ad files. "
        "Equivalent purpose to PI raw_spe_N24.RDS."
    )

    adata_all.uns["counts_layer"] = "counts"
    adata_all.uns["spatial_key"] = "spatial"

    # ---------------------------------------------------------
    # Final validation
    # ---------------------------------------------------------

    section("VALIDATE MERGED OBJECT")

    n_donors = (
        adata_all.obs["BrNum"]
        .astype(str)
        .nunique()
    )

    print(f"Cells:    {adata_all.n_obs:,}")
    print(f"Features: {adata_all.n_vars:,}")
    print(f"Donors:   {n_donors}")

    if adata_all.n_obs != total_cells:

        raise RuntimeError(
            "Merged cell count does not equal "
            "sum of donor cell counts."
        )

    if n_donors != 24:

        raise RuntimeError(
            f"Expected 24 donors, found {n_donors}"
        )

    if not adata_all.obs_names.is_unique:

        raise RuntimeError(
            "Merged object has duplicate cell IDs."
        )

    observed = set(
        adata_all.obs["BrNum"]
        .astype(str)
        .unique()
    )

    expected = set(DONORS)

    if observed != expected:

        raise RuntimeError(
            f"Donor mismatch.\n"
            f"Missing={expected-observed}\n"
            f"Extra={observed-expected}"
        )

    # ---------------------------------------------------------
    # Save
    # ---------------------------------------------------------

    section("SAVE FINAL N24")

    tmp_h5ad = OUT_H5AD.with_suffix(
        ".h5ad.tmp"
    )

    # h5py depends on recognizable extension,
    # therefore use a normal temporary .h5ad name.
    tmp_h5ad = (
        OUT_H5AD.parent
        / "xenium_N24_raw.tmp.h5ad"
    )

    adata_all.write_h5ad(
        tmp_h5ad,
        compression="gzip",
    )

    tmp_h5ad.replace(OUT_H5AD)

    print(f"Saved: {OUT_H5AD}")

    # ---------------------------------------------------------
    # Tables
    # ---------------------------------------------------------

    sample_summary = pd.DataFrame(
        summaries
    ).sort_values("BrNum")

    sample_summary.to_csv(
        TABLE_DIR
        / "xenium_N24_sample_summary.csv",
        index=False,
    )

    adata_all.obs.to_csv(
        TABLE_DIR
        / "xenium_N24_cell_metadata.csv.gz",
        compression="gzip",
    )

    adata_all.var.to_csv(
        TABLE_DIR
        / "xenium_N24_feature_metadata.csv",
    )

    build_status = sample_summary.copy()

    build_status["status"] = "ok"
    build_status["error"] = ""

    build_status.to_csv(
        TABLE_DIR
        / "xenium_N24_build_status.csv",
        index=False,
    )

    # ---------------------------------------------------------
    # Reuse original plotting function
    # ---------------------------------------------------------

    mod = load_original_module()

    mod.save_plots(
        adata_all,
        sample_summary,
    )

    section("DONE")

    print(f"Final h5ad: {OUT_H5AD}")
    print(f"Cells:      {adata_all.n_obs:,}")
    print(f"Features:   {adata_all.n_vars:,}")
    print(f"Donors:     {n_donors}")


if __name__ == "__main__":
    main()
