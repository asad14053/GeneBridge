#!/usr/bin/env python3

from pathlib import Path
import re
import h5py
import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parents[3]

SAMPLE = "Br6432"

SAMPLE_DIR = PROJECT_ROOT / "data" / "raw" / "xenium" / SAMPLE

OUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "xenium_branch"
    / "01_build_datastructure"
    / "debug_br6432"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)


def section(title: str):
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def read_cells_file(path: Path) -> pd.DataFrame:
    path_str = str(path)

    if path_str.endswith(".parquet") or path_str.endswith(".parquet.gz"):
        return pd.read_parquet(path)

    if path_str.endswith(".csv") or path_str.endswith(".csv.gz"):
        return pd.read_csv(path)

    raise ValueError(f"Unsupported file type: {path}")


def normalize_id(x: str) -> str:
    x = str(x).strip()

    # common transformations
    x = re.sub(r"^Br[0-9]+_", "", x)
    x = re.sub(r"-1$", "", x)
    x = re.sub(r"\.0$", "", x)

    return x


def id_report(h5_ids, meta_ids, label: str) -> dict:
    h5_set = set(h5_ids)
    meta_set = set(meta_ids)

    exact_common = h5_set & meta_set

    h5_norm = pd.Index([normalize_id(x) for x in h5_ids])
    meta_norm = pd.Index([normalize_id(x) for x in meta_ids])

    norm_common = set(h5_norm) & set(meta_norm)

    return {
        "metadata_source": label,
        "h5_cells": len(h5_ids),
        "metadata_rows": len(meta_ids),
        "metadata_unique_ids": len(meta_set),
        "exact_matches": len(exact_common),
        "normalized_matches": len(norm_common),
        "h5_first_10": ";".join(map(str, list(h5_ids[:10]))),
        "metadata_first_10": ";".join(map(str, list(meta_ids[:10]))),
    }


section("Find files")

matrix_files = sorted(SAMPLE_DIR.glob("*cell_feature_matrix.h5"))
metadata_files = sorted(
    list(SAMPLE_DIR.glob("*cells*.parquet*"))
    + list(SAMPLE_DIR.glob("*cells*.csv*"))
)

print("Sample dir:", SAMPLE_DIR)

print("\nMatrix files:")
for f in matrix_files:
    print(f)

print("\nMetadata candidate files:")
for f in metadata_files:
    print(f)

if len(matrix_files) == 0:
    raise FileNotFoundError("No cell_feature_matrix.h5 found.")

if len(metadata_files) == 0:
    raise FileNotFoundError("No cells metadata file found.")

matrix_path = matrix_files[0]

section("Read H5 matrix using scanpy")

adata = sc.read_10x_h5(matrix_path, gex_only=False)
adata.obs_names = adata.obs_names.astype(str)

h5_ids = pd.Index(adata.obs_names.astype(str))

print(adata)
print("H5 matrix cells:", adata.n_obs)
print("H5 features:", adata.n_vars)
print("First 20 h5 obs_names:")
print(list(h5_ids[:20]))

section("Read raw H5 barcodes directly")

with h5py.File(matrix_path, "r") as f:
    raw_barcodes = f["matrix"]["barcodes"][:20]
    raw_barcodes = [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in raw_barcodes
    ]

print("First 20 raw H5 barcodes:")
print(raw_barcodes)

section("Compare against all metadata candidates")

reports = []

for meta_path in metadata_files:
    print("\n" + "-" * 100)
    print("Metadata file:", meta_path)

    try:
        cells = read_cells_file(meta_path)
    except Exception as e:
        print("Could not read:", e)
        continue

    print("Shape:", cells.shape)
    print("Columns:", cells.columns.tolist())

    # candidate ID columns
    candidate_cols = [
        c for c in cells.columns
        if ("cell" in c.lower()) or ("index" in c.lower()) or ("barcode" in c.lower())
    ]

    if "cell_id" in cells.columns and "cell_id" not in candidate_cols:
        candidate_cols.insert(0, "cell_id")

    print("Candidate ID columns:", candidate_cols)

    for col in candidate_cols:
        meta_ids = pd.Index(cells[col].astype(str))
        rep = id_report(h5_ids, meta_ids, f"{meta_path.name}::{col}")
        reports.append(rep)

        print(f"\nColumn: {col}")
        print("Metadata rows:", len(meta_ids))
        print("Unique metadata IDs:", meta_ids.nunique())
        print("Exact matches:", rep["exact_matches"])
        print("Normalized matches:", rep["normalized_matches"])
        print("First 10 metadata IDs:")
        print(list(meta_ids[:10]))

reports_df = pd.DataFrame(reports)
reports_df = reports_df.sort_values(
    ["exact_matches", "normalized_matches"],
    ascending=False,
)

out_csv = OUT_DIR / "br6432_cell_id_match_report.csv"
reports_df.to_csv(out_csv, index=False)

section("Best matching reports")
print(reports_df.head(24).to_string(index=False))

print("\nSaved report:")
print(out_csv)