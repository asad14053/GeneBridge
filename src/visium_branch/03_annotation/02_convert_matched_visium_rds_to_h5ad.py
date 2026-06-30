#!/usr/bin/env python3

"""
02_convert_matched_visium_rds_to_h5ad.py

Convert matched Visium SpatialExperiment RDS to AnnData h5ad.

Input:
    data/processed/visium/visium_N24_matched_layer_annotated.rds

Output:
    data/processed/visium/visium_N24_matched_layer_annotated.h5ad

Why this uses R internally:
    The input .rds is a Bioconductor SpatialExperiment object.
    Pure Python cannot reliably decode this object.
    This script calls R/zellkonverter from Python and then validates the h5ad in Python.
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]

DEFAULT_INPUT_RDS = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "visium"
    / "visium_N24_matched_layer_annotated.rds"
)

DEFAULT_OUTPUT_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "visium"
    / "visium_N24_matched_layer_annotated.h5ad"
)

DEFAULT_SUMMARY_CSV = (
    PROJECT_ROOT
    / "outputs"
    / "visium_branch"
    / "03_annotation"
    / "layerwise_matched_N24"
    / "tables"
    / "visium_N24_h5ad_conversion_summary.csv"
)


def section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def check_rscript_available() -> None:
    result = subprocess.run(
        ["which", "Rscript"],
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Rscript not found. Load R or activate an environment with R installed."
        )

    print("Rscript found:", result.stdout.strip())


def build_r_conversion_script(input_rds: Path, output_h5ad: Path) -> str:
    """
    R code embedded inside Python.

    It:
      1. reads SpatialExperiment RDS
      2. ensures logcounts exist
      3. copies spatialCoords into reducedDim(spatial)
      4. writes h5ad using zellkonverter
    """

    r_code = f"""
suppressPackageStartupMessages({{
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scuttle)
  library(zellkonverter)
}})

input_rds <- "{input_rds}"
output_h5ad <- "{output_h5ad}"

cat("\\nReading RDS:\\n")
cat(input_rds, "\\n")

if (!file.exists(input_rds)) {{
  stop("Input RDS not found: ", input_rds)
}}

spe <- readRDS(input_rds)

cat("\\nLoaded object class:\\n")
print(class(spe))

cat("\\nObject summary:\\n")
print(spe)

cat("\\nAssays before conversion:\\n")
print(assayNames(spe))

# Make sure column names and row names are safe.
if (is.null(colnames(spe))) {{
  colnames(spe) <- paste0("spot_", seq_len(ncol(spe)))
}}

if (is.null(rownames(spe))) {{
  rownames(spe) <- paste0("gene_", seq_len(nrow(spe)))
}}

colnames(spe) <- make.unique(as.character(colnames(spe)))
rownames(spe) <- make.unique(as.character(rownames(spe)))

# Make sure logcounts exists for AnnData X.
if (!("logcounts" %in% assayNames(spe))) {{
  cat("\\nlogcounts not found. Creating logcounts with scuttle::logNormCounts...\\n")
  spe <- scuttle::logNormCounts(spe)
}} else {{
  cat("\\nlogcounts already exists.\\n")
}}

# Put spatial coordinates into reducedDim(spe, "spatial")
# zellkonverter writes reducedDims to AnnData obsm.
coords <- tryCatch(
  spatialCoords(spe),
  error = function(e) NULL
)

if (!is.null(coords) && ncol(coords) >= 2) {{
  reducedDim(spe, "spatial") <- as.matrix(coords[, 1:2])
  cat("\\nAdded spatial coordinates to reducedDim(spe, 'spatial').\\n")
}} else {{
  cat("\\nWARNING: spatialCoords(spe) not found or has fewer than 2 columns.\\n")
}}

# Convert common problematic S4/DataFrame columns to simple vectors if needed.
# This prevents some h5ad-writing issues caused by complex metadata columns.
cd <- as.data.frame(colData(spe))
colData(spe) <- S4Vectors::DataFrame(cd)

rd <- as.data.frame(rowData(spe))
rowData(spe) <- S4Vectors::DataFrame(rd)

dir.create(dirname(output_h5ad), recursive = TRUE, showWarnings = FALSE)

cat("\\nWriting h5ad:\\n")
cat(output_h5ad, "\\n")

zellkonverter::writeH5AD(
  spe,
  file = output_h5ad,
  X_name = "logcounts"
)

cat("\\nDone writing h5ad.\\n")
"""
    return r_code


def run_r_conversion(input_rds: Path, output_h5ad: Path) -> None:
    r_code = build_r_conversion_script(input_rds, output_h5ad)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".R",
        delete=False,
    ) as tmp:
        tmp.write(r_code)
        tmp_r_path = Path(tmp.name)

    print("Temporary R conversion script:", tmp_r_path)

    try:
        cmd = ["Rscript", str(tmp_r_path)]
        result = subprocess.run(cmd, text=True)

        if result.returncode != 0:
            raise RuntimeError(
                f"R conversion failed with exit code {result.returncode}"
            )

    finally:
        try:
            tmp_r_path.unlink()
        except Exception:
            pass


def validate_h5ad(output_h5ad: Path, summary_csv: Path) -> None:
    section("Validating output h5ad in Python")

    if not output_h5ad.exists():
        raise FileNotFoundError(f"Output h5ad was not created: {output_h5ad}")

    adata = ad.read_h5ad(output_h5ad)

    print(adata)

    print("\nobs columns:")
    print(adata.obs.columns.tolist())

    print("\nvar columns:")
    print(adata.var.columns.tolist())

    print("\nobsm keys:")
    print(list(adata.obsm.keys()))

    required_obs = [
        "BrNum_matched",
        "visium_layer_annotation",
    ]

    rows = []

    rows.append({"metric": "n_spots", "value": adata.n_obs})
    rows.append({"metric": "n_genes", "value": adata.n_vars})

    for col in required_obs:
        exists = col in adata.obs.columns
        rows.append({"metric": f"has_{col}", "value": exists})
        print(f"{col}: {exists}")

    if "BrNum_matched" in adata.obs.columns:
        n_donors = adata.obs["BrNum_matched"].astype(str).nunique()
        rows.append({"metric": "n_BrNum_matched", "value": n_donors})

        print("\nDonors:")
        print(n_donors)
        print(sorted(adata.obs["BrNum_matched"].astype(str).unique()))

    if "visium_layer_annotation" in adata.obs.columns:
        print("\nLayer counts:")
        print(adata.obs["visium_layer_annotation"].astype(str).value_counts())

    rows.append({"metric": "has_spatial_obsm", "value": "spatial" in adata.obsm})

    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(summary_csv, index=False)

    print("\nSaved summary:")
    print(summary_csv)


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-rds",
        default=str(DEFAULT_INPUT_RDS),
        help="Input matched Visium RDS.",
    )

    parser.add_argument(
        "--output-h5ad",
        default=str(DEFAULT_OUTPUT_H5AD),
        help="Output matched Visium h5ad.",
    )

    parser.add_argument(
        "--summary-csv",
        default=str(DEFAULT_SUMMARY_CSV),
        help="Conversion validation summary CSV.",
    )

    args = parser.parse_args()

    input_rds = Path(args.input_rds)
    output_h5ad = Path(args.output_h5ad)
    summary_csv = Path(args.summary_csv)

    section("Convert matched Visium RDS to h5ad")

    print("Input RDS:", input_rds)
    print("Output h5ad:", output_h5ad)

    if not input_rds.exists():
        raise FileNotFoundError(f"Input RDS not found: {input_rds}")

    check_rscript_available()
    run_r_conversion(input_rds, output_h5ad)
    validate_h5ad(output_h5ad, summary_csv)

    section("Done")
    print("Final Visium h5ad:")
    print(output_h5ad)


if __name__ == "__main__":
    main()