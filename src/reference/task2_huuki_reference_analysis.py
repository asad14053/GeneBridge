"""
Task 2: Huuki-Myers spatialDLPFC reference analysis in Python

Purpose
-------
Use the Huuki-Myers spatialDLPFC dataset as an accessible DLPFC reference for:
    1. Visium cortical layer / spatial-domain reference inspection
    2. snRNA-seq cell-type reference inspection
    3. marker overlap with your Xenium panel

Important
---------
The official processed Huuki-Myers data are distributed through the R/Bioconductor
package spatialLIBD. Therefore, this Python script can optionally create and run
a small R script to download the data and convert the R objects to .h5ad using
zellkonverter. After conversion, all inspection/summary analysis is done in Python.

Run from project root:
    python src/reference/task2_huuki_reference_analysis.py --download-convert

If you already have .h5ad files:
    python src/reference/task2_huuki_reference_analysis.py
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

REF_DIR = PROJECT_ROOT / "data" / "reference" / "huuki_spatialDLPFC"
OUT_DIR = PROJECT_ROOT / "outputs" / "task2_huuki_reference"
TABLE_DIR = OUT_DIR / "tables"

REF_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)


VISIUM_H5AD = REF_DIR / "huuki_spatialDLPFC_Visium.h5ad"
SNRNA_H5AD = REF_DIR / "huuki_spatialDLPFC_snRNAseq.h5ad"


# =============================================================================
# Markers used for your project-level checking
# =============================================================================

CELLTYPE_MARKERS = {
    "Excitatory_neuron": ["SLC17A7", "SATB2", "RORB", "PCP4"],
    "Inhibitory_neuron": ["GAD1", "GAD2", "PVALB", "SST", "VIP"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG"],
    "OPC": ["PDGFRA", "CSPG4"],
    "Microglia": ["CX3CR1", "P2RY12", "C3"],
    "Endothelial": ["CLDN5", "FLT1", "VWF"],
    "Mural": ["PDGFRB", "RGS5"],
}

LAYER_MARKERS = {
    "L1": ["RELN", "LAMP5", "CXCL14"],
    "L2_3": ["CUX2", "CALB1", "LINC00507", "SATB2"],
    "L4": ["RORB", "SCNN1A"],
    "L5": ["BCL11B", "FEZF2", "ETV1"],
    "L6": ["FOXP2", "TLE4", "THEMIS", "PCP4"],
    "WM": ["MBP", "PLP1", "MOG", "MOBP"],
}


# =============================================================================
# Utilities
# =============================================================================

def print_section(title: str) -> None:
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def require_python_package(package_name: str, import_name: Optional[str] = None):
    if import_name is None:
        import_name = package_name

    try:
        return __import__(import_name)
    except ImportError as exc:
        raise ImportError(
            f"Missing Python package: {package_name}\n"
            f"Install it with:\n"
            f"    pip install {package_name}\n"
            f"or:\n"
            f"    conda install -c conda-forge {package_name}"
        ) from exc


def safe_to_csv(df: pd.DataFrame, path: Path, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=index)
    print(f"Saved: {path}")


def find_existing_xenium_h5() -> Optional[Path]:
    preferred = PROJECT_ROOT / "data" / "raw" / "xenium" / "Br2039"
    if preferred.exists():
        files = sorted(preferred.glob("*cell_feature_matrix.h5"))
        if files:
            return files[0]

    xenium_root = PROJECT_ROOT / "data" / "raw" / "xenium"
    if xenium_root.exists():
        files = sorted(xenium_root.glob("**/*cell_feature_matrix.h5"))
        if files:
            return files[0]

    return None


def get_gene_symbols_from_anndata(adata) -> set[str]:
    genes = set(map(str, adata.var_names))

    for col in ["gene_name", "gene_search", "symbol", "gene_symbol", "name"]:
        if col in adata.var.columns:
            genes.update(adata.var[col].dropna().astype(str).tolist())

    return genes


def get_gene_symbols_from_xenium_h5(h5_path: Path) -> set[str]:
    sc = require_python_package("scanpy", "scanpy")
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    return set(map(str, adata.var_names))


# =============================================================================
# Step 1: Download and convert with R, controlled by Python
# =============================================================================

def write_r_download_convert_script(r_script_path: Path) -> None:
    ref_dir_r = str(REF_DIR).replace("\\", "/")

    r_code = f'''
message("Task 2 Huuki-Myers download + conversion")

ref_dir <- normalizePath("{ref_dir_r}", mustWork = FALSE)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)

message("Reference directory: ", ref_dir)

if (!requireNamespace("BiocManager", quietly = TRUE)) {{
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}}

packages <- c(
    "spatialLIBD",
    "SpatialExperiment",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "zellkonverter"
)

for (pkg in packages) {{
    if (!requireNamespace(pkg, quietly = TRUE)) {{
        message("Installing missing Bioconductor package: ", pkg)
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }}
}}

suppressPackageStartupMessages(library(spatialLIBD))
suppressPackageStartupMessages(library(SpatialExperiment))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(zellkonverter))

message("Downloading spatialDLPFC Visium...")
spe <- spatialLIBD::fetch_data(
    type = "spatialDLPFC_Visium",
    destdir = ref_dir
)

message("Downloading spatialDLPFC snRNA-seq...")
sce <- spatialLIBD::fetch_data(
    type = "spatialDLPFC_snRNAseq",
    destdir = ref_dir
)

message("Saving RDS copies...")
saveRDS(spe, file.path(ref_dir, "huuki_spatialDLPFC_Visium.rds"))
saveRDS(sce, file.path(ref_dir, "huuki_spatialDLPFC_snRNAseq.rds"))

message("Writing metadata CSV files...")
write.csv(
    as.data.frame(SummarizedExperiment::colData(spe)),
    file.path(ref_dir, "huuki_visium_coldata.csv"),
    row.names = TRUE
)

write.csv(
    as.data.frame(SummarizedExperiment::rowData(spe)),
    file.path(ref_dir, "huuki_visium_rowdata.csv"),
    row.names = TRUE
)

write.csv(
    as.data.frame(SummarizedExperiment::colData(sce)),
    file.path(ref_dir, "huuki_snRNAseq_coldata.csv"),
    row.names = TRUE
)

write.csv(
    as.data.frame(SummarizedExperiment::rowData(sce)),
    file.path(ref_dir, "huuki_snRNAseq_rowdata.csv"),
    row.names = TRUE
)

spatial_coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
write.csv(
    spatial_coords,
    file.path(ref_dir, "huuki_visium_spatial_coords.csv"),
    row.names = TRUE
)

message("Converting Visium to h5ad...")
zellkonverter::writeH5AD(
    spe,
    file = file.path(ref_dir, "huuki_spatialDLPFC_Visium.h5ad"),
    X_name = "logcounts"
)

message("Converting snRNA-seq to h5ad...")
zellkonverter::writeH5AD(
    sce,
    file = file.path(ref_dir, "huuki_spatialDLPFC_snRNAseq.h5ad"),
    X_name = "logcounts"
)

message("Done.")
'''

    r_script_path.write_text(r_code, encoding="utf-8")
    print(f"Saved R helper script: {r_script_path}")


def run_r_download_convert() -> None:
    print_section("Downloading and converting Huuki-Myers data through spatialLIBD")

    r_script_path = REF_DIR / "task2_download_convert_huuki_spatialDLPFC.R"
    write_r_download_convert_script(r_script_path)

    cmd = ["Rscript", str(r_script_path)]
    print("Running:")
    print(" ".join(cmd))

    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Could not find Rscript on your system PATH.\n\n"
            "Fix options:\n"
            "1. Install R.\n"
            "2. Add R/bin to PATH.\n"
            "3. Run this on HPC/cluster where R is available.\n\n"
            "Then rerun:\n"
            "    python src/reference/task2_huuki_reference_analysis.py --download-convert"
        ) from exc


# =============================================================================
# Step 2: Inspect converted AnnData objects
# =============================================================================

def load_h5ad(path: Path):
    sc = require_python_package("scanpy", "scanpy")

    if not path.exists():
        raise FileNotFoundError(
            f"Missing file: {path}\n\n"
            "Run first:\n"
            "    python src/reference/task2_huuki_reference_analysis.py --download-convert"
        )

    print(f"Loading: {path}")
    return sc.read_h5ad(path)


def save_metadata_columns(adata, name: str) -> pd.DataFrame:
    obs_cols = pd.DataFrame({
        "object": name,
        "axis": "obs",
        "column": list(adata.obs.columns),
    })

    var_cols = pd.DataFrame({
        "object": name,
        "axis": "var",
        "column": list(adata.var.columns),
    })

    out = pd.concat([obs_cols, var_cols], ignore_index=True)
    out_path = TABLE_DIR / f"huuki_{name}_metadata_columns.csv"
    safe_to_csv(out, out_path)

    return out


def summarize_adata(adata, name: str) -> dict:
    genes = get_gene_symbols_from_anndata(adata)

    return {
        "object": name,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "n_obs_metadata_columns": adata.obs.shape[1],
        "n_var_metadata_columns": adata.var.shape[1],
        "n_gene_symbols_detected": len(genes),
    }


def candidate_columns(adata, candidates: list[str], contains_any: list[str]) -> list[str]:
    cols = list(adata.obs.columns)
    lower_map = {c.lower(): c for c in cols}

    found = []

    for c in candidates:
        if c.lower() in lower_map:
            found.append(lower_map[c.lower()])

    for c in cols:
        cl = c.lower()
        if any(k.lower() in cl for k in contains_any):
            if c not in found:
                found.append(c)

    return found


def save_candidate_annotation_counts(adata, name: str) -> None:
    if name == "visium":
        candidates = [
            "sample_id",
            "BrNum",
            "position",
            "spatialLIBD",
            "layer_guess",
            "layer_guess_reordered",
            "layer",
            "ManualAnnotation",
            "BayesSpace_harmony_09",
            "BayesSpace_harmony_16",
            "BayesSpace_harmony_28",
        ]
        contains_any = [
            "sample",
            "brnum",
            "layer",
            "bayesspace",
            "domain",
            "spatial",
            "position",
        ]
    else:
        candidates = [
            "sample_id",
            "BrNum",
            "cellType",
            "cell_type",
            "celltype",
            "broad_cell_type",
            "cluster",
            "label",
            "annotation",
            "layer",
        ]
        contains_any = [
            "sample",
            "brnum",
            "cell",
            "type",
            "cluster",
            "label",
            "annot",
            "layer",
        ]

    cols = candidate_columns(adata, candidates, contains_any)

    rows = []

    for col in cols:
        try:
            counts = adata.obs[col].astype(str).value_counts(dropna=False).head(50)
            for value, n in counts.items():
                rows.append({
                    "object": name,
                    "column": col,
                    "value": value,
                    "n": int(n),
                })
        except Exception as exc:
            rows.append({
                "object": name,
                "column": col,
                "value": f"ERROR: {exc}",
                "n": np.nan,
            })

    out = pd.DataFrame(rows)

    if out.empty:
        out = pd.DataFrame([{
            "object": name,
            "column": "NO_CANDIDATE_COLUMN_FOUND",
            "value": "Inspect metadata columns manually",
            "n": np.nan,
        }])

    out_path = TABLE_DIR / f"huuki_{name}_candidate_annotation_counts.csv"
    safe_to_csv(out, out_path)


# =============================================================================
# Step 3: Marker overlap analysis
# =============================================================================

def marker_overlap_table(reference_genes: set[str], xenium_genes: Optional[set[str]]) -> pd.DataFrame:
    rows = []

    for marker_type, marker_dict in [
        ("celltype", CELLTYPE_MARKERS),
        ("layer", LAYER_MARKERS),
    ]:
        for group, genes in marker_dict.items():
            for gene in genes:
                rows.append({
                    "marker_type": marker_type,
                    "group": group,
                    "gene": gene,
                    "present_in_huuki_reference": gene in reference_genes,
                    "present_in_xenium_panel": (gene in xenium_genes) if xenium_genes is not None else np.nan,
                })

    return pd.DataFrame(rows)


def run_marker_overlap(visium_adata, sn_adata) -> None:
    print_section("Marker overlap: Huuki reference vs Xenium panel")

    huuki_genes = set()
    huuki_genes.update(get_gene_symbols_from_anndata(visium_adata))
    huuki_genes.update(get_gene_symbols_from_anndata(sn_adata))

    xenium_h5 = find_existing_xenium_h5()

    xenium_genes = None
    if xenium_h5 is not None:
        print(f"Found Xenium H5: {xenium_h5}")
        xenium_genes = get_gene_symbols_from_xenium_h5(xenium_h5)
    else:
        print("No local Xenium cell_feature_matrix.h5 found. Skipping Xenium panel overlap.")

    overlap = marker_overlap_table(huuki_genes, xenium_genes)

    out_path = TABLE_DIR / "xenium_panel_overlap_predefined_markers.csv"
    safe_to_csv(overlap, out_path)

    summary = (
        overlap
        .groupby(["marker_type", "group"], dropna=False)
        .agg(
            n_markers=("gene", "count"),
            n_present_in_huuki=("present_in_huuki_reference", "sum"),
            n_present_in_xenium=("present_in_xenium_panel", "sum"),
        )
        .reset_index()
    )

    out_summary = TABLE_DIR / "xenium_panel_overlap_predefined_markers_summary.csv"
    safe_to_csv(summary, out_summary)


# =============================================================================
# Optional: derive snRNA markers using Scanpy
# =============================================================================

def find_best_celltype_column(sn_adata) -> Optional[str]:
    preferred = [
        "cellType",
        "cell_type",
        "celltype",
        "broad_cell_type",
        "annotation",
        "cell_annotation",
        "cluster",
    ]

    cols = candidate_columns(
        sn_adata,
        candidates=preferred,
        contains_any=["celltype", "cell_type", "cell type", "annotation", "cluster"],
    )

    for col in cols:
        n_unique = sn_adata.obs[col].astype(str).nunique()
        if 2 <= n_unique <= 100:
            return col

    return None


def optional_sn_marker_analysis(sn_adata, max_cells_per_group: int = 1000) -> None:
    print_section("Optional snRNA-seq marker analysis")

    sc = require_python_package("scanpy", "scanpy")

    groupby = find_best_celltype_column(sn_adata)

    if groupby is None:
        print("No usable cell-type/group column detected. Skipping marker analysis.")
        return

    print(f"Using groupby column: {groupby}")

    adata = sn_adata.copy()

    selected_indices = []

    for group, idx in adata.obs.groupby(groupby).indices.items():
        idx = np.array(idx)
        if len(idx) > max_cells_per_group:
            idx = np.random.choice(idx, size=max_cells_per_group, replace=False)
        selected_indices.extend(idx.tolist())

    adata = adata[selected_indices].copy()
    print(f"Downsampled object: {adata.n_obs} cells/nuclei x {adata.n_vars} genes")

    try:
        max_val = adata.X.max()
        if hasattr(max_val, "compute"):
            max_val = max_val.compute()
        max_val = float(max_val)
    except Exception:
        max_val = 100.0

    if max_val > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method="wilcoxon",
        n_genes=50,
    )

    marker_df = sc.get.rank_genes_groups_df(adata, group=None)

    out_path = TABLE_DIR / "huuki_snRNAseq_scanpy_top_markers.csv"
    safe_to_csv(marker_df, out_path)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Task 2 Huuki-Myers spatialDLPFC reference analysis."
    )

    parser.add_argument(
        "--download-convert",
        action="store_true",
        help=(
            "Use Python to generate/run an R script that downloads spatialDLPFC "
            "with spatialLIBD and converts objects to h5ad."
        ),
    )

    parser.add_argument(
        "--make-sn-markers",
        action="store_true",
        help="Optionally compute snRNA-seq top markers with Scanpy. This may be slow/heavy.",
    )

    parser.add_argument(
        "--max-cells-per-group",
        type=int,
        default=1000,
        help="Maximum nuclei/cells per group for optional marker analysis.",
    )

    args = parser.parse_args()

    print_section("Task 2: Huuki-Myers spatialDLPFC reference analysis")

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Reference dir: {REF_DIR}")
    print(f"Output table dir: {TABLE_DIR}")

    if args.download_convert:
        run_r_download_convert()

    print_section("Loading converted h5ad files")
    visium = load_h5ad(VISIUM_H5AD)
    sn = load_h5ad(SNRNA_H5AD)

    print(visium)
    print(sn)

    print_section("Saving metadata inspection tables")
    save_metadata_columns(visium, "visium")
    save_metadata_columns(sn, "snRNAseq")

    save_candidate_annotation_counts(visium, "visium")
    save_candidate_annotation_counts(sn, "snRNAseq")

    summary_rows = [
        summarize_adata(visium, "Huuki spatialDLPFC Visium"),
        summarize_adata(sn, "Huuki spatialDLPFC snRNAseq"),
    ]

    summary = pd.DataFrame(summary_rows)
    safe_to_csv(summary, TABLE_DIR / "task2_huuki_reference_summary.csv")

    run_marker_overlap(visium, sn)

    if args.make_sn_markers:
        optional_sn_marker_analysis(sn, max_cells_per_group=args.max_cells_per_group)

    print_section("Done")

    print("Main outputs:")
    print(TABLE_DIR / "task2_huuki_reference_summary.csv")
    print(TABLE_DIR / "huuki_visium_metadata_columns.csv")
    print(TABLE_DIR / "huuki_snRNAseq_metadata_columns.csv")
    print(TABLE_DIR / "huuki_visium_candidate_annotation_counts.csv")
    print(TABLE_DIR / "huuki_snRNAseq_candidate_annotation_counts.csv")
    print(TABLE_DIR / "xenium_panel_overlap_predefined_markers.csv")

    print("\nMeeting sentence:")
    print(
        "I used the accessible Huuki-Myers spatialDLPFC reference as the first Task 2 reference. "
        "The workflow converts the official spatialLIBD Visium and snRNA-seq objects to h5ad, "
        "inspects metadata columns for layer/domain/cell-type labels, and checks which DLPFC "
        "cell-type/layer markers overlap with our Xenium panel."
    )


if __name__ == "__main__":
    main()
