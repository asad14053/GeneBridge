#!/usr/bin/env python3

"""
00_prepare_imputation_inputs.py

Task:
    Validate Xenium and Visium matched N24 imputation inputs.

Inputs:
    data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad
    data/processed/visium/visium_N24_matched_layer_annotated.h5ad

Outputs:
    data/processed/imputation_inputs/shared_genes_xenium_visium.csv
    data/processed/imputation_inputs/xenium_only_genes.csv
    data/processed/imputation_inputs/visium_only_candidate_imputation_genes.csv
    data/processed/imputation_inputs/matched_N24_donor_list.csv
    data/processed/imputation_inputs/imputation_input_summary.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_XENIUM_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "xenium"
    / "xenium_N24_layer_celltype_annotated.h5ad"
)

DEFAULT_VISIUM_H5AD = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "visium"
    / "visium_N24_matched_layer_annotated.h5ad"
)

DEFAULT_OUTDIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "imputation_inputs"
)


def section(title: str) -> None:
    print("\n" + "=" * 100)
    print(title)
    print("=" * 100)


def clean_gene_index(values: pd.Index | pd.Series) -> pd.Index:
    genes = pd.Index(values.astype(str))
    genes = pd.Index(
        [
            g.strip()
            for g in genes
            if g is not None
            and str(g).strip() != ""
            and str(g).strip().lower() not in {"nan", "none", "null", "na"}
        ]
    )
    return genes


def get_gene_symbols(adata: ad.AnnData, preferred_col: str | None) -> pd.Index:
    if preferred_col is not None and preferred_col in adata.var.columns:
        return clean_gene_index(adata.var[preferred_col])
    return clean_gene_index(pd.Index(adata.var_names.astype(str)))


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--xenium-h5ad",
        default=str(DEFAULT_XENIUM_H5AD),
        help="Annotated Xenium h5ad input.",
    )

    parser.add_argument(
        "--visium-h5ad",
        default=str(DEFAULT_VISIUM_H5AD),
        help="Matched annotated Visium h5ad input.",
    )

    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTDIR),
        help="Output directory for imputation input tables.",
    )

    args = parser.parse_args()

    xenium_h5ad = Path(args.xenium_h5ad)
    visium_h5ad = Path(args.visium_h5ad)
    outdir = Path(args.outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    section("Input files")

    print("Xenium:", xenium_h5ad)
    print("Visium:", visium_h5ad)
    print("Output directory:", outdir)

    if not xenium_h5ad.exists():
        raise FileNotFoundError(f"Xenium h5ad not found: {xenium_h5ad}")

    if not visium_h5ad.exists():
        raise FileNotFoundError(f"Visium h5ad not found: {visium_h5ad}")

    section("Loading AnnData objects")

    # backed="r" avoids loading full expression matrices into memory.
    xen = ad.read_h5ad(xenium_h5ad, backed="r")
    vis = ad.read_h5ad(visium_h5ad, backed="r")

    section("Xenium")

    print(xen)
    print("Xenium obs columns:")
    print(xen.obs.columns.tolist())
    print("\nXenium var columns:")
    print(xen.var.columns.tolist())

    section("Visium")

    print(vis)
    print("Visium obs columns:")
    print(vis.obs.columns.tolist())
    print("\nVisium var columns:")
    print(vis.var.columns.tolist())

    section("Required metadata check")

    required_xen_obs = [
        "BrNum",
        "layer_annotation",
        "cell_type_annotation",
    ]

    required_vis_obs = [
        "BrNum_matched",
        "visium_layer_annotation",
    ]

    missing_xen = [c for c in required_xen_obs if c not in xen.obs.columns]
    missing_vis = [c for c in required_vis_obs if c not in vis.obs.columns]

    print("Missing Xenium columns:", missing_xen)
    print("Missing Visium columns:", missing_vis)

    if missing_xen:
        raise ValueError(f"Missing required Xenium obs columns: {missing_xen}")

    if missing_vis:
        raise ValueError(f"Missing required Visium obs columns: {missing_vis}")

    section("Donor matching")

    xen_donors = set(xen.obs["BrNum"].astype(str))
    vis_donors = set(vis.obs["BrNum_matched"].astype(str))

    shared_donors = sorted(xen_donors & vis_donors)
    xen_missing_vis = sorted(xen_donors - vis_donors)
    vis_missing_xen = sorted(vis_donors - xen_donors)

    print("Xenium donors:", len(xen_donors))
    print(sorted(xen_donors))

    print("\nVisium donors:", len(vis_donors))
    print(sorted(vis_donors))

    print("\nShared donors:", len(shared_donors))
    print(shared_donors)

    print("\nXenium missing in Visium:")
    print(xen_missing_vis)

    print("\nVisium missing in Xenium:")
    print(vis_missing_xen)

    if xen_missing_vis or vis_missing_xen:
        raise ValueError("Xenium and Visium donor sets do not match exactly.")

    pd.DataFrame({"BrNum": shared_donors}).to_csv(
        outdir / "matched_N24_donor_list.csv",
        index=False,
    )

    section("Gene matching")

    xen_gene_col = "Symbol" if "Symbol" in xen.var.columns else None
    vis_gene_col = "gene_name" if "gene_name" in vis.var.columns else None

    print("Xenium gene column used:", xen_gene_col if xen_gene_col else "var_names")
    print("Visium gene column used:", vis_gene_col if vis_gene_col else "var_names")

    xen_genes = get_gene_symbols(xen, xen_gene_col)
    vis_genes = get_gene_symbols(vis, vis_gene_col)

    xen_gene_set = set(xen_genes)
    vis_gene_set = set(vis_genes)

    shared_genes = sorted(xen_gene_set & vis_gene_set)
    xen_only_genes = sorted(xen_gene_set - vis_gene_set)
    vis_only_genes = sorted(vis_gene_set - xen_gene_set)

    print("Xenium genes:", len(xen_gene_set))
    print("Visium genes:", len(vis_gene_set))
    print("Shared genes:", len(shared_genes))
    print("Xenium-only genes:", len(xen_only_genes))
    print("Visium-only candidate imputation genes:", len(vis_only_genes))

    print("\nFirst 20 shared genes:")
    print(shared_genes[:20])

    pd.DataFrame({"gene": shared_genes}).to_csv(
        outdir / "shared_genes_xenium_visium.csv",
        index=False,
    )

    pd.DataFrame({"gene": xen_only_genes}).to_csv(
        outdir / "xenium_only_genes.csv",
        index=False,
    )

    pd.DataFrame({"gene": vis_only_genes}).to_csv(
        outdir / "visium_only_candidate_imputation_genes.csv",
        index=False,
    )

    section("Layer and cell-type summaries")

    xen_layer_counts = (
        xen.obs["layer_annotation"]
        .astype(str)
        .value_counts()
        .rename_axis("layer_annotation")
        .reset_index(name="n_xenium_cells")
    )

    xen_celltype_counts = (
        xen.obs["cell_type_annotation"]
        .astype(str)
        .value_counts()
        .rename_axis("cell_type_annotation")
        .reset_index(name="n_xenium_cells")
    )

    vis_layer_counts = (
        vis.obs["visium_layer_annotation"]
        .astype(str)
        .value_counts()
        .rename_axis("visium_layer_annotation")
        .reset_index(name="n_visium_spots")
    )

    xen_layer_counts.to_csv(outdir / "xenium_layer_counts.csv", index=False)
    xen_celltype_counts.to_csv(outdir / "xenium_celltype_counts.csv", index=False)
    vis_layer_counts.to_csv(outdir / "visium_layer_counts.csv", index=False)

    print("\nXenium layer counts:")
    print(xen_layer_counts)

    print("\nXenium cell-type counts:")
    print(xen_celltype_counts)

    print("\nVisium layer counts:")
    print(vis_layer_counts)

    section("Saving summary")

    summary = pd.DataFrame(
        {
            "metric": [
                "xenium_cells",
                "xenium_genes",
                "visium_spots",
                "visium_genes",
                "xenium_donors",
                "visium_donors",
                "shared_donors",
                "shared_genes",
                "xenium_only_genes",
                "visium_only_candidate_imputation_genes",
                "xenium_gene_column_used",
                "visium_gene_column_used",
            ],
            "value": [
                xen.n_obs,
                xen.n_vars,
                vis.n_obs,
                vis.n_vars,
                len(xen_donors),
                len(vis_donors),
                len(shared_donors),
                len(shared_genes),
                len(xen_only_genes),
                len(vis_only_genes),
                xen_gene_col if xen_gene_col else "var_names",
                vis_gene_col if vis_gene_col else "var_names",
            ],
        }
    )

    summary.to_csv(outdir / "imputation_input_summary.csv", index=False)

    print(summary)

    section("Done")

    print("Saved:")
    print(outdir / "shared_genes_xenium_visium.csv")
    print(outdir / "xenium_only_genes.csv")
    print(outdir / "visium_only_candidate_imputation_genes.csv")
    print(outdir / "matched_N24_donor_list.csv")
    print(outdir / "xenium_layer_counts.csv")
    print(outdir / "xenium_celltype_counts.csv")
    print(outdir / "visium_layer_counts.csv")
    print(outdir / "imputation_input_summary.csv")


if __name__ == "__main__":
    main()