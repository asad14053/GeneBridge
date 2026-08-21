#!/usr/bin/env python

"""
03_align_vista_Br8667_genes.py

Purpose:
    Align Huuki snRNA-seq Br8667 and Huuki Visium Br8667 by shared genes
    for VISTA beta imputation.

Inputs:
    seq_data:
        data/processed/imputation_beta/Br8667/seq_data_huuki_snrna_Br8667_full.h5ad

    spatial_data:
        data/processed/imputation_beta/Br8667/spatial_data_huuki_visium_Br8667.h5ad

Outputs:
    data/processed/imputation_beta/Br8667/seq_data_huuki_snrna_Br8667_shared.h5ad
    data/processed/imputation_beta/Br8667/spatial_data_huuki_visium_Br8667_shared.h5ad
    data/processed/imputation_beta/Br8667/Br8667_vista_shared_genes.csv
    data/processed/imputation_beta/Br8667/Br8667_vista_alignment_summary.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

DEFAULT_SEQ_PATH = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/seq_data_huuki_snrna_Br8667_full.h5ad"
)

DEFAULT_SPATIAL_PATH = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667/spatial_data_huuki_visium_Br8667.h5ad"
)

DEFAULT_OUT_DIR = PROJECT_ROOT / "data/processed/imputation_beta/Br8667"


def print_file_created(label: str, path: Path) -> None:
    """
    Print file path and size after creating an output file.
    """
    path = Path(path)

    if path.exists():
        size_mb = path.stat().st_size / (1024 ** 2)
        print(f"\n✅ Created {label}:")
        print(f"   {path}")
        print(f"   Size: {size_mb:.2f} MB")
    else:
        print(f"\n❌ Expected {label}, but file was not found:")
        print(f"   {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align Huuki snRNA-seq and Huuki Visium Br8667 genes for VISTA."
    )

    parser.add_argument(
        "--seq_path",
        type=Path,
        default=DEFAULT_SEQ_PATH,
        help="Input snRNA-seq h5ad path.",
    )

    parser.add_argument(
        "--spatial_path",
        type=Path,
        default=DEFAULT_SPATIAL_PATH,
        help="Input Visium spatial h5ad path.",
    )

    parser.add_argument(
        "--out_dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory.",
    )

    parser.add_argument(
        "--max_genes",
        default="all",
        help="Use all shared genes by default. Optionally provide a number, e.g. --max_genes 5000.",
    )

    parser.add_argument(
        "--compression",
        default="gzip",
        choices=["gzip", "none"],
        help="Compression for output h5ad files. Default: gzip.",
    )

    return parser.parse_args()


def get_gene_symbols(adata, dataset_name: str) -> pd.Series:
    """
    Return one gene-symbol-like value per variable.

    Priority:
        1. gene_symbol
        2. gene_name
        3. Symbol
        4. symbol
        5. var_names
    """

    preferred_cols = ["gene_symbol", "gene_name", "Symbol", "symbol"]

    for col in preferred_cols:
        if col in adata.var.columns:
            genes = adata.var[col].astype(str).str.strip()
            print(f"{dataset_name}: using adata.var['{col}'] as gene symbols")
            return genes

    print(f"{dataset_name}: using adata.var_names as gene symbols")
    return pd.Series(adata.var_names.astype(str), index=adata.var_names).str.strip()


def make_gene_position_table(adata, dataset_name: str) -> pd.DataFrame:
    """
    Create table:
        gene_symbol, var_position, original_var_name

    Removes invalid gene names and duplicate gene symbols.
    For duplicate gene symbols, keeps the first occurrence.
    """

    genes = get_gene_symbols(adata, dataset_name).reset_index(drop=True)

    table = pd.DataFrame(
        {
            "gene_symbol": genes.astype(str).str.strip(),
            "var_position": range(adata.n_vars),
            "original_var_name": adata.var_names.astype(str),
        }
    )

    invalid = table["gene_symbol"].isin(["", "nan", "None", "NA", "NaN"])
    table = table.loc[~invalid].copy()

    n_before = len(table)
    table = table.drop_duplicates(subset="gene_symbol", keep="first").copy()
    n_after = len(table)

    print(f"\n{dataset_name} gene table:")
    print(f"  original n_vars: {adata.n_vars}")
    print(f"  valid genes: {n_before}")
    print(f"  duplicate symbols removed: {n_before - n_after}")
    print(f"  unique genes kept: {n_after}")

    return table


def select_shared_genes(
    seq_table: pd.DataFrame,
    spatial_table: pd.DataFrame,
    seq_adata,
    max_genes: str,
) -> list[str]:
    shared = sorted(
        set(seq_table["gene_symbol"])
        & set(spatial_table["gene_symbol"])
    )

    print(f"\nShared genes before optional filtering: {len(shared)}")

    if len(shared) == 0:
        raise ValueError("No shared genes found between seq_data and spatial_data.")

    if max_genes.lower() in ["all", "full", "none"]:
        return shared

    try:
        max_n = int(max_genes)
    except ValueError as exc:
        raise ValueError(
            f"Invalid --max_genes value: {max_genes}. Use 'all' or an integer."
        ) from exc

    if max_n <= 0:
        raise ValueError("--max_genes must be positive or 'all'.")

    if len(shared) <= max_n:
        print(f"Shared genes <= max_genes. Keeping all {len(shared)} genes.")
        return shared

    seq_gene_col = get_gene_symbols(seq_adata, "seq_data_for_ranking").reset_index(drop=True)

    rank_df = pd.DataFrame(
        {
            "gene_symbol": seq_gene_col.astype(str).str.strip(),
            "var_position": range(seq_adata.n_vars),
        }
    )

    if "binomial_deviance" in seq_adata.var.columns:
        print(f"Selecting top {max_n} shared genes by seq_data binomial_deviance.")

        rank_df["binomial_deviance"] = pd.to_numeric(
            seq_adata.var["binomial_deviance"].values,
            errors="coerce",
        )

        rank_df = rank_df[rank_df["gene_symbol"].isin(shared)].copy()
        rank_df = rank_df.sort_values("binomial_deviance", ascending=False)

        selected = (
            rank_df["gene_symbol"]
            .drop_duplicates()
            .head(max_n)
            .tolist()
        )
    else:
        print(f"binomial_deviance not found. Selecting first {max_n} alphabetically.")
        selected = shared[:max_n]

    selected = sorted(selected)

    print(f"Shared genes after --max_genes {max_genes}: {len(selected)}")

    return selected


def write_h5ad_with_optional_compression(adata, path: Path, compression: str) -> None:
    """
    Write h5ad with optional gzip compression.
    """
    if path.exists():
        print("\nRemoving old output before rewriting:")
        print(path)
        path.unlink()

    if compression == "gzip":
        adata.write_h5ad(path, compression="gzip")
    else:
        adata.write_h5ad(path)


def main() -> None:
    args = parse_args()

    seq_path = args.seq_path
    spatial_path = args.spatial_path
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    seq_out = out_dir / "seq_data_huuki_snrna_Br8667_shared.h5ad"
    spatial_out = out_dir / "spatial_data_huuki_visium_Br8667_shared.h5ad"
    genes_out = out_dir / "Br8667_vista_shared_genes.csv"
    summary_out = out_dir / "Br8667_vista_alignment_summary.csv"

    print("=" * 100)
    print("VISTA Br8667 gene alignment")
    print("=" * 100)

    print("\nInput files:")
    print(f"seq_data:     {seq_path}")
    print(f"spatial_data: {spatial_path}")

    print("\nOutput directory:")
    print(out_dir)

    print("\nPlanned output files:")
    print(f"1. aligned seq_data h5ad:     {seq_out}")
    print(f"2. aligned spatial_data h5ad: {spatial_out}")
    print(f"3. shared genes CSV:          {genes_out}")
    print(f"4. alignment summary CSV:     {summary_out}")

    if not seq_path.exists():
        raise FileNotFoundError(f"Missing seq_data file: {seq_path}")

    if not spatial_path.exists():
        raise FileNotFoundError(f"Missing spatial_data file: {spatial_path}")

    print("=" * 100)
    print("Loading input files")
    print("=" * 100)

    seq = ad.read_h5ad(seq_path)
    spatial = ad.read_h5ad(spatial_path)

    print("\nseq_data:")
    print(seq)

    print("\nspatial_data:")
    print(spatial)

    print("\nseq_data layers:")
    print(list(seq.layers.keys()))

    print("\nspatial_data layers:")
    print(list(spatial.layers.keys()))

    if "spatial" not in spatial.obsm:
        print("\nWARNING: spatial_data does not contain obsm['spatial'].")
    else:
        print("\nspatial_data has obsm['spatial']: True")

    print("=" * 100)
    print("Building gene maps")
    print("=" * 100)

    seq_table = make_gene_position_table(seq, "seq_data")
    spatial_table = make_gene_position_table(spatial, "spatial_data")

    shared_genes = select_shared_genes(
        seq_table=seq_table,
        spatial_table=spatial_table,
        seq_adata=seq,
        max_genes=args.max_genes,
    )

    if len(shared_genes) < 500:
        print("\nFirst 20 seq genes:")
        print(seq_table["gene_symbol"].head(20).tolist())

        print("\nFirst 20 spatial genes:")
        print(spatial_table["gene_symbol"].head(20).tolist())

        raise ValueError(
            f"Too few shared genes: {len(shared_genes)}. Check gene symbol columns."
        )

    seq_pos_map = dict(zip(seq_table["gene_symbol"], seq_table["var_position"]))
    spatial_pos_map = dict(zip(spatial_table["gene_symbol"], spatial_table["var_position"]))

    seq_positions = [seq_pos_map[g] for g in shared_genes]
    spatial_positions = [spatial_pos_map[g] for g in shared_genes]

    print("=" * 100)
    print("Subsetting and aligning")
    print("=" * 100)

    seq_shared = seq[:, seq_positions].copy()
    spatial_shared = spatial[:, spatial_positions].copy()

    seq_shared.var_names = pd.Index(shared_genes, dtype="str")
    spatial_shared.var_names = pd.Index(shared_genes, dtype="str")

    seq_shared.var["matched_gene_symbol"] = shared_genes
    spatial_shared.var["matched_gene_symbol"] = shared_genes

    seq_shared.var.index.name = None
    seq_shared.obs.index.name = None
    spatial_shared.var.index.name = None
    spatial_shared.obs.index.name = None

    print("\nAligned seq_data:")
    print(seq_shared)

    print("\nAligned spatial_data:")
    print(spatial_shared)

    same_genes = seq_shared.var_names.equals(spatial_shared.var_names)

    print("\nSame genes:")
    print(same_genes)

    if not same_genes:
        raise ValueError("Gene alignment failed: seq and spatial var_names differ.")

    print("=" * 100)
    print("Saving outputs")
    print("=" * 100)

    print("\nWriting shared genes CSV...")
    shared_gene_df = pd.DataFrame(
        {
            "gene": shared_genes,
            "seq_var_position": seq_positions,
            "spatial_var_position": spatial_positions,
        }
    )
    shared_gene_df.to_csv(genes_out, index=False)
    print_file_created("shared genes CSV", genes_out)

    print("\nWriting alignment summary CSV...")
    summary = pd.DataFrame(
        [
            {"metric": "seq_input_path", "value": str(seq_path)},
            {"metric": "spatial_input_path", "value": str(spatial_path)},
            {"metric": "seq_input_cells", "value": seq.n_obs},
            {"metric": "seq_input_genes", "value": seq.n_vars},
            {"metric": "spatial_input_spots", "value": spatial.n_obs},
            {"metric": "spatial_input_genes", "value": spatial.n_vars},
            {"metric": "shared_genes", "value": len(shared_genes)},
            {"metric": "max_genes_arg", "value": args.max_genes},
            {"metric": "compression", "value": args.compression},
            {"metric": "seq_output_path", "value": str(seq_out)},
            {"metric": "spatial_output_path", "value": str(spatial_out)},
            {"metric": "shared_gene_csv", "value": str(genes_out)},
            {"metric": "alignment_summary_csv", "value": str(summary_out)},
            {"metric": "spatial_has_obsm_spatial", "value": "spatial" in spatial_shared.obsm},
        ]
    )
    summary.to_csv(summary_out, index=False)
    print_file_created("alignment summary CSV", summary_out)

    print("\nWriting aligned seq_data h5ad...")
    print(seq_out)
    write_h5ad_with_optional_compression(seq_shared, seq_out, args.compression)
    print_file_created("aligned seq_data h5ad", seq_out)

    print("\nWriting aligned spatial_data h5ad...")
    print(spatial_out)
    write_h5ad_with_optional_compression(spatial_shared, spatial_out, args.compression)
    print_file_created("aligned spatial_data h5ad", spatial_out)

    print("=" * 100)
    print("Reload validation")
    print("=" * 100)

    seq_check = ad.read_h5ad(seq_out)
    spatial_check = ad.read_h5ad(spatial_out)

    print("\nReloaded seq_data:")
    print(seq_check)

    print("\nReloaded spatial_data:")
    print(spatial_check)

    same_genes_after_reload = seq_check.var_names.equals(spatial_check.var_names)

    print("\nSame genes after reload:")
    print(same_genes_after_reload)

    print("\nNumber of shared genes:")
    print(seq_check.n_vars)

    print("\nspatial obsm keys:")
    print(list(spatial_check.obsm.keys()))

    if not same_genes_after_reload:
        raise ValueError("Reload validation failed: genes differ after writing.")

    print("=" * 100)
    print("Final output files created")
    print("=" * 100)

    print_file_created("aligned seq_data h5ad", seq_out)
    print_file_created("aligned spatial_data h5ad", spatial_out)
    print_file_created("shared genes CSV", genes_out)
    print_file_created("alignment summary CSV", summary_out)

    print("\nDONE. VISTA-ready Br8667 files created.")


if __name__ == "__main__":
    main()
