#!/usr/bin/env python

"""
04_validate_vista_Br8667_inputs.py

Purpose:
    Validate final Br8667 VISTA-ready inputs before imputation.

Checks:
    1. Both h5ad files exist
    2. seq_data and spatial_data have same genes
    3. gene order is identical
    4. spatial_data has obsm["spatial"]
    5. no duplicated genes
    6. print final paths for VISTA
"""

from pathlib import Path
import pandas as pd
import anndata as ad


PROJECT_ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BASE = PROJECT_ROOT / "data/processed/imputation_beta/Br8667"

SEQ_PATH = BASE / "seq_data_huuki_snrna_Br8667_shared.h5ad"
SPATIAL_PATH = BASE / "spatial_data_huuki_visium_Br8667_shared.h5ad"

OUT_DIR = PROJECT_ROOT / "outputs/imputation_beta/Br8667"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_OUT = OUT_DIR / "04_validate_vista_Br8667_inputs_summary.csv"


def main():
    print("=" * 100)
    print("Step 6: Validate VISTA-ready Br8667 inputs")
    print("=" * 100)

    print("\nInput files:")
    print("seq_data:    ", SEQ_PATH)
    print("spatial_data:", SPATIAL_PATH)

    if not SEQ_PATH.exists():
        raise FileNotFoundError(f"Missing seq_data: {SEQ_PATH}")

    if not SPATIAL_PATH.exists():
        raise FileNotFoundError(f"Missing spatial_data: {SPATIAL_PATH}")

    print("\nLoading seq_data...")
    seq = ad.read_h5ad(SEQ_PATH)

    print("Loading spatial_data...")
    spatial = ad.read_h5ad(SPATIAL_PATH)

    print("\nseq_data:")
    print(seq)

    print("\nspatial_data:")
    print(spatial)

    same_gene_count = seq.n_vars == spatial.n_vars
    same_gene_order = seq.var_names.equals(spatial.var_names)
    seq_unique_genes = seq.var_names.is_unique
    spatial_unique_genes = spatial.var_names.is_unique
    spatial_has_coordinates = "spatial" in spatial.obsm

    print("\nValidation checks:")
    print("Same number of genes:       ", same_gene_count)
    print("Same gene order:            ", same_gene_order)
    print("seq_data unique genes:      ", seq_unique_genes)
    print("spatial_data unique genes:  ", spatial_unique_genes)
    print("spatial has obsm['spatial']:", spatial_has_coordinates)

    print("\nFinal dimensions:")
    print("seq cells:      ", seq.n_obs)
    print("spatial spots:  ", spatial.n_obs)
    print("shared genes:   ", seq.n_vars)

    print("\nFirst 10 shared genes:")
    print(seq.var_names[:10].tolist())

    print("\nSpatial coordinate shape:")
    if spatial_has_coordinates:
        print(spatial.obsm["spatial"].shape)
    else:
        print("MISSING")

    summary = pd.DataFrame(
        [
            {"metric": "seq_path", "value": str(SEQ_PATH)},
            {"metric": "spatial_path", "value": str(SPATIAL_PATH)},
            {"metric": "seq_cells", "value": seq.n_obs},
            {"metric": "spatial_spots", "value": spatial.n_obs},
            {"metric": "shared_genes", "value": seq.n_vars},
            {"metric": "same_gene_count", "value": same_gene_count},
            {"metric": "same_gene_order", "value": same_gene_order},
            {"metric": "seq_unique_genes", "value": seq_unique_genes},
            {"metric": "spatial_unique_genes", "value": spatial_unique_genes},
            {"metric": "spatial_has_obsm_spatial", "value": spatial_has_coordinates},
        ]
    )

    summary.to_csv(SUMMARY_OUT, index=False)

    print("\nSaved validation summary:")
    print(SUMMARY_OUT)

    if not same_gene_count:
        raise ValueError("FAILED: seq_data and spatial_data have different gene counts.")

    if not same_gene_order:
        raise ValueError("FAILED: seq_data and spatial_data genes are not in same order.")

    if not seq_unique_genes:
        raise ValueError("FAILED: seq_data has duplicated gene names.")

    if not spatial_unique_genes:
        raise ValueError("FAILED: spatial_data has duplicated gene names.")

    if not spatial_has_coordinates:
        raise ValueError("FAILED: spatial_data does not have obsm['spatial'].")

    print("=" * 100)
    print("PASS: Br8667 VISTA inputs are ready")
    print("=" * 100)

    print("\nUse these paths in VISTA:")
    print(f'seq_data = ad.read_h5ad("{SEQ_PATH}")')
    print(f'spatial_data = ad.read_h5ad("{SPATIAL_PATH}")')


if __name__ == "__main__":
    main()
