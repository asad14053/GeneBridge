#!/usr/bin/env python

from pathlib import Path
import re
import anndata as ad

PROJECT_ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

in_h5ad = PROJECT_ROOT / "data/processed/visium/spatialDLPFC_Visium_sce.h5ad"

out_dir = PROJECT_ROOT / "data/processed/imputation_beta/Br8667"
out_dir.mkdir(parents=True, exist_ok=True)

out_h5ad = out_dir / "spatial_data_huuki_visium_Br8667.h5ad"


def normalize_brain_id(x):
    s = str(x)
    m = re.search(r"(\d{4})", s)
    if m is None:
        return None
    return "Br" + m.group(1)


def find_brain_column(obs):
    candidates = [
        "BrNum",
        "BrNum_matched",
        "brnum",
        "sample_id",
        "sample",
        "Sample",
        "SAMPLE_ID",
        "brain",
        "Brain",
        "donor",
        "Donor",
        "donor_id",
        "subject",
        "Subject",
    ]

    for c in candidates:
        if c in obs.columns:
            vals = obs[c].map(normalize_brain_id)
            if vals.notna().sum() > 0:
                return c

    possible = [
        c for c in obs.columns
        if any(k in c.lower() for k in ["br", "brain", "sample", "donor", "subject"])
    ]

    print("Could not find brain column.")
    print("Possible columns:")
    print(possible)
    print("All columns:")
    print(obs.columns.tolist())

    raise ValueError("No usable brain column found.")


print("Loading Huuki Visium:")
print(in_h5ad)

adata = ad.read_h5ad(in_h5ad)

print("\nOriginal Huuki Visium:")
print(adata)

br_col = find_brain_column(adata.obs)

br = adata.obs[br_col].map(normalize_brain_id)
keep = br == "Br8667"

print("\nBrain column used:", br_col)
print("Br8667 Visium spots:", int(keep.sum()))

if keep.sum() == 0:
    raise ValueError("No Br8667 spots found in Huuki Visium.")

adata_b = adata[keep].copy()

# Use gene symbols as var_names if available.
if "gene_name" in adata_b.var.columns:
    adata_b.var["original_var_names"] = adata_b.var_names.astype(str)
    adata_b.var["gene_symbol"] = adata_b.var["gene_name"].astype(str)

    adata_b.var_names = adata_b.var["gene_symbol"].astype(str)
    adata_b.var_names_make_unique()
else:
    adata_b.var_names_make_unique()

adata_b.obs_names_make_unique()

# Important: avoid AnnData write error from index-name conflicts.
adata_b.var.index.name = None
adata_b.obs.index.name = None

print("\nBr8667 Huuki Visium subset:")
print(adata_b)

print("\nOutput h5ad:")
print(out_h5ad)

if out_h5ad.exists():
    out_h5ad.unlink()

adata_b.write_h5ad(out_h5ad)

print("\nDONE. Saved:")
print(out_h5ad)
