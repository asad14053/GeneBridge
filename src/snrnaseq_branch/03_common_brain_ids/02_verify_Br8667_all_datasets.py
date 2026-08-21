#!/usr/bin/env python

from pathlib import Path
import re
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

paths = {
    "xenium_N24": PROJECT_ROOT / "data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad",
    "visium_N24": PROJECT_ROOT / "data/processed/visium/visium_N24_matched_layer_annotated.h5ad",
    "huuki_visium": PROJECT_ROOT / "data/processed/visium/spatialDLPFC_Visium_sce.h5ad",
    "huuki_snrna_meta": PROJECT_ROOT / "outputs/huuki_myers/tables/huuki_snrna_metadata.csv",
}

out_dir = PROJECT_ROOT / "outputs/huuki_myers/tables"
out_dir.mkdir(parents=True, exist_ok=True)

def normalize_brain_id(x):
    if pd.isna(x):
        return None
    s = str(x)
    m = re.search(r"(\d{4})", s)
    if m is None:
        return None
    return "Br" + m.group(1)

def find_br_col(obs, dataset_name):
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

    print(f"\nCould not find brain column for {dataset_name}")
    print("Possible columns:")
    print(possible)
    raise ValueError(f"No brain ID column found for {dataset_name}")

brain_sets = {}
brain_counts = []
used_columns = []

# Xenium
xen = sc.read_h5ad(paths["xenium_N24"])
c = find_br_col(xen.obs, "xenium_N24")
xen_br = xen.obs[c].map(normalize_brain_id).dropna()
brain_sets["xenium_N24"] = set(xen_br)
used_columns.append({"dataset": "xenium_N24", "brain_column": c})
for br, n in xen_br.value_counts().items():
    brain_counts.append({"dataset": "xenium_N24", "BrNum": br, "count": int(n)})

# matched Visium N24
vis_n24 = sc.read_h5ad(paths["visium_N24"])
c = find_br_col(vis_n24.obs, "visium_N24")
vis_n24_br = vis_n24.obs[c].map(normalize_brain_id).dropna()
brain_sets["visium_N24"] = set(vis_n24_br)
used_columns.append({"dataset": "visium_N24", "brain_column": c})
for br, n in vis_n24_br.value_counts().items():
    brain_counts.append({"dataset": "visium_N24", "BrNum": br, "count": int(n)})

# Huuki Visium
hvis = sc.read_h5ad(paths["huuki_visium"])
c = find_br_col(hvis.obs, "huuki_visium")
hvis_br = hvis.obs[c].map(normalize_brain_id).dropna()
brain_sets["huuki_visium"] = set(hvis_br)
used_columns.append({"dataset": "huuki_visium", "brain_column": c})
for br, n in hvis_br.value_counts().items():
    brain_counts.append({"dataset": "huuki_visium", "BrNum": br, "count": int(n)})

# Huuki snRNA metadata
sn_meta = pd.read_csv(paths["huuki_snrna_meta"])
if "BrNum" not in sn_meta.columns:
    raise ValueError("Huuki snRNA metadata does not contain BrNum column.")

sn_br = sn_meta["BrNum"].map(normalize_brain_id).dropna()
brain_sets["huuki_snrna"] = set(sn_br)
used_columns.append({"dataset": "huuki_snrna", "brain_column": "BrNum"})
for br, n in sn_br.value_counts().items():
    brain_counts.append({"dataset": "huuki_snrna", "BrNum": br, "count": int(n)})

# Intersections
common_xenium_visium_N24 = sorted(brain_sets["xenium_N24"] & brain_sets["visium_N24"])
common_huuki = sorted(brain_sets["huuki_visium"] & brain_sets["huuki_snrna"])
common_all = sorted(
    brain_sets["xenium_N24"]
    & brain_sets["visium_N24"]
    & brain_sets["huuki_visium"]
    & brain_sets["huuki_snrna"]
)

print("\nBrain columns used:")
print(pd.DataFrame(used_columns))

print("\nXenium + matched Visium N24 common:")
print(len(common_xenium_visium_N24), common_xenium_visium_N24)

print("\nHuuki Visium + Huuki snRNA common:")
print(len(common_huuki), common_huuki)

print("\nCommon across all four datasets:")
print(len(common_all), common_all)

print("\nBr8667 check:")
for name, s in brain_sets.items():
    print(name, "Br8667" in s)

# Save outputs
pd.DataFrame(used_columns).to_csv(out_dir / "all_dataset_brain_columns_used.csv", index=False)

pd.DataFrame(brain_counts).to_csv(
    out_dir / "all_dataset_brain_counts_long.csv",
    index=False
)

pd.DataFrame({"BrNum": common_xenium_visium_N24}).to_csv(
    out_dir / "common_xenium_visium_N24_brains.csv",
    index=False
)

pd.DataFrame({"BrNum": common_huuki}).to_csv(
    out_dir / "common_huuki_visium_snrna_brains.csv",
    index=False
)

pd.DataFrame({"BrNum": common_all}).to_csv(
    out_dir / "common_all_four_datasets_brains.csv",
    index=False
)

br8667_check = pd.DataFrame([
    {"dataset": name, "Br8667_present": "Br8667" in s}
    for name, s in brain_sets.items()
])
br8667_check.to_csv(out_dir / "Br8667_all_dataset_check.csv", index=False)

print("\nSaved:")
print(out_dir / "common_all_four_datasets_brains.csv")
print(out_dir / "Br8667_all_dataset_check.csv")
print(out_dir / "all_dataset_brain_counts_long.csv")
