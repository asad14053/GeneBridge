#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

H5AD = (
    ROOT / "data/processed/xenium/"
    "xenium_N24_layer_celltype_annotated.h5ad"
)

COHORT = (
    ROOT / "data/metadata/"
    "patient_xenium_visium_24_common_with_dx.csv"
)

RUN_META = (
    ROOT / "data/metadata/"
    "xenium_PI_run_metadata_24.csv"
)

OUT = (
    ROOT / "outputs/imputation_full/DE/"
    "task2_3_preimputation_baseline/domain_pseudobulk"
)

OUT.mkdir(parents=True, exist_ok=True)

COUNTS_OUT = OUT / "original_xenium_N24_donor_domain_pseudobulk_300genes.csv.gz"
META_OUT   = OUT / "original_xenium_N24_donor_domain_pseudobulk_metadata.csv"
VALID_OUT  = OUT / "original_xenium_N24_donor_domain_validation.csv"


def find_col(df, candidates, label):
    for c in candidates:
        if c in df.columns:
            return c
    raise RuntimeError(
        f"Could not resolve {label}. Available: {list(df.columns)}"
    )


print("=" * 110)
print("PRE-IMPUTATION ORIGINAL XENIUM — DONOR × DOMAIN PSEUDOBULK")
print("=" * 110)

# =============================================================================
# LOAD ORIGINAL MEASURED XENIUM
# =============================================================================

adata = sc.read_h5ad(H5AD)

print("\nInput:", H5AD)
print("Shape:", adata.shape)
print("Layers:", list(adata.layers.keys()))

if adata.shape[1] != 300:
    raise RuntimeError(f"Expected 300 genes, found {adata.shape[1]}")

if "counts" not in adata.layers:
    raise RuntimeError("Required measured layer 'counts' is absent.")

for c in ["BrNum", "predictions_smooth", "layer_annotation"]:
    if c not in adata.obs.columns:
        raise RuntimeError(f"Missing obs column: {c}")


# =============================================================================
# CANONICAL DONOR METADATA
# =============================================================================

cohort = pd.read_csv(COHORT)

br_col  = find_col(cohort, ["BrNum", "Brnumbr", "brnum"], "donor")
dx_col  = find_col(cohort, ["Dx", "DX", "dx"], "Dx")
age_col = find_col(cohort, ["Age", "AGE", "age"], "Age")
sex_col = find_col(cohort, ["Sex", "SEX", "sex"], "Sex")

donor = cohort[[br_col, dx_col, age_col, sex_col]].copy()
donor.columns = ["BrNum", "Dx", "Age", "Sex"]

donor["BrNum"] = donor["BrNum"].astype(str).str.strip()

donor["Dx"] = (
    donor["Dx"].astype(str).str.strip().str.upper()
    .replace({
        "CONTROL": "NTC",
        "CON": "NTC",
        "CTRL": "NTC",
        "SCHIZOPHRENIA": "SCZ",
    })
)

donor["Age"] = pd.to_numeric(donor["Age"], errors="coerce")

donor["Sex"] = (
    donor["Sex"].astype(str).str.strip().str.upper()
    .replace({"MALE": "M", "FEMALE": "F"})
)

donor = donor.drop_duplicates()

if len(donor) != 24:
    raise RuntimeError(f"Expected 24 donors, got {len(donor)}")

if (donor["Dx"] == "NTC").sum() != 12:
    raise RuntimeError("Expected 12 NTC donors")

if (donor["Dx"] == "SCZ").sum() != 12:
    raise RuntimeError("Expected 12 SCZ donors")

if donor["Age"].isna().any():
    raise RuntimeError("Missing Age in canonical metadata")

if donor["Sex"].nunique() < 2:
    raise RuntimeError("Sex metadata invalid")


# =============================================================================
# TECHNICAL METADATA
# =============================================================================

run = pd.read_csv(RUN_META)

run_br = find_col(run, ["BrNum", "Brnumbr", "brnum"], "run-metadata donor")

for c in ["slide_id", "run_date"]:
    if c not in run.columns:
        raise RuntimeError(f"Missing {c} from {RUN_META}")

run = (
    run[[run_br, "slide_id", "run_date"]]
    .rename(columns={run_br: "BrNum"})
    .drop_duplicates()
)

run["BrNum"] = run["BrNum"].astype(str).str.strip()

donor = donor.merge(
    run,
    on="BrNum",
    how="left",
    validate="one_to_one"
)

if donor[["slide_id", "run_date"]].isna().any().any():
    raise RuntimeError("Missing slide_id/run_date after donor mapping")


# =============================================================================
# PREPARE CELL METADATA
# =============================================================================

obs = adata.obs[
    ["BrNum", "predictions_smooth", "layer_annotation"]
].copy()

obs["BrNum"] = obs["BrNum"].astype(str)
obs["predictions_smooth"] = obs["predictions_smooth"].astype(str)
obs["layer_annotation"] = obs["layer_annotation"].astype(str)

# Use original measured raw counts, NOT log1p_norm and NOT imputed count_scale.
X = adata.layers["counts"]

if sparse.issparse(X):
    X = X.tocsr()
else:
    X = np.asarray(X)

genes = adata.var_names.astype(str).tolist()

print("\nMeasured matrix used: layers['counts']")
print("Cells:", adata.n_obs)
print("Genes:", len(genes))


# =============================================================================
# DONOR × DOMAIN SUMMATION
#
# GitHub conceptual equivalent:
# var_registration = domain
# var_sample_id = BrNum
# min_ncells = 10
# =============================================================================

groups = (
    obs.groupby(
        ["BrNum", "predictions_smooth", "layer_annotation"],
        observed=True,
        sort=True
    )
    .indices
)

pb_rows = []
meta_rows = []

donor_lookup = donor.set_index("BrNum")

for (br, domain, layer), idx in groups.items():

    idx = np.asarray(idx)

    n_cells = len(idx)

    # Match GitHub min_ncells = 10
    if n_cells < 10:
        continue

    if sparse.issparse(X):
        summed = np.asarray(X[idx, :].sum(axis=0)).ravel()
    else:
        summed = X[idx, :].sum(axis=0)

    pb_id = f"{br}__{domain}"

    pb_rows.append(
        pd.Series(
            summed,
            index=genes,
            name=pb_id
        )
    )

    d = donor_lookup.loc[br]

    meta_rows.append({
        "pseudobulk_id": pb_id,
        "BrNum": br,
        "Dx": d["Dx"],
        "Age": d["Age"],
        "Sex": d["Sex"],
        "predictions_smooth": domain,
        "layer_annotation": layer,
        "slide_id": d["slide_id"],
        "run_date": d["run_date"],
        "n_cells": n_cells,
    })


pb = pd.DataFrame(pb_rows)
meta = pd.DataFrame(meta_rows)

pb.index.name = "pseudobulk_id"


# =============================================================================
# VALIDATION
# =============================================================================

print("\n" + "=" * 110)
print("VALIDATION")
print("=" * 110)

print("Pseudobulks :", len(meta))
print("Donors      :", meta["BrNum"].nunique())
print("Domains     :", meta["predictions_smooth"].nunique())
print("Genes       :", pb.shape[1])
print("NTC donors  :", meta.loc[meta.Dx == "NTC", "BrNum"].nunique())
print("SCZ donors  :", meta.loc[meta.Dx == "SCZ", "BrNum"].nunique())
print("Slides      :", meta["slide_id"].nunique())
print("Min cells/PB:", meta["n_cells"].min())
print("Max cells/PB:", meta["n_cells"].max())

if meta["BrNum"].nunique() != 24:
    raise RuntimeError("Expected 24 donors")

if meta["predictions_smooth"].nunique() != 7:
    raise RuntimeError("Expected 7 spatial domains")

if len(meta) != 168:
    raise RuntimeError(
        f"Expected 24 × 7 = 168 donor-domain PBs; found {len(meta)}"
    )

if pb.shape != (168, 300):
    raise RuntimeError(f"Expected 168 × 300 PB matrix; found {pb.shape}")

if meta["n_cells"].min() < 10:
    raise RuntimeError("Found PB with <10 cells")


# =============================================================================
# SAVE
# =============================================================================

pb.to_csv(COUNTS_OUT, compression="gzip")
meta.to_csv(META_OUT, index=False)

validation = (
    meta.groupby(
        ["BrNum", "predictions_smooth"],
        observed=True
    )
    .agg(
        n_cells=("n_cells", "first"),
        Dx=("Dx", "first"),
        Age=("Age", "first"),
        Sex=("Sex", "first"),
        slide_id=("slide_id", "first"),
        run_date=("run_date", "first"),
    )
    .reset_index()
)

validation.to_csv(VALID_OUT, index=False)

print("\nSaved counts:")
print(COUNTS_OUT)

print("\nSaved metadata:")
print(META_OUT)

print("\nSaved validation:")
print(VALID_OUT)

print("\nFINAL STATUS: PASS — ORIGINAL MEASURED XENIUM DONOR × DOMAIN PSEUDOBULK")
