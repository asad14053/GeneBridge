#!/usr/bin/env python3

from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

# ---------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------

METADATA = (
    ROOT
    / "data/metadata/patient_xenium_visium_24_common_with_dx.csv"
)

N24_H5AD = (
    ROOT
    / "data/processed/xenium/xenium_N24_imputation_ready.h5ad"
)

# ---------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------

META_DIR = ROOT / "data/metadata"

EX1_MANIFEST = META_DIR / "imputation_ex1_ntc_samples.csv"
EX2_MANIFEST = META_DIR / "imputation_ex2_scz_samples.csv"

EX1_DIR = (
    ROOT
    / "data/processed/imputation_full/ex1_ntc/targets"
)

EX2_DIR = (
    ROOT
    / "data/processed/imputation_full/ex2_scz/targets"
)

SUMMARY_CSV = (
    META_DIR
    / "imputation_ex1_ex2_target_summary.csv"
)

EX1_DIR.mkdir(parents=True, exist_ok=True)
EX2_DIR.mkdir(parents=True, exist_ok=True)

# Protect against broken donor extraction.
MIN_EXPECTED_CELLS = 10000


def patient_number(x):
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else 999999


def matrix_stats(x):

    if sparse.issparse(x):

        total = x.shape[0] * x.shape[1]

        if x.nnz == 0:
            return 0.0, 0.0, 0.0, 1.0

        mn = min(0.0, float(x.data.min()))
        mx = float(x.data.max())
        mean = float(x.sum()) / total
        zero_fraction = 1.0 - (x.nnz / total)

        return mn, mx, mean, zero_fraction

    x = np.asarray(x)

    return (
        float(np.min(x)),
        float(np.max(x)),
        float(np.mean(x)),
        float(np.mean(x == 0)),
    )


def integer_like(x, atol=1e-6):

    vals = x.data if sparse.issparse(x) else np.asarray(x).ravel()

    if vals.size == 0:
        return True

    return bool(
        np.all(
            np.abs(vals - np.rint(vals)) <= atol
        )
    )


print("=" * 110)
print("PREPARE PRODUCTION TARGETS")
print("Ex1 = 12 NTC Xenium")
print("Ex2 = 12 SCZ Xenium")
print("Reference for both = all 10 Huuki snRNA")
print("=" * 110)


# =====================================================================
# 1. Metadata
# =====================================================================

meta = pd.read_csv(METADATA)

required = {
    "patient_index",
    "patient_id",
    "xenium_gsm",
    "Dx",
}

missing = required - set(meta.columns)

if missing:
    raise RuntimeError(
        f"Missing metadata columns: {sorted(missing)}"
    )

meta["Dx"] = (
    meta["Dx"]
    .astype(str)
    .str.strip()
    .str.upper()
)

meta["patient_id"] = (
    meta["patient_id"]
    .astype(str)
    .str.strip()
)

if set(meta["Dx"]) != {"NTC", "SCZ"}:
    raise RuntimeError(
        f"Unexpected Dx labels: {sorted(meta['Dx'].unique())}"
    )

ntc = meta.loc[meta["Dx"] == "NTC"].copy()
scz = meta.loc[meta["Dx"] == "SCZ"].copy()

ntc["_patient_number"] = ntc["patient_index"].map(patient_number)
scz["_patient_number"] = scz["patient_index"].map(patient_number)

ntc = (
    ntc.sort_values("_patient_number")
    .drop(columns="_patient_number")
)

scz = (
    scz.sort_values("_patient_number")
    .drop(columns="_patient_number")
)

if len(ntc) != 12:
    raise RuntimeError(
        f"Expected 12 NTC donors; found {len(ntc)}"
    )

if len(scz) != 12:
    raise RuntimeError(
        f"Expected 12 SCZ donors; found {len(scz)}"
    )

ntc.to_csv(EX1_MANIFEST, index=False)
scz.to_csv(EX2_MANIFEST, index=False)

print("\nManifests created:")
print("Ex1:", EX1_MANIFEST)
print("Ex2:", EX2_MANIFEST)

print("\nEx1 NTC donors:")
print(
    ntc[
        ["patient_index", "patient_id", "Dx", "xenium_gsm"]
    ].to_string(index=False)
)

print("\nEx2 SCZ donors:")
print(
    scz[
        ["patient_index", "patient_id", "Dx", "xenium_gsm"]
    ].to_string(index=False)
)


# =====================================================================
# 2. Read N24 ONCE
# =====================================================================

print("\n" + "=" * 110)
print("READING N24 XENIUM OBJECT ONCE")
print("=" * 110)

adata = ad.read_h5ad(N24_H5AD)

print("Shape :", adata.shape)
print("Layers:", list(adata.layers.keys()))
print("obsm  :", list(adata.obsm.keys()))

if "counts" not in adata.layers:
    raise RuntimeError("layers['counts'] missing")

if "spatial" not in adata.obsm:
    raise RuntimeError("obsm['spatial'] missing")

if "BrNum" not in adata.obs:
    raise RuntimeError("obs['BrNum'] missing")

if "Dx" not in adata.obs:
    raise RuntimeError("obs['Dx'] missing")

if adata.n_vars != 300:
    raise RuntimeError(
        f"Expected 300 Xenium genes; found {adata.n_vars}"
    )


# =====================================================================
# 3. Convert donor column ONCE
# =====================================================================

brnum = (
    adata.obs["BrNum"]
    .astype(str)
    .to_numpy()
)

dx_array = (
    adata.obs["Dx"]
    .astype(str)
    .str.upper()
    .to_numpy()
)

# This avoids repeatedly converting 1.22 million strings.
print("\nDonor index prepared.")


# =====================================================================
# 4. Extraction function
# =====================================================================

summary_rows = []


def extract_experiment(
    experiment,
    diagnosis,
    manifest_df,
    out_dir,
):

    print("\n" + "=" * 110)
    print(
        f"{experiment.upper()} — {diagnosis} TARGET EXTRACTION"
    )
    print("=" * 110)

    n_written = 0
    n_blocked = 0

    for i, row in enumerate(
        manifest_df.itertuples(index=False),
        start=1,
    ):

        donor = str(row.patient_id)
        gsm = str(row.xenium_gsm)

        print(
            f"\n[{i:02d}/12] "
            f"{experiment} | {donor} | {diagnosis}",
            flush=True,
        )

        mask = brnum == donor
        n_cells = int(mask.sum())

        print("Cells found:", n_cells, flush=True)

        # ----------------------------------------------------------
        # Donor QC guard
        # ----------------------------------------------------------

        if n_cells < MIN_EXPECTED_CELLS:

            print(
                "BLOCKED: implausibly low cell count.",
                flush=True,
            )

            summary_rows.append(
                {
                    "experiment": experiment,
                    "Dx": diagnosis,
                    "patient_id": donor,
                    "xenium_gsm": gsm,
                    "status": "BLOCKED_LOW_CELLS",
                    "n_cells": n_cells,
                    "n_genes": adata.n_vars,
                    "target_h5ad": "",
                }
            )

            n_blocked += 1
            continue

        donor_dx = sorted(
            set(dx_array[mask])
        )

        if donor_dx != [diagnosis]:

            raise RuntimeError(
                f"{donor}: expected Dx={diagnosis}, "
                f"found {donor_dx}"
            )

        # ----------------------------------------------------------
        # Build lightweight production input
        #
        # IMPORTANT:
        # X = raw Xenium counts.
        # No normalized matrix is used as model input.
        # ----------------------------------------------------------

        counts = adata.layers["counts"][mask].copy()

        obs = adata.obs.loc[mask].copy()
        var = adata.var.copy()

        coords = (
            np.asarray(
                adata.obsm["spatial"][mask]
            )
            .copy()
        )

        target = ad.AnnData(
            X=counts,
            obs=obs,
            var=var,
        )

        target.obsm["spatial"] = coords

        # Keep an explicit counts layer for downstream compatibility.
        target.layers["counts"] = counts.copy()

        # Do NOT copy log1p_norm:
        # production imputation input is deliberately count scale.

        target.uns["production_imputation"] = {
            "experiment": experiment,
            "diagnosis": diagnosis,
            "reference": "all_10_Huuki_snRNA",
            "target_donor": donor,
            "gene_holdout": False,
            "n_xenium_input_genes": int(target.n_vars),
            "X_scale": "raw_counts",
            "spatial_key": "spatial",
        }

        mn, mx, mean, zf = matrix_stats(target.X)
        int_like = integer_like(target.X)

        spatial_finite = bool(
            np.isfinite(coords).all()
        )

        print("Shape           :", target.shape)
        print("X scale         : raw counts")
        print("Count min       :", mn)
        print("Count max       :", mx)
        print("Count mean      :", mean)
        print("Zero fraction   :", zf)
        print("Integer-like    :", int_like)
        print("Spatial shape   :", coords.shape)
        print("Spatial finite  :", spatial_finite)

        if mn < 0:
            raise RuntimeError(
                f"{donor}: negative counts detected"
            )

        if not int_like:
            raise RuntimeError(
                f"{donor}: counts are not integer-like"
            )

        if coords.shape != (n_cells, 2):
            raise RuntimeError(
                f"{donor}: bad spatial shape {coords.shape}"
            )

        if not spatial_finite:
            raise RuntimeError(
                f"{donor}: non-finite spatial coordinates"
            )

        out_file = (
            out_dir
            / f"spatial_data_xenium_{donor}_{experiment}.h5ad"
        )

        print("Writing:", out_file, flush=True)

        # LZF is much faster than gzip and still compressed.
        target.write_h5ad(
            out_file,
            compression="lzf",
        )

        print("DONE:", donor, flush=True)

        summary_rows.append(
            {
                "experiment": experiment,
                "Dx": diagnosis,
                "patient_id": donor,
                "xenium_gsm": gsm,
                "status": "READY",
                "n_cells": n_cells,
                "n_genes": target.n_vars,
                "x_scale": "raw_counts",
                "integer_like": int_like,
                "spatial_available": True,
                "target_h5ad": str(out_file),
            }
        )

        n_written += 1

        del target
        del counts
        del obs
        del var
        del coords

    print("\n" + "-" * 80)
    print(
        f"{experiment}: written={n_written}, "
        f"blocked={n_blocked}"
    )

    if n_written == 12:
        ready_flag = out_dir.parent / "TARGETS_READY.flag"

        ready_flag.write_text(
            f"{experiment}: 12/12 targets READY\n"
        )

        print("READY FLAG:", ready_flag)

    return n_written, n_blocked


# =====================================================================
# 5. Experiment 1 — NTC
# =====================================================================

ex1_written, ex1_blocked = extract_experiment(
    experiment="ex1_ntc",
    diagnosis="NTC",
    manifest_df=ntc,
    out_dir=EX1_DIR,
)


# =====================================================================
# 6. Experiment 2 — SCZ
# =====================================================================

ex2_written, ex2_blocked = extract_experiment(
    experiment="ex2_scz",
    diagnosis="SCZ",
    manifest_df=scz,
    out_dir=EX2_DIR,
)


# =====================================================================
# 7. Final summary
# =====================================================================

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    SUMMARY_CSV,
    index=False,
)

print("\n" + "=" * 110)
print("FINAL EX1 + EX2 TARGET SUMMARY")
print("=" * 110)

print(
    summary[
        [
            "experiment",
            "patient_id",
            "Dx",
            "status",
            "n_cells",
            "n_genes",
        ]
    ].to_string(index=False)
)

print("\nExperiment summary:")

print(
    summary.groupby(
        ["experiment", "status"]
    ).size()
)

print("\nEx1 NTC:")
print(f"  READY   = {ex1_written}/12")
print(f"  BLOCKED = {ex1_blocked}/12")

print("\nEx2 SCZ:")
print(f"  READY   = {ex2_written}/12")
print(f"  BLOCKED = {ex2_blocked}/12")

print("\nCombined summary:")
print(SUMMARY_CSV)

print("\nEx1 targets:")
print(EX1_DIR)

print("\nEx2 targets:")
print(EX2_DIR)

print("\n" + "=" * 110)

if ex1_written == 12 and ex2_written == 12:
    print("SUCCESS: BOTH EXPERIMENTS ARE 12/12 READY")
elif ex2_written == 12:
    print(
        "EX2 SCZ IS READY. "
        "EX1 HAS AT LEAST ONE BLOCKED DONOR."
    )
else:
    print(
        "WARNING: one or both experiments contain blocked donors."
    )

print("=" * 110)
