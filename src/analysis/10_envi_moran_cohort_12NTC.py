#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE = (
    ROOT
    / "results/analysis/envi_moran_visium_vs_imputed"
)

OUT = (
    BASE
    / "cohort_12NTC"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

DONORS = [
    "Br1113",
    "Br2719",
    "Br5314",
    "Br5400",
    "Br5436",
    "Br5622",
    "Br5639",
    "Br5931",
    "Br6389",
    "Br6432",
    "Br8433",
    "Br8667",
]

rows = []

for br in DONORS:

    f = (
        BASE
        / br
        / f"{br}_Visium_vs_ENVI_MoranI_summary.csv"
    )

    if not f.exists():
        raise FileNotFoundError(f)

    x = pd.read_csv(f)

    if len(x) != 1:
        raise RuntimeError(
            f"{br}: expected one summary row."
        )

    rows.append(x.iloc[0])

df = pd.DataFrame(rows)

df["BrNum"] = DONORS

# ------------------------------------------
# Save cohort table
# ------------------------------------------

TABLE = (
    OUT
    / "ENVI_MoranI_12NTC_summary.csv"
)

df.to_csv(
    TABLE,
    index=False
)

# ------------------------------------------
# Print
# ------------------------------------------

print("=" * 90)
print("12-NTC MORAN-I COHORT SUMMARY")
print("=" * 90)

cols = [
    "BrNum",
    "n_imputed_shared_genes",
    "n_imputed_valid_moran_both",
    "spearman_moran",
    "pearson_moran",
    "regression_slope",
    "median_delta_moran",
]

available = [
    c for c in cols
    if c in df.columns
]

print(
    df[available]
    .to_string(index=False)
)

print()
print(
    "median Spearman:",
    np.nanmedian(
        df["spearman_moran"]
    )
)

print(
    "median Pearson:",
    np.nanmedian(
        df["pearson_moran"]
    )
)

print(
    "median slope:",
    np.nanmedian(
        df["regression_slope"]
    )
)

print(
    "median delta Moran:",
    np.nanmedian(
        df["median_delta_moran"]
    )
)


# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(13, 8.5)
)

gs = fig.add_gridspec(
    2,
    2,
    hspace=0.38,
    wspace=0.30
)


# ------------------------------------------------------------
# A — Moran correlation by donor
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[0, 0]
)

d = df.sort_values(
    "spearman_moran",
    ascending=False
)

ax.bar(
    d["BrNum"],
    d["spearman_moran"]
)

ax.axhline(
    np.nanmedian(
        d["spearman_moran"]
    ),
    linestyle="--",
    linewidth=1
)

ax.set_ylabel(
    "Spearman correlation\n"
    "(Visium Moran's I vs ENVI Moran's I)"
)

ax.tick_params(
    axis="x",
    rotation=60
)

ax.set_title(
    "a   Gene-level Moran's I correspondence",
    loc="left",
    fontweight="bold"
)


# ------------------------------------------------------------
# B — regression slope
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[0, 1]
)

d = df.sort_values(
    "regression_slope",
    ascending=False
)

ax.bar(
    d["BrNum"],
    d["regression_slope"]
)

ax.axhline(
    1,
    linestyle="--",
    linewidth=1
)

ax.axhline(
    np.nanmedian(
        d["regression_slope"]
    ),
    linestyle=":",
    linewidth=1
)

ax.set_ylabel(
    "Regression slope\n"
    "ENVI Moran's I ~ Visium Moran's I"
)

ax.tick_params(
    axis="x",
    rotation=60
)

ax.set_title(
    "b   Preservation of spatial-autocorrelation magnitude",
    loc="left",
    fontweight="bold"
)


# ------------------------------------------------------------
# C — median delta Moran
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[1, 0]
)

d = df.sort_values(
    "median_delta_moran",
    ascending=True
)

ax.barh(
    d["BrNum"],
    d["median_delta_moran"]
)

ax.axvline(
    0,
    linewidth=1
)

ax.axvline(
    np.nanmedian(
        d["median_delta_moran"]
    ),
    linestyle="--",
    linewidth=1
)

ax.set_xlabel(
    "Median Δ Moran's I\n"
    "(ENVI − Visium)"
)

ax.set_title(
    "c   Direction of spatial-autocorrelation shift",
    loc="left",
    fontweight="bold"
)


# ------------------------------------------------------------
# D — valid genes
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[1, 1]
)

if "n_imputed_valid_moran_both" in df.columns:

    valid_col = (
        "n_imputed_valid_moran_both"
    )

else:

    valid_col = (
        "n_imputed_evaluable_genes"
    )

d = df.sort_values(
    valid_col,
    ascending=False
)

ax.bar(
    d["BrNum"],
    d[valid_col]
)

ax.tick_params(
    axis="x",
    rotation=60
)

ax.set_ylabel(
    "Imputed-only genes with\n"
    "valid Moran's I on both platforms"
)

ax.set_title(
    "d   Gene coverage",
    loc="left",
    fontweight="bold"
)


fig.suptitle(
    "Preservation of spatial autocorrelation in "
    "ENVI-imputed expression across 12 NTC donors",
    fontsize=14,
    fontweight="bold"
)

PNG = (
    OUT
    / "fig_ENVI_MoranI_12NTC_cohort.png"
)

PDF = (
    OUT
    / "fig_ENVI_MoranI_12NTC_cohort.pdf"
)

fig.savefig(
    PNG,
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    PDF,
    bbox_inches="tight"
)

plt.close(fig)

print()
print("Table:", TABLE)
print("PNG  :", PNG)
print("PDF  :", PDF)
