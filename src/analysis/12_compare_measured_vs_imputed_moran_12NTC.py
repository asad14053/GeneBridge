#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MEAS_BASE = (
    ROOT
    / "results/analysis/measured300_moran_visium_vs_xenium"
)

IMP_BASE = (
    ROOT
    / "results/analysis/envi_moran_visium_vs_imputed"
)

OUT = (
    ROOT
    / "results/analysis/moran_measured_vs_imputed_12NTC"
)

OUT.mkdir(parents=True, exist_ok=True)

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

    fm = (
        MEAS_BASE
        / br
        / f"{br}_measured300_MoranI_summary.csv"
    )

    fi = (
        IMP_BASE
        / br
        / f"{br}_Visium_vs_ENVI_MoranI_summary.csv"
    )

    if not fm.exists():
        raise FileNotFoundError(fm)

    if not fi.exists():
        raise FileNotFoundError(fi)

    m = pd.read_csv(fm).iloc[0]
    i = pd.read_csv(fi).iloc[0]

    rows.append({
        "BrNum": br,

        "measured_spearman":
            m["spearman_moran"],

        "imputed_spearman":
            i["spearman_moran"],

        "measured_pearson":
            m["pearson_moran"],

        "imputed_pearson":
            i["pearson_moran"],

        "measured_slope":
            m["regression_slope"],

        "imputed_slope":
            i["regression_slope"],

        "measured_delta_moran":
            m["median_delta_moran"],

        "imputed_delta_moran":
            i["median_delta_moran"],

        "n_measured_valid":
            m["n_measured_valid_moran_both"],

        "n_imputed_valid":
            i.get(
                "n_imputed_valid_moran_both",
                np.nan
            ),
    })


df = pd.DataFrame(rows)


# ============================================================
# DERIVED DIFFERENCES
# ============================================================

df["diff_spearman"] = (
    df["measured_spearman"]
    -
    df["imputed_spearman"]
)

df["diff_slope"] = (
    df["measured_slope"]
    -
    df["imputed_slope"]
)

df["diff_abs_delta"] = (
    np.abs(df["imputed_delta_moran"])
    -
    np.abs(df["measured_delta_moran"])
)


# ============================================================
# PAIRED TESTS
# ============================================================

w_rho = wilcoxon(
    df["measured_spearman"],
    df["imputed_spearman"]
)

w_slope = wilcoxon(
    df["measured_slope"],
    df["imputed_slope"]
)

w_delta = wilcoxon(
    np.abs(df["measured_delta_moran"]),
    np.abs(df["imputed_delta_moran"])
)


# ============================================================
# SAVE TABLE
# ============================================================

TABLE = OUT / "measured_vs_imputed_MoranI_12NTC.csv"

df.to_csv(
    TABLE,
    index=False
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print("=" * 90)
print("MEASURED-300 vs ENVI-IMPUTED MORAN-I CONTROL")
print("=" * 90)

print(df.to_string(index=False))

print()
print("Median measured Spearman:",
      np.nanmedian(df["measured_spearman"]))

print("Median imputed Spearman:",
      np.nanmedian(df["imputed_spearman"]))

print("Wilcoxon Spearman p:",
      w_rho.pvalue)

print()

print("Median measured slope:",
      np.nanmedian(df["measured_slope"]))

print("Median imputed slope:",
      np.nanmedian(df["imputed_slope"]))

print("Wilcoxon slope p:",
      w_slope.pvalue)

print()

print("Median measured |delta Moran|:",
      np.nanmedian(np.abs(df["measured_delta_moran"])))

print("Median imputed |delta Moran|:",
      np.nanmedian(np.abs(df["imputed_delta_moran"])))

print("Wilcoxon |delta Moran| p:",
      w_delta.pvalue)


# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(figsize=(13, 9))

gs = fig.add_gridspec(
    2,
    2,
    hspace=0.38,
    wspace=0.30
)


# ------------------------------------------------------------
# A — Spearman
# ------------------------------------------------------------

ax = fig.add_subplot(gs[0, 0])

for _, r in df.iterrows():

    ax.plot(
        [0, 1],
        [
            r["measured_spearman"],
            r["imputed_spearman"]
        ],
        alpha=0.5,
        linewidth=1
    )

    ax.scatter(
        [0, 1],
        [
            r["measured_spearman"],
            r["imputed_spearman"]
        ],
        s=30
    )

ax.axhline(0, linewidth=0.8)

ax.set_xticks([0, 1])

ax.set_xticklabels([
    "Measured 300",
    "ENVI imputed-only"
])

ax.set_ylabel(
    "Gene-level Spearman correlation\n"
    "(Visium Moran's I vs Xenium Moran's I)"
)

ax.set_title(
    "a   Preservation of gene ranking",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"Measured median = "
    f"{np.nanmedian(df['measured_spearman']):.3f}\n"
    f"Imputed median = "
    f"{np.nanmedian(df['imputed_spearman']):.3f}\n"
    f"paired Wilcoxon p={w_rho.pvalue:.3g}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ------------------------------------------------------------
# B — Slope
# ------------------------------------------------------------

ax = fig.add_subplot(gs[0, 1])

for _, r in df.iterrows():

    ax.plot(
        [0, 1],
        [
            r["measured_slope"],
            r["imputed_slope"]
        ],
        alpha=0.5,
        linewidth=1
    )

    ax.scatter(
        [0, 1],
        [
            r["measured_slope"],
            r["imputed_slope"]
        ],
        s=30
    )

ax.axhline(
    1,
    linestyle="--",
    linewidth=1
)

ax.axhline(
    0,
    linewidth=0.8
)

ax.set_xticks([0, 1])

ax.set_xticklabels([
    "Measured 300",
    "ENVI imputed-only"
])

ax.set_ylabel(
    "Regression slope\n"
    "(ideal = 1)"
)

ax.set_title(
    "b   Preservation of Moran-I magnitude",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"Measured median = "
    f"{np.nanmedian(df['measured_slope']):.3f}\n"
    f"Imputed median = "
    f"{np.nanmedian(df['imputed_slope']):.3f}\n"
    f"paired Wilcoxon p={w_slope.pvalue:.3g}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ------------------------------------------------------------
# C — median delta
# ------------------------------------------------------------

ax = fig.add_subplot(gs[1, 0])

for _, r in df.iterrows():

    ax.plot(
        [0, 1],
        [
            r["measured_delta_moran"],
            r["imputed_delta_moran"]
        ],
        alpha=0.5,
        linewidth=1
    )

    ax.scatter(
        [0, 1],
        [
            r["measured_delta_moran"],
            r["imputed_delta_moran"]
        ],
        s=30
    )

ax.axhline(
    0,
    linestyle="--",
    linewidth=1
)

ax.set_xticks([0, 1])

ax.set_xticklabels([
    "Measured 300",
    "ENVI imputed-only"
])

ax.set_ylabel(
    "Median Δ Moran's I\n"
    "(Xenium − Visium)"
)

ax.set_title(
    "c   Systematic spatial-autocorrelation shift",
    loc="left",
    fontweight="bold"
)


# ------------------------------------------------------------
# D — donor-level diagnostic plane
# ------------------------------------------------------------

ax = fig.add_subplot(gs[1, 1])

ax.scatter(
    df["imputed_spearman"],
    df["imputed_delta_moran"],
    s=55
)

for _, r in df.iterrows():

    ax.annotate(
        r["BrNum"],
        (
            r["imputed_spearman"],
            r["imputed_delta_moran"]
        ),
        fontsize=8,
        xytext=(3, 3),
        textcoords="offset points"
    )

ax.axvline(
    0,
    linewidth=0.8
)

ax.axhline(
    0,
    linewidth=0.8
)

ax.set_xlabel(
    "ENVI gene-level Moran Spearman"
)

ax.set_ylabel(
    "ENVI median Δ Moran's I"
)

ax.set_title(
    "d   ENVI spatial-fidelity diagnostic",
    loc="left",
    fontweight="bold"
)


fig.suptitle(
    "Measured-gene positive control versus "
    "ENVI-imputed spatial autocorrelation",
    fontsize=14,
    fontweight="bold"
)

fig.text(
    0.5,
    0.015,
    "Measured and imputed groups contain different gene sets; "
    "this comparison is a cross-platform positive control, "
    "not a same-gene imputation accuracy test.",
    ha="center",
    fontsize=8.5
)

PNG = OUT / "fig_measured_vs_imputed_MoranI_12NTC.png"
PDF = OUT / "fig_measured_vs_imputed_MoranI_12NTC.pdf"

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
