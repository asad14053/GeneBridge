#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon, spearmanr


# ============================================================
# CONFIG
# ============================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

INFILE = (
    ROOT
    / "results/analysis/crosstech_measured300"
    / "crosstech_summary_all24.csv"
)

OUT = (
    ROOT
    / "results/analysis/crosstech_measured300"
    / "cohort_summary"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INFILE)

required = [
    "BrNum",
    "n_visium_spots",
    "n_occupied_pseudospots",
    "n_retained_genes",
    "n_patterned_genes",
    "median_cosine_patterned",
    "median_cosine_nonpatterned",
    "median_spearman_patterned",
]

missing = [
    c for c in required
    if c not in df.columns
]

if missing:
    raise RuntimeError(
        f"Missing columns: {missing}\n"
        f"Available: {list(df.columns)}"
    )

df = df.copy()

df["BrNum"] = df["BrNum"].astype(str)

df["pseudospot_coverage"] = (
    df["n_occupied_pseudospots"]
    /
    df["n_visium_spots"]
)

df["delta_cosine"] = (
    df["median_cosine_patterned"]
    -
    df["median_cosine_nonpatterned"]
)


# ============================================================
# COHORT STATISTICS
# ============================================================

n = len(df)

pat = df[
    "median_cosine_patterned"
].to_numpy(float)

non = df[
    "median_cosine_nonpatterned"
].to_numpy(float)

delta = df[
    "delta_cosine"
].to_numpy(float)

coverage = df[
    "pseudospot_coverage"
].to_numpy(float)

npat = df[
    "n_patterned_genes"
].to_numpy(float)


def med_iqr(x):
    x = np.asarray(x, float)
    return (
        np.nanmedian(x),
        np.nanpercentile(x, 25),
        np.nanpercentile(x, 75),
    )


pat_med, pat_q1, pat_q3 = med_iqr(pat)
non_med, non_q1, non_q3 = med_iqr(non)
delta_med, delta_q1, delta_q3 = med_iqr(delta)
cov_med, cov_q1, cov_q3 = med_iqr(coverage)

good = np.isfinite(pat) & np.isfinite(non)

if good.sum() >= 2:
    w = wilcoxon(
        pat[good],
        non[good],
        alternative="greater"
    )
    wilcoxon_stat = float(w.statistic)
    wilcoxon_p = float(w.pvalue)
else:
    wilcoxon_stat = np.nan
    wilcoxon_p = np.nan


good_cov = (
    np.isfinite(coverage)
    &
    np.isfinite(pat)
)

if good_cov.sum() >= 3:
    rho_cov, p_cov = spearmanr(
        coverage[good_cov],
        pat[good_cov]
    )
else:
    rho_cov = np.nan
    p_cov = np.nan


positive = int(
    np.sum(
        delta[
            np.isfinite(delta)
        ] > 0
    )
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print("=" * 80)
print("24-DONOR CROSS-TECHNOLOGY COHORT SUMMARY")
print("=" * 80)

print(f"donors                         : {n}")
print(
    "patterned cosine              : "
    f"{pat_med:.3f} "
    f"[IQR {pat_q1:.3f}-{pat_q3:.3f}]"
)
print(
    "non-patterned cosine          : "
    f"{non_med:.3f} "
    f"[IQR {non_q1:.3f}-{non_q3:.3f}]"
)
print(
    "patterned - non-patterned     : "
    f"{delta_med:.3f} "
    f"[IQR {delta_q1:.3f}-{delta_q3:.3f}]"
)
print(
    "donors with positive delta    : "
    f"{positive}/{n}"
)
print(
    "paired Wilcoxon, patterned > "
    f"non-patterned: p={wilcoxon_p:.3e}"
)
print(
    "pseudospot coverage           : "
    f"{cov_med*100:.1f}% "
    f"[IQR {cov_q1*100:.1f}-{cov_q3*100:.1f}%]"
)
print(
    "coverage vs patterned cosine  : "
    f"Spearman rho={rho_cov:+.3f}, "
    f"p={p_cov:.3g}"
)


# ============================================================
# SAVE DONOR-LEVEL TABLE
# ============================================================

table_cols = [
    "BrNum",
    "n_visium_spots",
    "n_occupied_pseudospots",
    "pseudospot_coverage",
    "n_retained_genes",
    "n_patterned_genes",
    "median_cosine_patterned",
    "median_cosine_nonpatterned",
    "delta_cosine",
    "median_spearman_patterned",
]

df[
    table_cols
].to_csv(
    OUT / "crosstech_cohort_donor_metrics.csv",
    index=False
)


# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(13.8, 10.5)
)

gs = fig.add_gridspec(
    2,
    2,
    hspace=0.38,
    wspace=0.30
)


# ------------------------------------------------------------
# A. PAIRED PATTERNED vs NON-PATTERNED
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[0, 0]
)

rng = np.random.default_rng(42)

for _, r in df.iterrows():

    x = np.array(
        [0.0, 1.0]
    )

    y = np.array([
        r["median_cosine_nonpatterned"],
        r["median_cosine_patterned"],
    ])

    jitter = rng.normal(
        0,
        0.018
    )

    xx = x + jitter

    ax.plot(
        xx,
        y,
        linewidth=0.8,
        alpha=0.45
    )

    ax.scatter(
        xx,
        y,
        s=22,
        zorder=3
    )


ax.scatter(
    [0, 1],
    [non_med, pat_med],
    s=110,
    marker="D",
    edgecolor="black",
    linewidth=0.8,
    zorder=5
)

ax.set_xlim(
    -0.35,
    1.35
)

ax.set_ylim(
    0,
    max(
        0.65,
        np.nanmax(
            np.r_[pat, non]
        ) + 0.05
    )
)

ax.set_xticks(
    [0, 1]
)

ax.set_xticklabels([
    "Not patterned",
    "Spatially patterned"
])

ax.set_ylabel(
    "Median cosine similarity"
)

ax.set_title(
    "a   Cross-technology correspondence",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"{positive}/{n} donors: patterned > non-patterned\n"
    f"paired Wilcoxon p={wilcoxon_p:.2e}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ------------------------------------------------------------
# B. DELTA COSINE BY DONOR
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[0, 1]
)

d = df.sort_values(
    "delta_cosine",
    ascending=True
).reset_index(drop=True)

y = np.arange(
    len(d)
)

ax.barh(
    y,
    d["delta_cosine"]
)

ax.axvline(
    0,
    linewidth=1
)

ax.axvline(
    delta_med,
    linestyle="--",
    linewidth=1.3
)

ax.set_yticks(
    y
)

ax.set_yticklabels(
    d["BrNum"],
    fontsize=7.5
)

ax.set_xlabel(
    "Δ cosine = patterned − non-patterned"
)

ax.set_title(
    "b   Patterned/non-patterned separation",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.98,
    0.03,
    f"median Δ = {delta_med:.3f}\n"
    f"IQR {delta_q1:.3f}-{delta_q3:.3f}",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=9
)


# ------------------------------------------------------------
# C. PSEUDOSPOT COVERAGE vs COSINE
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[1, 0]
)

ax.scatter(
    100 * df["pseudospot_coverage"],
    df["median_cosine_patterned"],
    s=45
)

# annotate all donors lightly
for _, r in df.iterrows():

    ax.annotate(
        r["BrNum"],
        (
            100 * r["pseudospot_coverage"],
            r["median_cosine_patterned"]
        ),
        xytext=(3, 3),
        textcoords="offset points",
        fontsize=6.5,
        alpha=0.75
    )


ax.set_xlabel(
    "Visium spots receiving ≥1 Xenium cell (%)"
)

ax.set_ylabel(
    "Median cosine, spatially patterned genes"
)

ax.set_title(
    "c   Coverage versus correspondence",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"Spearman ρ={rho_cov:+.2f}\n"
    f"p={p_cov:.3g}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ------------------------------------------------------------
# D. NUMBER OF PATTERNED GENES
# ------------------------------------------------------------

ax = fig.add_subplot(
    gs[1, 1]
)

d2 = df.sort_values(
    "BrNum"
).reset_index(drop=True)

x = np.arange(
    len(d2)
)

ax.bar(
    x,
    d2["n_patterned_genes"]
)

ax.axhline(
    np.nanmedian(
        d2["n_patterned_genes"]
    ),
    linestyle="--",
    linewidth=1.2
)

ax.set_xticks(
    x
)

ax.set_xticklabels(
    d2["BrNum"],
    rotation=70,
    ha="right",
    fontsize=7
)

ax.set_ylabel(
    "Spatially patterned genes"
)

ax.set_title(
    "d   Spatially patterned genes by donor",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.98,
    0.97,
    f"median = "
    f"{np.nanmedian(d2['n_patterned_genes']):.0f}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=9
)


# ============================================================
# TITLE / FOOTNOTE
# ============================================================

fig.suptitle(
    "Cross-technology spatial expression correspondence across 24 donors",
    fontsize=14,
    fontweight="bold",
    y=0.98
)

fig.text(
    0.5,
    0.015,
    "Per donor, measured Xenium genes were compared with the paired Visium section after "
    "STalign registration. Xenium cells passing the matching threshold were aggregated "
    "into Visium-centered pseudospots. Genes were classified as spatially patterned "
    "from Visium Moran's I with BH-adjusted q < 0.05. Diamond markers in panel a denote "
    "the cohort medians.",
    ha="center",
    va="bottom",
    fontsize=8.5,
    wrap=True
)

plt.subplots_adjust(
    bottom=0.10,
    top=0.92
)


# ============================================================
# SAVE
# ============================================================

png = (
    OUT
    / "fig_crosstech_cohort_24donors.png"
)

pdf = (
    OUT
    / "fig_crosstech_cohort_24donors.pdf"
)

fig.savefig(
    png,
    dpi=600,
    bbox_inches="tight"
)

fig.savefig(
    pdf,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# SAVE COHORT STATS
# ============================================================

stats = pd.DataFrame([
    {
        "n_donors": n,
        "median_cosine_patterned": pat_med,
        "patterned_q1": pat_q1,
        "patterned_q3": pat_q3,
        "median_cosine_nonpatterned": non_med,
        "nonpatterned_q1": non_q1,
        "nonpatterned_q3": non_q3,
        "median_delta_cosine": delta_med,
        "delta_q1": delta_q1,
        "delta_q3": delta_q3,
        "n_positive_delta": positive,
        "wilcoxon_statistic": wilcoxon_stat,
        "wilcoxon_p_greater": wilcoxon_p,
        "median_pseudospot_coverage": cov_med,
        "coverage_q1": cov_q1,
        "coverage_q3": cov_q3,
        "coverage_vs_patterned_cosine_spearman": rho_cov,
        "coverage_vs_patterned_cosine_p": p_cov,
        "median_n_patterned_genes": np.nanmedian(npat),
    }
])

stats.to_csv(
    OUT / "crosstech_cohort_statistics.csv",
    index=False
)

print()
print("=" * 80)
print("OUTPUT")
print("=" * 80)
print("PNG   :", png)
print("PDF   :", pdf)
print(
    "TABLE :",
    OUT / "crosstech_cohort_donor_metrics.csv"
)
print(
    "STATS :",
    OUT / "crosstech_cohort_statistics.csv"
)
