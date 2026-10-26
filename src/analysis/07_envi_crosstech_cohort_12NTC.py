#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import wilcoxon


# ============================================================
# CONFIG
# ============================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

INFILE = (
    ROOT
    / "results/analysis/envi_fulltranscriptome_crosstech"
    / "ENVI_fulltranscriptome_summary_12NTC.csv"
)

OUT = (
    ROOT
    / "results/analysis/envi_fulltranscriptome_crosstech"
    / "cohort_summary_12NTC"
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
    "n_measured_evaluable",
    "n_imputed_evaluable",
    "n_patterned_measured",
    "n_patterned_imputed",
    "median_cosine_patterned_measured",
    "median_cosine_patterned_imputed",
    "median_spearman_patterned_measured",
    "median_spearman_patterned_imputed",
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

if len(df) != 12:
    raise RuntimeError(
        f"Expected 12 donors, found {len(df)}"
    )

if df["BrNum"].nunique() != 12:
    raise RuntimeError(
        "Duplicate donor IDs found."
    )

df = (
    df
    .sort_values("BrNum")
    .reset_index(drop=True)
)


# ============================================================
# DERIVED METRICS
# ============================================================

df["delta_cosine"] = (
    df["median_cosine_patterned_imputed"]
    -
    df["median_cosine_patterned_measured"]
)

df["delta_spearman"] = (
    df["median_spearman_patterned_imputed"]
    -
    df["median_spearman_patterned_measured"]
)

df["patterned_fraction_measured"] = (
    df["n_patterned_measured"]
    /
    df["n_measured_evaluable"]
)

df["patterned_fraction_imputed"] = (
    df["n_patterned_imputed"]
    /
    df["n_imputed_evaluable"]
)


# ============================================================
# HELPERS
# ============================================================

def med_iqr(x):

    x = np.asarray(
        x,
        dtype=float
    )

    return (
        np.nanmedian(x),
        np.nanpercentile(x, 25),
        np.nanpercentile(x, 75)
    )


# ============================================================
# VALUES
# ============================================================

cm = df[
    "median_cosine_patterned_measured"
].to_numpy(float)

ci = df[
    "median_cosine_patterned_imputed"
].to_numpy(float)

dc = df[
    "delta_cosine"
].to_numpy(float)

sm = df[
    "median_spearman_patterned_measured"
].to_numpy(float)

si = df[
    "median_spearman_patterned_imputed"
].to_numpy(float)

fm = df[
    "patterned_fraction_measured"
].to_numpy(float)

fi = df[
    "patterned_fraction_imputed"
].to_numpy(float)


cm_med, cm_q1, cm_q3 = med_iqr(cm)
ci_med, ci_q1, ci_q3 = med_iqr(ci)
dc_med, dc_q1, dc_q3 = med_iqr(dc)

sm_med, sm_q1, sm_q3 = med_iqr(sm)
si_med, si_q1, si_q3 = med_iqr(si)


# ============================================================
# TESTS
# ============================================================

w_cos = wilcoxon(
    ci,
    cm,
    alternative="two-sided"
)

w_sp = wilcoxon(
    si,
    sm,
    alternative="two-sided"
)

w_frac = wilcoxon(
    fi,
    fm,
    alternative="two-sided"
)


# ============================================================
# PRINT SUMMARY
# ============================================================

print("=" * 80)
print("12-NTC ENVI FULL-TRANSCRIPTOME CROSS-TECH SUMMARY")
print("=" * 80)

print(
    f"patterned measured cosine    : "
    f"{cm_med:.3f} "
    f"[IQR {cm_q1:.3f}-{cm_q3:.3f}]"
)

print(
    f"patterned imputed cosine     : "
    f"{ci_med:.3f} "
    f"[IQR {ci_q1:.3f}-{ci_q3:.3f}]"
)

print(
    f"imputed - measured cosine    : "
    f"{dc_med:+.3f} "
    f"[IQR {dc_q1:+.3f}-{dc_q3:+.3f}]"
)

print(
    f"paired Wilcoxon cosine       : "
    f"p={w_cos.pvalue:.4g}"
)

print()

print(
    f"imputed > measured cosine    : "
    f"{int((dc > 0).sum())}/12"
)

print(
    f"imputed < measured cosine    : "
    f"{int((dc < 0).sum())}/12"
)

print()

print(
    f"patterned measured Spearman  : "
    f"{sm_med:.3f} "
    f"[IQR {sm_q1:.3f}-{sm_q3:.3f}]"
)

print(
    f"patterned imputed Spearman   : "
    f"{si_med:.3f} "
    f"[IQR {si_q1:.3f}-{si_q3:.3f}]"
)

print(
    f"paired Wilcoxon Spearman     : "
    f"p={w_sp.pvalue:.4g}"
)


# ============================================================
# FIGURE
# ============================================================

fig = plt.figure(
    figsize=(13.5, 10.0)
)

gs = fig.add_gridspec(
    2,
    2,
    hspace=0.38,
    wspace=0.32
)


# ============================================================
# PANEL A
# PAIRED COSINE
# ============================================================

ax = fig.add_subplot(
    gs[0, 0]
)

rng = np.random.default_rng(42)

for _, r in df.iterrows():

    jitter = rng.normal(
        0,
        0.015
    )

    xx = np.array([
        0 + jitter,
        1 + jitter
    ])

    yy = np.array([
        r[
            "median_cosine_patterned_measured"
        ],
        r[
            "median_cosine_patterned_imputed"
        ],
    ])

    ax.plot(
        xx,
        yy,
        alpha=0.45,
        linewidth=0.9
    )

    ax.scatter(
        xx,
        yy,
        s=28,
        zorder=3
    )


ax.scatter(
    [0, 1],
    [cm_med, ci_med],
    marker="D",
    s=115,
    edgecolor="black",
    linewidth=0.8,
    zorder=5
)


ax.set_xticks(
    [0, 1]
)

ax.set_xticklabels([
    "Measured Xenium",
    "ENVI imputed-only"
])

ax.set_ylabel(
    "Median cosine similarity\n"
    "(spatially patterned genes)"
)

ax.set_xlim(
    -0.35,
    1.35
)

ax.set_ylim(
    0.30,
    0.65
)

ax.set_title(
    "a   Spatial correspondence",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"Measured median = {cm_med:.3f}\n"
    f"Imputed median = {ci_med:.3f}\n"
    f"paired Wilcoxon p={w_cos.pvalue:.3g}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ============================================================
# PANEL B
# DELTA COSINE
# ============================================================

ax = fig.add_subplot(
    gs[0, 1]
)

d = (
    df
    .sort_values(
        "delta_cosine",
        ascending=True
    )
    .reset_index(drop=True)
)

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
    dc_med,
    linestyle="--",
    linewidth=1.2
)

ax.set_yticks(
    y
)

ax.set_yticklabels(
    d["BrNum"],
    fontsize=8
)

ax.set_xlabel(
    "Δ cosine = imputed − measured"
)

ax.set_title(
    "b   Donor-level difference",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.98,
    0.03,
    f"median Δ = {dc_med:+.3f}\n"
    f"IQR {dc_q1:+.3f} to {dc_q3:+.3f}\n"
    f"imputed higher: {(dc > 0).sum()}/12",
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=9
)


# ============================================================
# PANEL C
# FRACTION PATTERNED
# ============================================================

ax = fig.add_subplot(
    gs[1, 0]
)

for _, r in df.iterrows():

    ax.plot(
        [0, 1],
        [
            100
            *
            r[
                "patterned_fraction_measured"
            ],
            100
            *
            r[
                "patterned_fraction_imputed"
            ]
        ],
        alpha=0.40,
        linewidth=0.8
    )

    ax.scatter(
        [0, 1],
        [
            100
            *
            r[
                "patterned_fraction_measured"
            ],
            100
            *
            r[
                "patterned_fraction_imputed"
            ]
        ],
        s=25
    )


fm_med = np.nanmedian(
    100 * fm
)

fi_med = np.nanmedian(
    100 * fi
)


ax.scatter(
    [0, 1],
    [fm_med, fi_med],
    marker="D",
    s=110,
    edgecolor="black",
    linewidth=0.8,
    zorder=5
)

ax.set_xticks(
    [0, 1]
)

ax.set_xticklabels([
    "Measured Xenium",
    "ENVI imputed-only"
])

ax.set_ylabel(
    "Spatially patterned genes\n"
    "(% of evaluable genes)"
)

ax.set_xlim(
    -0.35,
    1.35
)

ax.set_title(
    "c   Visium-defined spatial patterning",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"Measured median = {fm_med:.1f}%\n"
    f"Imputed median = {fi_med:.1f}%\n"
    f"paired Wilcoxon p={w_frac.pvalue:.3g}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ============================================================
# PANEL D
# PAIRED SPEARMAN
# ============================================================

ax = fig.add_subplot(
    gs[1, 1]
)

for _, r in df.iterrows():

    ax.plot(
        [0, 1],
        [
            r[
                "median_spearman_patterned_measured"
            ],
            r[
                "median_spearman_patterned_imputed"
            ]
        ],
        alpha=0.45,
        linewidth=0.9
    )

    ax.scatter(
        [0, 1],
        [
            r[
                "median_spearman_patterned_measured"
            ],
            r[
                "median_spearman_patterned_imputed"
            ]
        ],
        s=28
    )


ax.scatter(
    [0, 1],
    [sm_med, si_med],
    marker="D",
    s=115,
    edgecolor="black",
    linewidth=0.8,
    zorder=5
)


ax.axhline(
    0,
    linewidth=0.7
)

ax.set_xticks(
    [0, 1]
)

ax.set_xticklabels([
    "Measured Xenium",
    "ENVI imputed-only"
])

ax.set_ylabel(
    "Median Spearman correlation\n"
    "(spatially patterned genes)"
)

ax.set_xlim(
    -0.35,
    1.35
)

ax.set_title(
    "d   Spot-wise rank correspondence",
    loc="left",
    fontweight="bold"
)

ax.text(
    0.03,
    0.97,
    f"Measured median = {sm_med:.3f}\n"
    f"Imputed median = {si_med:.3f}\n"
    f"paired Wilcoxon p={w_sp.pvalue:.3g}",
    transform=ax.transAxes,
    va="top",
    fontsize=9
)


# ============================================================
# MAIN TITLE / FOOTNOTE
# ============================================================

fig.suptitle(
    "Spatial validation of ENVI full-transcriptome imputation across 12 NTC donors",
    fontsize=14,
    fontweight="bold",
    y=0.98
)

fig.text(
    0.5,
    0.015,
    "ENVI count-scale expression was aggregated from aligned Xenium cells whose centroids "
    "fell inside physical Visium spots. Correspondence was evaluated only for genes shared "
    "with paired Visium and sufficiently observed for evaluation. Spatial patterning was "
    "defined independently from Visium Moran's I (BH q < 0.05). Measured and imputed-only "
    "groups contain different genes; therefore the paired donor comparison reflects the "
    "magnitude of spatial correspondence, not a head-to-head accuracy test on identical genes.",
    ha="center",
    va="bottom",
    fontsize=8.5,
    wrap=True
)

plt.subplots_adjust(
    bottom=0.11,
    top=0.92
)


# ============================================================
# SAVE
# ============================================================

PNG = (
    OUT
    / "fig_ENVI_fulltranscriptome_12NTC_summary.png"
)

PDF = (
    OUT
    / "fig_ENVI_fulltranscriptome_12NTC_summary.pdf"
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


# ============================================================
# SAVE DONOR TABLE
# ============================================================

OUT_TABLE = (
    OUT
    / "ENVI_fulltranscriptome_12NTC_donor_metrics.csv"
)

df.to_csv(
    OUT_TABLE,
    index=False
)


# ============================================================
# SAVE COHORT STATS
# ============================================================

stats = pd.DataFrame([
    {
        "n_donors":
            12,

        "median_cosine_patterned_measured":
            cm_med,

        "measured_cosine_q1":
            cm_q1,

        "measured_cosine_q3":
            cm_q3,

        "median_cosine_patterned_imputed":
            ci_med,

        "imputed_cosine_q1":
            ci_q1,

        "imputed_cosine_q3":
            ci_q3,

        "median_delta_cosine":
            dc_med,

        "delta_cosine_q1":
            dc_q1,

        "delta_cosine_q3":
            dc_q3,

        "n_imputed_higher":
            int(
                (dc > 0).sum()
            ),

        "n_measured_higher":
            int(
                (dc < 0).sum()
            ),

        "wilcoxon_cosine_p":
            w_cos.pvalue,

        "median_spearman_measured":
            sm_med,

        "median_spearman_imputed":
            si_med,

        "wilcoxon_spearman_p":
            w_sp.pvalue,

        "median_patterned_fraction_measured":
            np.nanmedian(fm),

        "median_patterned_fraction_imputed":
            np.nanmedian(fi),

        "wilcoxon_patterned_fraction_p":
            w_frac.pvalue,
    }
])

STATS = (
    OUT
    / "ENVI_fulltranscriptome_12NTC_cohort_statistics.csv"
)

stats.to_csv(
    STATS,
    index=False
)


print()
print("=" * 80)
print("OUTPUT")
print("=" * 80)

print(
    "PNG   :",
    PNG
)

print(
    "PDF   :",
    PDF
)

print(
    "TABLE :",
    OUT_TABLE
)

print(
    "STATS :",
    STATS
)
