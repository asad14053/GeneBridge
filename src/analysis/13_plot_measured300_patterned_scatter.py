#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")
BR = os.environ.get("BR", "Br8667")

MORAN_CSV = (
    ROOT
    / "results/analysis/measured300_moran_visium_vs_xenium"
    / BR
    / f"{BR}_measured300_Visium_vs_Xenium_MoranI.csv"
)

GENE_METRICS = (
    ROOT
    / "results/analysis/crosstech_measured300"
    / BR
    / f"{BR}_measured300_gene_metrics.csv"
)

OUTDIR = (
    ROOT
    / "results/analysis/measured300_moran_visium_vs_xenium"
    / BR
)
OUTDIR.mkdir(parents=True, exist_ok=True)


def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def as_bool(s):
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    x = s.astype(str).str.strip().str.lower()
    return x.isin(["true", "1", "t", "yes"])


m = pd.read_csv(MORAN_CSV)
g = pd.read_csv(GENE_METRICS)

gene_col_m = find_col(m, ["gene", "Gene"])
gene_col_g = find_col(g, ["gene", "Gene"])

if gene_col_m is None or gene_col_g is None:
    raise RuntimeError("Could not find gene column.")

m = m.rename(columns={gene_col_m: "gene"})
g = g.rename(columns={gene_col_g: "gene"})

pattern_col = find_col(
    g,
    [
        "visium_patterned",
        "patterned",
        "spatially_patterned",
        "patterned_visium",
    ],
)

if pattern_col is not None:
    g["visium_patterned"] = as_bool(g[pattern_col])
else:
    q_col = find_col(
        g,
        [
            "visium_moran_q",
            "moran_q_visium",
            "q_value_visium",
            "bh_q_visium",
            "padj_moran_visium",
        ],
    )
    if q_col is None:
        raise RuntimeError(
            "Could not find a Visium-patterned flag or Visium Moran q-value column "
            f"in {GENE_METRICS}"
        )
    g["visium_patterned"] = pd.to_numeric(g[q_col], errors="coerce") < 0.05

df = m.merge(
    g[["gene", "visium_patterned"]],
    on="gene",
    how="left"
)

df["visium_patterned"] = df["visium_patterned"].fillna(False)

x = pd.to_numeric(df["moran_visium"], errors="coerce")
y = pd.to_numeric(df["moran_xenium_measured"], errors="coerce")

ok = np.isfinite(x) & np.isfinite(y)
df = df.loc[ok].copy()

x = df["moran_visium"].to_numpy(float)
y = df["moran_xenium_measured"].to_numpy(float)
patterned = df["visium_patterned"].to_numpy(bool)

rho, _ = spearmanr(x, y)
r, _ = pearsonr(x, y)
slope, intercept = np.polyfit(x, y, 1)
delta = y - x

fig, ax = plt.subplots(figsize=(7.8, 7.8))

# non-patterned
mask0 = ~patterned
ax.scatter(
    x[mask0],
    y[mask0],
    s=20,
    alpha=0.30,
    linewidths=0,
    label="Not Visium-patterned",
)

# patterned
mask1 = patterned
ax.scatter(
    x[mask1],
    y[mask1],
    s=22,
    alpha=0.70,
    linewidths=0,
    label="Visium-patterned",
)

lo = min(np.nanmin(x), np.nanmin(y))
hi = max(np.nanmax(x), np.nanmax(y))
pad = max(0.02, 0.05 * (hi - lo))
lo -= pad
hi += pad

# y=x
ax.plot(
    [lo, hi],
    [lo, hi],
    linestyle="--",
    linewidth=1.5,
    label="y = x",
)

# fit
xx = np.linspace(lo, hi, 200)
ax.plot(
    xx,
    intercept + slope * xx,
    linewidth=1.5,
    label="Gene-level fit",
)

ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)

ax.set_xlabel("Visium Moran's I", fontsize=14)
ax.set_ylabel("Measured Xenium pseudospot Moran's I", fontsize=14)

ax.set_title(
    f"{BR}: preservation of spatial autocorrelation",
    fontsize=18,
    fontweight="bold",
    pad=10,
)

txt = (
    f"Shared measured genes: {len(df):,}\n"
    f"Valid Moran's I in both: {len(df):,}\n"
    f"Spearman ρ={rho:.3f}\n"
    f"Pearson r={r:.3f}\n"
    f"median ΔI={np.median(delta):+.3f}\n"
    f"slope={slope:.3f}"
)

ax.text(
    0.03,
    0.97,
    txt,
    transform=ax.transAxes,
    va="top",
    ha="left",
    fontsize=12,
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.9),
)

ax.legend(loc="lower right", frameon=False, fontsize=12)

fig.tight_layout()

png = OUTDIR / f"fig_{BR}_measured300_patterned_scatter.png"
pdf = OUTDIR / f"fig_{BR}_measured300_patterned_scatter.pdf"

fig.savefig(png, dpi=600, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

print("Saved:")
print(png)
print(pdf)
