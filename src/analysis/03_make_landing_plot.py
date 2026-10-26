#!/usr/bin/env python3

from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROJECT_ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

DEFAULT_CSV = (
    PROJECT_ROOT
    / "results"
    / "analysis"
    / "stalign_final_aligned_qc"
    / "final_aligned_qc_per_layer.csv"
)

DEFAULT_OUT = (
    PROJECT_ROOT
    / "results"
    / "analysis"
    / "stalign_final_aligned_qc"
    / "landing_plots"
)

LAY = ["L1/M", "L2/3", "L3/4", "L5", "L6", "WM", "WMtz"]


def make_one_plot(df, br, outdir):
    sub = df[df["BrNum"] == br].copy()

    if sub.empty:
        raise ValueError(f"No rows found for {br}")

    # keep layer order consistent
    sub["layer"] = pd.Categorical(sub["layer"], categories=LAY, ordered=True)
    sub = sub.sort_values("layer")

    x = np.arange(len(LAY))
    y = sub["mean_landed_visium_spd"].to_numpy(dtype=float)
    sd = sub["sd_landed_visium_spd"].to_numpy(dtype=float)
    n = sub["n"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    # perfect landing diagonal
    ax.plot([0, 6], [0, 6], ls="--", lw=1, color="0.6", label="perfect")

    # GM shading
    ax.axvspan(0.5, 4.5, alpha=0.08)
    ax.text(2.5, 6.35, "GM", ha="center", fontsize=10)

    # observed
    ax.errorbar(
        x, y, yerr=sd,
        fmt="-o",
        lw=1.5,
        ms=5,
        capsize=3,
        label="observed"
    )

    # annotate n
    for xi, yi, ni in zip(x, y, n):
        if np.isfinite(yi):
            ax.text(xi, yi + 0.18, f"n={int(ni)}", ha="center", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(LAY, rotation=45)
    ax.set_yticks(np.arange(7))
    ax.set_yticklabels(LAY)

    ax.set_xlim(-0.3, 6.3)
    ax.set_ylim(-0.3, 6.7)

    ax.set_xlabel("Xenium Lieber layer")
    ax.set_ylabel("Mean landed Visium domain")
    ax.set_title(f"{br} — landing plot")

    ax.legend(fontsize=9)
    plt.tight_layout()

    out_png = outdir / f"{br}_landing_plot.png"
    out_pdf = outdir / f"{br}_landing_plot.pdf"

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", out_png)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--br", type=str, default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    args.out.mkdir(parents=True, exist_ok=True)

    donors = sorted(df["BrNum"].dropna().unique())

    if args.br is not None:
        donors = [args.br]

    for br in donors:
        make_one_plot(df, br, args.out)


if __name__ == "__main__":
    main()
