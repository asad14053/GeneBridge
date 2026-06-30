
#!/usr/bin/env python

from pathlib import Path
import re
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

visium_path = PROJECT_ROOT / "data/processed/visium/spatialDLPFC_Visium_sce.h5ad"
snrna_meta_path = PROJECT_ROOT / "outputs/huuki_myers/tables/huuki_snrna_metadata.csv"

out_table_dir = PROJECT_ROOT / "outputs/huuki_myers/tables"
out_fig_dir = PROJECT_ROOT / "outputs/huuki_myers/figures"

out_table_dir.mkdir(parents=True, exist_ok=True)
out_fig_dir.mkdir(parents=True, exist_ok=True)

def normalize_brain_id(x):
    """
    Converts values like Br8667, 8667, br_8667, etc. into Br8667.
    Returns None if no 4-digit brain number is found.
    """
    if pd.isna(x):
        return None
    s = str(x).strip()
    m = re.search(r"(\d{4})", s)
    if m is None:
        return None
    return "Br" + m.group(1)

def find_brain_column(obs):
    candidates = [
        "BrNum",
        "brnum",
        "BrNum_matched",
        "brain",
        "Brain",
        "brain_id",
        "BrainID",
        "sample_id",
        "sample",
        "Sample",
        "subject",
        "Subject",
        "donor",
        "Donor",
        "donor_id",
        "SAMPLE_ID",
    ]

    for c in candidates:
        if c in obs.columns:
            vals = obs[c].astype(str).map(normalize_brain_id)
            if vals.notna().sum() > 0:
                return c

    possible = [
        c for c in obs.columns
        if any(k in c.lower() for k in ["br", "brain", "sample", "donor", "subject"])
    ]

    print("\nCould not automatically find brain ID column.")
    print("\nPossible columns:")
    print(possible)
    print("\nAll obs columns:")
    print(obs.columns.tolist())

    raise ValueError("No usable brain ID column found in Visium obs.")

print("Loading Huuki Visium h5ad:")
print(visium_path)
vis = sc.read_h5ad(visium_path)

print("\nHuuki Visium object:")
print(vis)

print("\nLoading Huuki snRNA metadata:")
print(snrna_meta_path)
sn_meta = pd.read_csv(snrna_meta_path)

print("\nHuuki snRNA metadata shape:")
print(sn_meta.shape)

# -----------------------------
# Brain IDs
# -----------------------------
vis_br_col = find_brain_column(vis.obs)

if "BrNum" not in sn_meta.columns:
    raise ValueError("snRNA metadata has no BrNum column.")

vis_br = vis.obs[vis_br_col].map(normalize_brain_id)
sn_br = sn_meta["BrNum"].map(normalize_brain_id)

vis.obs["BrNum_clean"] = vis_br.values
sn_meta["BrNum_clean"] = sn_br.values

vis_br_set = set(vis_br.dropna())
sn_br_set = set(sn_br.dropna())
common_brains = sorted(vis_br_set & sn_br_set)

print("\nVisium brain column used:", vis_br_col)
print("Visium unique brains:", len(vis_br_set))
print("snRNA unique brains:", len(sn_br_set))
print("Common brains:", len(common_brains))
print(common_brains)

# -----------------------------
# Tables
# -----------------------------
common_df = pd.DataFrame({"BrNum": common_brains})
common_df.to_csv(out_table_dir / "huuki_visium_snrna_common_brains.csv", index=False)

vis_counts = (
    vis.obs["BrNum_clean"]
    .dropna()
    .value_counts()
    .rename_axis("BrNum")
    .reset_index(name="visium_spots")
)

sn_counts = (
    sn_meta["BrNum_clean"]
    .dropna()
    .value_counts()
    .rename_axis("BrNum")
    .reset_index(name="snrna_cells")
)

summary = pd.merge(vis_counts, sn_counts, on="BrNum", how="outer").fillna(0)
summary["visium_spots"] = summary["visium_spots"].astype(int)
summary["snrna_cells"] = summary["snrna_cells"].astype(int)
summary["is_common"] = summary["BrNum"].isin(common_brains)
summary = summary.sort_values(["is_common", "BrNum"], ascending=[False, True])

summary.to_csv(out_table_dir / "huuki_visium_snrna_brain_summary.csv", index=False)
vis_counts.to_csv(out_table_dir / "huuki_visium_brain_counts.csv", index=False)
sn_counts.to_csv(out_table_dir / "huuki_snrna_brain_counts_from_metadata.csv", index=False)

# -----------------------------
# Figure 1: Visium spots per brain
# -----------------------------
plot_df = vis_counts.sort_values("visium_spots", ascending=False)

plt.figure(figsize=(12, 5))
plt.bar(plot_df["BrNum"], plot_df["visium_spots"])
plt.xticks(rotation=90)
plt.ylabel("Visium spots")
plt.title("Huuki-Myers Visium spots per brain")
plt.tight_layout()
plt.savefig(out_fig_dir / "task1b_huuki_visium_spots_per_brain.png", dpi=300)
plt.close()

# -----------------------------
# Figure 2: snRNA cells per brain
# -----------------------------
plot_df = sn_counts.sort_values("snrna_cells", ascending=False)

plt.figure(figsize=(12, 5))
plt.bar(plot_df["BrNum"], plot_df["snrna_cells"])
plt.xticks(rotation=90)
plt.ylabel("snRNA-seq cells")
plt.title("Huuki-Myers snRNA-seq cells per brain")
plt.tight_layout()
plt.savefig(out_fig_dir / "task1b_huuki_snrna_cells_per_brain.png", dpi=300)
plt.close()

# -----------------------------
# Figure 3: Common brain comparison
# -----------------------------
common_summary = summary[summary["is_common"]].copy()
common_summary = common_summary.sort_values("BrNum")

plt.figure(figsize=(10, 5))
x = list(range(len(common_summary)))
plt.bar([i - 0.2 for i in x], common_summary["visium_spots"], width=0.4, label="Visium spots")
plt.bar([i + 0.2 for i in x], common_summary["snrna_cells"], width=0.4, label="snRNA cells")
plt.xticks(x, common_summary["BrNum"], rotation=45)
plt.ylabel("Count")
plt.title("Huuki-Myers common brain IDs")
plt.legend()
plt.tight_layout()
plt.savefig(out_fig_dir / "task1c_huuki_common_brains_visium_snrna_counts.png", dpi=300)
plt.close()

# -----------------------------
# Figure 4: Visium spatial plot, common vs non-common
# -----------------------------
if "spatial" in vis.obsm:
    coords = pd.DataFrame(vis.obsm["spatial"], columns=["x", "y"], index=vis.obs_names)
    coords["is_common"] = vis.obs["BrNum_clean"].isin(common_brains).values

    plt.figure(figsize=(7, 7))
    plt.scatter(
        coords.loc[~coords["is_common"], "x"],
        coords.loc[~coords["is_common"], "y"],
        s=1,
        alpha=0.15,
        label="Other brains"
    )
    plt.scatter(
        coords.loc[coords["is_common"], "x"],
        coords.loc[coords["is_common"], "y"],
        s=1,
        alpha=0.4,
        label="Common brains"
    )
    plt.gca().invert_yaxis()
    plt.xlabel("spatial x")
    plt.ylabel("spatial y")
    plt.title("Huuki Visium spatial coordinates: common vs other brains")
    plt.legend(markerscale=8)
    plt.tight_layout()
    plt.savefig(out_fig_dir / "task1b_huuki_visium_spatial_common_vs_other.png", dpi=300)
    plt.close()
else:
    print("\nWARNING: vis.obsm['spatial'] not found. Skipping Visium spatial plot.")

# -----------------------------
# Figure 5: snRNA UMAP, common vs non-common
# -----------------------------
if {"UMAP1", "UMAP2"}.issubset(sn_meta.columns):
    sn_meta["is_common"] = sn_meta["BrNum_clean"].isin(common_brains)

    plt.figure(figsize=(7, 6))
    plt.scatter(
        sn_meta.loc[~sn_meta["is_common"], "UMAP1"],
        sn_meta.loc[~sn_meta["is_common"], "UMAP2"],
        s=1,
        alpha=0.15,
        label="Other brains"
    )
    plt.scatter(
        sn_meta.loc[sn_meta["is_common"], "UMAP1"],
        sn_meta.loc[sn_meta["is_common"], "UMAP2"],
        s=1,
        alpha=0.4,
        label="Common brains"
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Huuki snRNA-seq UMAP: common vs other brains")
    plt.legend(markerscale=8)
    plt.tight_layout()
    plt.savefig(out_fig_dir / "task1b_huuki_snrna_umap_common_vs_other.png", dpi=300)
    plt.close()
else:
    print("\nWARNING: UMAP1/UMAP2 not found in snRNA metadata. Skipping snRNA UMAP plot.")

# -----------------------------
# Save short text report
# -----------------------------
report_path = out_table_dir / "task1_huuki_common_brains_report.txt"

report = [
    "Task-1 Huuki-Myers common brain report",
    "",
    f"Visium path: {visium_path}",
    f"snRNA metadata path: {snrna_meta_path}",
    "",
    f"Visium brain column used: {vis_br_col}",
    f"Visium unique brains: {len(vis_br_set)}",
    f"snRNA unique brains: {len(sn_br_set)}",
    f"Common brains: {len(common_brains)}",
    "",
    "Common brain IDs:",
    ", ".join(common_brains),
    "",
    "Output tables:",
    str(out_table_dir / "huuki_visium_snrna_common_brains.csv"),
    str(out_table_dir / "huuki_visium_snrna_brain_summary.csv"),
    "",
    "Output figures:",
    str(out_fig_dir / "task1b_huuki_visium_spots_per_brain.png"),
    str(out_fig_dir / "task1b_huuki_snrna_cells_per_brain.png"),
    str(out_fig_dir / "task1c_huuki_common_brains_visium_snrna_counts.png"),
    str(out_fig_dir / "task1b_huuki_visium_spatial_common_vs_other.png"),
    str(out_fig_dir / "task1b_huuki_snrna_umap_common_vs_other.png"),
]

report_path.write_text("\n".join(report))

print("\nSaved tables:")
print(out_table_dir / "huuki_visium_snrna_common_brains.csv")
print(out_table_dir / "huuki_visium_snrna_brain_summary.csv")
print(out_table_dir / "task1_huuki_common_brains_report.txt")

print("\nSaved figures to:")
print(out_fig_dir)

print("\nDONE Task-1.")