from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


# -------------------------
# Paths
# -------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

XENIUM_DIR = PROJECT_ROOT / "data" / "raw" / "xenium"
TABLE_DIR = PROJECT_ROOT / "outputs" / "tables"
FIGURE_DIR = PROJECT_ROOT / "outputs" / "figures"

TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


# -------------------------
# Load Xenium expression matrix
# -------------------------

xenium_h5 = list(XENIUM_DIR.rglob("*cell_feature_matrix.h5"))[0]

xenium = sc.read_10x_h5(xenium_h5)
xenium.var_names_make_unique()

print("Loaded Xenium:")
print(xenium)


# -------------------------
# Load QC metrics from previous step
# -------------------------

qc_file = TABLE_DIR / "xenium_cell_qc_metrics.csv"

if not qc_file.exists():
    raise FileNotFoundError(
        "xenium_cell_qc_metrics.csv not found. "
        "Run src/analysis/xenium_qc.py first."
    )

qc = pd.read_csv(qc_file)

if "cell_id" not in qc.columns:
    qc = qc.rename(columns={qc.columns[0]: "cell_id"})

qc["cell_id"] = qc["cell_id"].astype(str)
xenium.obs_names = xenium.obs_names.astype(str)


# -------------------------
# Match cells
# -------------------------

common_cells = sorted(
    set(qc["cell_id"])
    &
    set(xenium.obs_names)
)

print("Cells in matrix:", xenium.n_obs)
print("Cells in QC file:", qc.shape[0])
print("Matched cells:", len(common_cells))

if len(common_cells) == 0:
    raise ValueError("No matching cell IDs between Xenium matrix and QC file.")

xenium = xenium[common_cells, :].copy()

qc = (
    qc
    .set_index("cell_id")
    .loc[common_cells]
    .copy()
)

for col in qc.columns:
    xenium.obs[col] = qc[col].values


# -------------------------
# Remove outlier cells
# -------------------------

if "is_outlier_lieber_like" in xenium.obs.columns:
    keep_cells = ~xenium.obs["is_outlier_lieber_like"].astype(bool)
    xenium = xenium[keep_cells, :].copy()
    print("After removing Lieber-like outliers:")
    print(xenium)
else:
    print("WARNING: is_outlier_lieber_like not found. PCA will use all cells.")


# -------------------------
# Basic preprocessing
# -------------------------

sc.pp.filter_genes(xenium, min_cells=3)

sc.pp.normalize_total(
    xenium,
    target_sum=1e4
)

sc.pp.log1p(xenium)

sc.pp.highly_variable_genes(
    xenium,
    n_top_genes=200,
    flavor="seurat"
)

xenium = xenium[:, xenium.var["highly_variable"]].copy()

sc.pp.scale(
    xenium,
    max_value=10
)


# -------------------------
# PCA
# -------------------------

sc.tl.pca(
    xenium,
    n_comps=10,
    svd_solver="arpack"
)

pca_coords = pd.DataFrame(
    xenium.obsm["X_pca"],
    index=xenium.obs_names,
    columns=[f"PC{i+1}" for i in range(xenium.obsm["X_pca"].shape[1])]
)

pca_coords.to_csv(TABLE_DIR / "xenium_pca_coordinates.csv")


# -------------------------
# Save variance explained
# -------------------------

variance = pd.DataFrame({
    "PC": [f"PC{i+1}" for i in range(len(xenium.uns["pca"]["variance_ratio"]))],
    "variance_ratio": xenium.uns["pca"]["variance_ratio"]
})

variance.to_csv(
    TABLE_DIR / "xenium_pca_variance_explained.csv",
    index=False
)


# -------------------------
# Plot PCA colored by QC metrics
# -------------------------

def plot_pca(color_col, filename):
    if color_col not in xenium.obs.columns:
        print(f"Skipping {color_col}: column not found.")
        return

    plt.figure(figsize=(7, 6))

    values = xenium.obs[color_col]

    if values.dtype == bool:
        values = values.astype(int)

    scatter = plt.scatter(
        xenium.obsm["X_pca"][:, 0],
        xenium.obsm["X_pca"][:, 1],
        c=values,
        s=2,
        alpha=0.8,
    )

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"PCA colored by {color_col}")
    plt.colorbar(scatter, label=color_col)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / filename, dpi=300)
    plt.close()


plot_pca("matrix_total_counts", "xenium_pca_total_counts.png")
plot_pca("n_genes_by_counts", "xenium_pca_detected_genes.png")
plot_pca("pct_control_probe", "xenium_pca_control_probe.png")
plot_pca("pct_control_codeword", "xenium_pca_control_codeword.png")
plot_pca("pct_unassigned", "xenium_pca_unassigned.png")
plot_pca("pct_counts_mt", "xenium_pca_mito.png")
plot_pca("cell_area", "xenium_pca_cell_area.png")
plot_pca("nucleus_area", "xenium_pca_nucleus_area.png")


# -------------------------
# Spatial PCA plots
# -------------------------

def plot_spatial(color_col, filename):
    required = ["x_centroid", "y_centroid"]

    for col in required:
        if col not in xenium.obs.columns:
            print(f"Skipping spatial {color_col}: {col} missing.")
            return

    if color_col not in xenium.obs.columns:
        print(f"Skipping spatial {color_col}: column not found.")
        return

    plt.figure(figsize=(7, 7))

    values = xenium.obs[color_col]

    if values.dtype == bool:
        values = values.astype(int)

    scatter = plt.scatter(
        xenium.obs["x_centroid"],
        xenium.obs["y_centroid"],
        c=values,
        s=1,
        alpha=0.8,
    )

    plt.gca().invert_yaxis()
    plt.xlabel("x_centroid")
    plt.ylabel("y_centroid")
    plt.title(f"Spatial map colored by {color_col}")
    plt.colorbar(scatter, label=color_col)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / filename, dpi=300)
    plt.close()


plot_spatial("matrix_total_counts", "xenium_spatial_total_counts.png")
plot_spatial("n_genes_by_counts", "xenium_spatial_detected_genes.png")
plot_spatial("pct_control_probe", "xenium_spatial_control_probe.png")
plot_spatial("pct_control_codeword", "xenium_spatial_control_codeword.png")
plot_spatial("pct_unassigned", "xenium_spatial_unassigned.png")
plot_spatial("pct_counts_mt", "xenium_spatial_mito.png")
plot_spatial("PC1", "xenium_spatial_PC1.png")
plot_spatial("PC2", "xenium_spatial_PC2.png")


# Add PC1/PC2 into obs and replot spatial PCs
xenium.obs["PC1"] = xenium.obsm["X_pca"][:, 0]
xenium.obs["PC2"] = xenium.obsm["X_pca"][:, 1]

plot_spatial("PC1", "xenium_spatial_PC1.png")
plot_spatial("PC2", "xenium_spatial_PC2.png")


# -------------------------
# Save processed PCA object
# -------------------------

out_h5ad = PROJECT_ROOT / "data" / "processed" / "xenium_pca_qc.h5ad"
out_h5ad.parent.mkdir(parents=True, exist_ok=True)

xenium.write_h5ad(out_h5ad)


print("\nSaved tables:")
print(TABLE_DIR / "xenium_pca_coordinates.csv")
print(TABLE_DIR / "xenium_pca_variance_explained.csv")

print("\nSaved figures to:")
print(FIGURE_DIR)

print("\nSaved AnnData:")
print(out_h5ad)