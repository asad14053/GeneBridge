from pathlib import Path
import gzip
import shutil

import numpy as np
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
# Find files
# -------------------------

xenium_h5_files = list(XENIUM_DIR.rglob("*cell_feature_matrix.h5"))
cell_parquet_gz_files = list(XENIUM_DIR.rglob("*cells.parquet.gz"))

if len(xenium_h5_files) == 0:
    raise FileNotFoundError("No Xenium cell_feature_matrix.h5 file found.")

if len(cell_parquet_gz_files) == 0:
    raise FileNotFoundError("No Xenium cells.parquet.gz file found.")

xenium_h5 = xenium_h5_files[0]
cell_parquet_gz = cell_parquet_gz_files[0]
cell_parquet = cell_parquet_gz.with_suffix("")

print("Project root:", PROJECT_ROOT)
print("Xenium H5:", xenium_h5)
print("Cells parquet gz:", cell_parquet_gz)


# -------------------------
# Read Xenium expression matrix
# -------------------------

xenium = sc.read_10x_h5(xenium_h5)
xenium.var_names_make_unique()

print("\n===== Xenium expression matrix =====")
print(xenium)
print("Cells:", xenium.n_obs)
print("Features/genes:", xenium.n_vars)


# -------------------------
# Read cells.parquet.gz
# -------------------------

if not cell_parquet.exists():
    print("\nUnzipping cells.parquet.gz...")
    with gzip.open(cell_parquet_gz, "rb") as f_in:
        with open(cell_parquet, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

cells = pd.read_parquet(cell_parquet)

print("\n===== Xenium cell metadata =====")
print("Cells parquet shape:", cells.shape)
print("Columns:", cells.columns.tolist())


# -------------------------
# Check required columns
# -------------------------

required_columns = [
    "cell_id",
    "x_centroid",
    "y_centroid",
    "transcript_counts",
    "control_probe_counts",
    "control_codeword_counts",
    "unassigned_codeword_counts",
    "total_counts",
    "cell_area",
    "nucleus_area",
]

missing_columns = [
    col for col in required_columns
    if col not in cells.columns
]

if len(missing_columns) > 0:
    raise ValueError(f"Missing required columns in cells.parquet: {missing_columns}")


# -------------------------
# Match cell IDs between H5 and parquet
# -------------------------

cells["cell_id"] = cells["cell_id"].astype(str)
xenium.obs_names = xenium.obs_names.astype(str)

common_cell_ids = sorted(
    set(cells["cell_id"])
    &
    set(xenium.obs_names)
)

print("\n===== Cell ID matching =====")
print("Cells in H5:", xenium.n_obs)
print("Cells in parquet:", cells.shape[0])
print("Matched cells:", len(common_cell_ids))

if len(common_cell_ids) == 0:
    print("\nWARNING: No matching cell IDs found.")
    print("First 5 H5 cell IDs:")
    print(list(xenium.obs_names[:5]))
    print("First 5 parquet cell IDs:")
    print(cells["cell_id"].head().tolist())
    raise ValueError("Cell IDs do not match between H5 and cells.parquet.")


# Keep only matched cells
xenium = xenium[common_cell_ids, :].copy()
cells = (
    cells
    .set_index("cell_id")
    .loc[common_cell_ids]
    .copy()
)


# -------------------------
# Expression-based QC
# -------------------------

xenium.var["mt"] = xenium.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    xenium,
    qc_vars=["mt"],
    inplace=True,
    percent_top=None,
)

cells["n_genes_by_counts"] = xenium.obs["n_genes_by_counts"].values
cells["matrix_total_counts"] = xenium.obs["total_counts"].values
cells["pct_counts_mt"] = xenium.obs["pct_counts_mt"].values


# -------------------------
# Metadata-based QC percentages
# -------------------------

safe_total = cells["total_counts"].replace(0, np.nan)

cells["pct_control_probe"] = (
    100 * cells["control_probe_counts"] / safe_total
)

cells["pct_control_codeword"] = (
    100 * cells["control_codeword_counts"] / safe_total
)

cells["pct_unassigned"] = (
    100 * cells["unassigned_codeword_counts"] / safe_total
)


# -------------------------
# Outlier thresholds
# Lieber-inspired QC logic
# -------------------------

cells["empty_cell_out"] = cells["matrix_total_counts"] == 0

probe_thresh = cells["pct_control_probe"].quantile(0.99)
codeword_thresh = cells["pct_control_codeword"].quantile(0.99)
unassigned_thresh = cells["pct_unassigned"].quantile(0.99)

mito_thresh = cells["pct_counts_mt"].quantile(0.99)
low_detected_thresh = cells["n_genes_by_counts"].quantile(0.01)

low_total_counts_thresh = cells["matrix_total_counts"].quantile(0.01)
high_total_counts_thresh = cells["matrix_total_counts"].quantile(0.99)

cells["neg_probe_out"] = cells["pct_control_probe"] >= probe_thresh
cells["neg_codeword_out"] = cells["pct_control_codeword"] >= codeword_thresh
cells["unassigned_out"] = cells["pct_unassigned"] >= unassigned_thresh

cells["mito_out"] = cells["pct_counts_mt"] >= mito_thresh
cells["detected_out"] = cells["n_genes_by_counts"] <= low_detected_thresh

cells["total_counts_out"] = (
    (cells["matrix_total_counts"] <= low_total_counts_thresh)
    |
    (cells["matrix_total_counts"] >= high_total_counts_thresh)
)

# Lieber final union did NOT include mito_out.
# We keep both versions for transparency.

cells["is_outlier_lieber_like"] = (
    cells["empty_cell_out"]
    |
    cells["neg_probe_out"]
    |
    cells["neg_codeword_out"]
    |
    cells["unassigned_out"]
    |
    cells["detected_out"]
    |
    cells["total_counts_out"]
)

cells["is_outlier_with_mito"] = (
    cells["is_outlier_lieber_like"]
    |
    cells["mito_out"]
)


# -------------------------
# Save QC metrics
# -------------------------

qc_metrics_file = TABLE_DIR / "xenium_cell_qc_metrics.csv"
outlier_file = TABLE_DIR / "xenium_outlier_cells.csv"
summary_file = TABLE_DIR / "xenium_qc_summary.csv"

cells.to_csv(qc_metrics_file)

outlier_cells = (
    cells[cells["is_outlier_lieber_like"]]
    .reset_index()[["cell_id"]]
)

outlier_cells.to_csv(outlier_file, index=False)

summary = pd.DataFrame({
    "metric": [
        "total_cells_h5",
        "total_cells_parquet",
        "matched_cells",
        "total_features",
        "empty_cell_out",
        "neg_probe_threshold_99pct",
        "neg_probe_out",
        "neg_codeword_threshold_99pct",
        "neg_codeword_out",
        "unassigned_threshold_99pct",
        "unassigned_out",
        "mito_threshold_99pct",
        "mito_out",
        "low_detected_threshold_1pct",
        "detected_out",
        "low_total_counts_threshold_1pct",
        "high_total_counts_threshold_99pct",
        "total_counts_out",
        "final_outliers_lieber_like",
        "final_outliers_with_mito",
        "kept_cells_lieber_like",
    ],
    "value": [
        xenium.n_obs,
        cells.shape[0],
        len(common_cell_ids),
        xenium.n_vars,
        int(cells["empty_cell_out"].sum()),
        float(probe_thresh),
        int(cells["neg_probe_out"].sum()),
        float(codeword_thresh),
        int(cells["neg_codeword_out"].sum()),
        float(unassigned_thresh),
        int(cells["unassigned_out"].sum()),
        float(mito_thresh),
        int(cells["mito_out"].sum()),
        float(low_detected_thresh),
        int(cells["detected_out"].sum()),
        float(low_total_counts_thresh),
        float(high_total_counts_thresh),
        int(cells["total_counts_out"].sum()),
        int(cells["is_outlier_lieber_like"].sum()),
        int(cells["is_outlier_with_mito"].sum()),
        int((~cells["is_outlier_lieber_like"]).sum()),
    ]
})

summary.to_csv(summary_file, index=False)


# -------------------------
# Basic plots
# -------------------------

def save_hist(column, filename, threshold=None):
    plt.figure(figsize=(7, 5))
    plt.hist(cells[column].dropna(), bins=60)
    if threshold is not None:
        plt.axvline(threshold, linestyle="--")
    plt.xlabel(column)
    plt.ylabel("Number of cells")
    plt.title(column)
    plt.tight_layout()
    plt.savefig(FIGURE_DIR / filename, dpi=300)
    plt.close()


save_hist("pct_control_probe", "xenium_pct_control_probe_hist.png", probe_thresh)
save_hist("pct_control_codeword", "xenium_pct_control_codeword_hist.png", codeword_thresh)
save_hist("pct_unassigned", "xenium_pct_unassigned_hist.png", unassigned_thresh)
save_hist("pct_counts_mt", "xenium_pct_mito_hist.png", mito_thresh)
save_hist("n_genes_by_counts", "xenium_detected_genes_hist.png", low_detected_thresh)
save_hist("matrix_total_counts", "xenium_total_counts_hist.png", high_total_counts_thresh)


# Spatial plot: all cells
plt.figure(figsize=(7, 7))
plt.scatter(
    cells["x_centroid"],
    cells["y_centroid"],
    s=0.2,
)
plt.gca().invert_yaxis()
plt.xlabel("x_centroid")
plt.ylabel("y_centroid")
plt.title("Xenium cell locations")
plt.tight_layout()
plt.savefig(FIGURE_DIR / "xenium_cell_locations.png", dpi=300)
plt.close()


# Spatial plot: outliers
plt.figure(figsize=(7, 7))

normal = cells[~cells["is_outlier_lieber_like"]]
outliers = cells[cells["is_outlier_lieber_like"]]

plt.scatter(
    normal["x_centroid"],
    normal["y_centroid"],
    s=0.2,
    alpha=0.4,
)

plt.scatter(
    outliers["x_centroid"],
    outliers["y_centroid"],
    s=1,
    alpha=0.8,
)

plt.gca().invert_yaxis()
plt.xlabel("x_centroid")
plt.ylabel("y_centroid")
plt.title("Xenium QC outliers")
plt.tight_layout()
plt.savefig(FIGURE_DIR / "xenium_spatial_outliers.png", dpi=300)
plt.close()


# -------------------------
# Print final summary
# -------------------------

print("\n===== Xenium QC summary =====")
print(summary)

print("\nSaved tables:")
print(qc_metrics_file)
print(outlier_file)
print(summary_file)

print("\nSaved figures:")
print(FIGURE_DIR / "xenium_pct_control_probe_hist.png")
print(FIGURE_DIR / "xenium_pct_control_codeword_hist.png")
print(FIGURE_DIR / "xenium_pct_unassigned_hist.png")
print(FIGURE_DIR / "xenium_pct_mito_hist.png")
print(FIGURE_DIR / "xenium_detected_genes_hist.png")
print(FIGURE_DIR / "xenium_total_counts_hist.png")
print(FIGURE_DIR / "xenium_cell_locations.png")
print(FIGURE_DIR / "xenium_spatial_outliers.png")