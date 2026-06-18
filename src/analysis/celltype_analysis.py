# """
# Task 1: Patient-level spatial map for Br2039

# Goal:
#     Generate 1 Xenium cell-level map and 1 Visium spot-level map
#     with first-pass marker-based broad cell-type annotation.

# This version includes a robust reader for Xenium cells.parquet.gz files.
# Some GEO files are gzip-wrapped Parquet files, so normal pd.read_parquet()
# can fail with:
#     "Parquet magic bytes not found in footer"

# Expected input folders:

#     data/raw/xenium/Br2039/
#         *cell_feature_matrix.h5
#         *cells.parquet.gz

#     data/raw/visium/Br2039/
#         *matrix.mtx.gz
#         *features.tsv.gz
#         *barcodes.tsv.gz
#         *tissue_positions*.csv.gz

# Generated outputs:

#     outputs/task1_patient_map/Br2039/
#         figures/
#             Br2039_xenium_celltype_map.png
#             Br2039_visium_celltype_map.png
#             Br2039_combined_xenium_visium_map.png

#         tables/
#             Br2039_xenium_celltype_annotations.csv
#             Br2039_visium_spot_annotations.csv
#             Br2039_xenium_celltype_counts.csv
#             Br2039_visium_spot_celltype_counts.csv
#             Br2039_xenium_marker_availability.csv
#             Br2039_visium_marker_availability.csv
#             celltype_color_legend.csv

#         docs/
#             task1_meeting_notes.md
# """

# from pathlib import Path
# import gzip
# import shutil
# import tempfile
# import os

# import numpy as np
# import pandas as pd
# import scanpy as sc
# import scipy.io
# import scipy.sparse as sp
# import matplotlib.pyplot as plt
# from anndata import AnnData


# # =============================================================================
# # Project paths
# # =============================================================================

# PROJECT_ROOT = Path(__file__).resolve().parents[2]

# SAMPLE_ID = "Br2039"

# XENIUM_DIR = PROJECT_ROOT / "data" / "raw" / "xenium" / SAMPLE_ID
# VISIUM_DIR = PROJECT_ROOT / "data" / "raw" / "visium" / SAMPLE_ID

# OUT_BASE = PROJECT_ROOT / "outputs" / "task1_patient_map" / SAMPLE_ID
# FIG_DIR = OUT_BASE / "figures"
# TABLE_DIR = OUT_BASE / "tables"
# DOC_DIR = OUT_BASE / "docs"

# FIG_DIR.mkdir(parents=True, exist_ok=True)
# TABLE_DIR.mkdir(parents=True, exist_ok=True)
# DOC_DIR.mkdir(parents=True, exist_ok=True)


# # =============================================================================
# # PI-provided categorical color map
# # =============================================================================

# cat_color = [
#     "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
#     "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
#     "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
#     "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785"
# ]


# # =============================================================================
# # Broad cell-type markers for first-pass annotation
# # =============================================================================

# MARKERS = {
#     "Excitatory_neuron": ["SLC17A7", "SATB2", "RORB", "PCP4"],
#     "Inhibitory_neuron": ["GAD1", "GAD2", "PVALB", "SST", "VIP"],
#     "Astrocyte": ["AQP4", "GFAP", "SLC1A2"],
#     "Oligodendrocyte": ["MBP", "PLP1", "MOG"],
#     "OPC": ["PDGFRA", "CSPG4"],
#     "Microglia": ["CX3CR1", "P2RY12", "C3"],
#     "Endothelial": ["CLDN5", "FLT1", "VWF"],
#     "Mural": ["PDGFRB", "RGS5"],
# }

# CELLTYPE_ORDER = [
#     "Excitatory_neuron",
#     "Inhibitory_neuron",
#     "Astrocyte",
#     "Oligodendrocyte",
#     "OPC",
#     "Microglia",
#     "Endothelial",
#     "Mural",
#     "Unknown",
# ]

# CELLTYPE_COLORS = {
#     celltype: cat_color[i]
#     for i, celltype in enumerate(CELLTYPE_ORDER)
# }


# # =============================================================================
# # Utility functions
# # =============================================================================

# def print_section(title: str):
#     print("\n" + "=" * 80)
#     print(title)
#     print("=" * 80)


# def find_one_file(folder: Path, pattern: str) -> Path:
#     """
#     Find one matching file.
#     If multiple files match, use the first sorted file and print a warning.
#     """
#     files = sorted(folder.glob(pattern))

#     if len(files) == 0:
#         raise FileNotFoundError(f"No file matching '{pattern}' found in: {folder}")

#     if len(files) > 1:
#         print(f"WARNING: Multiple files found for pattern '{pattern}'. Using: {files[0]}")

#     return files[0]


# def sparse_mean_by_rows(X):
#     """
#     Calculate row means for dense or sparse matrices.
#     """
#     if sp.issparse(X):
#         return np.asarray(X.mean(axis=1)).ravel()

#     return np.mean(X, axis=1)


# def read_parquet_maybe_gzip(path: Path) -> pd.DataFrame:
#     """
#     Read a Parquet file that may be either:
#         1. normal Parquet file
#         2. gzip-wrapped Parquet file, often named .parquet.gz

#     Why this is needed:
#         Some GEO files are externally gzip-compressed. In that case,
#         pd.read_parquet("file.parquet.gz") may fail with:
#         "Parquet magic bytes not found in footer."

#     This function:
#         - First tries normal pd.read_parquet().
#         - If that fails, checks whether the file starts with gzip bytes.
#         - If gzip, decompresses to a temporary .parquet file.
#         - Reads the temporary parquet file.
#         - Deletes the temporary file.
#     """

#     path = Path(path)

#     print(f"Reading parquet-like file: {path}")

#     # First try the normal reader.
#     try:
#         return pd.read_parquet(path)
#     except Exception as normal_error:
#         print("Normal pd.read_parquet() failed.")
#         print(f"Reason: {normal_error}")

#     # Check first bytes.
#     with open(path, "rb") as f:
#         first_bytes = f.read(4)

#     # Normal parquet usually starts with b"PAR1".
#     if first_bytes == b"PAR1":
#         return pd.read_parquet(path)

#     # Gzip files start with 1f 8b.
#     if first_bytes[:2] == b"\x1f\x8b":
#         print("Detected gzip-compressed parquet. Decompressing temporarily...")

#         temp_path = None

#         try:
#             with gzip.open(path, "rb") as gz:
#                 with tempfile.NamedTemporaryFile(
#                     suffix=".parquet",
#                     delete=False
#                 ) as tmp:
#                     shutil.copyfileobj(gz, tmp)
#                     temp_path = Path(tmp.name)

#             print(f"Temporary decompressed parquet: {temp_path}")

#             df = pd.read_parquet(temp_path)
#             return df

#         finally:
#             if temp_path is not None and temp_path.exists():
#                 try:
#                     os.remove(temp_path)
#                     print("Temporary decompressed parquet removed.")
#                 except Exception as cleanup_error:
#                     print(f"WARNING: Could not remove temporary file: {temp_path}")
#                     print(cleanup_error)

#     # If it is neither parquet nor gzip, it may be a bad download.
#     raise ValueError(
#         f"Could not read file as Parquet or gzip-compressed Parquet: {path}\n"
#         f"First bytes were: {first_bytes}\n\n"
#         "Possible causes:\n"
#         "1. The file download is incomplete/corrupted.\n"
#         "2. The file is actually an HTML download page, not the data file.\n"
#         "3. The file has the wrong extension.\n\n"
#         "To inspect first bytes, run:\n"
#         f"python -c \"from pathlib import Path; p=Path(r'{str(path)}'); "
#         "print(open(p,'rb').read(50))\""
#     )


# # =============================================================================
# # Load Xenium
# # =============================================================================

# def load_xenium_br2039() -> AnnData:
#     """
#     Load Xenium cell_feature_matrix.h5 and cells.parquet.gz.
#     Adds cell centroid coordinates to adata.obs.
#     """

#     print_section("Loading Xenium data")

#     h5_path = find_one_file(XENIUM_DIR, "*cell_feature_matrix.h5")
#     cells_path = find_one_file(XENIUM_DIR, "*cells.parquet.gz")

#     print(f"Xenium matrix file: {h5_path}")
#     print(f"Xenium cells metadata file: {cells_path}")

#     adata = sc.read_10x_h5(h5_path)
#     adata.var_names_make_unique()

#     cells = read_parquet_maybe_gzip(cells_path)

#     print("Xenium cells.parquet columns:")
#     print(list(cells.columns))

#     if "cell_id" not in cells.columns:
#         raise ValueError("Xenium cells.parquet.gz must contain column: cell_id")

#     if "x_centroid" not in cells.columns:
#         raise ValueError("Xenium cells.parquet.gz must contain column: x_centroid")

#     if "y_centroid" not in cells.columns:
#         raise ValueError("Xenium cells.parquet.gz must contain column: y_centroid")

#     cells["cell_id"] = cells["cell_id"].astype(str)
#     cells = cells.set_index("cell_id")

#     adata.obs_names = adata.obs_names.astype(str)

#     common_cells = adata.obs_names.intersection(cells.index)

#     print(f"Xenium matrix cells: {adata.n_obs}")
#     print(f"Xenium metadata cells: {cells.shape[0]}")
#     print(f"Matched cells: {len(common_cells)}")
#     print(f"Xenium genes: {adata.n_vars}")

#     if len(common_cells) == 0:
#         print("\nFirst 5 matrix cell IDs:")
#         print(list(adata.obs_names[:5]))
#         print("\nFirst 5 metadata cell IDs:")
#         print(list(cells.index[:5]))

#         raise ValueError(
#             "No matching cell IDs between Xenium matrix and cells.parquet.gz. "
#             "Check whether cell_id format differs."
#         )

#     adata = adata[common_cells].copy()
#     adata.obs = adata.obs.join(cells.loc[common_cells], how="left")

#     adata.obs["sample_id"] = SAMPLE_ID
#     adata.obs["technology"] = "Xenium"

#     return adata


# # =============================================================================
# # Load Visium
# # =============================================================================

# def load_visium_br2039() -> AnnData:
#     """
#     Load GEO-style Visium matrix, features, barcodes, and tissue positions.
#     This does not require a full Space Ranger folder.
#     """

#     print_section("Loading Visium data")

#     matrix_path = find_one_file(VISIUM_DIR, "*matrix.mtx.gz")
#     features_path = find_one_file(VISIUM_DIR, "*features.tsv.gz")
#     barcodes_path = find_one_file(VISIUM_DIR, "*barcodes.tsv.gz")
#     positions_path = find_one_file(VISIUM_DIR, "*tissue_positions*.csv.gz")

#     print(f"Visium matrix file: {matrix_path}")
#     print(f"Visium features file: {features_path}")
#     print(f"Visium barcodes file: {barcodes_path}")
#     print(f"Visium positions file: {positions_path}")

#     X = scipy.io.mmread(matrix_path).tocsr()

#     features = pd.read_csv(
#         features_path,
#         sep="\t",
#         header=None,
#     )

#     if features.shape[1] >= 3:
#         features = features.iloc[:, :3]
#         features.columns = ["gene_id", "gene_name", "feature_type"]
#     elif features.shape[1] == 2:
#         features.columns = ["gene_id", "gene_name"]
#         features["feature_type"] = "Gene Expression"
#     else:
#         raise ValueError("Unexpected features.tsv.gz format")

#     barcodes = pd.read_csv(
#         barcodes_path,
#         sep="\t",
#         header=None,
#         names=["barcode"],
#     )

#     # 10x matrix is genes x spots.
#     # AnnData expects observations x variables = spots x genes.
#     adata = AnnData(X=X.T)

#     adata.obs_names = barcodes["barcode"].astype(str).values
#     adata.var_names = features["gene_name"].astype(str).values
#     adata.var["gene_id"] = features["gene_id"].astype(str).values
#     adata.var["feature_type"] = features["feature_type"].astype(str).values
#     adata.var_names_make_unique()

#     positions_raw = pd.read_csv(positions_path, header=None)

#     # Some files include a header row; some do not.
#     # Expected columns:
#     # barcode, in_tissue, array_row, array_col, pxl_col_in_fullres, pxl_row_in_fullres
#     first_cell = str(positions_raw.iloc[0, 0]).lower()

#     if first_cell in ["barcode", "barcodes"]:
#         positions = pd.read_csv(positions_path)
#     else:
#         positions = positions_raw

#     if positions.shape[1] < 6:
#         raise ValueError("Unexpected tissue_positions.csv.gz format: fewer than 6 columns")

#     positions = positions.iloc[:, :6]
#     positions.columns = [
#         "barcode",
#         "in_tissue",
#         "array_row",
#         "array_col",
#         "pxl_col_in_fullres",
#         "pxl_row_in_fullres",
#     ]

#     positions["barcode"] = positions["barcode"].astype(str)
#     positions = positions.set_index("barcode")

#     common_spots = adata.obs_names.intersection(positions.index)

#     print(f"Visium matrix spots: {adata.n_obs}")
#     print(f"Visium position spots: {positions.shape[0]}")
#     print(f"Matched spots: {len(common_spots)}")
#     print(f"Visium genes: {adata.n_vars}")

#     if len(common_spots) == 0:
#         print("\nFirst 5 matrix barcodes:")
#         print(list(adata.obs_names[:5]))
#         print("\nFirst 5 position barcodes:")
#         print(list(positions.index[:5]))

#         raise ValueError(
#             "No matching barcodes between Visium matrix and tissue_positions file."
#         )

#     adata = adata[common_spots].copy()
#     adata.obs = adata.obs.join(positions.loc[common_spots], how="left")

#     # Keep tissue-covered spots only.
#     if "in_tissue" in adata.obs.columns:
#         before = adata.n_obs
#         adata = adata[adata.obs["in_tissue"].astype(int) == 1].copy()
#         after = adata.n_obs
#         print(f"Kept in-tissue Visium spots: {after} / {before}")

#     adata.obs["sample_id"] = SAMPLE_ID
#     adata.obs["technology"] = "Visium"

#     return adata


# # =============================================================================
# # Marker-based annotation
# # =============================================================================

# def normalize_log1p(adata: AnnData) -> AnnData:
#     """
#     Normalize total counts and log-transform.
#     """
#     adata = adata.copy()

#     # Keep genes detected in at least one cell/spot.
#     sc.pp.filter_genes(adata, min_cells=1)

#     # Normalize and log-transform.
#     sc.pp.normalize_total(adata, target_sum=1e4)
#     sc.pp.log1p(adata)

#     return adata


# def marker_score_annotation(adata: AnnData, label: str) -> AnnData:
#     """
#     First-pass broad cell-type annotation using mean expression of marker genes.

#     Method:
#         - Normalize/log1p data.
#         - For each cell type, find available marker genes.
#         - Compute average expression of those markers per cell/spot.
#         - Assign cell type with highest marker score.
#         - If all marker scores are <= 0, assign Unknown.
#     """

#     print_section(f"Marker-based annotation: {label}")

#     adata = normalize_log1p(adata)

#     marker_summary = []

#     for celltype, marker_genes in MARKERS.items():
#         available_genes = [g for g in marker_genes if g in adata.var_names]
#         score_col = f"score_{celltype}"

#         if len(available_genes) == 0:
#             adata.obs[score_col] = 0.0
#             print(f"{celltype}: no markers found")
#         else:
#             X_sub = adata[:, available_genes].X
#             adata.obs[score_col] = sparse_mean_by_rows(X_sub)
#             print(f"{celltype}: using {len(available_genes)} markers: {available_genes}")

#         marker_summary.append(
#             {
#                 "technology": label,
#                 "celltype": celltype,
#                 "requested_markers": ",".join(marker_genes),
#                 "available_markers": ",".join(available_genes),
#                 "n_available_markers": len(available_genes),
#             }
#         )

#     score_cols = [f"score_{ct}" for ct in MARKERS.keys()]
#     score_df = adata.obs[score_cols].copy()

#     adata.obs["predicted_celltype"] = (
#         score_df.idxmax(axis=1).str.replace("score_", "", regex=False)
#     )
#     adata.obs["max_marker_score"] = score_df.max(axis=1)

#     adata.obs.loc[
#         adata.obs["max_marker_score"] <= 0,
#         "predicted_celltype"
#     ] = "Unknown"

#     marker_summary_df = pd.DataFrame(marker_summary)
#     marker_summary_out = TABLE_DIR / f"{SAMPLE_ID}_{label.lower()}_marker_availability.csv"
#     marker_summary_df.to_csv(marker_summary_out, index=False)

#     print(f"Saved marker availability: {marker_summary_out}")

#     print(f"\n{label} predicted cell-type counts:")
#     print(adata.obs["predicted_celltype"].value_counts())

#     return adata


# # =============================================================================
# # Plotting
# # =============================================================================

# def plot_xenium_map(adata: AnnData):
#     """
#     Plot Xenium cell-level spatial map.
#     Each dot is one segmented Xenium cell.
#     """

#     print_section("Plotting Xenium map")

#     plot_df = adata.obs.copy()

#     plt.figure(figsize=(8, 8))

#     for celltype in CELLTYPE_ORDER:
#         sub = plot_df[plot_df["predicted_celltype"] == celltype]

#         if sub.shape[0] == 0:
#             continue

#         plt.scatter(
#             sub["x_centroid"],
#             sub["y_centroid"],
#             s=1,
#             alpha=0.75,
#             color=CELLTYPE_COLORS.get(celltype, "#787878"),
#             label=f"{celltype} ({sub.shape[0]})",
#             linewidths=0,
#         )

#     plt.gca().invert_yaxis()
#     plt.xlabel("X centroid")
#     plt.ylabel("Y centroid")
#     plt.title(f"{SAMPLE_ID} Xenium cell-level map\nMarker-based broad cell type")

#     plt.legend(
#         bbox_to_anchor=(1.05, 1),
#         loc="upper left",
#         frameon=False,
#         markerscale=6,
#         fontsize=8,
#     )

#     plt.tight_layout()

#     out_path = FIG_DIR / f"{SAMPLE_ID}_xenium_celltype_map.png"
#     plt.savefig(out_path, dpi=300)
#     plt.close()

#     print(f"Saved: {out_path}")


# def plot_visium_map(adata: AnnData):
#     """
#     Plot Visium spot-level spatial map.
#     Each dot is one Visium spot.
#     """

#     print_section("Plotting Visium map")

#     plot_df = adata.obs.copy()

#     plt.figure(figsize=(8, 8))

#     for celltype in CELLTYPE_ORDER:
#         sub = plot_df[plot_df["predicted_celltype"] == celltype]

#         if sub.shape[0] == 0:
#             continue

#         plt.scatter(
#             sub["pxl_col_in_fullres"],
#             sub["pxl_row_in_fullres"],
#             s=22,
#             alpha=0.85,
#             color=CELLTYPE_COLORS.get(celltype, "#787878"),
#             label=f"{celltype} ({sub.shape[0]})",
#             linewidths=0,
#         )

#     plt.gca().invert_yaxis()
#     plt.xlabel("Pixel column in full-resolution image")
#     plt.ylabel("Pixel row in full-resolution image")
#     plt.title(f"{SAMPLE_ID} Visium spot-level map\nDominant marker-based broad cell type")

#     plt.legend(
#         bbox_to_anchor=(1.05, 1),
#         loc="upper left",
#         frameon=False,
#         markerscale=1.5,
#         fontsize=8,
#     )

#     plt.tight_layout()

#     out_path = FIG_DIR / f"{SAMPLE_ID}_visium_celltype_map.png"
#     plt.savefig(out_path, dpi=300)
#     plt.close()

#     print(f"Saved: {out_path}")


# def plot_combined_map(xenium: AnnData, visium: AnnData):
#     """
#     Plot side-by-side Xenium and Visium maps for meeting slide.
#     """

#     print_section("Plotting combined Xenium + Visium map")

#     xenium_df = xenium.obs.copy()
#     visium_df = visium.obs.copy()

#     fig, axes = plt.subplots(1, 2, figsize=(16, 7))

#     ax = axes[0]

#     for celltype in CELLTYPE_ORDER:
#         sub = xenium_df[xenium_df["predicted_celltype"] == celltype]

#         if sub.shape[0] == 0:
#             continue

#         ax.scatter(
#             sub["x_centroid"],
#             sub["y_centroid"],
#             s=1,
#             alpha=0.75,
#             color=CELLTYPE_COLORS.get(celltype, "#787878"),
#             label=celltype,
#             linewidths=0,
#         )

#     ax.invert_yaxis()
#     ax.set_xlabel("X centroid")
#     ax.set_ylabel("Y centroid")
#     ax.set_title("Xenium: cell-level")

#     ax = axes[1]

#     for celltype in CELLTYPE_ORDER:
#         sub = visium_df[visium_df["predicted_celltype"] == celltype]

#         if sub.shape[0] == 0:
#             continue

#         ax.scatter(
#             sub["pxl_col_in_fullres"],
#             sub["pxl_row_in_fullres"],
#             s=22,
#             alpha=0.85,
#             color=CELLTYPE_COLORS.get(celltype, "#787878"),
#             label=celltype,
#             linewidths=0,
#         )

#     ax.invert_yaxis()
#     ax.set_xlabel("Pixel column")
#     ax.set_ylabel("Pixel row")
#     ax.set_title("Visium: spot-level")

#     handles = []
#     labels = []

#     for celltype in CELLTYPE_ORDER:
#         if (
#             celltype in xenium_df["predicted_celltype"].values
#             or celltype in visium_df["predicted_celltype"].values
#         ):
#             handles.append(
#                 plt.Line2D(
#                     [0],
#                     [0],
#                     marker="o",
#                     color="w",
#                     markerfacecolor=CELLTYPE_COLORS.get(celltype, "#787878"),
#                     markersize=8,
#                     label=celltype,
#                 )
#             )
#             labels.append(celltype)

#     fig.legend(
#         handles,
#         labels,
#         loc="center right",
#         frameon=False,
#         fontsize=9,
#     )

#     fig.suptitle(
#         f"{SAMPLE_ID} patient-level spatial map: Xenium and Visium\n"
#         "First-pass marker-based broad cell-type annotation",
#         fontsize=14,
#     )

#     plt.tight_layout(rect=[0, 0, 0.88, 0.92])

#     out_path = FIG_DIR / f"{SAMPLE_ID}_combined_xenium_visium_map.png"
#     plt.savefig(out_path, dpi=300)
#     plt.close()

#     print(f"Saved: {out_path}")


# # =============================================================================
# # Save outputs
# # =============================================================================

# def save_tables(xenium: AnnData, visium: AnnData):
#     """
#     Save annotation tables, count tables, and color legend.
#     """

#     print_section("Saving tables")

#     xenium_cols = [
#         "sample_id",
#         "technology",
#         "x_centroid",
#         "y_centroid",
#         "predicted_celltype",
#         "max_marker_score",
#     ]

#     xenium_cols = [c for c in xenium_cols if c in xenium.obs.columns]

#     xenium_table = xenium.obs[xenium_cols].copy()
#     xenium_table.index.name = "cell_id"

#     xenium_out = TABLE_DIR / f"{SAMPLE_ID}_xenium_celltype_annotations.csv"
#     xenium_table.to_csv(xenium_out)

#     visium_cols = [
#         "sample_id",
#         "technology",
#         "in_tissue",
#         "array_row",
#         "array_col",
#         "pxl_col_in_fullres",
#         "pxl_row_in_fullres",
#         "predicted_celltype",
#         "max_marker_score",
#     ]

#     visium_cols = [c for c in visium_cols if c in visium.obs.columns]

#     visium_table = visium.obs[visium_cols].copy()
#     visium_table.index.name = "barcode"

#     visium_out = TABLE_DIR / f"{SAMPLE_ID}_visium_spot_annotations.csv"
#     visium_table.to_csv(visium_out)

#     xenium_counts = (
#         xenium.obs["predicted_celltype"]
#         .value_counts()
#         .rename_axis("predicted_celltype")
#         .reset_index(name="n_cells")
#     )

#     xenium_counts_out = TABLE_DIR / f"{SAMPLE_ID}_xenium_celltype_counts.csv"
#     xenium_counts.to_csv(xenium_counts_out, index=False)

#     visium_counts = (
#         visium.obs["predicted_celltype"]
#         .value_counts()
#         .rename_axis("predicted_celltype")
#         .reset_index(name="n_spots")
#     )

#     visium_counts_out = TABLE_DIR / f"{SAMPLE_ID}_visium_spot_celltype_counts.csv"
#     visium_counts.to_csv(visium_counts_out, index=False)

#     color_legend = pd.DataFrame(
#         {
#             "celltype": CELLTYPE_ORDER,
#             "color": [CELLTYPE_COLORS[ct] for ct in CELLTYPE_ORDER],
#         }
#     )

#     color_legend_out = TABLE_DIR / "celltype_color_legend.csv"
#     color_legend.to_csv(color_legend_out, index=False)

#     print(f"Saved: {xenium_out}")
#     print(f"Saved: {visium_out}")
#     print(f"Saved: {xenium_counts_out}")
#     print(f"Saved: {visium_counts_out}")
#     print(f"Saved: {color_legend_out}")


# def save_meeting_notes():
#     """
#     Save a short meeting note file you can copy into slides or use while presenting.
#     """

#     print_section("Saving meeting notes")

#     notes = f"""# Task 1 Meeting Notes: {SAMPLE_ID} Patient Map

# ## Goal

# Generate one patient-level spatial map using the same donor:

# - Xenium: cell-level spatial map
# - Visium: spot-level spatial map
# - Annotation: first-pass broad cell-type assignment using canonical marker genes
# - Color map: PI-provided categorical color palette

# ## Key statement for advisor

# This is the same donor/patient, {SAMPLE_ID}, but Xenium and Visium are not necessarily the exact same physical tissue section.

# Xenium provides cell-level resolution. Each point in the Xenium map represents one segmented cell.

# Visium provides spot-level resolution. Each point in the Visium map represents one Visium spot, which may contain multiple cells.

# ## Annotation method

# This version uses first-pass marker-based scoring.

# For each cell or spot:
# 1. Normalize total counts.
# 2. Log-transform expression.
# 3. Compute mean expression score for marker genes of each broad cell type.
# 4. Assign the cell type with the highest marker score.
# 5. If all marker scores are zero, assign Unknown.

# ## Broad cell types

# - Excitatory neuron
# - Inhibitory neuron
# - Astrocyte
# - Oligodendrocyte
# - OPC
# - Microglia
# - Endothelial
# - Mural
# - Unknown

# ## Important limitation

# This is not final reference-based annotation.

# The next step is to use a DLPFC single-cell or single-nucleus reference, such as Huuki-Myers or PsychENCODE, for more rigorous label transfer or spot deconvolution.

# ## Generated outputs

# Figures:

# - outputs/task1_patient_map/{SAMPLE_ID}/figures/{SAMPLE_ID}_xenium_celltype_map.png
# - outputs/task1_patient_map/{SAMPLE_ID}/figures/{SAMPLE_ID}_visium_celltype_map.png
# - outputs/task1_patient_map/{SAMPLE_ID}/figures/{SAMPLE_ID}_combined_xenium_visium_map.png

# Tables:

# - outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_xenium_celltype_annotations.csv
# - outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_visium_spot_annotations.csv
# - outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_xenium_celltype_counts.csv
# - outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_visium_spot_celltype_counts.csv
# - outputs/task1_patient_map/{SAMPLE_ID}/tables/celltype_color_legend.csv
# """

#     notes_out = DOC_DIR / "task1_meeting_notes.md"

#     with open(notes_out, "w", encoding="utf-8") as f:
#         f.write(notes)

#     print(f"Saved: {notes_out}")


# # =============================================================================
# # Main
# # =============================================================================

# def main():
#     print_section(f"Task 1 patient map for {SAMPLE_ID}")

#     print(f"Project root: {PROJECT_ROOT}")
#     print(f"Xenium dir: {XENIUM_DIR}")
#     print(f"Visium dir: {VISIUM_DIR}")
#     print(f"Output dir: {OUT_BASE}")

#     if not XENIUM_DIR.exists():
#         raise FileNotFoundError(f"Xenium folder does not exist: {XENIUM_DIR}")

#     if not VISIUM_DIR.exists():
#         raise FileNotFoundError(f"Visium folder does not exist: {VISIUM_DIR}")

#     # Load data
#     xenium = load_xenium_br2039()
#     visium = load_visium_br2039()

#     # Annotate
#     xenium = marker_score_annotation(xenium, label="Xenium")
#     visium = marker_score_annotation(visium, label="Visium")

#     # Plot maps
#     plot_xenium_map(xenium)
#     plot_visium_map(visium)
#     plot_combined_map(xenium, visium)

#     # Save tables and notes
#     save_tables(xenium, visium)
#     save_meeting_notes()

#     print_section("Done")

#     print("Generated main figures:")
#     print(FIG_DIR / f"{SAMPLE_ID}_xenium_celltype_map.png")
#     print(FIG_DIR / f"{SAMPLE_ID}_visium_celltype_map.png")
#     print(FIG_DIR / f"{SAMPLE_ID}_combined_xenium_visium_map.png")

#     print("\nUse this sentence in the meeting:")
#     print(
#         f"This is a first-pass patient-level spatial map for {SAMPLE_ID}. "
#         "Xenium is shown at cell-level resolution and Visium at spot-level resolution. "
#         "Both are colored using the PI-provided categorical color palette and broad "
#         "marker-based cell-type annotation."
#     )


# if __name__ == "__main__":
#     main()
"""
Task 1 corrected: Patient-level spatial map for Br2039

Advisor correction:
    - Xenium: cell-level map with broad cell-type annotation.
    - Visium: spot-level map with cortical layer-by-layer annotation.

This script generates:
    1. Br2039 Xenium cell-type map
    2. Br2039 Visium cortical layer map
    3. Combined Xenium + Visium meeting figure
    4. Annotation/count tables
    5. Meeting notes

Important:
    Visium layer annotation here is a first-pass marker-based approach.
    For publication-quality layer annotation, use manual annotation / reference-based
    spatial registration from spatialLIBD / LIBD DLPFC references.
"""

from pathlib import Path
import gzip
import shutil
import tempfile
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import matplotlib.pyplot as plt
from anndata import AnnData


# =============================================================================
# Project paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SAMPLE_ID = "Br2039"

XENIUM_DIR = PROJECT_ROOT / "data" / "raw" / "xenium" / SAMPLE_ID
VISIUM_DIR = PROJECT_ROOT / "data" / "raw" / "visium" / SAMPLE_ID

OUT_BASE = PROJECT_ROOT / "outputs" / "task1_patient_map" / SAMPLE_ID
FIG_DIR = OUT_BASE / "figures"
TABLE_DIR = OUT_BASE / "tables"
DOC_DIR = OUT_BASE / "docs"

FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
DOC_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# PI-provided categorical color map
# =============================================================================

cat_color = [
    "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
    "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
    "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
    "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785"
]


# =============================================================================
# Xenium broad cell-type markers
# =============================================================================

CELLTYPE_MARKERS = {
    "Excitatory_neuron": ["SLC17A7", "SATB2", "RORB", "PCP4"],
    "Inhibitory_neuron": ["GAD1", "GAD2", "PVALB", "SST", "VIP"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG"],
    "OPC": ["PDGFRA", "CSPG4"],
    "Microglia": ["CX3CR1", "P2RY12", "C3"],
    "Endothelial": ["CLDN5", "FLT1", "VWF"],
    "Mural": ["PDGFRB", "RGS5"],
}

CELLTYPE_ORDER = [
    "Excitatory_neuron",
    "Inhibitory_neuron",
    "Astrocyte",
    "Oligodendrocyte",
    "OPC",
    "Microglia",
    "Endothelial",
    "Mural",
    "Unknown",
]

CELLTYPE_COLORS = {
    celltype: cat_color[i]
    for i, celltype in enumerate(CELLTYPE_ORDER)
}


# =============================================================================
# Visium cortical layer markers
# =============================================================================
# First-pass marker set for human DLPFC cortical layer annotation.
# This is intentionally broad and conservative.
#
# Interpretation:
#     L1    = superficial / layer 1-associated markers
#     L2_3  = upper-layer excitatory neuron markers
#     L4    = granular / RORB-enriched layer
#     L5    = deep-layer excitatory neuron markers
#     L6    = deeper corticothalamic/corticocortical markers
#     WM    = white matter / oligodendrocyte-rich region

LAYER_MARKERS = {
    "L1": ["RELN", "LAMP5", "CXCL14"],
    "L2_3": ["CUX2", "CALB1", "LINC00507", "SATB2"],
    "L4": ["RORB", "SCNN1A"],
    "L5": ["BCL11B", "FEZF2", "ETV1"],
    "L6": ["FOXP2", "TLE4", "THEMIS", "PCP4"],
    "WM": ["MBP", "PLP1", "MOG", "MOBP"],
}

LAYER_ORDER = ["L1", "L2_3", "L4", "L5", "L6", "WM", "Unknown"]

LAYER_COLORS = {
    "L1": cat_color[0],
    "L2_3": cat_color[1],
    "L4": cat_color[2],
    "L5": cat_color[3],
    "L6": cat_color[4],
    "WM": cat_color[5],
    "Unknown": cat_color[10],
}


# =============================================================================
# Utility functions
# =============================================================================

def print_section(title: str):
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def find_one_file(folder: Path, pattern: str) -> Path:
    files = sorted(folder.glob(pattern))
    if len(files) == 0:
        raise FileNotFoundError(f"No file matching '{pattern}' found in: {folder}")
    if len(files) > 1:
        print(f"WARNING: Multiple files found for pattern '{pattern}'. Using: {files[0]}")
    return files[0]


def sparse_mean_by_rows(X):
    if sp.issparse(X):
        return np.asarray(X.mean(axis=1)).ravel()
    return np.mean(X, axis=1)


def read_parquet_maybe_gzip(path: Path) -> pd.DataFrame:
    """
    Read a Parquet file that may be normal Parquet or gzip-wrapped Parquet.
    This fixes the common GEO error: Parquet magic bytes not found in footer.
    """
    path = Path(path)
    print(f"Reading parquet-like file: {path}")

    try:
        return pd.read_parquet(path)
    except Exception as normal_error:
        print("Normal pd.read_parquet() failed.")
        print(f"Reason: {normal_error}")

    with open(path, "rb") as f:
        first_bytes = f.read(4)

    if first_bytes == b"PAR1":
        return pd.read_parquet(path)

    if first_bytes[:2] == b"\x1f\x8b":
        print("Detected gzip-compressed parquet. Decompressing temporarily...")
        temp_path = None
        try:
            with gzip.open(path, "rb") as gz:
                with tempfile.NamedTemporaryFile(suffix=".parquet", delete=False) as tmp:
                    shutil.copyfileobj(gz, tmp)
                    temp_path = Path(tmp.name)
            print(f"Temporary decompressed parquet: {temp_path}")
            df = pd.read_parquet(temp_path)
            return df
        finally:
            if temp_path is not None and temp_path.exists():
                try:
                    os.remove(temp_path)
                    print("Temporary decompressed parquet removed.")
                except Exception as cleanup_error:
                    print(f"WARNING: Could not remove temporary file: {temp_path}")
                    print(cleanup_error)

    raise ValueError(
        f"Could not read file as Parquet or gzip-compressed Parquet: {path}\n"
        f"First bytes were: {first_bytes}\n\n"
        "Possible causes: incomplete/corrupted download, HTML file, or wrong extension."
    )


# =============================================================================
# Load Xenium
# =============================================================================

def load_xenium_br2039() -> AnnData:
    print_section("Loading Xenium data")

    h5_path = find_one_file(XENIUM_DIR, "*cell_feature_matrix.h5")
    cells_path = find_one_file(XENIUM_DIR, "*cells.parquet.gz")

    print(f"Xenium matrix file: {h5_path}")
    print(f"Xenium cells metadata file: {cells_path}")

    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    cells = read_parquet_maybe_gzip(cells_path)

    print("Xenium cells.parquet columns:")
    print(list(cells.columns))

    for col in ["cell_id", "x_centroid", "y_centroid"]:
        if col not in cells.columns:
            raise ValueError(f"Xenium cells.parquet.gz must contain column: {col}")

    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id")
    adata.obs_names = adata.obs_names.astype(str)

    common_cells = adata.obs_names.intersection(cells.index)

    print(f"Xenium matrix cells: {adata.n_obs}")
    print(f"Xenium metadata cells: {cells.shape[0]}")
    print(f"Matched cells: {len(common_cells)}")
    print(f"Xenium genes: {adata.n_vars}")

    if len(common_cells) == 0:
        print("\nFirst 5 matrix cell IDs:")
        print(list(adata.obs_names[:5]))
        print("\nFirst 5 metadata cell IDs:")
        print(list(cells.index[:5]))
        raise ValueError("No matching cell IDs between Xenium matrix and cells.parquet.gz.")

    adata = adata[common_cells].copy()
    adata.obs = adata.obs.join(cells.loc[common_cells], how="left")
    adata.obs["sample_id"] = SAMPLE_ID
    adata.obs["technology"] = "Xenium"
    return adata


# =============================================================================
# Load Visium
# =============================================================================

def load_visium_br2039() -> AnnData:
    print_section("Loading Visium data")

    matrix_path = find_one_file(VISIUM_DIR, "*matrix.mtx.gz")
    features_path = find_one_file(VISIUM_DIR, "*features.tsv.gz")
    barcodes_path = find_one_file(VISIUM_DIR, "*barcodes.tsv.gz")
    positions_path = find_one_file(VISIUM_DIR, "*tissue_positions*.csv.gz")

    print(f"Visium matrix file: {matrix_path}")
    print(f"Visium features file: {features_path}")
    print(f"Visium barcodes file: {barcodes_path}")
    print(f"Visium positions file: {positions_path}")

    X = scipy.io.mmread(matrix_path).tocsr()

    features = pd.read_csv(features_path, sep="\t", header=None)
    if features.shape[1] >= 3:
        features = features.iloc[:, :3]
        features.columns = ["gene_id", "gene_name", "feature_type"]
    elif features.shape[1] == 2:
        features.columns = ["gene_id", "gene_name"]
        features["feature_type"] = "Gene Expression"
    else:
        raise ValueError("Unexpected features.tsv.gz format")

    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, names=["barcode"])

    adata = AnnData(X=X.T)
    adata.obs_names = barcodes["barcode"].astype(str).values
    adata.var_names = features["gene_name"].astype(str).values
    adata.var["gene_id"] = features["gene_id"].astype(str).values
    adata.var["feature_type"] = features["feature_type"].astype(str).values
    adata.var_names_make_unique()

    positions_raw = pd.read_csv(positions_path, header=None)
    first_cell = str(positions_raw.iloc[0, 0]).lower()

    if first_cell in ["barcode", "barcodes"]:
        positions = pd.read_csv(positions_path)
    else:
        positions = positions_raw

    if positions.shape[1] < 6:
        raise ValueError("Unexpected tissue_positions.csv.gz format: fewer than 6 columns")

    positions = positions.iloc[:, :6]
    positions.columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]

    positions["barcode"] = positions["barcode"].astype(str)
    positions = positions.set_index("barcode")
    common_spots = adata.obs_names.intersection(positions.index)

    print(f"Visium matrix spots: {adata.n_obs}")
    print(f"Visium position spots: {positions.shape[0]}")
    print(f"Matched spots: {len(common_spots)}")
    print(f"Visium genes: {adata.n_vars}")

    if len(common_spots) == 0:
        print("\nFirst 5 matrix barcodes:")
        print(list(adata.obs_names[:5]))
        print("\nFirst 5 position barcodes:")
        print(list(positions.index[:5]))
        raise ValueError("No matching barcodes between Visium matrix and tissue_positions file.")

    adata = adata[common_spots].copy()
    adata.obs = adata.obs.join(positions.loc[common_spots], how="left")

    if "in_tissue" in adata.obs.columns:
        before = adata.n_obs
        adata = adata[adata.obs["in_tissue"].astype(int) == 1].copy()
        after = adata.n_obs
        print(f"Kept in-tissue Visium spots: {after} / {before}")

    adata.obs["sample_id"] = SAMPLE_ID
    adata.obs["technology"] = "Visium"
    return adata


# =============================================================================
# Normalization and marker scoring
# =============================================================================

def normalize_log1p(adata: AnnData) -> AnnData:
    adata = adata.copy()
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def marker_score_annotation(
    adata: AnnData,
    markers: dict,
    ordered_labels: list,
    output_column: str,
    score_prefix: str,
    label: str,
    marker_out_name: str,
) -> AnnData:
    print_section(f"Marker-based annotation: {label}")

    adata = normalize_log1p(adata)
    marker_summary = []

    for class_name, marker_genes in markers.items():
        available_genes = [g for g in marker_genes if g in adata.var_names]
        score_col = f"{score_prefix}_{class_name}"

        if len(available_genes) == 0:
            adata.obs[score_col] = 0.0
            print(f"{class_name}: no markers found")
        else:
            X_sub = adata[:, available_genes].X
            adata.obs[score_col] = sparse_mean_by_rows(X_sub)
            print(f"{class_name}: using {len(available_genes)} markers: {available_genes}")

        marker_summary.append(
            {
                "technology": label,
                "class": class_name,
                "requested_markers": ",".join(marker_genes),
                "available_markers": ",".join(available_genes),
                "n_available_markers": len(available_genes),
            }
        )

    score_cols = [f"{score_prefix}_{ct}" for ct in markers.keys()]
    score_df = adata.obs[score_cols].copy()

    adata.obs[output_column] = score_df.idxmax(axis=1).str.replace(f"{score_prefix}_", "", regex=False)
    adata.obs[f"max_{score_prefix}_score"] = score_df.max(axis=1)
    adata.obs.loc[adata.obs[f"max_{score_prefix}_score"] <= 0, output_column] = "Unknown"

    final_order = [x for x in ordered_labels if x in list(adata.obs[output_column].unique())]
    adata.obs[output_column] = pd.Categorical(
        adata.obs[output_column],
        categories=final_order,
        ordered=True,
    )

    marker_summary_df = pd.DataFrame(marker_summary)
    marker_summary_out = TABLE_DIR / marker_out_name
    marker_summary_df.to_csv(marker_summary_out, index=False)
    print(f"Saved marker availability: {marker_summary_out}")

    print(f"\n{label} annotation counts:")
    print(adata.obs[output_column].value_counts(dropna=False))

    return adata


# =============================================================================
# Plotting
# =============================================================================

def plot_xenium_celltype_map(adata: AnnData):
    print_section("Plotting Xenium cell-type map")
    plot_df = adata.obs.copy()

    plt.figure(figsize=(8, 8))

    for celltype in CELLTYPE_ORDER:
        sub = plot_df[plot_df["predicted_celltype"] == celltype]
        if sub.shape[0] == 0:
            continue
        plt.scatter(
            sub["x_centroid"],
            sub["y_centroid"],
            s=5,
            alpha=0.75,
            color=CELLTYPE_COLORS.get(celltype, "#787878"),
            label=f"{celltype} ({sub.shape[0]})",
            linewidths=0,
        )

    plt.gca().invert_yaxis()
    plt.xlabel("X centroid")
    plt.ylabel("Y centroid")
    plt.title(f"{SAMPLE_ID} Xenium cell-level map\nBroad cell-type annotation")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, markerscale=6, fontsize=8)
    plt.tight_layout()

    out_path = FIG_DIR / f"{SAMPLE_ID}_xenium_celltype_map.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def plot_visium_layer_map(adata: AnnData):
    print_section("Plotting Visium cortical layer map")
    plot_df = adata.obs.copy()

    plt.figure(figsize=(8, 8))

    for layer in LAYER_ORDER:
        sub = plot_df[plot_df["predicted_layer"] == layer]
        if sub.shape[0] == 0:
            continue
        plt.scatter(
            sub["pxl_col_in_fullres"],
            sub["pxl_row_in_fullres"],
            s=60,
            alpha=0.85,
            color=LAYER_COLORS.get(layer, "#787878"),
            label=f"{layer} ({sub.shape[0]})",
            linewidths=0,
        )

    plt.gca().invert_yaxis()
    plt.xlabel("Pixel column in full-resolution image")
    plt.ylabel("Pixel row in full-resolution image")
    plt.title(f"{SAMPLE_ID} Visium spot-level map\nCortical layer annotation")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, markerscale=1.5, fontsize=8)
    plt.tight_layout()

    out_path = FIG_DIR / f"{SAMPLE_ID}_visium_layer_map.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def plot_combined_map(xenium: AnnData, visium: AnnData):
    print_section("Plotting combined Xenium cell-type + Visium layer map")
    xenium_df = xenium.obs.copy()
    visium_df = visium.obs.copy()

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    ax = axes[0]
    for celltype in CELLTYPE_ORDER:
        sub = xenium_df[xenium_df["predicted_celltype"] == celltype]
        if sub.shape[0] == 0:
            continue
        ax.scatter(
            sub["x_centroid"],
            sub["y_centroid"],
            s=5,
            alpha=0.75,
            color=CELLTYPE_COLORS.get(celltype, "#787878"),
            label=celltype,
            linewidths=0,
        )
    ax.invert_yaxis()
    ax.set_xlabel("X centroid")
    ax.set_ylabel("Y centroid")
    ax.set_title("Xenium: cell-level cell types")

    ax = axes[1]
    for layer in LAYER_ORDER:
        sub = visium_df[visium_df["predicted_layer"] == layer]
        if sub.shape[0] == 0:
            continue
        ax.scatter(
            sub["pxl_col_in_fullres"],
            sub["pxl_row_in_fullres"],
            s=22,
            alpha=0.85,
            color=LAYER_COLORS.get(layer, "#787878"),
            label=layer,
            linewidths=0,
        )
    ax.invert_yaxis()
    ax.set_xlabel("Pixel column")
    ax.set_ylabel("Pixel row")
    ax.set_title("Visium: spot-level cortical layers")

    celltype_handles = []
    celltype_labels = []
    for celltype in CELLTYPE_ORDER:
        if celltype in xenium_df["predicted_celltype"].astype(str).values:
            celltype_handles.append(
                plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=CELLTYPE_COLORS.get(celltype, "#787878"), markersize=7)
            )
            celltype_labels.append(celltype)

    layer_handles = []
    layer_labels = []
    for layer in LAYER_ORDER:
        if layer in visium_df["predicted_layer"].astype(str).values:
            layer_handles.append(
                plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=LAYER_COLORS.get(layer, "#787878"), markersize=8)
            )
            layer_labels.append(layer)

    legend1 = fig.legend(
        celltype_handles,
        celltype_labels,
        title="Xenium cell types",
        loc="center right",
        bbox_to_anchor=(1.02, 0.62),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    fig.add_artist(legend1)

    fig.legend(
        layer_handles,
        layer_labels,
        title="Visium layers",
        loc="center right",
        bbox_to_anchor=(1.02, 0.28),
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )

    fig.suptitle(
        f"{SAMPLE_ID} patient-level spatial map\n"
        "Xenium cell-type annotation and Visium cortical layer annotation",
        fontsize=14,
    )

    plt.tight_layout(rect=[0, 0, 0.86, 0.92])

    out_path = FIG_DIR / f"{SAMPLE_ID}_combined_xenium_celltype_visium_layer_map.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


# =============================================================================
# Save outputs
# =============================================================================

def save_tables(xenium: AnnData, visium: AnnData):
    print_section("Saving tables")

    xenium_cols = [
        "sample_id", "technology", "x_centroid", "y_centroid",
        "predicted_celltype", "max_celltype_score",
    ]
    xenium_cols = [c for c in xenium_cols if c in xenium.obs.columns]
    xenium_table = xenium.obs[xenium_cols].copy()
    xenium_table.index.name = "cell_id"
    xenium_out = TABLE_DIR / f"{SAMPLE_ID}_xenium_celltype_annotations.csv"
    xenium_table.to_csv(xenium_out)

    visium_cols = [
        "sample_id", "technology", "in_tissue", "array_row", "array_col",
        "pxl_col_in_fullres", "pxl_row_in_fullres", "predicted_layer", "max_layer_score",
    ]
    visium_cols = [c for c in visium_cols if c in visium.obs.columns]
    visium_table = visium.obs[visium_cols].copy()
    visium_table.index.name = "barcode"
    visium_out = TABLE_DIR / f"{SAMPLE_ID}_visium_layer_annotations.csv"
    visium_table.to_csv(visium_out)

    xenium_counts = (
        xenium.obs["predicted_celltype"]
        .value_counts(dropna=False)
        .rename_axis("predicted_celltype")
        .reset_index(name="n_cells")
    )
    xenium_counts_out = TABLE_DIR / f"{SAMPLE_ID}_xenium_celltype_counts.csv"
    xenium_counts.to_csv(xenium_counts_out, index=False)

    visium_counts = (
        visium.obs["predicted_layer"]
        .value_counts(dropna=False)
        .rename_axis("predicted_layer")
        .reset_index(name="n_spots")
    )
    visium_counts_out = TABLE_DIR / f"{SAMPLE_ID}_visium_layer_counts.csv"
    visium_counts.to_csv(visium_counts_out, index=False)

    celltype_legend = pd.DataFrame({"celltype": CELLTYPE_ORDER, "color": [CELLTYPE_COLORS[ct] for ct in CELLTYPE_ORDER]})
    celltype_legend_out = TABLE_DIR / "xenium_celltype_color_legend.csv"
    celltype_legend.to_csv(celltype_legend_out, index=False)

    layer_legend = pd.DataFrame({"layer": LAYER_ORDER, "color": [LAYER_COLORS[ly] for ly in LAYER_ORDER]})
    layer_legend_out = TABLE_DIR / "visium_layer_color_legend.csv"
    layer_legend.to_csv(layer_legend_out, index=False)

    print(f"Saved: {xenium_out}")
    print(f"Saved: {visium_out}")
    print(f"Saved: {xenium_counts_out}")
    print(f"Saved: {visium_counts_out}")
    print(f"Saved: {celltype_legend_out}")
    print(f"Saved: {layer_legend_out}")


def save_meeting_notes():
    print_section("Saving meeting notes")

    notes = f"""# Task 1 Meeting Notes: {SAMPLE_ID} Patient Map

## Advisor-requested output

The corrected Task 1 output contains:

1. Xenium cell-level spatial map with broad cell-type annotation.
2. Visium spot-level spatial map with cortical layer-by-layer annotation.

## Key statement

This is the same donor/patient, {SAMPLE_ID}, but Xenium and Visium are not necessarily the exact same physical tissue section.

Xenium is cell-level. Each point in the Xenium map is one segmented cell.

Visium is spot-level. Each point in the Visium map is one spot, annotated as a predicted cortical layer.

## Xenium annotation

Xenium is annotated using broad marker-based cell-type scoring.

## Visium layer annotation

Visium is annotated layer-by-layer using first-pass marker-based cortical layer scoring:

- L1
- L2/3
- L4
- L5
- L6
- WM
- Unknown

## Important limitation

This is a first-pass marker-based layer annotation.

For a stronger version, the next step should be to use reference-based cortical layer annotation from LIBD/spatialLIBD, where DLPFC Visium spots have manual layer labels for the six cortical layers plus white matter.

## Generated figures

- outputs/task1_patient_map/{SAMPLE_ID}/figures/{SAMPLE_ID}_xenium_celltype_map.png
- outputs/task1_patient_map/{SAMPLE_ID}/figures/{SAMPLE_ID}_visium_layer_map.png
- outputs/task1_patient_map/{SAMPLE_ID}/figures/{SAMPLE_ID}_combined_xenium_celltype_visium_layer_map.png

## Generated tables

- outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_xenium_celltype_annotations.csv
- outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_visium_layer_annotations.csv
- outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_xenium_celltype_counts.csv
- outputs/task1_patient_map/{SAMPLE_ID}/tables/{SAMPLE_ID}_visium_layer_counts.csv
"""

    notes_out = DOC_DIR / "task1_meeting_notes.md"
    with open(notes_out, "w", encoding="utf-8") as f:
        f.write(notes)
    print(f"Saved: {notes_out}")


# =============================================================================
# Main
# =============================================================================

def main():
    print_section(f"Corrected Task 1 patient map for {SAMPLE_ID}")

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Xenium dir: {XENIUM_DIR}")
    print(f"Visium dir: {VISIUM_DIR}")
    print(f"Output dir: {OUT_BASE}")

    if not XENIUM_DIR.exists():
        raise FileNotFoundError(f"Xenium folder does not exist: {XENIUM_DIR}")
    if not VISIUM_DIR.exists():
        raise FileNotFoundError(f"Visium folder does not exist: {VISIUM_DIR}")

    xenium = load_xenium_br2039()
    visium = load_visium_br2039()

    xenium = marker_score_annotation(
        adata=xenium,
        markers=CELLTYPE_MARKERS,
        ordered_labels=CELLTYPE_ORDER,
        output_column="predicted_celltype",
        score_prefix="celltype",
        label="Xenium cell type",
        marker_out_name=f"{SAMPLE_ID}_xenium_marker_availability.csv",
    )

    visium = marker_score_annotation(
        adata=visium,
        markers=LAYER_MARKERS,
        ordered_labels=LAYER_ORDER,
        output_column="predicted_layer",
        score_prefix="layer",
        label="Visium cortical layer",
        marker_out_name=f"{SAMPLE_ID}_visium_layer_marker_availability.csv",
    )

    plot_xenium_celltype_map(xenium)
    plot_visium_layer_map(visium)
    plot_combined_map(xenium, visium)

    save_tables(xenium, visium)
    save_meeting_notes()

    print_section("Done")
    print("Generated main figures:")
    print(FIG_DIR / f"{SAMPLE_ID}_xenium_celltype_map.png")
    print(FIG_DIR / f"{SAMPLE_ID}_visium_layer_map.png")
    print(FIG_DIR / f"{SAMPLE_ID}_combined_xenium_celltype_visium_layer_map.png")

    print("\nUse this sentence in the meeting:")
    print(
        f"This is a first-pass patient-level spatial map for {SAMPLE_ID}. "
        "Xenium is shown with broad cell-type annotation at cell-level resolution. "
        "Visium is shown with cortical layer-by-layer annotation at spot-level resolution. "
        "Both use the PI-provided categorical color palette."
    )


if __name__ == "__main__":
    main()
