# import scanpy as sc
# import pandas as pd
# from pathlib import Path

# PROJECT_ROOT = Path(__file__).resolve().parents[2]

# OUT_DIR = PROJECT_ROOT / "outputs" / "tables"
# OUT_DIR.mkdir(parents=True, exist_ok=True)

# # -------------------------
# # Xenium
# # -------------------------

# xenium_dir = (
#     PROJECT_ROOT
#     / "data"
#     / "raw"
#     / "xenium"
#     / "Br2039"
# )

# xenium_h5 = list(
#     xenium_dir.glob("*cell_feature_matrix.h5")
# )[0]

# xenium = sc.read_10x_h5(xenium_h5)

# # -------------------------
# # Visium
# # -------------------------

# visium_dir = (
#     PROJECT_ROOT
#     / "data"
#     / "raw"
#     / "visium"
#     / "Br2039"
# )

# matrix_file = list(
#     visium_dir.glob("*matrix.mtx.gz")
# )[0]

# features_file = list(
#     visium_dir.glob("*features.tsv.gz")
# )[0]

# visium = sc.read_mtx(matrix_file).T

# features = pd.read_csv(
#     features_file,
#     sep="\t",
#     header=None
# )

# visium.var_names = (
#     features.iloc[:, 1]
#     .astype(str)
#     .values
# )

# visium.var_names_make_unique()

# # -------------------------
# # Shared genes
# # -------------------------

# shared_genes = sorted(
#     set(visium.var_names)
#     &
#     set(xenium.var_names)
# )

# print("\n===== Shared Genes =====")
# print(f"Visium genes : {visium.n_vars:,}")
# print(f"Xenium genes : {xenium.n_vars:,}")
# print(f"Shared genes : {len(shared_genes):,}")

# shared_df = pd.DataFrame({
#     "gene": shared_genes
# })

# output_file = (
#     OUT_DIR
#     / "shared_genes.csv"
# )

# shared_df.to_csv(
#     output_file,
#     index=False
# )

# print(f"\nSaved to:")
# print(output_file)

from pathlib import Path
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parents[2]

VISIUM_DIR = PROJECT_ROOT / "data" / "raw" / "visium"
XENIUM_DIR = PROJECT_ROOT / "data" / "raw" / "xenium"
OUT_DIR = PROJECT_ROOT / "outputs" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("Project root:", PROJECT_ROOT)
print("Visium dir:", VISIUM_DIR)
print("Xenium dir:", XENIUM_DIR)
print("Visium exists:", VISIUM_DIR.exists())
print("Xenium exists:", XENIUM_DIR.exists())

features_files = list(VISIUM_DIR.rglob("*features.tsv.gz"))
xenium_h5_files = list(XENIUM_DIR.rglob("*cell_feature_matrix.h5"))

print("Features files found:", features_files)
print("Xenium h5 files found:", xenium_h5_files)

if len(features_files) == 0:
    raise FileNotFoundError("No Visium features.tsv.gz file found.")

if len(xenium_h5_files) == 0:
    raise FileNotFoundError("No Xenium cell_feature_matrix.h5 file found.")

features_file = features_files[0]
xenium_h5 = xenium_h5_files[0]

visium_features = pd.read_csv(features_file, sep="\t", header=None)
visium_genes = set(visium_features[1].astype(str))

xenium = sc.read_10x_h5(xenium_h5)
xenium_genes = set(xenium.var_names.astype(str))

shared = sorted(visium_genes & xenium_genes)


visium_only = sorted(visium_genes - xenium_genes)
xenium_only = sorted(xenium_genes - visium_genes)

print("\n===== Gene overlap summary =====")
print("Visium genes:", len(visium_genes))
print("Xenium genes:", len(xenium_genes))
print("Shared genes:", len(shared))
print("Visium-only genes:", len(visium_only))
print("Xenium-only genes:", len(xenium_only))

print("\nFirst 20 shared genes:")
print(shared[:20])

pd.Series(shared, name="gene").to_csv(OUT_DIR / "shared_genes.csv", index=False)
#pd.Series(visium_only, name="gene").to_csv(OUT_DIR / "visium_only_genes.csv", index=False)
#pd.Series(xenium_only, name="gene").to_csv(OUT_DIR / "xenium_only_genes.csv", index=False)

summary = pd.DataFrame({
    "category": ["visium_genes", "xenium_genes", "shared_genes", "visium_only_genes", "xenium_only_genes"],
    "count": [len(visium_genes), len(xenium_genes), len(shared), len(visium_only), len(xenium_only)]
})

pd.Series(shared, name="gene").to_csv(
    OUT_DIR / "shared_genes.csv",
    index=False
)

summary.to_csv(OUT_DIR / "gene_overlap_summary.csv", index=False)

print("\nSaved files:")
print(OUT_DIR / "shared_genes.csv")
#print(OUT_DIR / "visium_only_genes.csv")
#print(OUT_DIR / "xenium_only_genes.csv")
print(OUT_DIR / "gene_overlap_summary.csv")