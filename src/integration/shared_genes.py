import scanpy as sc
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUT_DIR = PROJECT_ROOT / "outputs" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------
# Xenium
# -------------------------

xenium_dir = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "xenium"
    / "Br2039"
)

xenium_h5 = list(
    xenium_dir.glob("*cell_feature_matrix.h5")
)[0]

xenium = sc.read_10x_h5(xenium_h5)

# -------------------------
# Visium
# -------------------------

visium_dir = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "visium"
    / "Br2039"
)

matrix_file = list(
    visium_dir.glob("*matrix.mtx.gz")
)[0]

features_file = list(
    visium_dir.glob("*features.tsv.gz")
)[0]

visium = sc.read_mtx(matrix_file).T

features = pd.read_csv(
    features_file,
    sep="\t",
    header=None
)

visium.var_names = (
    features.iloc[:, 1]
    .astype(str)
    .values
)

visium.var_names_make_unique()

# -------------------------
# Shared genes
# -------------------------

shared_genes = sorted(
    set(visium.var_names)
    &
    set(xenium.var_names)
)

print("\n===== Shared Genes =====")
print(f"Visium genes : {visium.n_vars:,}")
print(f"Xenium genes : {xenium.n_vars:,}")
print(f"Shared genes : {len(shared_genes):,}")

shared_df = pd.DataFrame({
    "gene": shared_genes
})

output_file = (
    OUT_DIR
    / "shared_genes.csv"
)

shared_df.to_csv(
    output_file,
    index=False
)

print(f"\nSaved to:")
print(output_file)