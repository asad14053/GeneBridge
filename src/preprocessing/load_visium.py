import scanpy as sc
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

sample_dir = PROJECT_ROOT / "data" / "raw" / "visium" / "Br2039"
OUT_DIR = PROJECT_ROOT / "outputs" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("Project root:", PROJECT_ROOT)
print("Sample dir:", sample_dir)
print("Exists:", sample_dir.exists())

matrix_files = list(sample_dir.glob("*matrix.mtx.gz"))
features_files = list(sample_dir.glob("*features.tsv.gz"))
barcodes_files = list(sample_dir.glob("*barcodes.tsv.gz"))

print("Matrix files found:", matrix_files)
print("Features files found:", features_files)
print("Barcodes files found:", barcodes_files)

matrix_file = matrix_files[0]
features_file = features_files[0]
barcodes_file = barcodes_files[0]

adata = sc.read_mtx(matrix_file).T

features = pd.read_csv(
    features_file,
    sep="\t",
    header=None
)

barcodes = pd.read_csv(
    barcodes_file,
    sep="\t",
    header=None
)

adata.var_names = features.iloc[:, 1].astype(str).values
adata.obs_names = barcodes.iloc[:, 0].astype(str).values

adata.var_names_make_unique()

print("\n===== Visium Summary =====")
print(f"Spots: {adata.n_obs:,}")
print(f"Genes: {adata.n_vars:,}")

gene_counts = adata.X.sum(axis=0)

if hasattr(gene_counts, "A1"):
    gene_counts = gene_counts.A1

top_genes = pd.DataFrame({
    "gene": adata.var_names,
    "counts": gene_counts
})

top_genes = top_genes.sort_values(
    "counts",
    ascending=False
)

print("\nTop 20 genes")
print(top_genes.head(20))

output_file = OUT_DIR / "visium_top_genes.csv"

top_genes.to_csv(
    output_file,
    index=False
)

print(f"Saved to: {output_file}")