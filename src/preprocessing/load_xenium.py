import scanpy as sc
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

sample_dir = PROJECT_ROOT / "data" / "raw" / "xenium" / "Br2039"
OUT_DIR = PROJECT_ROOT / "outputs" / "tables"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("Project root:", PROJECT_ROOT)
print("Sample dir:", sample_dir)
print("Exists:", sample_dir.exists())

h5_files = list(sample_dir.glob("*cell_feature_matrix.h5"))
print("H5 files found:", h5_files)

h5_file = h5_files[0]

adata = sc.read_10x_h5(h5_file)

print("\n===== Xenium Summary =====")
print(f"Cells: {adata.n_obs:,}")
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

output_file = OUT_DIR / "xenium_top_genes.csv"

top_genes.to_csv(
    output_file,
    index=False
)

print(f"Saved to: {output_file}")