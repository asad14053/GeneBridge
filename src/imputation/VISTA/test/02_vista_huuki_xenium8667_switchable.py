from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import torch

# ============================================================
# VISTA switchable script
# Huuki-Myers full snRNA reference + Xenium Br8667
#
# smoke = small subset, quick training
# full  = all cells, heavier training
#
# No donor/reference matching:
#   Huuki is treated as one pooled reference atlas.
# ============================================================

parser = argparse.ArgumentParser()

parser.add_argument(
    "--mode",
    choices=["smoke", "full"],
    default="smoke",
    help="smoke = quick subset test; full = full Huuki + full Xenium training",
)

parser.add_argument(
    "--epochs",
    type=int,
    default=None,
    help="Override number of training epochs",
)

parser.add_argument(
    "--extra-reference-genes",
    type=int,
    default=None,
    help="Number of Huuki-only genes to include beyond shared Huuki-Xenium genes",
)

parser.add_argument(
    "--ref-cells",
    type=int,
    default=None,
    help="Override number of Huuki reference cells for smoke/debug",
)

parser.add_argument(
    "--xenium-cells",
    type=int,
    default=None,
    help="Override number of Xenium cells for smoke/debug",
)

args = parser.parse_args()

# -----------------------------
# Paths
# -----------------------------
PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

VISTA_ROOT = PROJECT_ROOT / "src/imputation/VISTA"
sys.path.insert(0, str(VISTA_ROOT))

from vista import GIMVI_GCN

# ============================================================
# Patch VISTA spatial KNN graph to avoid FAISS -1 neighbor bug
# ============================================================
import scipy.sparse as sp_sparse
from scipy.spatial import cKDTree
import vista.faiss_neig as faiss_neig


def safe_knn_graph(coords, n_neighbors=20):
    coords = np.asarray(coords, dtype=np.float32)

    if coords.ndim != 2:
        raise ValueError(f"Spatial coordinates must be 2D, got shape {coords.shape}")

    if coords.shape[0] == 0:
        raise ValueError("No spatial cells available for KNN graph.")

    if not np.isfinite(coords).all():
        bad = np.where(~np.isfinite(coords).all(axis=1))[0][:10]
        raise ValueError(f"Spatial coordinates contain NaN/Inf. First bad rows: {bad}")

    n_samples = coords.shape[0]

    if n_samples < 2:
        return sp_sparse.csr_matrix((n_samples, n_samples), dtype=np.float32)

    # +1 because nearest neighbor may be the point itself
    k = min(int(n_neighbors) + 1, n_samples)

    tree = cKDTree(coords)
    distances, indices = tree.query(coords, k=k)

    if k == 1:
        distances = distances[:, None]
        indices = indices[:, None]

    row_indices = np.repeat(np.arange(n_samples), k)
    col_indices = indices.reshape(-1)
    distance_data = distances.reshape(-1)

    # Remove self edges and invalid neighbors
    valid = (
        (col_indices >= 0)
        & (col_indices < n_samples)
        & np.isfinite(distance_data)
        & (row_indices != col_indices)
    )

    row_indices = row_indices[valid]
    col_indices = col_indices[valid]
    distance_data = distance_data[valid]

    return sp_sparse.csr_matrix(
        (distance_data, (row_indices, col_indices)),
        shape=(n_samples, n_samples),
        dtype=np.float32,
    )


faiss_neig.knn_graph = safe_knn_graph
print("Patched VISTA faiss_neig.knn_graph with safe scipy cKDTree version.")

HUUKI_REF_H5AD = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

XENIUM_H5AD = (
    PROJECT_ROOT
    / "data/processed/imputation_beta/Br8667"
    / "spatial_data_xenium_Br8667_vista.h5ad"
)

OUT_BASE = PROJECT_ROOT / "outputs/imputation_beta/Br8667/vista_full_huuki_pooled"
PRED_DIR = PROJECT_ROOT / "data/processed/imputation_beta/Br8667/benchmark/predictions"

OUT_BASE.mkdir(parents=True, exist_ok=True)
PRED_DIR.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Mode settings
# -----------------------------
if args.mode == "smoke":
    MAX_REF_CELLS = 5000
    MAX_XENIUM_CELLS = 5000
    MAX_EPOCHS = 2
    N_EXTRA_REFERENCE_GENES = 100
    BATCH_SIZE = 128

else:
    MAX_REF_CELLS = None
    MAX_XENIUM_CELLS = None
    MAX_EPOCHS = 50
    N_EXTRA_REFERENCE_GENES = 500
    BATCH_SIZE = 128

if args.epochs is not None:
    MAX_EPOCHS = args.epochs

if args.extra_reference_genes is not None:
    N_EXTRA_REFERENCE_GENES = args.extra_reference_genes

if args.ref_cells is not None:
    MAX_REF_CELLS = args.ref_cells

if args.xenium_cells is not None:
    MAX_XENIUM_CELLS = args.xenium_cells

N_LATENT = 32
NEIGHBOR_SIZE = 20
SEED = 42
USE_GPU = torch.cuda.is_available()

print("=" * 80)
print("VISTA Huuki pooled reference + Xenium Br8667")
print("=" * 80)
print("Mode:", args.mode)
print("Huuki reference:", HUUKI_REF_H5AD)
print("Xenium target:", XENIUM_H5AD)
print("Max Huuki cells:", MAX_REF_CELLS)
print("Max Xenium cells:", MAX_XENIUM_CELLS)
print("Epochs:", MAX_EPOCHS)
print("Extra Huuki-only genes:", N_EXTRA_REFERENCE_GENES)
print("Batch size:", BATCH_SIZE)
print("GPU available:", USE_GPU)
if USE_GPU:
    print("GPU:", torch.cuda.get_device_name(0))

if not HUUKI_REF_H5AD.exists():
    raise FileNotFoundError(HUUKI_REF_H5AD)

if not XENIUM_H5AD.exists():
    raise FileNotFoundError(XENIUM_H5AD)

# -----------------------------
# Load backed first, so smoke does not load full Huuki into RAM
# -----------------------------
print("\nOpening Huuki reference in backed mode...")
huuki_backed = sc.read_h5ad(HUUKI_REF_H5AD, backed="r")

print("\nOpening Xenium in backed mode...")
xenium_backed = sc.read_h5ad(XENIUM_H5AD, backed="r")

print("\nHuuki backed:", huuki_backed)
print("Xenium backed:", xenium_backed)

huuki_genes = pd.Index(huuki_backed.var_names.astype(str))
xenium_genes = pd.Index(xenium_backed.var_names.astype(str))

shared_genes = sorted(set(huuki_genes).intersection(set(xenium_genes)))
huuki_only_genes = sorted(set(huuki_genes) - set(xenium_genes))

print("\nShared Huuki-Xenium genes:", len(shared_genes))
print("Huuki-only genes:", len(huuki_only_genes))
print("First shared genes:", shared_genes[:20])

if len(shared_genes) < 50:
    raise ValueError("Too few shared genes. Check gene symbols / var_names.")

info_genes = shared_genes
if N_EXTRA_REFERENCE_GENES == -1:
    extra_genes = huuki_only_genes
else:
    extra_genes = huuki_only_genes[:N_EXTRA_REFERENCE_GENES]
    
gene_for_impute = info_genes + extra_genes

print("\nGenes used:")
print("Xenium info genes:", len(info_genes))
print("Extra Huuki-only genes:", len(extra_genes))
print("Total model genes:", len(gene_for_impute))

# -----------------------------
# Select cells
# -----------------------------
rng = np.random.default_rng(SEED)

if MAX_REF_CELLS is not None and huuki_backed.n_obs > MAX_REF_CELLS:
    ref_idx = rng.choice(huuki_backed.n_obs, size=MAX_REF_CELLS, replace=False)
else:
    ref_idx = np.arange(huuki_backed.n_obs)

if MAX_XENIUM_CELLS is not None and xenium_backed.n_obs > MAX_XENIUM_CELLS:
    xenium_idx = rng.choice(xenium_backed.n_obs, size=MAX_XENIUM_CELLS, replace=False)
else:
    xenium_idx = np.arange(xenium_backed.n_obs)

print("\nLoading selected data into memory...")
huuki_model = huuki_backed[ref_idx, gene_for_impute].to_memory()
xenium_model = xenium_backed[xenium_idx, info_genes].to_memory()

huuki_backed.file.close()
xenium_backed.file.close()

huuki_model.var_names = huuki_model.var_names.astype(str)
xenium_model.var_names = xenium_model.var_names.astype(str)
huuki_model.obs_names = huuki_model.obs_names.astype(str)
xenium_model.obs_names = xenium_model.obs_names.astype(str)

huuki_model.var_names_make_unique()
xenium_model.var_names_make_unique()
huuki_model.obs_names_make_unique()
xenium_model.obs_names_make_unique()

print("\nHuuki model:", huuki_model)
print("Xenium model:", xenium_model)

if "spatial" not in xenium_model.obsm:
    raise ValueError("Xenium AnnData does not contain obsm['spatial'].")

print("Xenium spatial shape:", xenium_model.obsm["spatial"].shape)

# -----------------------------
# Required obs columns
# -----------------------------
huuki_model.obs["names"] = huuki_model.obs_names.astype(str)
xenium_model.obs["names"] = xenium_model.obs_names.astype(str)

# No donor/reference matching
huuki_model.obs["reference_batch"] = "Huuki_pooled"
xenium_model.obs["batch"] = "Br8667_xenium"

# -----------------------------
# VISTA setup
# -----------------------------
print("\nSetting up VISTA AnnData...")

GIMVI_GCN.setup_anndata(
    huuki_model,
)

GIMVI_GCN.setup_anndata(
    xenium_model,
    batch_key="batch",
    obs_names="names",
)

# -----------------------------
# Train
# -----------------------------
print("\nTraining VISTA...")

model = GIMVI_GCN(
    huuki_model,
    xenium_model,
    n_latent=N_LATENT,
    neighbor_size=NEIGHBOR_SIZE,
)

model.train(
    max_epochs=MAX_EPOCHS,
    train_size=1.0,
    validation_size=0.0,
    use_gpu=USE_GPU,
    batch_size=BATCH_SIZE,
)

print("\nTraining finished.")

# -----------------------------
# Predict
# -----------------------------
print("\nGenerating imputed values...")

pred = model.get_imputed_values(
    normalized=False,
    batch_size=BATCH_SIZE,
)[0]

print("Prediction shape:", pred.shape)

pred_adata = sc.AnnData(
    X=pred,
    obs=xenium_model.obs.copy(),
    var=pd.DataFrame(index=pd.Index(gene_for_impute, name="gene")),
)

pred_adata.obsm["spatial"] = xenium_model.obsm["spatial"].copy()

out_tag = f"{args.mode}_epochs{MAX_EPOCHS}_extra{N_EXTRA_REFERENCE_GENES}"

OUT_H5AD = PRED_DIR / f"vista_huuki_pooled_xenium8667_{out_tag}_predictions.h5ad"
OUT_SUMMARY = OUT_BASE / f"vista_huuki_pooled_xenium8667_{out_tag}_summary.csv"

pred_adata.write_h5ad(OUT_H5AD, compression="gzip")

print("\nSaved prediction:")
print(OUT_H5AD)
print(pred_adata)

summary = pd.DataFrame(
    {
        "item": [
            "mode",
            "huuki_input_shape",
            "xenium_input_shape",
            "huuki_model_shape",
            "xenium_model_shape",
            "shared_genes",
            "extra_reference_genes",
            "total_model_genes",
            "epochs",
            "batch_size",
            "n_latent",
            "neighbor_size",
            "use_gpu",
            "output_prediction_h5ad",
        ],
        "value": [
            args.mode,
            str((77604, 36601)),
            str((66164, 300)),
            str(huuki_model.shape),
            str(xenium_model.shape),
            len(shared_genes),
            len(extra_genes),
            len(gene_for_impute),
            MAX_EPOCHS,
            BATCH_SIZE,
            N_LATENT,
            NEIGHBOR_SIZE,
            USE_GPU,
            str(OUT_H5AD),
        ],
    }
)

summary.to_csv(OUT_SUMMARY, index=False)

print("\nSaved summary:")
print(OUT_SUMMARY)
print("\nDONE.")
