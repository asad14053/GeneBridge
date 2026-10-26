#!/usr/bin/env python3

import argparse
import gc
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


parser = argparse.ArgumentParser()
parser.add_argument("--work-dir", required=True)
args = parser.parse_args()


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

R = ROOT / "results/analysis/ling_reference_build"
CHUNKS = R / "array_chunks"

OUTPUT = (
    ROOT
    / "data/processed/snrnaseq/ling_compatible/"
      "ling_snrna_reference_full_allgenes.h5ad"
)

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

WORK = Path(args.work_dir)
WORK.mkdir(parents=True, exist_ok=True)


cells = pd.read_csv(
    R / "array_selected_cells.tsv.gz",
    sep="\t",
)

genes = pd.read_csv(
    R / "array_gene_manifest.csv"
)

libs = pd.read_csv(
    R / "array_libraries.tsv",
    sep="\t",
)


N = len(cells)
G = len(genes)

assert N == 77604, N
assert G == 31491, G
assert len(libs) == 87


print("Cells:", N)
print("Genes:", G)
print("Libraries:", len(libs))


X_PATH = WORK / "ling_X_float64.dat"
L_PATH = WORK / "ling_logcounts_float64.dat"

X = np.memmap(
    X_PATH,
    dtype=np.float64,
    mode="w+",
    shape=(N, G),
)

X[:] = 0.0


libsize = np.zeros(N, dtype=np.float64)
detected = np.zeros(N, dtype=np.float64)

symbols = (
    genes["canonical_symbol"]
    .astype(str)
    .tolist()
)

mito_cols = np.asarray(
    [
        i for i, g in enumerate(symbols)
        if g.upper().startswith("MT-")
    ],
    dtype=int,
)

mito_sum = np.zeros(N, dtype=np.float64)
mito_detected = np.zeros(N, dtype=np.float64)


for task in range(87):

    p = CHUNKS / f"chunk_{task:03d}.npz"
    r = CHUNKS / f"rows_{task:03d}.npy"

    if not p.exists():
        raise FileNotFoundError(p)

    if not r.exists():
        raise FileNotFoundError(r)

    S = sparse.load_npz(p).tocsr()
    rows = np.load(r)

    if S.shape != (
        len(rows),
        G,
    ):
        raise RuntimeError(
            f"Task {task}: wrong shape {S.shape}"
        )

    block = S.toarray().astype(
        np.float64,
        copy=False,
    )

    X[rows, :] = block

    libsize[rows] = block.sum(axis=1)

    detected[rows] = np.count_nonzero(
        block,
        axis=1,
    )

    if len(mito_cols):
        M = block[:, mito_cols]

        mito_sum[rows] = M.sum(axis=1)

        mito_detected[rows] = np.count_nonzero(
            M,
            axis=1,
        )

    print(
        f"{task:02d}/86 "
        f"cells={len(rows):,} "
        f"nnz={S.nnz:,}"
    )

    del S, block
    gc.collect()


X.flush()


if np.any(libsize <= 0):
    raise RuntimeError(
        f"Zero-count cells: {(libsize <= 0).sum()}"
    )


# ======================================================================
# Huuki-identical normalization
# ======================================================================

mean_library = float(
    libsize.mean()
)

sizeFactor = (
    libsize
    / mean_library
)

print()
print("Mean library size:", mean_library)
print("Mean sizeFactor:", sizeFactor.mean())


L = np.memmap(
    L_PATH,
    dtype=np.float64,
    mode="w+",
    shape=(N, G),
)


for start in range(0, N, 512):

    end = min(
        start + 512,
        N,
    )

    L[start:end, :] = np.log2(
        X[start:end, :]
        / sizeFactor[start:end, None]
        + 1.0
    )

    if start % 5120 == 0:
        print(
            f"logcounts {end:,}/{N:,}"
        )

L.flush()


# ======================================================================
# Huuki-compatible obs = exactly 34 columns
# ======================================================================

cells = cells.sort_values(
    "reference_row"
).reset_index(drop=True)

donor = (
    cells["DONOR"]
    .astype(str)
    .to_numpy()
)

barcode = (
    cells["CELL_BARCODE"]
    .astype(str)
    .to_numpy()
)

prefix = (
    cells["PREFIX"]
    .astype(str)
    .to_numpy()
)

cell_id = (
    cells["cell_id"]
    .astype(str)
    .to_numpy()
)

broad = (
    cells["cellclass"]
    .fillna("NA")
    .astype(str)
    .to_numpy()
)

subclass = (
    cells["subclass"]
    .fillna("NA")
    .astype(str)
    .to_numpy()
)


obs = pd.DataFrame(
    index=pd.Index(cell_id)
)

obs["Sample"] = donor
obs["Barcode"] = barcode
obs["key"] = cell_id
obs["SAMPLE_ID"] = donor
obs["pos"] = np.nan
obs["BrNum"] = donor
obs["round"] = np.nan
obs["Position"] = np.nan

# Until donor-level Ling table is added:
obs["age"] = np.nan
obs["sex"] = "NA"
obs["diagnosis"] = "NA"

obs["sum"] = libsize
obs["detected"] = detected

obs["subsets_Mito_sum"] = mito_sum
obs["subsets_Mito_detected"] = mito_detected
obs["subsets_Mito_percent"] = (
    100.0 * mito_sum / libsize
)

obs["total"] = libsize

obs["high_mito"] = np.nan
obs["low_sum"] = np.nan
obs["low_detected"] = np.nan
obs["discard_auto"] = np.nan
obs["doubletScore"] = np.nan

obs["prelimCluster"] = "NA"
obs["collapsedCluster"] = "NA"
obs["kmeans"] = "NA"

obs["sizeFactor"] = sizeFactor

obs["cellType_broad_k"] = broad
obs["cellType_k"] = subclass
obs["cellType_broad_hc"] = broad
obs["cellType_hc"] = subclass

obs["cellType_layer"] = "NA"
obs["layer_annotation"] = "NA"

obs["source_chunk"] = prefix
obs["merge_batch"] = prefix


HUUKI = [
    "Sample",
    "Barcode",
    "key",
    "SAMPLE_ID",
    "pos",
    "BrNum",
    "round",
    "Position",
    "age",
    "sex",
    "diagnosis",
    "sum",
    "detected",
    "subsets_Mito_sum",
    "subsets_Mito_detected",
    "subsets_Mito_percent",
    "total",
    "high_mito",
    "low_sum",
    "low_detected",
    "discard_auto",
    "doubletScore",
    "prelimCluster",
    "collapsedCluster",
    "kmeans",
    "sizeFactor",
    "cellType_broad_k",
    "cellType_k",
    "cellType_broad_hc",
    "cellType_hc",
    "cellType_layer",
    "layer_annotation",
    "source_chunk",
    "merge_batch",
]

obs = obs[HUUKI]

assert obs.shape == (77604, 34)


# ======================================================================
# AnnData
# ======================================================================

var = pd.DataFrame(
    index=pd.Index(symbols)
)

assert var.index.is_unique

adata = ad.AnnData(
    X=X,
    obs=obs,
    var=var,
)

adata.layers["logcounts"] = L

print()
print("Writing:")
print(OUTPUT)

adata.write_h5ad(
    OUTPUT,
    compression="gzip",
    compression_opts=4,
)

del adata
X.flush()
L.flush()
del X, L

gc.collect()


# ======================================================================
# Validation
# ======================================================================

a = ad.read_h5ad(
    OUTPUT,
    backed="r",
)

print()
print("FINAL")
print("Shape:", a.shape)
print("OBS:", a.obs.shape)
print("VAR:", a.var.shape)
print("Layers:", list(a.layers.keys()))
print("OBSM:", list(a.obsm.keys()))
print("X dtype:", a.X.dtype)
print("Unique var_names:", a.var_names.is_unique)

assert a.shape == (77604, 31491)
assert a.obs.shape == (77604, 34)
assert "logcounts" in a.layers
assert len(a.obsm) == 0
assert a.var_names.is_unique

a.file.close()


with h5py.File(OUTPUT, "r") as h:
    print("HDF5 X:", h["X"])
    print(
        "HDF5 logcounts:",
        h["layers/logcounts"]
    )


print()
print("REFERENCE BUILD: PASS")
print(OUTPUT)
