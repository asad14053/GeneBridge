#!/usr/bin/env python3

import argparse
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


parser = argparse.ArgumentParser()
parser.add_argument("--task-id", type=int, required=True)
args = parser.parse_args()

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")
OUT = ROOT / "results/analysis/ling_reference_build"
CHUNKS = OUT / "array_chunks"
HARM = ROOT / "results/analysis/ling_huuki_kwon_gene_harmonization"

CHUNKS.mkdir(parents=True, exist_ok=True)

BUILDER = ROOT / "scripts/ling_harmonization/05_build_ling_huuki_compatible_reference.py"

spec = importlib.util.spec_from_file_location("lingbuilder", BUILDER)
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)


libs = pd.read_csv(
    OUT / "array_libraries.tsv",
    sep="\t",
)

row = libs.loc[
    libs["task_id"] == args.task_id
]

if len(row) != 1:
    raise RuntimeError(
        f"Task {args.task_id} not uniquely found"
    )

row = row.iloc[0]

prefix = row["prefix"]
h5_path = Path(row["h5_path"])

print("Task:", args.task_id)
print("Prefix:", prefix)
print("H5:", h5_path)


selected = pd.read_csv(
    OUT / "array_selected_cells.tsv.gz",
    sep="\t",
)

cells = (
    selected[
        selected["PREFIX"] == prefix
    ]
    .sort_values("reference_row")
    .copy()
)

if len(cells) != int(row["n_selected"]):
    raise RuntimeError(
        "Selected-cell count mismatch"
    )


genes = pd.read_csv(
    OUT / "array_gene_manifest.csv"
)

target_ens = genes["ensembl_id"].astype(str).tolist()

target_map = {
    ens: j
    for j, ens in enumerate(target_ens)
}


ling_path = B.find_file(
    [
        HARM / "Ling_union_symbol_to_ensembl.csv",
        HARM / "ling_union_symbol_to_ensembl.csv",
    ],
    search_root=HARM,
    contains="ling_union_symbol",
)

_, symbol_to_ens, _ = B.load_crosswalk(
    ling_path,
    "Ling",
)


with h5py.File(h5_path, "r") as h:

    barcodes = B.decode_array(
        h["matrix/barcodes"][:]
    )

    source_genes = B.decode_array(
        h["matrix/features/name"][:]
    )

    shape = tuple(
        int(v)
        for v in h["matrix/shape"][:]
    )

    data = h["matrix/data"][:]
    indices = h["matrix/indices"][:]
    indptr = h["matrix/indptr"][:]


if shape != (
    len(source_genes),
    len(barcodes),
):
    raise RuntimeError(
        f"Unexpected H5 matrix shape: {shape}"
    )


barcode_map = {
    b: i
    for i, b in enumerate(barcodes)
}

wanted = (
    cells["CELL_BARCODE"]
    .astype(str)
    .tolist()
)

missing = [
    b for b in wanted
    if b not in barcode_map
]

if missing:
    raise RuntimeError(
        f"{len(missing)} selected barcodes absent from H5"
    )

h5_cols = np.asarray(
    [barcode_map[b] for b in wanted],
    dtype=np.int64,
)


X_gene_cell = sparse.csc_matrix(
    (data, indices, indptr),
    shape=shape,
)

X_gene_cell = X_gene_cell[
    :,
    h5_cols,
]


source_rows = []
target_cols = []

for src_j, symbol in enumerate(source_genes):

    ens = symbol_to_ens.get(symbol)

    if ens is None:
        continue

    target_j = target_map.get(ens)

    if target_j is None:
        continue

    source_rows.append(src_j)
    target_cols.append(target_j)


M = sparse.csr_matrix(
    (
        np.ones(
            len(source_rows),
            dtype=np.float64,
        ),
        (
            np.asarray(source_rows),
            np.asarray(target_cols),
        ),
    ),
    shape=(
        len(source_genes),
        len(target_ens),
    ),
)


X = (
    X_gene_cell.T.tocsr()
    @ M
).tocsr()


if X.shape != (
    len(cells),
    len(target_ens),
):
    raise RuntimeError(
        f"Wrong output shape: {X.shape}"
    )


chunk_path = (
    CHUNKS
    / f"chunk_{args.task_id:03d}.npz"
)

row_path = (
    CHUNKS
    / f"rows_{args.task_id:03d}.npy"
)

sparse.save_npz(
    chunk_path,
    X,
    compressed=True,
)

np.save(
    row_path,
    cells["reference_row"].to_numpy(
        dtype=np.int64
    ),
)


print()
print("Output:", chunk_path)
print("Shape:", X.shape)
print("nnz:", X.nnz)
print("Rows:", row_path)
print("PASS")
