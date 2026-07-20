from pathlib import Path
import re
import shutil
import scanpy as sc

try:
    from anndata.experimental import concat_on_disk
except Exception as e:
    raise ImportError(
        "Your anndata does not support concat_on_disk. "
        "Try: pip install -U anndata\n"
        f"Original error: {e}"
    )

PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

SCE_DIR = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
)

CHUNK_DIR = SCE_DIR / "huuki_h5ad_chunks_allgenes_full_maxcells0_chunk50"

OUT_H5AD = SCE_DIR / "huuki_snrna_reference_full_allgenes.h5ad"

TMP_DIR = SCE_DIR / "huuki_merge_batches_tmp"

EXPECTED_CHUNKS = 1553
EXPECTED_CELLS = 77604
EXPECTED_GENES = 36601

# Merge 50 chunks at a time first.
# 1553 chunks -> around 32 batch h5ad files.
BATCH_SIZE = 50


def chunk_number(path: Path) -> int:
    m = re.search(r"chunk_(\d+)\.h5ad$", path.name)
    if m is None:
        raise ValueError(f"Cannot parse chunk number from: {path.name}")
    return int(m.group(1))


def verify_h5ad(path: Path, expected_cells=None, expected_genes=None):
    a = sc.read_h5ad(path, backed="r")
    print(a)
    n_obs, n_vars = a.shape
    a.file.close()

    if expected_cells is not None and n_obs != expected_cells:
        raise ValueError(f"{path.name}: expected {expected_cells} cells, got {n_obs}")

    if expected_genes is not None and n_vars != expected_genes:
        raise ValueError(f"{path.name}: expected {expected_genes} genes, got {n_vars}")

    return n_obs, n_vars


print("Chunk dir:", CHUNK_DIR)
print("Output:", OUT_H5AD)
print("Temp dir:", TMP_DIR)

if not CHUNK_DIR.exists():
    raise FileNotFoundError(CHUNK_DIR)

chunk_files = sorted(
    CHUNK_DIR.glob("huuki_allgenes_chunk_*.h5ad"),
    key=chunk_number,
)

print("Found chunks:", len(chunk_files))
print("First:", chunk_files[0].name)
print("Last:", chunk_files[-1].name)

if len(chunk_files) != EXPECTED_CHUNKS:
    raise ValueError(f"Expected {EXPECTED_CHUNKS} chunks, found {len(chunk_files)}")

# ------------------------------------------------------------
# Step 1: verify chunks quickly
# ------------------------------------------------------------
print("\nStep 1: verifying source chunks...")

first = sc.read_h5ad(chunk_files[0], backed="r")
reference_genes = list(first.var_names.astype(str))
first.file.close()

total_cells = 0
seen_cells = set()

for i, f in enumerate(chunk_files, start=1):
    a = sc.read_h5ad(f, backed="r")

    if a.n_vars != EXPECTED_GENES:
        raise ValueError(f"{f.name}: expected {EXPECTED_GENES} genes, got {a.n_vars}")

    if list(a.var_names.astype(str)) != reference_genes:
        raise ValueError(f"{f.name}: gene order mismatch")

    cells = list(a.obs_names.astype(str))
    overlap = seen_cells.intersection(cells)

    if overlap:
        raise ValueError(f"Duplicated cells in {f.name}: {list(overlap)[:10]}")

    seen_cells.update(cells)
    total_cells += a.n_obs

    if i % 100 == 0 or i == len(chunk_files):
        print(f"Checked {i}/{len(chunk_files)} chunks | cells: {total_cells}")

    a.file.close()

if total_cells != EXPECTED_CELLS:
    raise ValueError(f"Expected {EXPECTED_CELLS} total cells, got {total_cells}")

print("PASS: source chunks are valid")

# ------------------------------------------------------------
# Step 2: remove broken final output and old temp batches
# ------------------------------------------------------------
if OUT_H5AD.exists():
    print("\nRemoving broken/old final output:", OUT_H5AD)
    OUT_H5AD.unlink()

if TMP_DIR.exists():
    print("Removing old temp dir:", TMP_DIR)
    shutil.rmtree(TMP_DIR)

TMP_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Step 3: merge chunks into batch h5ads
# ------------------------------------------------------------
print("\nStep 2: creating batch h5ad files...")

batch_files = []

for batch_start in range(0, len(chunk_files), BATCH_SIZE):
    batch_id = (batch_start // BATCH_SIZE) + 1
    batch_chunk_files = chunk_files[batch_start: batch_start + BATCH_SIZE]

    batch_out = TMP_DIR / f"huuki_batch_{batch_id:03d}.h5ad"
    batch_tmp = TMP_DIR / f"huuki_batch_{batch_id:03d}.tmp.h5ad"

    if batch_tmp.exists():
        batch_tmp.unlink()
    if batch_out.exists():
        batch_out.unlink()

    expected_batch_cells = 0
    for f in batch_chunk_files:
        a = sc.read_h5ad(f, backed="r")
        expected_batch_cells += a.n_obs
        a.file.close()

    print(
        f"\nBatch {batch_id}: "
        f"{len(batch_chunk_files)} chunks, "
        f"{expected_batch_cells} cells"
    )

    concat_on_disk(
        [str(f) for f in batch_chunk_files],
        str(batch_tmp),
        axis=0,
        join="inner",
        label="source_chunk",
        keys=[f"chunk_{chunk_number(f)}" for f in batch_chunk_files],
        index_unique=None,
        max_loaded_elems=20_000_000,
    )

    print("Verifying batch:", batch_tmp.name)
    verify_h5ad(
        batch_tmp,
        expected_cells=expected_batch_cells,
        expected_genes=EXPECTED_GENES,
    )

    batch_tmp.rename(batch_out)
    batch_files.append(batch_out)

print("\nCreated batch files:", len(batch_files))

# ------------------------------------------------------------
# Step 4: merge batch h5ads into final h5ad
# ------------------------------------------------------------
print("\nStep 3: merging batch h5ads into final h5ad...")

FINAL_TMP = SCE_DIR / "huuki_snrna_reference_full_allgenes.tmp.h5ad"

if FINAL_TMP.exists():
    FINAL_TMP.unlink()

concat_on_disk(
    [str(f) for f in batch_files],
    str(FINAL_TMP),
    axis=0,
    join="inner",
    label="merge_batch",
    keys=[f"batch_{i+1:03d}" for i in range(len(batch_files))],
    index_unique=None,
    max_loaded_elems=20_000_000,
)

print("\nVerifying final temp file...")
verify_h5ad(
    FINAL_TMP,
    expected_cells=EXPECTED_CELLS,
    expected_genes=EXPECTED_GENES,
)

FINAL_TMP.rename(OUT_H5AD)

print("\n============================================================")
print("MERGE SUCCESSFUL")
print("============================================================")
print("Final h5ad:", OUT_H5AD)
print("Cells:", EXPECTED_CELLS)
print("Genes:", EXPECTED_GENES)
print("Temporary batch folder:", TMP_DIR)
print("============================================================")