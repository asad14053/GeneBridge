from pathlib import Path
import re
import scanpy as sc

PROJECT_ROOT = Path("/users/mjabin/projects/GeneBridge")

CHUNK_DIR = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
    / "huuki_h5ad_chunks_allgenes_full_maxcells0_chunk50"
)

MERGED_H5AD = (
    PROJECT_ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated"
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

EXPECTED_CHUNKS = 1553
EXPECTED_CELLS = 77604
EXPECTED_GENES = 36601


def chunk_number(path: Path) -> int:
    m = re.search(r"chunk_(\d+)\.h5ad$", path.name)
    if m is None:
        raise ValueError(f"Cannot parse chunk number from: {path.name}")
    return int(m.group(1))


print("Chunk dir:", CHUNK_DIR)
print("Merged h5ad:", MERGED_H5AD)

if not CHUNK_DIR.exists():
    raise FileNotFoundError(CHUNK_DIR)

if not MERGED_H5AD.exists():
    raise FileNotFoundError(MERGED_H5AD)

chunk_files = sorted(
    CHUNK_DIR.glob("huuki_allgenes_chunk_*.h5ad"),
    key=chunk_number,
)

print("\nNumber of chunks found:", len(chunk_files))
print("First chunk:", chunk_files[0].name)
print("Last chunk:", chunk_files[-1].name)

if len(chunk_files) != EXPECTED_CHUNKS:
    raise ValueError(f"Expected {EXPECTED_CHUNKS} chunks, found {len(chunk_files)}")

# ------------------------------------------------------------
# Verify chunk files
# ------------------------------------------------------------
print("\nVerifying chunks...")

first = sc.read_h5ad(chunk_files[0], backed="r")
reference_genes = list(first.var_names.astype(str))
first.file.close()

all_chunk_cells = []
total_cells = 0

for i, f in enumerate(chunk_files, start=1):
    adata = sc.read_h5ad(f, backed="r")

    genes = list(adata.var_names.astype(str))
    cells = list(adata.obs_names.astype(str))

    if adata.n_vars != EXPECTED_GENES:
        raise ValueError(f"{f.name}: expected {EXPECTED_GENES} genes, got {adata.n_vars}")

    if genes != reference_genes:
        raise ValueError(f"{f.name}: gene names/order do not match first chunk")

    total_cells += adata.n_obs
    all_chunk_cells.extend(cells)

    if i % 100 == 0 or i == len(chunk_files):
        print(f"Checked {i}/{len(chunk_files)} chunks | cells so far: {total_cells}")

    adata.file.close()

print("\nTotal cells across chunks:", total_cells)

if total_cells != EXPECTED_CELLS:
    raise ValueError(f"Expected {EXPECTED_CELLS} cells across chunks, got {total_cells}")

if len(all_chunk_cells) != len(set(all_chunk_cells)):
    seen = set()
    duplicates = []
    for c in all_chunk_cells:
        if c in seen:
            duplicates.append(c)
            if len(duplicates) >= 10:
                break
        seen.add(c)
    raise ValueError(f"Duplicated cell IDs in chunks. Examples: {duplicates}")

print("PASS: chunk count is correct")
print("PASS: total chunk cell count is correct")
print("PASS: all chunks have same genes/order")
print("PASS: no duplicated cell IDs across chunks")

# ------------------------------------------------------------
# Verify merged file
# ------------------------------------------------------------
print("\nVerifying merged h5ad...")

merged = sc.read_h5ad(MERGED_H5AD, backed="r")

print(merged)
print("Merged shape:", merged.shape)

if merged.n_obs != EXPECTED_CELLS:
    raise ValueError(f"Merged cell count wrong: expected {EXPECTED_CELLS}, got {merged.n_obs}")

if merged.n_vars != EXPECTED_GENES:
    raise ValueError(f"Merged gene count wrong: expected {EXPECTED_GENES}, got {merged.n_vars}")

merged_genes = list(merged.var_names.astype(str))
merged_cells = list(merged.obs_names.astype(str))

if merged_genes != reference_genes:
    raise ValueError("Merged gene names/order do not match chunk genes")

if len(merged_cells) != len(set(merged_cells)):
    raise ValueError("Merged h5ad has duplicated cell IDs")

missing_cells = set(all_chunk_cells) - set(merged_cells)
extra_cells = set(merged_cells) - set(all_chunk_cells)

if missing_cells:
    raise ValueError(f"Merged h5ad is missing cells. Examples: {list(missing_cells)[:10]}")

if extra_cells:
    raise ValueError(f"Merged h5ad has extra cells. Examples: {list(extra_cells)[:10]}")

merged.file.close()

print("\n============================================================")
print("ALL CHECKS PASSED")
print("============================================================")
print("Chunks:", len(chunk_files))
print("Cells in chunks:", total_cells)
print("Cells in merged:", EXPECTED_CELLS)
print("Genes:", EXPECTED_GENES)
print("Output:", MERGED_H5AD)
print("============================================================")