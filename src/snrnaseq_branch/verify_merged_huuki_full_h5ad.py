from pathlib import Path
import re

import scanpy as sc


# =============================================================================
# Paths
# =============================================================================

PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

SCE_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "sce_DLPFC_annotated"
)

CHUNK_DIR = (
    SCE_DIR
    / "huuki_h5ad_chunks_allgenes_full_maxcells0_chunk200"
)

MERGED_H5AD = (
    SCE_DIR
    / "huuki_snrna_reference_full_allgenes.h5ad"
)


# =============================================================================
# Expected dataset dimensions
# =============================================================================

EXPECTED_CHUNKS = 389
EXPECTED_CELLS = 77604
EXPECTED_GENES = 36601


# =============================================================================
# Helper function
# =============================================================================

def chunk_number(path: Path) -> int:
    """
    Extract chunk number from file names such as:
    huuki_allgenes_chunk_0001.h5ad
    """

    match = re.search(
        r"chunk_(\d+)\.h5ad$",
        path.name,
    )

    if match is None:
        raise ValueError(
            f"Cannot parse chunk number from: {path.name}"
        )

    return int(match.group(1))


# =============================================================================
# Start
# =============================================================================

print("=" * 70)
print("VERIFY HUUKI CHUNKS AND MERGED H5AD")
print("=" * 70)

print("Chunk directory:", CHUNK_DIR)
print("Merged H5AD:", MERGED_H5AD)
print("Expected chunks:", EXPECTED_CHUNKS)
print("Expected cells:", EXPECTED_CELLS)
print("Expected genes:", EXPECTED_GENES)


# =============================================================================
# Check chunk directory
# =============================================================================

if not CHUNK_DIR.exists():
    raise FileNotFoundError(
        f"Chunk directory does not exist:\n{CHUNK_DIR}"
    )


# =============================================================================
# Find chunk files
# =============================================================================

chunk_files = sorted(
    CHUNK_DIR.glob(
        "huuki_allgenes_chunk_*.h5ad"
    ),
    key=chunk_number,
)

print(
    "\nNumber of chunks found:",
    len(chunk_files),
)


# IMPORTANT:
# Check for zero chunks BEFORE using chunk_files[0]
if len(chunk_files) == 0:
    raise FileNotFoundError(
        "No H5AD chunks found in:\n"
        f"{CHUNK_DIR}"
    )


print(
    "First chunk:",
    chunk_files[0].name,
)

print(
    "Last chunk:",
    chunk_files[-1].name,
)


if len(chunk_files) != EXPECTED_CHUNKS:
    raise ValueError(
        f"Expected {EXPECTED_CHUNKS} chunks, "
        f"found {len(chunk_files)}"
    )


# =============================================================================
# Step 1: Verify chunk files
# =============================================================================

print("\n" + "=" * 70)
print("STEP 1: VERIFYING SOURCE CHUNKS")
print("=" * 70)


# Read first chunk and establish reference gene names/order
first = sc.read_h5ad(
    chunk_files[0],
    backed="r",
)

reference_genes = list(
    first.var_names.astype(str)
)

print(
    "First chunk shape:",
    first.shape,
)

first.file.close()


if len(reference_genes) != EXPECTED_GENES:
    raise ValueError(
        f"First chunk expected {EXPECTED_GENES} genes, "
        f"got {len(reference_genes)}"
    )


# =============================================================================
# Verify every chunk
# =============================================================================

all_chunk_cells = []

total_cells = 0


for i, chunk_file in enumerate(
    chunk_files,
    start=1,
):

    adata = sc.read_h5ad(
        chunk_file,
        backed="r",
    )


    # -------------------------------------------------------------------------
    # Check number of genes
    # -------------------------------------------------------------------------

    if adata.n_vars != EXPECTED_GENES:

        actual_genes = adata.n_vars

        adata.file.close()

        raise ValueError(
            f"{chunk_file.name}: "
            f"expected {EXPECTED_GENES} genes, "
            f"got {actual_genes}"
        )


    # -------------------------------------------------------------------------
    # Check gene names and order
    # -------------------------------------------------------------------------

    genes = list(
        adata.var_names.astype(str)
    )


    if genes != reference_genes:

        adata.file.close()

        raise ValueError(
            f"{chunk_file.name}: "
            "gene names/order do not match first chunk"
        )


    # -------------------------------------------------------------------------
    # Collect cell IDs
    # -------------------------------------------------------------------------

    cells = list(
        adata.obs_names.astype(str)
    )


    total_cells += adata.n_obs

    all_chunk_cells.extend(
        cells
    )


    adata.file.close()


    # -------------------------------------------------------------------------
    # Progress
    # -------------------------------------------------------------------------

    if (
        i % 50 == 0
        or i == len(chunk_files)
    ):

        print(
            f"Checked "
            f"{i}/{len(chunk_files)} chunks "
            f"| cells so far: "
            f"{total_cells}"
        )


# =============================================================================
# Check total cells
# =============================================================================

print(
    "\nTotal cells across chunks:",
    total_cells,
)


if total_cells != EXPECTED_CELLS:
    raise ValueError(
        f"Expected {EXPECTED_CELLS} cells across chunks, "
        f"got {total_cells}"
    )


# =============================================================================
# Check duplicate cell IDs
# =============================================================================

if len(all_chunk_cells) != len(
    set(all_chunk_cells)
):

    seen = set()
    duplicates = []


    for cell in all_chunk_cells:

        if cell in seen:

            duplicates.append(
                cell
            )


            if len(duplicates) >= 10:
                break


        seen.add(
            cell
        )


    raise ValueError(
        "Duplicated cell IDs in chunks. "
        f"Examples: {duplicates}"
    )


print(
    "\nPASS: chunk count is correct"
)

print(
    "PASS: total chunk cell count is correct"
)

print(
    "PASS: all chunks have the same genes/order"
)

print(
    "PASS: no duplicated cell IDs across chunks"
)


# =============================================================================
# Step 2: Check merged H5AD exists
# =============================================================================

print("\n" + "=" * 70)
print("STEP 2: VERIFYING MERGED H5AD")
print("=" * 70)


if not MERGED_H5AD.exists():
    raise FileNotFoundError(
        "Merged H5AD does not exist yet:\n"
        f"{MERGED_H5AD}\n\n"
        "Run the merge script successfully before "
        "running this verification script."
    )


# =============================================================================
# Open merged H5AD
# =============================================================================

merged = sc.read_h5ad(
    MERGED_H5AD,
    backed="r",
)


print(
    merged
)

print(
    "Merged shape:",
    merged.shape,
)


# =============================================================================
# Check merged dimensions
# =============================================================================

if merged.n_obs != EXPECTED_CELLS:

    actual_cells = merged.n_obs

    merged.file.close()

    raise ValueError(
        f"Merged cell count wrong: "
        f"expected {EXPECTED_CELLS}, "
        f"got {actual_cells}"
    )


if merged.n_vars != EXPECTED_GENES:

    actual_genes = merged.n_vars

    merged.file.close()

    raise ValueError(
        f"Merged gene count wrong: "
        f"expected {EXPECTED_GENES}, "
        f"got {actual_genes}"
    )


# =============================================================================
# Read merged genes and cells
# =============================================================================

merged_genes = list(
    merged.var_names.astype(str)
)


merged_cells = list(
    merged.obs_names.astype(str)
)


# =============================================================================
# Verify gene order
# =============================================================================

if merged_genes != reference_genes:

    merged.file.close()

    raise ValueError(
        "Merged gene names/order "
        "do not match source chunk genes"
    )


# =============================================================================
# Verify no duplicated merged cell IDs
# =============================================================================

if len(merged_cells) != len(
    set(merged_cells)
):

    merged.file.close()

    raise ValueError(
        "Merged H5AD has duplicated cell IDs"
    )


# =============================================================================
# Verify same cells exist before and after merge
# =============================================================================

chunk_cell_set = set(
    all_chunk_cells
)


merged_cell_set = set(
    merged_cells
)


missing_cells = (
    chunk_cell_set
    - merged_cell_set
)


extra_cells = (
    merged_cell_set
    - chunk_cell_set
)


if missing_cells:

    merged.file.close()

    raise ValueError(
        "Merged H5AD is missing cells. "
        f"Examples: {list(missing_cells)[:10]}"
    )


if extra_cells:

    merged.file.close()

    raise ValueError(
        "Merged H5AD has extra cells. "
        f"Examples: {list(extra_cells)[:10]}"
    )


merged.file.close()


# =============================================================================
# Success
# =============================================================================

print("\n" + "=" * 70)

print(
    "ALL CHECKS PASSED"
)

print(
    "=" * 70
)

print(
    "Chunks:",
    len(chunk_files),
)

print(
    "Cells in chunks:",
    total_cells,
)

print(
    "Cells in merged:",
    EXPECTED_CELLS,
)

print(
    "Genes:",
    EXPECTED_GENES,
)

print(
    "Output:",
    MERGED_H5AD,
)

print(
    "=" * 70
)