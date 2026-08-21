from pathlib import Path
import re
import shutil

import scanpy as sc

try:
    from anndata.experimental import concat_on_disk
except Exception as e:
    raise ImportError(
        "Your anndata does not support concat_on_disk. "
        "Try upgrading anndata.\n"
        f"Original error: {e}"
    )


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

# IMPORTANT:
# The R conversion used:
#
#   full 200 0
#
# Therefore:
#   CELLS_PER_CHUNK = 200
#   MAX_CELLS = 0
#
CHUNK_DIR = (
    SCE_DIR
    / "huuki_h5ad_chunks_allgenes_full_maxcells0_chunk200"
)

OUT_H5AD = (
    SCE_DIR
    / "huuki_snrna_reference_full_allgenes.h5ad"
)

TMP_DIR = (
    SCE_DIR
    / "huuki_merge_batches_tmp"
)


# =============================================================================
# Expected dimensions
# =============================================================================

EXPECTED_CHUNKS = 389
EXPECTED_CELLS = 77604
EXPECTED_GENES = 36601

# Merge 50 H5AD chunks at a time.
#
# 389 chunks / 50
# -> approximately 8 temporary batch H5AD files.
#
BATCH_SIZE = 50


# =============================================================================
# Helper: extract chunk number
# =============================================================================

def chunk_number(path: Path) -> int:
    """
    Extract numeric chunk ID from:

        huuki_allgenes_chunk_0001.h5ad
        huuki_allgenes_chunk_0002.h5ad
        ...

    Returns:
        1, 2, ...
    """

    match = re.search(
        r"chunk_(\d+)\.h5ad$",
        path.name,
    )

    if match is None:
        raise ValueError(
            f"Cannot parse chunk number from: {path.name}"
        )

    return int(
        match.group(1)
    )


# =============================================================================
# Helper: verify H5AD
# =============================================================================

def verify_h5ad(
    path: Path,
    expected_cells: int | None = None,
    expected_genes: int | None = None,
):
    """
    Open H5AD in backed mode and verify dimensions.
    """

    adata = sc.read_h5ad(
        path,
        backed="r",
    )

    n_obs = adata.n_obs
    n_vars = adata.n_vars

    print(
        f"{path.name}: "
        f"{n_obs} cells x {n_vars} genes"
    )

    adata.file.close()

    if (
        expected_cells is not None
        and n_obs != expected_cells
    ):
        raise ValueError(
            f"{path.name}: expected "
            f"{expected_cells} cells, "
            f"got {n_obs}"
        )

    if (
        expected_genes is not None
        and n_vars != expected_genes
    ):
        raise ValueError(
            f"{path.name}: expected "
            f"{expected_genes} genes, "
            f"got {n_vars}"
        )

    return n_obs, n_vars


# =============================================================================
# Start
# =============================================================================

print("=" * 70)
print("Huuki H5AD chunk merge")
print("=" * 70)

print(
    "Chunk directory:",
    CHUNK_DIR,
)

print(
    "Final output:",
    OUT_H5AD,
)

print(
    "Temporary directory:",
    TMP_DIR,
)

print(
    "Expected chunks:",
    EXPECTED_CHUNKS,
)

print(
    "Expected cells:",
    EXPECTED_CELLS,
)

print(
    "Expected genes:",
    EXPECTED_GENES,
)


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
    "\nFound chunks:",
    len(chunk_files),
)


# Important:
# Avoid chunk_files[0] crash when folder is empty.
if len(chunk_files) == 0:
    raise FileNotFoundError(
        "No H5AD chunk files were found in:\n"
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
# Step 1: Verify source chunks
# =============================================================================

print(
    "\n============================================================"
)

print(
    "Step 1: Verifying source chunks"
)

print(
    "============================================================"
)


# -----------------------------------------------------------------------------
# Use first chunk as reference for gene names/order
# -----------------------------------------------------------------------------

first = sc.read_h5ad(
    chunk_files[0],
    backed="r",
)

reference_genes = list(
    first.var_names.astype(str)
)

print(
    "Reference genes:",
    len(reference_genes),
)

print(
    "First chunk shape:",
    first.shape,
)

first.file.close()


if len(reference_genes) != EXPECTED_GENES:
    raise ValueError(
        f"First chunk expected "
        f"{EXPECTED_GENES} genes, "
        f"got {len(reference_genes)}"
    )


# -----------------------------------------------------------------------------
# Check all chunks
# -----------------------------------------------------------------------------

total_cells = 0

seen_cells = set()


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
    # Check identical gene order
    # -------------------------------------------------------------------------

    current_genes = list(
        adata.var_names.astype(str)
    )


    if current_genes != reference_genes:

        adata.file.close()

        raise ValueError(
            f"{chunk_file.name}: "
            "gene order mismatch"
        )


    # -------------------------------------------------------------------------
    # Check duplicate cell IDs
    # -------------------------------------------------------------------------

    cells = list(
        adata.obs_names.astype(str)
    )


    overlap = seen_cells.intersection(
        cells
    )


    if overlap:

        adata.file.close()

        raise ValueError(
            f"Duplicated cells found in "
            f"{chunk_file.name}: "
            f"{list(overlap)[:10]}"
        )


    seen_cells.update(
        cells
    )


    total_cells += adata.n_obs


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
            f"| cumulative cells: "
            f"{total_cells}"
        )


# -----------------------------------------------------------------------------
# Validate total cell count
# -----------------------------------------------------------------------------

if total_cells != EXPECTED_CELLS:

    raise ValueError(
        f"Expected "
        f"{EXPECTED_CELLS} total cells, "
        f"got {total_cells}"
    )


print(
    "\nPASS: All source chunks are valid."
)

print(
    "Total cells:",
    total_cells,
)

print(
    "Genes:",
    EXPECTED_GENES,
)


# =============================================================================
# Step 2: Prepare output directories
# =============================================================================

print(
    "\n============================================================"
)

print(
    "Step 2: Preparing merge directories"
)

print(
    "============================================================"
)


# Remove old final output
if OUT_H5AD.exists():

    print(
        "Removing old final output:",
        OUT_H5AD,
    )

    OUT_H5AD.unlink()


# Remove old temporary merge files
if TMP_DIR.exists():

    print(
        "Removing old temporary directory:",
        TMP_DIR,
    )

    shutil.rmtree(
        TMP_DIR
    )


TMP_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =============================================================================
# Step 3: Merge source chunks into batch H5AD files
# =============================================================================

print(
    "\n============================================================"
)

print(
    "Step 3: Creating batch H5AD files"
)

print(
    "============================================================"
)


batch_files = []


for batch_start in range(
    0,
    len(chunk_files),
    BATCH_SIZE,
):

    batch_id = (
        batch_start
        // BATCH_SIZE
    ) + 1


    batch_chunk_files = chunk_files[
        batch_start:
        batch_start + BATCH_SIZE
    ]


    batch_out = (
        TMP_DIR
        / f"huuki_batch_{batch_id:03d}.h5ad"
    )


    batch_tmp = (
        TMP_DIR
        / f"huuki_batch_{batch_id:03d}.tmp.h5ad"
    )


    # -------------------------------------------------------------------------
    # Clean existing batch outputs
    # -------------------------------------------------------------------------

    if batch_tmp.exists():
        batch_tmp.unlink()


    if batch_out.exists():
        batch_out.unlink()


    # -------------------------------------------------------------------------
    # Calculate expected cells in this batch
    # -------------------------------------------------------------------------

    expected_batch_cells = 0


    for chunk_file in batch_chunk_files:

        adata = sc.read_h5ad(
            chunk_file,
            backed="r",
        )

        expected_batch_cells += (
            adata.n_obs
        )

        adata.file.close()


    print(
        f"\nBatch {batch_id}"
    )

    print(
        "Chunks:",
        len(batch_chunk_files),
    )

    print(
        "Expected cells:",
        expected_batch_cells,
    )


    # -------------------------------------------------------------------------
    # Merge chunks on disk
    # -------------------------------------------------------------------------

    concat_on_disk(

        [
            str(f)
            for f in batch_chunk_files
        ],

        str(
            batch_tmp
        ),

        axis=0,

        join="inner",

        label="source_chunk",

        keys=[
            f"chunk_{chunk_number(f)}"
            for f in batch_chunk_files
        ],

        index_unique=None,

        max_loaded_elems=20_000_000,

    )


    # -------------------------------------------------------------------------
    # Verify temporary batch
    # -------------------------------------------------------------------------

    print(
        "Verifying batch:",
        batch_tmp.name,
    )


    verify_h5ad(

        batch_tmp,

        expected_cells=(
            expected_batch_cells
        ),

        expected_genes=(
            EXPECTED_GENES
        ),

    )


    # -------------------------------------------------------------------------
    # Promote temporary file to completed batch file
    # -------------------------------------------------------------------------

    batch_tmp.rename(
        batch_out
    )


    batch_files.append(
        batch_out
    )


print(
    "\nCreated batch files:",
    len(batch_files),
)


# =============================================================================
# Step 4: Merge batch files into final H5AD
# =============================================================================

print(
    "\n============================================================"
)

print(
    "Step 4: Merging batch H5AD files into final H5AD"
)

print(
    "============================================================"
)


FINAL_TMP = (
    SCE_DIR
    / "huuki_snrna_reference_full_allgenes.tmp.h5ad"
)


if FINAL_TMP.exists():

    FINAL_TMP.unlink()


concat_on_disk(

    [
        str(f)
        for f in batch_files
    ],

    str(
        FINAL_TMP
    ),

    axis=0,

    join="inner",

    label="merge_batch",

    keys=[
        f"batch_{i + 1:03d}"
        for i in range(
            len(batch_files)
        )
    ],

    index_unique=None,

    max_loaded_elems=20_000_000,

)


# =============================================================================
# Step 5: Verify final file
# =============================================================================

print(
    "\n============================================================"
)

print(
    "Step 5: Verifying final merged H5AD"
)

print(
    "============================================================"
)


verify_h5ad(

    FINAL_TMP,

    expected_cells=(
        EXPECTED_CELLS
    ),

    expected_genes=(
        EXPECTED_GENES
    ),

)


# =============================================================================
# Rename final temporary file
# =============================================================================

FINAL_TMP.rename(
    OUT_H5AD
)


# =============================================================================
# Success
# =============================================================================

print(
    "\n============================================================"
)

print(
    "MERGE SUCCESSFUL"
)

print(
    "============================================================"
)

print(
    "Final H5AD:",
    OUT_H5AD,
)

print(
    "Cells:",
    EXPECTED_CELLS,
)

print(
    "Genes:",
    EXPECTED_GENES,
)

print(
    "Source chunks:",
    EXPECTED_CHUNKS,
)

print(
    "Batch files:",
    len(batch_files),
)

print(
    "Temporary batch folder:",
    TMP_DIR,
)

print(
    "============================================================"
)