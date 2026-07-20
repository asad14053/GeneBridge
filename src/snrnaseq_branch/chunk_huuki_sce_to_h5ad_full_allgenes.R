#!/usr/bin/env Rscript

################################################################################
# Full Huuki-Myers snRNA SCE/RDS -> many small H5AD chunks
#
# This writes ALL genes but splits CELLS into small chunks.
# Do NOT use zellkonverter::writeH5AD on the full object at once.
#
# Usage:
#   # tiny test: 200 cells total, 50 cells per chunk
#   Rscript chunk_huuki_sce_to_h5ad_full_allgenes.R test 50 200
#
#   # full: all 77,604 cells, 200 cells per chunk
#   Rscript chunk_huuki_sce_to_h5ad_full_allgenes.R full 200 0
#
# Args:
#   1 MODE             test or full
#   2 CELLS_PER_CHUNK  e.g. 50, 100, 200
#   3 MAX_CELLS        test max cells; 0 means all cells
################################################################################

options(stringsAsFactors = FALSE)
options(width = 120)

args <- commandArgs(trailingOnly = TRUE)

MODE <- ifelse(length(args) >= 1, args[[1]], "full")
CELLS_PER_CHUNK <- ifelse(length(args) >= 2, as.integer(args[[2]]), 50)
MAX_CELLS <- ifelse(length(args) >= 3, as.integer(args[[3]]), 0)

if (!MODE %in% c("test", "full")) {
    stop("MODE must be either 'test' or 'full'")
}

PROJECT_ROOT <- "/users/mjabin/projects/GeneBridge"
SCE_DIR <- file.path(PROJECT_ROOT, "data/processed/snrnaseq/sce_DLPFC_annotated")

SE_RDS <- file.path(SCE_DIR, "se.rds")
ASSAYS_H5 <- file.path(SCE_DIR, "assays.h5")

OUT_DIR <- file.path(
    SCE_DIR,
    paste0("huuki_h5ad_chunks_allgenes_", MODE, "_maxcells", MAX_CELLS, "_chunk", CELLS_PER_CHUNK)
)

MANIFEST <- file.path(OUT_DIR, "chunk_manifest.csv")
LOG_FILE <- file.path(OUT_DIR, "chunk_conversion_log.txt")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sink(LOG_FILE, split = TRUE)

message("============================================================")
message("Full Huuki SCE/RDS -> H5AD chunks, all genes")
message("============================================================")
message("MODE:            ", MODE)
message("CELLS_PER_CHUNK: ", CELLS_PER_CHUNK)
message("MAX_CELLS:       ", MAX_CELLS)
message("SE_RDS:          ", SE_RDS)
message("ASSAYS_H5:       ", ASSAYS_H5)
message("OUT_DIR:         ", OUT_DIR)
message("MANIFEST:        ", MANIFEST)
message("LOG_FILE:        ", LOG_FILE)

if (!file.exists(SE_RDS)) stop("Missing se.rds: ", SE_RDS)
if (!file.exists(ASSAYS_H5)) stop("Missing assays.h5: ", ASSAYS_H5)

setwd(SCE_DIR)
message("Working directory: ", getwd())

required_pkgs <- c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "HDF5Array",
    "DelayedArray",
    "zellkonverter"
)

missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
    stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(DelayedArray)
    library(zellkonverter)
})

# Smaller blocks reduce memory spikes during delayed operations.
DelayedArray::setAutoBlockSize(50e6)

message("\nReading SCE...")
sce <- readRDS(SE_RDS)

message("\nSCE summary:")
print(sce)

message("\nOriginal dimensions genes x cells:")
print(dim(sce))

message("\nAssay names:")
print(assayNames(sce))

assay_names <- assayNames(sce)
assay_to_use <- if ("counts" %in% assay_names) {
    "counts"
} else if ("logcounts" %in% assay_names) {
    "logcounts"
} else {
    assay_names[1]
}

message("Using assay for AnnData X: ", assay_to_use)

# Use gene symbols as var_names if possible.
rowdata <- as.data.frame(rowData(sce))
symbol_candidates <- c("gene_name", "gene_symbol", "symbol", "gene", "external_gene_name")
symbol_col <- intersect(symbol_candidates, colnames(rowdata))[1]

if (!is.na(symbol_col)) {
    symbols <- as.character(rowdata[[symbol_col]])
    valid <- !is.na(symbols) & symbols != ""
    message("Gene symbol column: ", symbol_col)
    message("Fraction non-empty symbols: ", round(mean(valid), 3))
    if (mean(valid) > 0.8) {
        rownames(sce) <- make.unique(symbols)
        message("Using gene symbols from rowData column: ", symbol_col)
    }
}

if (is.null(rownames(sce)) || any(is.na(rownames(sce))) || any(rownames(sce) == "")) {
    rownames(sce) <- paste0("gene_", seq_len(nrow(sce)))
}
if (is.null(colnames(sce)) || any(is.na(colnames(sce))) || any(colnames(sce) == "")) {
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
}

rownames(sce) <- make.unique(as.character(rownames(sce)))
colnames(sce) <- make.unique(as.character(colnames(sce)))

all_cells <- colnames(sce)

if (MAX_CELLS > 0) {
    set.seed(42)
    n_keep <- min(MAX_CELLS, length(all_cells))
    selected_cells <- sample(all_cells, size = n_keep, replace = FALSE)
} else {
    selected_cells <- all_cells
}

message("Selected cells: ", length(selected_cells))
message("All genes kept: ", nrow(sce))

n_chunks <- ceiling(length(selected_cells) / CELLS_PER_CHUNK)
message("Number of chunks: ", n_chunks)

manifest <- data.frame(
    chunk_id = integer(),
    file = character(),
    n_genes = integer(),
    n_cells = integer(),
    stringsAsFactors = FALSE
)

for (i in seq_len(n_chunks)) {
    start_idx <- ((i - 1) * CELLS_PER_CHUNK) + 1
    end_idx <- min(i * CELLS_PER_CHUNK, length(selected_cells))
    cells_i <- selected_cells[start_idx:end_idx]

    chunk_file <- file.path(OUT_DIR, sprintf("huuki_allgenes_chunk_%04d.h5ad", i))

    message("\n------------------------------------------------------------")
    message("Writing chunk ", i, "/", n_chunks)
    message("Cells: ", length(cells_i))
    message("Genes: ", nrow(sce))
    message("File: ", chunk_file)
    message("------------------------------------------------------------")

    sce_i <- sce[, cells_i]

    if (file.exists(chunk_file)) {
        file.remove(chunk_file)
    }

    gc()

    zellkonverter::writeH5AD(
        sce_i,
        file = chunk_file,
        X_name = assay_to_use
    )

    if (!file.exists(chunk_file)) {
        stop("Chunk write failed: ", chunk_file)
    }

    manifest <- rbind(
        manifest,
        data.frame(
            chunk_id = i,
            file = chunk_file,
            n_genes = nrow(sce),
            n_cells = length(cells_i),
            stringsAsFactors = FALSE
        )
    )

    write.csv(manifest, MANIFEST, row.names = FALSE)

    rm(sce_i)
    gc()
}

message("\nDONE writing all chunks.")
message("Manifest: ", MANIFEST)
message("Output dir: ", OUT_DIR)

sink()
