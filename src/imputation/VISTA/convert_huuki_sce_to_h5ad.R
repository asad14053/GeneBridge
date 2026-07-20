
# Convert Huuki-Myers SingleCellExperiment to H5AD

args <- commandArgs(trailingOnly = TRUE)
se_rds <- args[[1]]
out_h5ad <- args[[2]]

message("Input se.rds: ", se_rds)
message("Output h5ad: ", out_h5ad)

# Use the SCE directory as working directory so relative HDF5Array links to assays.h5 resolve.
setwd(dirname(se_rds))
message("Working directory: ", getwd())

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

pkgs <- c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "HDF5Array",
    "DelayedArray",
    "zellkonverter"
)

for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(HDF5Array)
    library(DelayedArray)
    library(zellkonverter)
})

message("Reading SCE...")
sce <- readRDS(se_rds)

message("SCE object:")
print(sce)

message("Dimensions genes x cells:")
print(dim(sce))

message("Assay names:")
print(assayNames(sce))

message("colData columns:")
print(colnames(colData(sce)))

message("rowData columns:")
print(colnames(rowData(sce)))

# Pick assay for AnnData.X
assay_names <- assayNames(sce)
preferred_assays <- c("counts", "logcounts", "X")
assay_to_use <- intersect(preferred_assays, assay_names)[1]

if (is.na(assay_to_use)) {
    assay_to_use <- assay_names[1]
}

message("Using assay for AnnData X: ", assay_to_use)

# Prefer gene symbols as rownames if available, because Xenium usually uses gene symbols.
rowdata <- as.data.frame(rowData(sce))
symbol_candidates <- c("gene_name", "gene_symbol", "symbol", "gene", "external_gene_name")
symbol_col <- intersect(symbol_candidates, colnames(rowdata))[1]

if (!is.na(symbol_col)) {
    symbols <- as.character(rowdata[[symbol_col]])
    valid <- !is.na(symbols) & symbols != "" & !duplicated(symbols)
    frac_valid <- mean(valid)

    message("Candidate gene symbol column: ", symbol_col)
    message("Fraction valid unique symbols: ", round(frac_valid, 3))

    if (frac_valid > 0.8) {
        rownames(sce) <- make.unique(symbols)
        message("Set rownames(sce) to unique gene symbols from: ", symbol_col)
    } else {
        message("Not using symbol column because too many missing/duplicated values.")
    }
} else {
    message("No obvious gene symbol column found; keeping existing rownames(sce).")
}

# Make sure cell names are unique and non-empty.
if (is.null(colnames(sce)) || any(is.na(colnames(sce))) || any(colnames(sce) == "")) {
    colnames(sce) <- paste0("cell_", seq_len(ncol(sce)))
}
colnames(sce) <- make.unique(as.character(colnames(sce)))

# Make sure gene names are unique and non-empty.
if (is.null(rownames(sce)) || any(is.na(rownames(sce))) || any(rownames(sce) == "")) {
    rownames(sce) <- paste0("gene_", seq_len(nrow(sce)))
}
rownames(sce) <- make.unique(as.character(rownames(sce)))

message("Final dimensions genes x cells:")
print(dim(sce))

message("First genes:")
print(head(rownames(sce)))

message("First cells:")
print(head(colnames(sce)))

message("Writing H5AD...")
zellkonverter::writeH5AD(
    sce,
    file = out_h5ad,
    X_name = assay_to_use
)

message("Done writing: ", out_h5ad)
