#!/usr/bin/env Rscript

####################################################################################################
# download_huuki_myers_direct_h5ad.R
#
# Task-01:
#   Download Huuki-Myers / spatialDLPFC data from spatialLIBD.
#
# Important fixes:
#   1. Do NOT save full Visium SpatialExperiment as RDS.
#   2. Do NOT save Visium HDF5-SE here; it can fail with "error writing to connection".
#   3. Save Visium only as h5ad for Python.
#   4. Do NOT use R unzip() for snRNA-seq zip; use Linux unzip.
#   5. Load snRNA-seq using readRDS("se.rds") from inside sce_DLPFC_annotated/.
#
# Main outputs:
#   Huuki Visium:
#     data/processed/visium/spatialDLPFC_Visium_sce.h5ad
#
#   Huuki snRNA-seq:
#     data/processed/snrnaseq/sce_DLPFC_annotated/se.rds
#     data/processed/snrnaseq/sce_DLPFC_annotated/assays.h5
#
#   Journal summary:
#     outputs/huuki_myers/tables/huuki_download_summary.txt
####################################################################################################

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

project_root <- "/users/mjabin/projects/GeneBridge"

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

pkgs <- c(
    "spatialLIBD",
    "SpatialExperiment",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "zellkonverter",
    "HDF5Array",
    "rhdf5"
)

for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
}

suppressPackageStartupMessages({
    library(spatialLIBD)
    library(SpatialExperiment)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(zellkonverter)
    library(HDF5Array)
    library(rhdf5)
})

visium_dir <- file.path(project_root, "data", "processed", "visium")
snrna_dir  <- file.path(project_root, "data", "processed", "snrnaseq")
out_table_dir <- file.path(project_root, "outputs", "huuki_myers", "tables")

dir.create(visium_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(snrna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

section <- function(title) {
    cat("\n", paste(rep("=", 100), collapse = ""), "\n", sep = "")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n", sep = "")
}

####################################################################################################
# 1. Download Visium SpatialExperiment
####################################################################################################

section("1. Downloading Huuki-Myers Visium SpatialExperiment")

spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium")

cat("\nVisium object:\n")
print(spe)

cat("\nAssays:\n")
print(assayNames(spe))

cat("\nDimensions:\n")
print(dim(spe))

cat("\nFirst 50 colData columns:\n")
print(colnames(colData(spe))[1:min(50, ncol(colData(spe)))])

cat("\nrowData columns:\n")
print(colnames(rowData(spe)))

cat("\nspatialCoords columns:\n")
print(colnames(spatialCoords(spe)))

####################################################################################################
# 2. Prepare Visium for h5ad export
####################################################################################################

section("2. Preparing Visium object for h5ad export")

coords <- spatialCoords(spe)

if (!is.null(coords) && ncol(coords) >= 2) {
    reducedDim(spe, "spatial") <- as.matrix(coords[, 1:2])
    cat("Added spatial coordinates to reducedDim(spe, 'spatial').\n")
} else {
    warning("No valid spatialCoords found.")
}

sce_visium <- as(spe, "SingleCellExperiment")

# Remove heavy metadata. This helps avoid serialization/write problems.
metadata(sce_visium) <- list()

# Keep only useful reducedDims for now.
if ("spatial" %in% reducedDimNames(sce_visium)) {
    spatial_mat <- reducedDim(sce_visium, "spatial")
    reducedDims(sce_visium) <- SimpleList(spatial = spatial_mat)
} else {
    reducedDims(sce_visium) <- SimpleList()
}

rownames(sce_visium) <- make.unique(rownames(sce_visium))
colnames(sce_visium) <- make.unique(colnames(sce_visium))

cat("\nPrepared Visium SingleCellExperiment:\n")
print(sce_visium)

####################################################################################################
# 3. Save Visium as h5ad for Python
####################################################################################################

section("3. Saving Huuki Visium as h5ad")

visium_h5ad <- file.path(visium_dir, "spatialDLPFC_Visium_sce.h5ad")

if (file.exists(visium_h5ad)) {
    file.remove(visium_h5ad)
}

zellkonverter::writeH5AD(
    sce_visium,
    file = visium_h5ad,
    X_name = "logcounts"
)

cat("Saved Visium h5ad to:\n")
cat(visium_h5ad, "\n")

####################################################################################################
# 4. Download snRNA-seq zip
####################################################################################################

section("4. Downloading Huuki-Myers snRNA-seq zip")

sce_path_zip <- spatialLIBD::fetch_data(type = "spatialDLPFC_snRNAseq")

cat("snRNA-seq zip path:\n")
cat(sce_path_zip, "\n")

if (!file.exists(sce_path_zip)) {
    stop("snRNA-seq zip file does not exist: ", sce_path_zip)
}

cat("\nZip file size in GB:\n")
print(file.info(sce_path_zip)$size / 1e9)

####################################################################################################
# 5. Extract snRNA-seq using Linux unzip
####################################################################################################

section("5. Extracting snRNA-seq using Linux unzip")

sce_extract_dir <- file.path(snrna_dir, "sce_DLPFC_annotated")
macosx_dir <- file.path(snrna_dir, "__MACOSX")

if (dir.exists(sce_extract_dir)) {
    cat("Removing old extracted snRNA-seq directory:\n")
    cat(sce_extract_dir, "\n")
    unlink(sce_extract_dir, recursive = TRUE)
}

if (dir.exists(macosx_dir)) {
    unlink(macosx_dir, recursive = TRUE)
}

unzip_log_path <- file.path(snrna_dir, "unzip_snrna.log")

unzip_log <- system2(
    command = "unzip",
    args = c("-o", sce_path_zip, "-d", snrna_dir),
    stdout = TRUE,
    stderr = TRUE
)

writeLines(unzip_log, unzip_log_path)

unzip_status <- attr(unzip_log, "status")
if (is.null(unzip_status)) {
    unzip_status <- 0
}

if (unzip_status != 0) {
    stop("Linux unzip failed. Check log: ", unzip_log_path)
}

cat("Unzip completed.\n")
cat("Unzip log saved to:\n")
cat(unzip_log_path, "\n")

####################################################################################################
# 6. Validate extracted snRNA-seq files
####################################################################################################

section("6. Validating extracted snRNA-seq files")

se_rds <- file.path(sce_extract_dir, "se.rds")
assays_h5 <- file.path(sce_extract_dir, "assays.h5")

cat("Expected snRNA-seq files:\n")
cat(se_rds, "\n")
cat(assays_h5, "\n")

if (!file.exists(se_rds)) {
    stop("Missing se.rds: ", se_rds)
}

if (!file.exists(assays_h5)) {
    stop("Missing assays.h5: ", assays_h5)
}

assays_h5_size <- file.info(assays_h5)$size

cat("\nassays.h5 size in GB:\n")
print(assays_h5_size / 1e9)

if (is.na(assays_h5_size) || assays_h5_size < 3e9) {
    stop(
        "assays.h5 looks truncated. Expected about 3.9 GB, but got ",
        assays_h5_size / 1e9,
        " GB."
    )
}

cat("\nTesting assays.h5 with rhdf5::h5ls():\n")
print(rhdf5::h5ls(assays_h5, recursive = FALSE))

####################################################################################################
# 7. Load snRNA-seq object correctly
####################################################################################################

section("7. Loading snRNA-seq object using readRDS from inside extracted directory")

old_wd <- getwd()
setwd(sce_extract_dir)

sce_snrna <- readRDS("se.rds")

cat("\nsnRNA-seq object:\n")
print(sce_snrna)

cat("\nAssays:\n")
print(assayNames(sce_snrna))

cat("\nDimensions:\n")
print(dim(sce_snrna))

cat("\nFirst 50 colData columns:\n")
print(colnames(colData(sce_snrna))[1:min(50, ncol(colData(sce_snrna)))])

cat("\nrowData columns:\n")
print(colnames(rowData(sce_snrna)))

cat("\nTesting counts assay access:\n")
print(assay(sce_snrna, "counts", withDimnames = FALSE)[1:5, 1:5])

cat("\nTesting logcounts assay access:\n")
print(assay(sce_snrna, "logcounts", withDimnames = FALSE)[1:5, 1:5])

setwd(old_wd)

####################################################################################################
# 8. Save download summary
####################################################################################################

section("8. Saving download summary")

summary_path <- file.path(out_table_dir, "huuki_download_summary.txt")

summary_lines <- c(
    "Huuki-Myers / spatialDLPFC download summary",
    "",
    "Project root:",
    project_root,
    "",
    "Huuki Visium output:",
    visium_h5ad,
    "",
    "Huuki snRNA-seq outputs:",
    se_rds,
    assays_h5,
    "",
    paste0("snRNA-seq assays.h5 size GB: ", round(assays_h5_size / 1e9, 3)),
    "",
    "Important notes:",
    "1. Visium SpatialExperiment was converted to SingleCellExperiment before h5ad writing.",
    "2. Visium spatial coordinates were stored in reducedDim(spe, 'spatial').",
    "3. Visium h5ad uses X_name = 'logcounts'. Counts are stored as a layer.",
    "4. Visium HDF5-SE was intentionally skipped because it failed with error writing to connection.",
    "5. snRNA-seq zip was extracted using Linux unzip, not R unzip().",
    "6. snRNA-seq should be loaded with readRDS('se.rds') from inside sce_DLPFC_annotated/.",
    "7. Full snRNA-seq h5ad conversion was not done here because it may exceed memory.",
    "",
    "Visium object:",
    capture.output(print(sce_visium)),
    "",
    "snRNA-seq object:",
    capture.output(print(sce_snrna))
)

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
writeLines(summary_lines, summary_path)

cat("Saved summary to:\n")
cat(summary_path, "\n")

section("DONE")

cat("\nMain outputs:\n")

cat("\nHuuki Visium h5ad:\n")
cat(visium_h5ad, "\n")

cat("\nHuuki snRNA-seq se.rds:\n")
cat(se_rds, "\n")

cat("\nHuuki snRNA-seq assays.h5:\n")
cat(assays_h5, "\n")

cat("\nDownload summary:\n")
cat(summary_path, "\n")