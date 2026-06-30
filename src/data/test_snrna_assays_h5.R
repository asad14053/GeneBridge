#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(HDF5Array)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

project_root <- "/users/mjabin/projects/GeneBridge"
sce_dir <- file.path(project_root, "data", "processed", "snrnaseq", "sce_DLPFC_annotated")

assays_h5 <- file.path(sce_dir, "assays.h5")
se_rds <- file.path(sce_dir, "se.rds")

cat("assays.h5:\n")
cat(assays_h5, "\n\n")

cat("File exists:\n")
print(file.exists(assays_h5))

cat("\nFile info:\n")
print(file.info(assays_h5))

cat("\nTesting h5ls:\n")
print(rhdf5::h5ls(assays_h5))

cat("\nReading se.rds after setwd:\n")
old_wd <- getwd()
setwd(sce_dir)
on.exit(setwd(old_wd), add = TRUE)

sce <- readRDS("se.rds")

print(sce)
print(assayNames(sce))

cat("\nTesting assay access:\n")
print(assay(sce, "counts", withDimnames = FALSE)[1:5, 1:5])

cat("\nDONE\n")
