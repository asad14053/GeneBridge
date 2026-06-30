#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(zellkonverter)
})

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

project_root <- "/users/mjabin/projects/GeneBridge"

sce_dir <- file.path(
  project_root,
  "data",
  "processed",
  "snrnaseq",
  "sce_DLPFC_annotated"
)

out_h5ad <- file.path(
  project_root,
  "data",
  "processed",
  "snrnaseq",
  "spatialDLPFC_snRNAseq_sce.h5ad"
)

old_wd <- getwd()
setwd(sce_dir)
on.exit(setwd(old_wd), add = TRUE)

cat("Reading snRNA-seq object from:\n")
cat(sce_dir, "\n")

sce <- readRDS("se.rds")

cat("\nLoaded object:\n")
print(sce)

cat("\nAssays:\n")
print(assayNames(sce))

cat("\nTesting counts access:\n")
print(assay(sce, "counts", withDimnames = FALSE)[1:5, 1:5])

cat("\nMaking names unique...\n")
colnames(sce) <- make.unique(colnames(sce))
rownames(sce) <- make.unique(rownames(sce))

cat("\nWriting h5ad:\n")
cat(out_h5ad, "\n")

zellkonverter::writeH5AD(
  sce,
  file = out_h5ad,
  X_name = "logcounts"
)

cat("\nDONE. Saved:\n")
cat(out_h5ad, "\n")
