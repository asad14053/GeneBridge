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
  "spatialDLPFC_snRNAseq_top5000_hvg.h5ad"
)

old_wd <- getwd()
setwd(sce_dir)
on.exit(setwd(old_wd), add = TRUE)

cat("Reading snRNA-seq object from:\n")
cat(sce_dir, "\n")

sce <- readRDS("se.rds")

cat("\nOriginal object:\n")
print(sce)

cat("\nSelecting top 5000 genes by binomial_deviance...\n")

rd <- as.data.frame(rowData(sce))

if (!"binomial_deviance" %in% colnames(rd)) {
  stop("binomial_deviance column not found in rowData.")
}

gene_rank <- order(rd$binomial_deviance, decreasing = TRUE, na.last = NA)
top_n <- min(5000, length(gene_rank))
top_genes <- gene_rank[seq_len(top_n)]

sce_sub <- sce[top_genes, ]

cat("\nSubset object:\n")
print(sce_sub)

cat("\nMaking names unique...\n")
rownames(sce_sub) <- make.unique(rownames(sce_sub))
colnames(sce_sub) <- make.unique(colnames(sce_sub))

cat("\nTesting counts access:\n")
print(assay(sce_sub, "counts", withDimnames = FALSE)[1:5, 1:5])

cat("\nWriting h5ad:\n")
cat(out_h5ad, "\n")

zellkonverter::writeH5AD(
  sce_sub,
  file = out_h5ad,
  X_name = "logcounts"
)

cat("\nDONE. Saved:\n")
cat(out_h5ad, "\n")
