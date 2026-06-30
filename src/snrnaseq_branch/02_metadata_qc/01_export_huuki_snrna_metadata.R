#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(DelayedArray)
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

out_dir <- file.path(
  project_root,
  "outputs",
  "huuki_myers",
  "tables"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old_wd <- getwd()
setwd(sce_dir)
on.exit(setwd(old_wd), add = TRUE)

cat("Reading snRNA-seq object from:\n")
cat(sce_dir, "\n")

sce <- readRDS("se.rds")

cat("\nsnRNA-seq object:\n")
print(sce)

meta <- as.data.frame(colData(sce))
meta$cell_id <- colnames(sce)

if ("UMAP" %in% reducedDimNames(sce)) {
  umap <- as.data.frame(reducedDim(sce, "UMAP"))
  colnames(umap) <- c("UMAP1", "UMAP2")
  meta$UMAP1 <- umap$UMAP1
  meta$UMAP2 <- umap$UMAP2
}

out_meta <- file.path(out_dir, "huuki_snrna_metadata.csv")
write.csv(meta, out_meta, row.names = FALSE)

br_counts <- as.data.frame(table(meta$BrNum))
colnames(br_counts) <- c("BrNum", "snrna_cells")

out_counts <- file.path(out_dir, "huuki_snrna_brain_counts.csv")
write.csv(br_counts, out_counts, row.names = FALSE)

cat("\nSaved:\n")
cat(out_meta, "\n")
cat(out_counts, "\n")

cat("\nUnique snRNA BrNum:\n")
print(sort(unique(as.character(meta$BrNum))))