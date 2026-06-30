#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(DelayedArray)
  library(zellkonverter)
})

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

####################################################################################################
# Command-line argument parser
####################################################################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  # Supports both:
  #   --top_genes 5000
  #   --top_genes=5000

  eq_pattern <- paste0("^", flag, "=")
  eq_hit <- grep(eq_pattern, args, value = TRUE)

  if (length(eq_hit) > 0) {
    return(sub(eq_pattern, "", eq_hit[1]))
  }

  idx <- which(args == flag)

  if (length(idx) > 0 && idx[1] < length(args)) {
    return(args[idx[1] + 1])
  }

  return(default)
}

top_genes_arg <- get_arg("--top_genes", "5000")
brain_id_arg  <- get_arg("--brain_id", "Br8667")

brain_id_clean <- gsub("^Br", "", brain_id_arg)

if (tolower(top_genes_arg) %in% c("full", "all", "none", "inf", "infinite")) {
  TOP_N_GENES <- Inf
  gene_mode <- "full"
} else {
  TOP_N_GENES <- suppressWarnings(as.integer(top_genes_arg))

  if (is.na(TOP_N_GENES) || TOP_N_GENES <= 0) {
    stop("Invalid --top_genes value: ", top_genes_arg, ". Use 5000, 10000, or full.")
  }

  gene_mode <- paste0("top", TOP_N_GENES)
}

cat("\nConfiguration:\n")
cat("Brain ID:   ", paste0("Br", brain_id_clean), "\n")
cat("Gene mode:  ", gene_mode, "\n")

####################################################################################################
# Paths
####################################################################################################

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
  "data",
  "processed",
  "imputation_beta",
  paste0("Br", brain_id_clean)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_h5ad <- file.path(
  out_dir,
  paste0("seq_data_huuki_snrna_Br", brain_id_clean, "_", gene_mode, ".h5ad")
)

out_gene_list <- file.path(
  out_dir,
  paste0("seq_data_huuki_snrna_Br", brain_id_clean, "_", gene_mode, "_genes.csv")
)

####################################################################################################
# Load snRNA-seq
####################################################################################################

old_wd <- getwd()
setwd(sce_dir)
on.exit(setwd(old_wd), add = TRUE)

cat("\nReading Huuki snRNA-seq object from:\n")
cat(sce_dir, "\n")

sce <- readRDS("se.rds")

cat("\nOriginal object:\n")
print(sce)

####################################################################################################
# Subset cells by brain ID
####################################################################################################

if (!"BrNum" %in% colnames(colData(sce))) {
  stop("BrNum column not found in snRNA-seq colData.")
}

br <- as.character(colData(sce)$BrNum)
br_clean <- gsub("^Br", "", br)

keep_cells <- br_clean == brain_id_clean

cat("\nCells from Br", brain_id_clean, ":\n", sep = "")
print(sum(keep_cells))

if (sum(keep_cells) == 0) {
  stop("No cells found for Br", brain_id_clean, " in Huuki snRNA-seq.")
}

sce_b <- sce[, keep_cells]

####################################################################################################
# Subset genes
####################################################################################################

rd <- as.data.frame(rowData(sce_b))

if (is.infinite(TOP_N_GENES)) {

  cat("\nKeeping FULL gene set.\n")
  cat("Genes kept:\n")
  print(nrow(sce_b))

} else {

  cat("\nSelecting top genes by binomial_deviance.\n")
  cat("Requested top genes:\n")
  print(TOP_N_GENES)

  if (!"binomial_deviance" %in% colnames(rd)) {
    stop("binomial_deviance column not found in rowData. Cannot select top genes.")
  }

  gene_rank <- order(rd$binomial_deviance, decreasing = TRUE, na.last = NA)
  top_n <- min(TOP_N_GENES, length(gene_rank))
  keep_genes <- gene_rank[seq_len(top_n)]

  sce_b <- sce_b[keep_genes, ]

  cat("Genes kept:\n")
  print(nrow(sce_b))
}

####################################################################################################
# Clean names
####################################################################################################

rownames(sce_b) <- make.unique(rownames(sce_b))
colnames(sce_b) <- make.unique(colnames(sce_b))

####################################################################################################
# Save selected gene list
####################################################################################################

gene_df <- data.frame(
  gene = rownames(sce_b),
  stringsAsFactors = FALSE
)

if ("gene_id" %in% colnames(rowData(sce_b))) {
  gene_df$gene_id <- rowData(sce_b)$gene_id
}

if ("gene_name" %in% colnames(rowData(sce_b))) {
  gene_df$gene_name <- rowData(sce_b)$gene_name
}

if ("binomial_deviance" %in% colnames(rowData(sce_b))) {
  gene_df$binomial_deviance <- rowData(sce_b)$binomial_deviance
}

write.csv(gene_df, out_gene_list, row.names = FALSE)

####################################################################################################
# Test assay access
####################################################################################################

cat("\nFinal Br", brain_id_clean, " snRNA subset:\n", sep = "")
print(sce_b)

cat("\nTesting counts access:\n")
print(assay(sce_b, "counts", withDimnames = FALSE)[1:5, 1:5])

####################################################################################################
# Write h5ad
####################################################################################################

cat("\nWriting h5ad:\n")
cat(out_h5ad, "\n")

if (file.exists(out_h5ad)) {
  file.remove(out_h5ad)
}

zellkonverter::writeH5AD(
  sce_b,
  file = out_h5ad,
  X_name = "logcounts"
)

cat("\nDONE. Saved h5ad:\n")
cat(out_h5ad, "\n")

cat("\nSaved gene list:\n")
cat(out_gene_list, "\n")