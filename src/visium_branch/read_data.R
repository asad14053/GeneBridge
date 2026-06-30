#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript export_celltype_rds_to_csv.R input.rds output.csv")
}

input_rds <- args[[1]]
output_csv <- args[[2]]

cat("Reading RDS:\n")
cat(input_rds, "\n")

res <- readRDS(input_rds)

cat("RDS class:\n")
print(class(res))

# Expected nmfLabelTransfer output:
# res$targets is a list of target SPE/SCE objects.
if (!("targets" %in% names(res))) {
  stop("This RDS does not contain res$targets. Please inspect object structure.")
}

all_preds <- list()

for (i in seq_along(res$targets)) {
  spe <- res$targets[[i]]

  cell_ids <- colnames(spe)

  if (!("nmf_preds" %in% colnames(colData(spe)))) {
    stop(paste0("nmf_preds not found in target ", i))
  }

  pred <- as.character(colData(spe)$nmf_preds)

  brnum <- NA
  if ("BrNum" %in% colnames(colData(spe))) {
    brnum <- as.character(colData(spe)$BrNum)
  } else {
    brnum <- sub("_.*", "", cell_ids)
  }

  dx <- NA
  if ("Dx" %in% colnames(colData(spe))) {
    dx <- as.character(colData(spe)$Dx)
  }

  df <- data.frame(
    cell_id = cell_ids,
    BrNum = brnum,
    Dx = dx,
    cell_type_prediction = pred,
    stringsAsFactors = FALSE
  )

  all_preds[[i]] <- df
}

out <- do.call(rbind, all_preds)

cat("Output rows:\n")
print(nrow(out))

cat("Head:\n")
print(head(out))

write.csv(out, output_csv, row.names = FALSE)

cat("Saved CSV:\n")
cat(output_csv, "\n")