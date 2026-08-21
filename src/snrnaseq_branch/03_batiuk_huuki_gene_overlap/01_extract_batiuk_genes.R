#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
  PROJECT_ROOT,
  "data",
  "raw",
  "Baituk",
  "snRNA-seq_raw_countmatrices.RDS"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "snrnaseq_branch",
  "03_batiuk_huuki_gene_overlap"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

SAMPLE_SUMMARY_CSV <- file.path(
  OUTPUT_DIR,
  "batiuk_sample_gene_summary.csv"
)

UNION_GENES_CSV <- file.path(
  OUTPUT_DIR,
  "batiuk_gene_union.csv"
)

INTERSECTION_GENES_CSV <- file.path(
  OUTPUT_DIR,
  "batiuk_gene_intersection.csv"
)

cat("Input Batiuk RDS:\n")
cat(INPUT_RDS, "\n")

if (!file.exists(INPUT_RDS)) {
  stop("Input Batiuk RDS not found: ", INPUT_RDS)
}

cat("\nReading Batiuk RDS...\n")

batiuk <- readRDS(INPUT_RDS)

if (!is.list(batiuk)) {
  stop("Expected the Batiuk RDS to contain a list.")
}

sample_names <- names(batiuk)

if (is.null(sample_names)) {
  sample_names <- paste0(
    "sample_",
    seq_along(batiuk)
  )
}

gene_sets <- list()
sample_records <- vector(
  "list",
  length(batiuk)
)

for (i in seq_along(batiuk)) {

  sample_name <- sample_names[i]
  object <- batiuk[[i]]

  dimensions <- tryCatch(
    dim(object),
    error = function(e) NULL
  )

  genes <- tryCatch(
    rownames(object),
    error = function(e) NULL
  )

  genes <- as.character(genes)
  genes <- trimws(genes)
  genes <- genes[
    !is.na(genes) &
    nzchar(genes)
  ]

  unique_genes <- unique(genes)

  n_rows <- if (
    !is.null(dimensions) &&
    length(dimensions) >= 1
  ) {
    dimensions[1]
  } else {
    NA_integer_
  }

  n_cells <- if (
    !is.null(dimensions) &&
    length(dimensions) >= 2
  ) {
    dimensions[2]
  } else {
    NA_integer_
  }

  status <- if (length(unique_genes) > 0) {
    "gene_names_found"
  } else {
    "no_gene_names"
  }

  if (length(unique_genes) > 0) {
    gene_sets[[sample_name]] <- unique_genes
  }

  sample_records[[i]] <- data.frame(
    sample = sample_name,
    object_class = paste(
      class(object),
      collapse = ";"
    ),
    n_rows = n_rows,
    n_cells = n_cells,
    n_gene_names = length(genes),
    n_unique_genes = length(unique_genes),
    duplicated_gene_names = length(genes) -
      length(unique_genes),
    status = status,
    stringsAsFactors = FALSE
  )
}

sample_summary <- do.call(
  rbind,
  sample_records
)

if (length(gene_sets) == 0) {
  stop("No gene names were found in the Batiuk RDS.")
}

reference_sample <- names(gene_sets)[1]
reference_genes <- gene_sets[[reference_sample]]

sample_summary$same_gene_set_as_reference <- vapply(
  sample_summary$sample,
  function(sample_name) {

    if (!(sample_name %in% names(gene_sets))) {
      return(FALSE)
    }

    setequal(
      gene_sets[[sample_name]],
      reference_genes
    )
  },
  logical(1)
)

sample_summary$same_gene_order_as_reference <- vapply(
  sample_summary$sample,
  function(sample_name) {

    if (!(sample_name %in% names(gene_sets))) {
      return(FALSE)
    }

    identical(
      gene_sets[[sample_name]],
      reference_genes
    )
  },
  logical(1)
)

batiuk_union <- Reduce(
  union,
  gene_sets
)

batiuk_intersection <- Reduce(
  intersect,
  gene_sets
)

write.csv(
  sample_summary,
  SAMPLE_SUMMARY_CSV,
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = sort(batiuk_union),
    stringsAsFactors = FALSE
  ),
  UNION_GENES_CSV,
  row.names = FALSE
)

write.csv(
  data.frame(
    gene = sort(batiuk_intersection),
    stringsAsFactors = FALSE
  ),
  INTERSECTION_GENES_CSV,
  row.names = FALSE
)

cat("\nBatiuk sample summary:\n")
print(sample_summary)

cat("\nNumber of Batiuk objects:\n")
print(length(batiuk))

cat("\nSamples containing gene names:\n")
print(length(gene_sets))

cat("\nReference sample:\n")
print(reference_sample)

cat("\nBatiuk union genes:\n")
print(length(batiuk_union))

cat("\nBatiuk intersection genes:\n")
print(length(batiuk_intersection))

cat("\nAll valid samples have the same gene set:\n")
print(
  all(
    sample_summary$same_gene_set_as_reference[
      sample_summary$status == "gene_names_found"
    ]
  )
)

cat("\nCreated:\n")
cat(SAMPLE_SUMMARY_CSV, "\n")
cat(UNION_GENES_CSV, "\n")
cat(INTERSECTION_GENES_CSV, "\n")

cat("\nDONE\n")
