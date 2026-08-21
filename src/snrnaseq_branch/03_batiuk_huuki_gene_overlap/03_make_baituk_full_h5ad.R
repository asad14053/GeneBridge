#!/usr/bin/env Rscript

################################################################################
# Create one complete Baituk snRNA-seq H5AD from the sample-level RDS matrices.
#
# Input:
#   data/raw/Baituk/snRNA-seq_raw_countmatrices.RDS
#
# Output:
#   data/processed/snrnaseq/Baituk/
#   baituk_snrna_reference_full_allgenes.h5ad
#
# AnnData:
#   X = raw counts
#   obs["sample"] = Baituk sample ID
#   obs["original_cell_id"] = original barcode
#   var["gene_name"] = Baituk gene symbol
#
# Uses anndataR, not zellkonverter.
################################################################################

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(anndataR)
})

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
    PROJECT_ROOT,
    "data",
    "raw",
    "Baituk",
    "snRNA-seq_raw_countmatrices.RDS"
)

OUTPUT_H5AD <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "snrnaseq",
    "Baituk",
    "baituk_snrna_reference_full_allgenes.h5ad"
)

SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "03_batiuk_huuki_gene_overlap",
    "baituk_full_h5ad_summary.csv"
)

SAMPLE_SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "03_batiuk_huuki_gene_overlap",
    "baituk_full_h5ad_sample_summary.csv"
)


section <- function(title) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


extract_count_matrix <- function(object, sample_name) {

    if (
        inherits(object, "Matrix") ||
        is.matrix(object)
    ) {
        matrix_object <- object

    } else if (
        inherits(object, "SummarizedExperiment")
    ) {
        assay_names <- assayNames(object)

        if ("counts" %in% assay_names) {
            matrix_object <- assay(
                object,
                "counts",
                withDimnames = TRUE
            )
        } else {
            matrix_object <- assay(
                object,
                1,
                withDimnames = TRUE
            )
        }

    } else if (is.list(object)) {

        matrix_candidates <- which(
            vapply(
                object,
                function(x) {
                    inherits(x, "Matrix") || is.matrix(x)
                },
                logical(1)
            )
        )

        if (length(matrix_candidates) != 1) {
            stop(
                "Could not identify exactly one matrix inside sample ",
                sample_name
            )
        }

        matrix_object <- object[[matrix_candidates[1]]]

    } else {
        stop(
            "Unsupported Baituk object for sample ",
            sample_name,
            ". Class: ",
            paste(class(object), collapse = ";")
        )
    }

    if (is.null(dim(matrix_object))) {
        stop(
            "The object for sample ",
            sample_name,
            " does not have matrix dimensions."
        )
    }

    if (is.null(rownames(matrix_object))) {
        stop(
            "The matrix for sample ",
            sample_name,
            " has no gene row names."
        )
    }

    if (anyDuplicated(rownames(matrix_object)) > 0) {
        stop(
            "The matrix for sample ",
            sample_name,
            " contains duplicated gene names."
        )
    }

    if (!inherits(matrix_object, "dgCMatrix")) {
        matrix_object <- as(
            matrix_object,
            "dgCMatrix"
        )
    }

    return(matrix_object)
}


################################################################################
# Read RDS
################################################################################

section("Create full Baituk snRNA-seq H5AD")

cat("Input RDS:\n")
cat(INPUT_RDS, "\n")

cat("\nOutput H5AD:\n")
cat(OUTPUT_H5AD, "\n")

if (!file.exists(INPUT_RDS)) {
    stop(
        "Input RDS not found:\n",
        INPUT_RDS
    )
}

section("Reading Baituk RDS")

baituk_objects <- readRDS(INPUT_RDS)

if (!is.list(baituk_objects)) {
    stop(
        "Expected the Baituk RDS to contain a list of sample matrices."
    )
}

sample_names <- names(baituk_objects)

if (is.null(sample_names)) {
    sample_names <- paste0(
        "sample_",
        seq_along(baituk_objects)
    )
}

cat("Number of Baituk objects:", length(baituk_objects), "\n")
cat("Sample names:\n")
print(sample_names)

################################################################################
# Extract and validate matrices
################################################################################

section("Extracting sample matrices")

matrices <- list()
obs_tables <- list()
sample_records <- list()

reference_genes <- NULL

for (i in seq_along(baituk_objects)) {

    sample_name <- sample_names[i]

    cat(
        "\nProcessing sample ",
        sample_name,
        "...\n",
        sep = ""
    )

    matrix_object <- extract_count_matrix(
        baituk_objects[[i]],
        sample_name
    )

    current_genes <- as.character(
        rownames(matrix_object)
    )

    if (is.null(reference_genes)) {

        reference_genes <- current_genes

    } else if (!identical(
        current_genes,
        reference_genes
    )) {

        if (!setequal(
            current_genes,
            reference_genes
        )) {
            stop(
                "Sample ",
                sample_name,
                " does not contain the same gene set as the reference sample."
            )
        }

        cat(
            "Gene set matches, but order differs. Reordering sample.\n"
        )

        matrix_object <- matrix_object[
            match(
                reference_genes,
                current_genes
            ),
            ,
            drop = FALSE
        ]

        rownames(matrix_object) <- reference_genes
    }

    original_cell_ids <- colnames(matrix_object)

    if (is.null(original_cell_ids)) {
        original_cell_ids <- paste0(
            "cell_",
            seq_len(ncol(matrix_object))
        )
    }

    original_cell_ids <- as.character(
        original_cell_ids
    )

    global_cell_ids <- paste0(
        sample_name,
        "::",
        original_cell_ids
    )

    global_cell_ids <- make.unique(
        global_cell_ids
    )

    colnames(matrix_object) <- global_cell_ids

    matrices[[sample_name]] <- matrix_object

    obs_tables[[sample_name]] <- data.frame(
        sample = rep(
            sample_name,
            ncol(matrix_object)
        ),
        original_cell_id = original_cell_ids,
        stringsAsFactors = FALSE,
        row.names = global_cell_ids
    )

    sample_records[[sample_name]] <- data.frame(
        sample = sample_name,
        n_genes = nrow(matrix_object),
        n_cells = ncol(matrix_object),
        nonzero_values = length(matrix_object@x),
        stringsAsFactors = FALSE
    )

    cat(
        "Shape [genes x cells]: ",
        nrow(matrix_object),
        " x ",
        ncol(matrix_object),
        "\n",
        sep = ""
    )
}

################################################################################
# Merge matrices
################################################################################

section("Combining Baituk samples")

combined_counts <- do.call(
    cbind,
    unname(matrices)
)

if (!inherits(combined_counts, "dgCMatrix")) {
    combined_counts <- as(
        combined_counts,
        "dgCMatrix"
    )
}

combined_obs <- do.call(
    rbind,
    unname(obs_tables)
)

sample_summary <- do.call(
    rbind,
    sample_records
)

if (nrow(combined_obs) != ncol(combined_counts)) {
    stop(
        "Combined observation metadata does not match combined cell count."
    )
}

if (!identical(
    rownames(combined_obs),
    colnames(combined_counts)
)) {
    stop(
        "Combined observation row names do not match matrix column names."
    )
}

cat("Combined shape [genes x cells]:\n")
print(dim(combined_counts))

cat("\nSample summary:\n")
print(sample_summary)

################################################################################
# Build SingleCellExperiment
################################################################################

section("Building SingleCellExperiment")

gene_metadata <- data.frame(
    gene_name = reference_genes,
    stringsAsFactors = FALSE,
    row.names = reference_genes
)

sce <- SingleCellExperiment(
    assays = list(
        counts = combined_counts
    ),
    rowData = S4Vectors::DataFrame(
        gene_metadata
    ),
    colData = S4Vectors::DataFrame(
        combined_obs
    )
)

rownames(sce) <- reference_genes
colnames(sce) <- colnames(combined_counts)

print(sce)

################################################################################
# Write H5AD
################################################################################

dir.create(
    dirname(OUTPUT_H5AD),
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    dirname(SUMMARY_CSV),
    recursive = TRUE,
    showWarnings = FALSE
)

if (file.exists(OUTPUT_H5AD)) {
    cat("\nRemoving existing output:\n")
    cat(OUTPUT_H5AD, "\n")
    unlink(OUTPUT_H5AD)
}

section("Writing full Baituk H5AD")

anndataR::write_h5ad(
    object = sce,
    path = OUTPUT_H5AD,
    compression = "lzf",
    x_mapping = "counts",
    layers_mapping = FALSE,
    obs_mapping = TRUE,
    var_mapping = TRUE,
    obsm_mapping = FALSE,
    varm_mapping = FALSE,
    obsp_mapping = FALSE,
    varp_mapping = FALSE,
    uns_mapping = FALSE
)

if (!file.exists(OUTPUT_H5AD)) {
    stop(
        "Baituk H5AD was not created."
    )
}

################################################################################
# Save summaries
################################################################################

summary_df <- data.frame(
    metric = c(
        "input_rds",
        "output_h5ad",
        "n_samples",
        "n_cells",
        "n_genes",
        "X_content",
        "output_size_GB"
    ),
    value = c(
        INPUT_RDS,
        OUTPUT_H5AD,
        length(matrices),
        ncol(sce),
        nrow(sce),
        "raw_counts",
        round(
            file.info(OUTPUT_H5AD)$size / 1024^3,
            4
        )
    ),
    stringsAsFactors = FALSE
)

write.csv(
    summary_df,
    SUMMARY_CSV,
    row.names = FALSE
)

write.csv(
    sample_summary,
    SAMPLE_SUMMARY_CSV,
    row.names = FALSE
)

section("DONE")

print(summary_df)

cat("\nCreated Baituk H5AD:\n")
cat(OUTPUT_H5AD, "\n")

cat("\nSample summary:\n")
cat(SAMPLE_SUMMARY_CSV, "\n")
