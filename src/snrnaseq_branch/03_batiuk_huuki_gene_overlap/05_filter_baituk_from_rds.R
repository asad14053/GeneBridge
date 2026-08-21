#!/usr/bin/env Rscript

################################################################################
# Create the Baituk H5AD restricted to Baituk-Huuki shared genes.
#
# This does NOT recreate or modify the full Baituk H5AD.
#
# Input:
#   data/raw/Baituk/snRNA-seq_raw_countmatrices.RDS
#
# Shared genes:
#   outputs/snrnaseq_branch/03_batiuk_huuki_gene_overlap/
#   batiuk_huuki_shared_genes.csv
#
# Output:
#   data/processed/snrnaseq/Baituk/
#   baituk_snrna_reference_huuki_overlap.h5ad
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

SHARED_GENE_CANDIDATES <- c(
    file.path(
        PROJECT_ROOT,
        "outputs",
        "snrnaseq_branch",
        "03_batiuk_huuki_gene_overlap",
        "batiuk_huuki_shared_genes.csv"
    ),
    file.path(
        PROJECT_ROOT,
        "outputs",
        "snrnaseq_branch",
        "03_batiuk_huuki_gene_overlap",
        "baituk_huuki_shared_genes_final.csv"
    )
)

OUTPUT_H5AD <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "snrnaseq",
    "Baituk",
    "baituk_snrna_reference_huuki_overlap.h5ad"
)

SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "03_batiuk_huuki_gene_overlap",
    "baituk_filtered_h5ad_summary.csv"
)

SAMPLE_SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "03_batiuk_huuki_gene_overlap",
    "baituk_filtered_h5ad_sample_summary.csv"
)

EXPECTED_SHARED_GENES <- 36119L


section <- function(title) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


extract_count_matrix <- function(object, sample_name) {

    if (inherits(object, "Matrix") || is.matrix(object)) {

        matrix_object <- object

    } else if (inherits(object, "SummarizedExperiment")) {

        if ("counts" %in% assayNames(object)) {
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

        if (length(matrix_candidates) != 1L) {
            stop(
                "Could not identify exactly one matrix for sample ",
                sample_name,
                ". Candidate matrices: ",
                length(matrix_candidates)
            )
        }

        matrix_object <- object[[matrix_candidates[1]]]

    } else {

        stop(
            "Unsupported object for sample ",
            sample_name,
            ". Class: ",
            paste(class(object), collapse = ", ")
        )
    }

    if (is.null(dim(matrix_object))) {
        stop(
            "Sample ",
            sample_name,
            " does not have matrix dimensions."
        )
    }

    if (is.null(rownames(matrix_object))) {
        stop(
            "Sample ",
            sample_name,
            " has no gene row names."
        )
    }

    if (anyDuplicated(rownames(matrix_object)) > 0L) {
        stop(
            "Sample ",
            sample_name,
            " contains duplicated gene row names."
        )
    }

    if (!inherits(matrix_object, "dgCMatrix")) {
        matrix_object <- as(
            matrix_object,
            "dgCMatrix"
        )
    }

    matrix_object
}


################################################################################
# Validate files
################################################################################

section("INPUT VALIDATION")

cat("Input Baituk RDS:\n")
cat(INPUT_RDS, "\n")

cat("\nOutput filtered Baituk H5AD:\n")
cat(OUTPUT_H5AD, "\n")

if (!file.exists(INPUT_RDS)) {
    stop(
        "Baituk RDS was not found:\n",
        INPUT_RDS
    )
}

available_shared_files <- SHARED_GENE_CANDIDATES[
    file.exists(SHARED_GENE_CANDIDATES)
]

if (length(available_shared_files) == 0L) {
    stop(
        "No shared-gene CSV was found. Checked:\n",
        paste(
            SHARED_GENE_CANDIDATES,
            collapse = "\n"
        )
    )
}

SHARED_GENES_CSV <- available_shared_files[1]

cat("\nUsing shared-gene CSV:\n")
cat(SHARED_GENES_CSV, "\n")


################################################################################
# Read shared genes
################################################################################

section("READING SHARED GENES")

shared_table <- read.csv(
    SHARED_GENES_CSV,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (ncol(shared_table) == 0L) {
    stop("Shared-gene CSV contains no columns.")
}

preferred_columns <- c(
    "gene",
    "shared_gene",
    "gene_name",
    "genes"
)

matching_columns <- preferred_columns[
    preferred_columns %in% colnames(shared_table)
]

if (length(matching_columns) > 0L) {
    gene_column <- matching_columns[1]
} else {
    excluded_columns <- c(
        "X",
        "index",
        "gene_order",
        "baituk_var_index",
        "huuki_var_index"
    )

    available_columns <- setdiff(
        colnames(shared_table),
        excluded_columns
    )

    if (length(available_columns) == 0L) {
        gene_column <- colnames(shared_table)[1]
    } else {
        gene_column <- available_columns[1]
    }
}

shared_genes <- trimws(
    as.character(
        shared_table[[gene_column]]
    )
)

shared_genes <- shared_genes[
    !is.na(shared_genes) &
    shared_genes != "" &
    !tolower(shared_genes) %in% c(
        "nan",
        "none",
        "na"
    )
]

shared_genes <- unique(shared_genes)

cat("Gene column: ", gene_column, "\n", sep = "")
cat("Shared genes: ", length(shared_genes), "\n", sep = "")

if (length(shared_genes) != EXPECTED_SHARED_GENES) {
    stop(
        "Expected ",
        EXPECTED_SHARED_GENES,
        " shared genes, but found ",
        length(shared_genes),
        "."
    )
}


################################################################################
# Read Baituk RDS
################################################################################

section("READING BAITUK RDS")

baituk_objects <- readRDS(INPUT_RDS)

if (!is.list(baituk_objects)) {
    stop(
        "Expected the Baituk RDS to contain a list."
    )
}

sample_names <- names(baituk_objects)

if (is.null(sample_names)) {
    sample_names <- paste0(
        "sample_",
        seq_along(baituk_objects)
    )
}

cat("Number of Baituk samples: ", length(baituk_objects), "\n", sep = "")
print(sample_names)


################################################################################
# Filter each sample before combining
################################################################################

section("FILTERING INDIVIDUAL BAITUK SAMPLES")

filtered_matrices <- vector(
    mode = "list",
    length = length(baituk_objects)
)

obs_tables <- vector(
    mode = "list",
    length = length(baituk_objects)
)

sample_records <- vector(
    mode = "list",
    length = length(baituk_objects)
)

for (i in seq_along(baituk_objects)) {

    sample_name <- sample_names[i]

    cat(
        "\nProcessing ",
        sample_name,
        " [",
        i,
        "/",
        length(baituk_objects),
        "]\n",
        sep = ""
    )

    matrix_object <- extract_count_matrix(
        baituk_objects[[i]],
        sample_name
    )

    gene_indices <- match(
        shared_genes,
        rownames(matrix_object)
    )

    if (anyNA(gene_indices)) {

        missing_genes <- shared_genes[
            is.na(gene_indices)
        ]

        stop(
            "Sample ",
            sample_name,
            " is missing ",
            length(missing_genes),
            " required genes. Examples: ",
            paste(
                head(missing_genes, 20),
                collapse = ", "
            )
        )
    }

    # Important optimization:
    # filter every sample before combining the 24 matrices.
    matrix_object <- matrix_object[
        gene_indices,
        ,
        drop = FALSE
    ]

    rownames(matrix_object) <- shared_genes

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

    filtered_matrices[[i]] <- matrix_object

    obs_tables[[i]] <- data.frame(
        sample = rep(
            sample_name,
            ncol(matrix_object)
        ),
        original_cell_id = original_cell_ids,
        stringsAsFactors = FALSE,
        row.names = global_cell_ids
    )

    sample_records[[i]] <- data.frame(
        sample = sample_name,
        n_cells = ncol(matrix_object),
        n_genes = nrow(matrix_object),
        nonzero_values = length(matrix_object@x),
        stringsAsFactors = FALSE
    )

    cat(
        "Filtered shape [genes x cells]: ",
        nrow(matrix_object),
        " x ",
        ncol(matrix_object),
        "\n",
        sep = ""
    )
}

rm(baituk_objects)
invisible(gc())


################################################################################
# Combine filtered matrices
################################################################################

section("COMBINING FILTERED BAITUK MATRICES")

combined_counts <- do.call(
    cbind,
    filtered_matrices
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

if (nrow(combined_obs) != ncol(combined_counts)) {
    stop(
        "Observation metadata does not match the combined matrix.\n",
        "Metadata rows: ",
        nrow(combined_obs),
        "\nMatrix columns: ",
        ncol(combined_counts)
    )
}

# Force metadata row names to exactly match matrix column names.
rownames(combined_obs) <- colnames(combined_counts)

if (!identical(
    rownames(combined_obs),
    colnames(combined_counts)
)) {
    stop(
        "Observation row names still do not match matrix columns."
    )
}

cat("Combined shape [genes x cells]:\n")
print(dim(combined_counts))

rm(filtered_matrices, obs_tables)
invisible(gc())


################################################################################
# Build SingleCellExperiment
################################################################################

section("BUILDING FILTERED SINGLECELLEXPERIMENT")

gene_metadata <- data.frame(
    gene_name = shared_genes,
    alignment_gene = shared_genes,
    stringsAsFactors = FALSE,
    row.names = shared_genes
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

rownames(sce) <- shared_genes
colnames(sce) <- colnames(combined_counts)

print(sce)

if (nrow(sce) != EXPECTED_SHARED_GENES) {
    stop(
        "Filtered SCE has an unexpected number of genes: ",
        nrow(sce)
    )
}


################################################################################
# Write filtered H5AD
################################################################################

section("WRITING FILTERED BAITUK H5AD")

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
    cat("Removing partial existing output:\n")
    cat(OUTPUT_H5AD, "\n")
    unlink(OUTPUT_H5AD)
}

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
        "Filtered Baituk H5AD was not created."
    )
}


################################################################################
# Save summaries
################################################################################

summary_df <- data.frame(
    metric = c(
        "input_rds",
        "shared_genes_csv",
        "output_h5ad",
        "n_samples",
        "n_cells",
        "n_genes",
        "output_size_GB"
    ),
    value = c(
        INPUT_RDS,
        SHARED_GENES_CSV,
        OUTPUT_H5AD,
        length(sample_names),
        ncol(sce),
        nrow(sce),
        round(
            file.info(OUTPUT_H5AD)$size / 1024^3,
            4
        )
    ),
    stringsAsFactors = FALSE
)

sample_summary <- do.call(
    rbind,
    sample_records
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

cat("\nCreated filtered Baituk H5AD:\n")
cat(OUTPUT_H5AD, "\n")

cat("\nThe existing full Baituk H5AD was not modified.\n")
