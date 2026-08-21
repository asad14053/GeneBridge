#!/usr/bin/env Rscript

################################################################################
# Create separate SCZ and NTC Batiuk snRNA-seq H5AD files directly from RDS.
#
# Input:
#   data/raw/Baituk/snRNA-seq_raw_countmatrices.RDS
#
# Outputs:
#   baituk_snrna_reference_SCZ_full_allgenes.h5ad
#   baituk_snrna_reference_NTC_full_allgenes.h5ad
#
# AnnData:
#   X = raw counts
#   obs["sample"] = Batiuk sample ID
#   obs["diagnosis"] = SCZ or NTC
#   obs["donor_id"] = donor ID; MB8-2 is assigned donor MB8
#
# Uses native anndataR. No zellkonverter and no Python installation.
################################################################################

options(stringsAsFactors = FALSE)
options(width = 120)

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

BAITUK_OUTPUT_DIR <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "snrnaseq",
    "Baituk"
)

SCZ_OUTPUT_H5AD <- file.path(
    BAITUK_OUTPUT_DIR,
    "baituk_snrna_reference_SCZ_full_allgenes.h5ad"
)

NTC_OUTPUT_H5AD <- file.path(
    BAITUK_OUTPUT_DIR,
    "baituk_snrna_reference_NTC_full_allgenes.h5ad"
)

SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "03_batiuk_huuki_gene_overlap",
    "baituk_scz_ntc_h5ad_summary.csv"
)

SAMPLE_SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "03_batiuk_huuki_gene_overlap",
    "baituk_scz_ntc_sample_summary.csv"
)

SCZ_SAMPLES <- c(
    "MB6",
    "MB8",
    "MB8-2",
    "MB10",
    "MB12",
    "MB14",
    "MB22",
    "MB23",
    "MB54",
    "MB56"
)

NTC_SAMPLES <- c(
    "MB7",
    "MB9",
    "MB11",
    "MB13",
    "MB15",
    "MB16",
    "MB17",
    "MB18-2",
    "MB19",
    "MB21",
    "MB51",
    "MB53",
    "MB55",
    "MB57"
)

ALL_SAMPLES <- c(SCZ_SAMPLES, NTC_SAMPLES)

EXPECTED_GENES <- 60617L
EXPECTED_SCZ_CELLS <- 81817L
EXPECTED_NTC_CELLS <- 127236L


section <- function(title) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    flush.console()
}


normalize_sample_id <- function(sample_name) {
    sample_name <- trimws(as.character(sample_name))

    if (sample_name %in% ALL_SAMPLES) {
        return(sample_name)
    }

    candidates <- ALL_SAMPLES[
        order(nchar(ALL_SAMPLES), decreasing = TRUE)
    ]

    matches <- candidates[
        vapply(
            candidates,
            function(candidate) {
                grepl(candidate, sample_name, fixed = TRUE)
            },
            logical(1)
        )
    ]

    if (length(matches) == 1L) {
        return(matches[[1]])
    }

    NA_character_
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
        available_assays <- assayNames(object)

        if (length(available_assays) == 0L) {
            stop("No assays found for sample ", sample_name)
        }

        assay_name <- if ("counts" %in% available_assays) {
            "counts"
        } else {
            available_assays[[1]]
        }

        matrix_object <- assay(
            object,
            assay_name,
            withDimnames = TRUE
        )

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
                "Could not identify exactly one matrix inside sample ",
                sample_name
            )
        }

        matrix_object <- object[[matrix_candidates[[1]]]]

    } else {
        stop(
            "Unsupported object for sample ",
            sample_name,
            ". Class: ",
            paste(class(object), collapse = "; ")
        )
    }

    if (is.null(dim(matrix_object))) {
        stop("No matrix dimensions for sample ", sample_name)
    }

    if (is.null(rownames(matrix_object))) {
        stop("No gene names for sample ", sample_name)
    }

    if (anyDuplicated(rownames(matrix_object)) > 0L) {
        stop("Duplicated gene names in sample ", sample_name)
    }

    if (!inherits(matrix_object, "dgCMatrix")) {
        matrix_object <- methods::as(
            matrix_object,
            "dgCMatrix"
        )
    }

    matrix_object
}


write_group_h5ad <- function(
    diagnosis,
    sample_ids,
    matrices,
    obs_tables,
    reference_genes,
    output_h5ad,
    expected_cells
) {
    section(
        paste0("Combining and writing ", diagnosis)
    )

    missing_group_samples <- setdiff(
        sample_ids,
        names(matrices)
    )

    if (length(missing_group_samples) > 0L) {
        stop(
            diagnosis,
            " samples missing from matrix list: ",
            paste(missing_group_samples, collapse = ", ")
        )
    }

    cat(
        diagnosis,
        " samples:\n",
        paste(sample_ids, collapse = ", "),
        "\n"
    )

    combined_counts <- do.call(
        cbind,
        unname(matrices[sample_ids])
    )

    if (!inherits(combined_counts, "dgCMatrix")) {
        combined_counts <- methods::as(
            combined_counts,
            "dgCMatrix"
        )
    }

    combined_obs <- do.call(
        rbind,
        unname(obs_tables[sample_ids])
    )

    if (nrow(combined_obs) != ncol(combined_counts)) {
        stop(
            diagnosis,
            " observation metadata does not match cell count."
        )
    }

    if (!identical(
        rownames(combined_obs),
        colnames(combined_counts)
    )) {
        stop(
            diagnosis,
            " observation names do not match matrix column names."
        )
    }

    if (nrow(combined_counts) != EXPECTED_GENES) {
        stop(
            diagnosis,
            " gene-count mismatch. Expected ",
            EXPECTED_GENES,
            ", found ",
            nrow(combined_counts)
        )
    }

    if (ncol(combined_counts) != expected_cells) {
        stop(
            diagnosis,
            " cell-count mismatch. Expected ",
            expected_cells,
            ", found ",
            ncol(combined_counts)
        )
    }

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

    cat(
        diagnosis,
        " dimensions [genes x cells]: ",
        nrow(sce),
        " x ",
        ncol(sce),
        "\n",
        sep = ""
    )

    temporary_h5ad <- sub(
        "\\.h5ad$",
        ".tmp.h5ad",
        output_h5ad
    )

    unlink(temporary_h5ad)

    cat(
        "Writing temporary H5AD without compression:\n",
        temporary_h5ad,
        "\n"
    )

    anndataR::write_h5ad(
        object = sce,
        path = temporary_h5ad,
        compression = NULL,
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

    if (
        !file.exists(temporary_h5ad) ||
        file.info(temporary_h5ad)$size == 0
    ) {
        stop(
            diagnosis,
            " temporary H5AD was not created correctly."
        )
    }

    if (file.exists(output_h5ad)) {
        cat(
            "Removing existing final output:\n",
            output_h5ad,
            "\n"
        )
        unlink(output_h5ad)
    }

    if (!file.rename(
        temporary_h5ad,
        output_h5ad
    )) {
        stop(
            "Could not rename temporary file to:\n",
            output_h5ad
        )
    }

    output_size_gb <- file.info(output_h5ad)$size / 1024^3

    cat(
        diagnosis,
        " completed successfully.\n",
        "Output: ",
        output_h5ad,
        "\n",
        "Output size: ",
        round(output_size_gb, 3),
        " GB\n",
        sep = ""
    )

    result <- data.frame(
        diagnosis = diagnosis,
        n_samples = length(sample_ids),
        n_cells = ncol(sce),
        n_genes = nrow(sce),
        output_h5ad = output_h5ad,
        output_size_gb = round(output_size_gb, 4),
        stringsAsFactors = FALSE
    )

    rm(
        combined_counts,
        combined_obs,
        gene_metadata,
        sce
    )

    invisible(gc())

    result
}


################################################################################
# Read and validate RDS
################################################################################

section("Create separate Batiuk SCZ and NTC H5AD files")

cat("Input RDS:\n", INPUT_RDS, "\n")
cat("SCZ output:\n", SCZ_OUTPUT_H5AD, "\n")
cat("NTC output:\n", NTC_OUTPUT_H5AD, "\n")

if (!file.exists(INPUT_RDS)) {
    stop("Input RDS not found:\n", INPUT_RDS)
}

dir.create(
    BAITUK_OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    dirname(SUMMARY_CSV),
    recursive = TRUE,
    showWarnings = FALSE
)

section("Reading Batiuk RDS")

baituk_objects <- readRDS(INPUT_RDS)

if (!is.list(baituk_objects)) {
    stop(
        "Expected the Batiuk RDS to contain a list of sample matrices."
    )
}

original_names <- names(baituk_objects)

if (is.null(original_names)) {
    stop(
        "The Batiuk RDS list has no sample names."
    )
}

resolved_names <- vapply(
    original_names,
    normalize_sample_id,
    character(1)
)

if (anyNA(resolved_names)) {
    stop(
        "Could not resolve these RDS object names:\n",
        paste(
            original_names[is.na(resolved_names)],
            collapse = ", "
        )
    )
}

if (anyDuplicated(resolved_names) > 0L) {
    duplicated_ids <- unique(
        resolved_names[
            duplicated(resolved_names)
        ]
    )

    stop(
        "Multiple RDS entries mapped to the same sample ID: ",
        paste(duplicated_ids, collapse = ", ")
    )
}

names(baituk_objects) <- resolved_names

missing_samples <- setdiff(
    ALL_SAMPLES,
    names(baituk_objects)
)

if (length(missing_samples) > 0L) {
    stop(
        "Expected samples missing from the RDS: ",
        paste(missing_samples, collapse = ", ")
    )
}

baituk_objects <- baituk_objects[ALL_SAMPLES]

cat(
    "Resolved samples:\n",
    paste(names(baituk_objects), collapse = ", "),
    "\n"
)


################################################################################
# Extract matrices and metadata
################################################################################

section("Extracting sample matrices")

matrices <- list()
obs_tables <- list()
sample_records <- list()

reference_genes <- NULL

for (sample_name in ALL_SAMPLES) {

    diagnosis <- if (
        sample_name %in% SCZ_SAMPLES
    ) {
        "SCZ"
    } else {
        "NTC"
    }

    donor_id <- if (
        sample_name == "MB8-2"
    ) {
        "MB8"
    } else {
        sample_name
    }

    cat(
        "\nProcessing ",
        sample_name,
        " [",
        diagnosis,
        "]...\n",
        sep = ""
    )

    matrix_object <- extract_count_matrix(
        baituk_objects[[sample_name]],
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
                " does not contain the reference gene set."
            )
        }

        cat(
            "Gene order differs; reordering sample.\n"
        )

        matrix_object <- matrix_object[
            match(reference_genes, current_genes),
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

    global_cell_ids <- make.unique(
        paste0(
            sample_name,
            "::",
            original_cell_ids
        )
    )

    colnames(matrix_object) <- global_cell_ids

    matrices[[sample_name]] <- matrix_object

    obs_tables[[sample_name]] <- data.frame(
        sample = rep(
            sample_name,
            ncol(matrix_object)
        ),
        baituk_sample_id = rep(
            sample_name,
            ncol(matrix_object)
        ),
        donor_id = rep(
            donor_id,
            ncol(matrix_object)
        ),
        diagnosis = rep(
            diagnosis,
            ncol(matrix_object)
        ),
        original_cell_id = original_cell_ids,
        stringsAsFactors = FALSE,
        row.names = global_cell_ids
    )

    sample_records[[sample_name]] <- data.frame(
        sample = sample_name,
        donor_id = donor_id,
        diagnosis = diagnosis,
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

rm(baituk_objects)
invisible(gc())

sample_summary <- do.call(
    rbind,
    sample_records
)

cat("\nSample summary:\n")
print(sample_summary)

if (length(reference_genes) != EXPECTED_GENES) {
    stop(
        "Expected ",
        EXPECTED_GENES,
        " genes, found ",
        length(reference_genes)
    )
}


################################################################################
# Write final SCZ and NTC H5AD files
################################################################################

scz_summary <- write_group_h5ad(
    diagnosis = "SCZ",
    sample_ids = SCZ_SAMPLES,
    matrices = matrices,
    obs_tables = obs_tables,
    reference_genes = reference_genes,
    output_h5ad = SCZ_OUTPUT_H5AD,
    expected_cells = EXPECTED_SCZ_CELLS
)

ntc_summary <- write_group_h5ad(
    diagnosis = "NTC",
    sample_ids = NTC_SAMPLES,
    matrices = matrices,
    obs_tables = obs_tables,
    reference_genes = reference_genes,
    output_h5ad = NTC_OUTPUT_H5AD,
    expected_cells = EXPECTED_NTC_CELLS
)


################################################################################
# Save summaries
################################################################################

summary_df <- rbind(
    scz_summary,
    ntc_summary
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

cat("\nSCZ H5AD:\n", SCZ_OUTPUT_H5AD, "\n")
cat("\nNTC H5AD:\n", NTC_OUTPUT_H5AD, "\n")
cat("\nSummary:\n", SUMMARY_CSV, "\n")
cat("\nSample summary:\n", SAMPLE_SUMMARY_CSV, "\n")
cat("\nConversion completed successfully.\n")
