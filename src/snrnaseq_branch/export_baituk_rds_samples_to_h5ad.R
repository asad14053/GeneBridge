#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Matrix)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(zellkonverter)
})

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
    PROJECT_ROOT,
    "data",
    "raw",
    "Baituk",
    "snRNA-seq_raw_countmatrices.RDS"
)

SAMPLE_OUTPUT_DIR <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "snrnaseq",
    "Baituk",
    "baituk_rds_sample_h5ad"
)

MANIFEST_CSV <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "snrnaseq",
    "Baituk",
    "baituk_rds_sample_h5ad_manifest.csv"
)

SCZ_SAMPLES <- c(
    "MB6", "MB8", "MB8-2", "MB10", "MB12",
    "MB14", "MB22", "MB23", "MB54", "MB56"
)

NTC_SAMPLES <- c(
    "MB7", "MB9", "MB11", "MB13", "MB15",
    "MB16", "MB17", "MB18-2", "MB19", "MB21",
    "MB51", "MB53", "MB55", "MB57"
)

KNOWN_SAMPLES <- c(SCZ_SAMPLES, NTC_SAMPLES)
KNOWN_SAMPLES_LONGEST_FIRST <- KNOWN_SAMPLES[
    order(nchar(KNOWN_SAMPLES), decreasing = TRUE)
]

EXPECTED_GENES <- 60617L

dir.create(SAMPLE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)

cat(paste0(rep("=", 100), collapse = ""), "\n")
cat("EXPORT BATIUK RDS AS SAMPLE-LEVEL H5AD FILES\n")
cat(paste0(rep("=", 100), collapse = ""), "\n")
cat("Input RDS:", INPUT_RDS, "\n")
cat("Sample output directory:", SAMPLE_OUTPUT_DIR, "\n")

if (!file.exists(INPUT_RDS)) {
    stop("Input RDS does not exist: ", INPUT_RDS)
}

resolve_sample_id <- function(label) {
    label <- as.character(label)

    hits <- KNOWN_SAMPLES_LONGEST_FIRST[
        vapply(
            KNOWN_SAMPLES_LONGEST_FIRST,
            function(sample_id) grepl(sample_id, label, fixed = TRUE),
            logical(1)
        )
    ]

    if (length(hits) == 0L) {
        return(NA_character_)
    }

    hits[[1]]
}

is_count_container <- function(x) {
    if (inherits(x, "Matrix") || is.matrix(x)) {
        return(TRUE)
    }

    if (inherits(x, "SingleCellExperiment") ||
        inherits(x, "SummarizedExperiment")) {
        return(TRUE)
    }

    if (inherits(x, "Seurat")) {
        return(TRUE)
    }

    if (is.list(x) && !is.null(names(x))) {
        possible_names <- c(
            "counts", "count", "raw_counts",
            "count_matrix", "matrix", "mat", "X"
        )
        return(any(possible_names %in% names(x)))
    }

    FALSE
}

extract_counts <- function(x) {
    if (inherits(x, "SingleCellExperiment") ||
        inherits(x, "SummarizedExperiment")) {

        available_assays <- SummarizedExperiment::assayNames(x)

        if (length(available_assays) == 0L) {
            stop("SummarizedExperiment has no assays.")
        }

        assay_name <- if ("counts" %in% available_assays) {
            "counts"
        } else {
            available_assays[[1]]
        }

        return(SummarizedExperiment::assay(x, assay_name))
    }

    if (inherits(x, "Seurat")) {
        if (!requireNamespace("SeuratObject", quietly = TRUE)) {
            stop(
                "A Seurat object was detected, but SeuratObject ",
                "is not installed."
            )
        }

        return(
            SeuratObject::GetAssayData(
                x,
                assay = SeuratObject::DefaultAssay(x),
                layer = "counts"
            )
        )
    }

    if (inherits(x, "Matrix") || is.matrix(x)) {
        return(x)
    }

    if (is.list(x) && !is.null(names(x))) {
        possible_names <- c(
            "counts", "count", "raw_counts",
            "count_matrix", "matrix", "mat", "X"
        )

        for (candidate_name in possible_names) {
            if (candidate_name %in% names(x)) {
                return(extract_counts(x[[candidate_name]]))
            }
        }
    }

    stop(
        "Could not extract a count matrix from object with class: ",
        paste(class(x), collapse = ", ")
    )
}

find_sample_list <- function(object) {
    if (!is.list(object)) {
        stop(
            "Expected the RDS to contain a list of sample count matrices. ",
            "Detected class: ",
            paste(class(object), collapse = ", ")
        )
    }

    candidates <- list(root = object)

    if (!is.null(names(object))) {
        for (object_name in names(object)) {
            child <- object[[object_name]]

            if (is.list(child) &&
                !is.data.frame(child) &&
                !inherits(child, "Matrix")) {
                candidates[[paste0("root$", object_name)]] <- child
            }
        }
    }

    candidate_scores <- vapply(
        candidates,
        function(candidate) {
            if (is.null(names(candidate))) {
                return(0L)
            }

            resolved <- vapply(
                names(candidate),
                resolve_sample_id,
                character(1)
            )

            count_like <- vapply(
                candidate,
                is_count_container,
                logical(1)
            )

            sum(!is.na(resolved) & count_like)
        },
        integer(1)
    )

    best_name <- names(which.max(candidate_scores))
    best_candidate <- candidates[[best_name]]
    best_score <- candidate_scores[[best_name]]

    cat(
        "Selected RDS container:", best_name,
        "with", best_score, "recognized samples\n"
    )

    if (best_score < length(KNOWN_SAMPLES)) {
        stop(
            "Only ", best_score, " of ", length(KNOWN_SAMPLES),
            " expected sample objects were detected.\n",
            "Available top-level names: ",
            paste(names(object), collapse = ", ")
        )
    }

    best_candidate
}

prepare_matrix <- function(matrix_object, sample_id) {
    counts <- extract_counts(matrix_object)

    if (!inherits(counts, "Matrix")) {
        counts <- Matrix::Matrix(counts, sparse = TRUE)
    }

    if (nrow(counts) == EXPECTED_GENES) {
        # Already genes × cells.
    } else if (ncol(counts) == EXPECTED_GENES) {
        cat("Transposing matrix for", sample_id, "\n")
        counts <- Matrix::t(counts)
    } else {
        stop(
            "Could not determine matrix orientation for ", sample_id,
            ". Shape: ", nrow(counts), " × ", ncol(counts),
            "; expected one dimension to equal ", EXPECTED_GENES
        )
    }

    counts <- methods::as(counts, "dgCMatrix")

    if (is.null(rownames(counts))) {
        stop("Gene names are missing for sample ", sample_id)
    }

    if (is.null(colnames(counts))) {
        colnames(counts) <- paste0(
            sample_id,
            "_cell_",
            seq_len(ncol(counts))
        )
    }

    original_cell_names <- colnames(counts)

    colnames(counts) <- make.unique(
        paste0(sample_id, "__", original_cell_names)
    )

    counts
}

cat("Reading RDS...\n")
raw_object <- readRDS(INPUT_RDS)
cat("RDS loaded. Class:", paste(class(raw_object), collapse = ", "), "\n")

sample_container <- find_sample_list(raw_object)

container_names <- names(sample_container)

resolved_ids <- vapply(
    container_names,
    resolve_sample_id,
    character(1)
)

valid_entries <- !is.na(resolved_ids) &
    vapply(sample_container, is_count_container, logical(1))

sample_objects <- sample_container[valid_entries]
names(sample_objects) <- resolved_ids[valid_entries]

if (anyDuplicated(names(sample_objects))) {
    duplicate_ids <- unique(
        names(sample_objects)[duplicated(names(sample_objects))]
    )

    stop(
        "Multiple RDS entries mapped to the same sample ID: ",
        paste(duplicate_ids, collapse = ", ")
    )
}

missing_samples <- setdiff(KNOWN_SAMPLES, names(sample_objects))

if (length(missing_samples) > 0L) {
    stop(
        "Missing expected samples: ",
        paste(missing_samples, collapse = ", ")
    )
}

manifest_rows <- vector("list", length(KNOWN_SAMPLES))

for (sample_index in seq_along(KNOWN_SAMPLES)) {
    sample_id <- KNOWN_SAMPLES[[sample_index]]

    diagnosis <- if (sample_id %in% SCZ_SAMPLES) {
        "SCZ"
    } else {
        "NTC"
    }

    donor_id <- if (sample_id == "MB8-2") {
        "MB8"
    } else {
        sample_id
    }

    output_h5ad <- file.path(
        SAMPLE_OUTPUT_DIR,
        paste0(sample_id, ".h5ad")
    )

    temporary_h5ad <- file.path(
        SAMPLE_OUTPUT_DIR,
        paste0(sample_id, ".tmp.h5ad")
    )

    cat("\n", paste0(rep("-", 100), collapse = ""), "\n", sep = "")
    cat(
        "Sample", sample_index, "/", length(KNOWN_SAMPLES),
        ":", sample_id, "(", diagnosis, ")\n"
    )

    counts <- prepare_matrix(
        sample_objects[[sample_id]],
        sample_id
    )

    original_gene_ids <- rownames(counts)
    unique_gene_ids <- make.unique(original_gene_ids)
    rownames(counts) <- unique_gene_ids

    cell_metadata <- S4Vectors::DataFrame(
        sample = rep(sample_id, ncol(counts)),
        baituk_sample_id = rep(sample_id, ncol(counts)),
        donor_id = rep(donor_id, ncol(counts)),
        diagnosis = rep(diagnosis, ncol(counts)),
        row.names = colnames(counts)
    )

    gene_metadata <- S4Vectors::DataFrame(
        original_gene_id = original_gene_ids,
        row.names = unique_gene_ids
    )

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = counts),
        colData = cell_metadata,
        rowData = gene_metadata
    )

    cat(
        "Dimensions:",
        nrow(sce), "genes ×", ncol(sce), "cells\n"
    )

    unlink(temporary_h5ad)

    cat("Writing:", output_h5ad, "\n")

    zellkonverter::writeH5AD(
        sce,
        temporary_h5ad,
        X_name = "counts"
    )

    if (!file.exists(temporary_h5ad) ||
        file.info(temporary_h5ad)$size == 0) {
        stop("Temporary H5AD was not written correctly: ", temporary_h5ad)
    }

    if (file.exists(output_h5ad)) {
        unlink(output_h5ad)
    }

    if (!file.rename(temporary_h5ad, output_h5ad)) {
        stop(
            "Could not rename temporary H5AD to final output: ",
            output_h5ad
        )
    }

    manifest_rows[[sample_index]] <- data.frame(
        sample_id = sample_id,
        donor_id = donor_id,
        diagnosis = diagnosis,
        n_cells = ncol(sce),
        n_genes = nrow(sce),
        h5ad = output_h5ad,
        stringsAsFactors = FALSE
    )

    rm(counts, sce, cell_metadata, gene_metadata)
    invisible(gc())
}

manifest <- do.call(rbind, manifest_rows)
write.csv(manifest, MANIFEST_CSV, row.names = FALSE)

cat("\n", paste0(rep("=", 100), collapse = ""), "\n", sep = "")
cat("SAMPLE EXPORT COMPLETE\n")
cat("Manifest:", MANIFEST_CSV, "\n")
print(manifest)
cat(paste0(rep("=", 100), collapse = ""), "\n")
