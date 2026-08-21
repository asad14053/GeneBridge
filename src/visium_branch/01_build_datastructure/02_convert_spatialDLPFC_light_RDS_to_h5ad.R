#!/usr/bin/env Rscript

################################################################################
# Convert the lightweight Huuki Visium RDS to H5AD.
#
# Input:
#   data/processed/visium/spatialDLPFC_Visium_sce_light.RDS
#
# Output:
#   data/processed/visium/spatialDLPFC_Visium_sce.h5ad
#
# AnnData structure:
#   X                   <- logcounts
#   layers["counts"]    <- raw counts
#   obs                 <- colData
#   var                 <- rowData
#   obsm["spatial"]     <- Visium spatial coordinates
#
# This script uses anndataR, not zellkonverter.
################################################################################

Sys.setenv(
    HDF5_USE_FILE_LOCKING = "FALSE"
)

options(
    anndataR.write_null = FALSE
)

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(scuttle)
    library(anndataR)
})

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "visium",
    "spatialDLPFC_Visium_sce_light.RDS"
)

OUTPUT_H5AD <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "visium",
    "spatialDLPFC_Visium_sce.h5ad"
)

SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "visium_branch",
    "01_build_datastructure",
    "spatialDLPFC_Visium_sce_h5ad_conversion_summary.csv"
)


section <- function(title) {

    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


clean_metadata <- function(metadata_object, metadata_name) {

    metadata_df <- as.data.frame(
        metadata_object,
        stringsAsFactors = FALSE,
        optional = TRUE
    )

    if (ncol(metadata_df) == 0) {
        return(
            list(
                data = metadata_df,
                dropped = character(0)
            )
        )
    }

    keep_columns <- logical(
        ncol(metadata_df)
    )

    names(keep_columns) <- colnames(
        metadata_df
    )

    dropped_columns <- character(0)

    for (column_name in colnames(metadata_df)) {

        values <- metadata_df[[column_name]]

        if (is.factor(values)) {

            metadata_df[[column_name]] <- as.character(values)

            keep_columns[
                column_name
            ] <- TRUE

        } else if (
            is.character(values) ||
            is.numeric(values) ||
            is.integer(values) ||
            is.logical(values)
        ) {

            keep_columns[
                column_name
            ] <- TRUE

        } else {

            converted_values <- tryCatch(
                as.character(values),
                error = function(e) NULL
            )

            if (
                !is.null(converted_values) &&
                length(converted_values) == nrow(metadata_df)
            ) {

                metadata_df[[column_name]] <- converted_values

                keep_columns[
                    column_name
                ] <- TRUE

            } else {

                keep_columns[
                    column_name
                ] <- FALSE

                dropped_columns <- c(
                    dropped_columns,
                    column_name
                )
            }
        }
    }

    if (length(dropped_columns) > 0) {

        cat(
            "\nDropping unsupported ",
            metadata_name,
            " columns:\n",
            sep = ""
        )

        print(
            dropped_columns
        )
    }

    metadata_df <- metadata_df[
        ,
        keep_columns,
        drop = FALSE
    ]

    return(
        list(
            data = metadata_df,
            dropped = dropped_columns
        )
    )
}


extract_spatial_coordinates <- function(sce) {

    coordinates <- NULL
    source_name <- NULL

    # First choice: reducedDim(sce, "spatial")
    if (
        "spatial"
        %in%
        reducedDimNames(sce)
    ) {

        coordinates <- reducedDim(
            sce,
            "spatial"
        )

        source_name <- "reducedDim_spatial"
    }

    # Second choice: SpatialExperiment::spatialCoords()
    if (
        is.null(coordinates) &&
        requireNamespace(
            "SpatialExperiment",
            quietly = TRUE
        ) &&
        inherits(
            sce,
            "SpatialExperiment"
        )
    ) {

        coordinates <- tryCatch(
            SpatialExperiment::spatialCoords(
                sce
            ),
            error = function(e) NULL
        )

        if (!is.null(coordinates)) {
            source_name <- "spatialCoords"
        }
    }

    # Third choice: full-resolution pixel coordinates in colData.
    if (
        is.null(coordinates) &&
        all(
            c(
                "pxl_col_in_fullres",
                "pxl_row_in_fullres"
            )
            %in%
            colnames(
                colData(sce)
            )
        )
    ) {

        coordinates <- cbind(
            x = as.numeric(
                colData(sce)$pxl_col_in_fullres
            ),
            y = as.numeric(
                colData(sce)$pxl_row_in_fullres
            )
        )

        source_name <- "colData_fullres_pixels"
    }

    if (is.null(coordinates)) {

        stop(
            "Could not find spatial coordinates in:\n",
            "  reducedDim(sce, 'spatial')\n",
            "  SpatialExperiment::spatialCoords(sce)\n",
            "  colData pixel-coordinate columns"
        )
    }

    coordinates <- as.matrix(
        coordinates
    )

    if (nrow(coordinates) != ncol(sce)) {

        stop(
            "Spatial-coordinate row count does not match spot count.\n",
            "Spatial rows: ",
            nrow(coordinates),
            "\nSpots: ",
            ncol(sce)
        )
    }

    if (ncol(coordinates) < 2) {

        stop(
            "Spatial coordinates have fewer than two columns."
        )
    }

    coordinates <- coordinates[
        ,
        1:2,
        drop = FALSE
    ]

    storage.mode(
        coordinates
    ) <- "double"

    return(
        list(
            coordinates = coordinates,
            source = source_name
        )
    )
}


################################################################################
# Read the RDS
################################################################################

section(
    "Convert spatialDLPFC Visium light RDS to H5AD"
)

cat(
    "Input RDS:\n",
    INPUT_RDS,
    "\n"
)

cat(
    "\nOutput H5AD:\n",
    OUTPUT_H5AD,
    "\n"
)

if (!file.exists(INPUT_RDS)) {

    stop(
        "Input RDS was not found:\n",
        INPUT_RDS
    )
}

section(
    "Reading RDS"
)

sce_original <- readRDS(
    INPUT_RDS
)

cat(
    "Object class:\n"
)

print(
    class(sce_original)
)

cat(
    "\nObject dimensions [genes x spots]:\n"
)

print(
    dim(sce_original)
)

if (!inherits(
    sce_original,
    "SingleCellExperiment"
)) {

    stop(
        "The RDS does not contain a SingleCellExperiment-compatible object."
    )
}

cat(
    "\nAvailable assays:\n"
)

print(
    assayNames(sce_original)
)

cat(
    "\nAvailable reduced dimensions:\n"
)

print(
    reducedDimNames(sce_original)
)

################################################################################
# Prepare names
################################################################################

section(
    "Preparing gene and spot names"
)

if (is.null(
    rownames(sce_original)
)) {

    gene_names <- paste0(
        "gene_",
        seq_len(
            nrow(sce_original)
        )
    )

} else {

    gene_names <- make.unique(
        as.character(
            rownames(sce_original)
        )
    )
}


if (is.null(
    colnames(sce_original)
)) {

    spot_names <- paste0(
        "spot_",
        seq_len(
            ncol(sce_original)
        )
    )

} else {

    spot_names <- make.unique(
        as.character(
            colnames(sce_original)
        )
    )
}


cat(
    "Genes: ",
    length(gene_names),
    "\n",
    sep = ""
)

cat(
    "Spots: ",
    length(spot_names),
    "\n",
    sep = ""
)

################################################################################
# Check counts and logcounts
################################################################################

section(
    "Preparing assays"
)

if (!(
    "counts"
    %in%
    assayNames(sce_original)
)) {

    stop(
        "The input object does not contain a 'counts' assay."
    )
}


if (!(
    "logcounts"
    %in%
    assayNames(sce_original)
)) {

    cat(
        "logcounts not found. Creating logcounts with scuttle::logNormCounts().\n"
    )

    sce_original <- scuttle::logNormCounts(
        sce_original
    )

} else {

    cat(
        "logcounts assay already exists.\n"
    )
}


# Remove assay dimnames during extraction.
# This avoids the previous SummarizedExperiment dimname mismatch error.

counts_matrix <- assay(
    sce_original,
    "counts",
    withDimnames = FALSE
)

logcounts_matrix <- assay(
    sce_original,
    "logcounts",
    withDimnames = FALSE
)


cat(
    "\nCounts dimensions:\n"
)

print(
    dim(counts_matrix)
)


cat(
    "\nLogcounts dimensions:\n"
)

print(
    dim(logcounts_matrix)
)


if (!identical(
    dim(counts_matrix),
    dim(logcounts_matrix)
)) {

    stop(
        "counts and logcounts have different dimensions."
    )
}

################################################################################
# Prepare metadata
################################################################################

section(
    "Preparing metadata"
)

coldata_result <- clean_metadata(
    colData(sce_original),
    "colData"
)

rowdata_result <- clean_metadata(
    rowData(sce_original),
    "rowData"
)

clean_coldata <- coldata_result$data
clean_rowdata <- rowdata_result$data


cat(
    "Retained colData columns: ",
    ncol(clean_coldata),
    "\n",
    sep = ""
)

cat(
    "Retained rowData columns: ",
    ncol(clean_rowdata),
    "\n",
    sep = ""
)

################################################################################
# Prepare spatial coordinates
################################################################################

section(
    "Preparing spatial coordinates"
)

spatial_result <- extract_spatial_coordinates(
    sce_original
)

spatial_coordinates <- spatial_result$coordinates
spatial_source <- spatial_result$source


cat(
    "Spatial-coordinate source: ",
    spatial_source,
    "\n",
    sep = ""
)

cat(
    "Spatial dimensions [spots x 2]:\n"
)

print(
    dim(spatial_coordinates)
)

################################################################################
# Build a clean SCE
#
# Build assays first with no dimnames, then assign names and metadata.
################################################################################

section(
    "Building clean SingleCellExperiment"
)

sce_clean <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
        counts = counts_matrix,
        logcounts = logcounts_matrix
    )
)


rownames(sce_clean) <- gene_names
colnames(sce_clean) <- spot_names


rownames(clean_rowdata) <- gene_names
rownames(clean_coldata) <- spot_names


rowData(sce_clean) <- S4Vectors::DataFrame(
    clean_rowdata
)

colData(sce_clean) <- S4Vectors::DataFrame(
    clean_coldata
)


rownames(
    spatial_coordinates
) <- spot_names


reducedDim(
    sce_clean,
    "spatial"
) <- spatial_coordinates


cat(
    "\nClean object:\n"
)

print(
    sce_clean
)


cat(
    "\nAssays:\n"
)

print(
    assayNames(sce_clean)
)


cat(
    "\nReduced dimensions:\n"
)

print(
    reducedDimNames(sce_clean)
)

################################################################################
# Validate before writing
################################################################################

section(
    "Pre-write validation"
)

if (nrow(sce_clean) != nrow(sce_original)) {

    stop(
        "Gene count changed while creating the clean object."
    )
}


if (ncol(sce_clean) != ncol(sce_original)) {

    stop(
        "Spot count changed while creating the clean object."
    )
}


if (
    nrow(
        reducedDim(
            sce_clean,
            "spatial"
        )
    )
    !=
    ncol(sce_clean)
) {

    stop(
        "Spatial-coordinate count does not match spot count."
    )
}


cat(
    "Pre-write validation passed.\n"
)

################################################################################
# Write H5AD
################################################################################

dir.create(
    dirname(
        OUTPUT_H5AD
    ),
    recursive = TRUE,
    showWarnings = FALSE
)


dir.create(
    dirname(
        SUMMARY_CSV
    ),
    recursive = TRUE,
    showWarnings = FALSE
)


if (file.exists(
    OUTPUT_H5AD
)) {

    cat(
        "\nRemoving existing output:\n",
        OUTPUT_H5AD,
        "\n"
    )

    unlink(
        OUTPUT_H5AD
    )
}


section(
    "Writing H5AD using anndataR"
)


anndataR::write_h5ad(
    object = sce_clean,
    path = OUTPUT_H5AD,
    compression = "lzf",
    mode = "w",
    x_mapping = "logcounts",
    layers_mapping = c(
        counts = "counts"
    ),
    obs_mapping = TRUE,
    var_mapping = TRUE,
    obsm_mapping = list(
        spatial = "spatial"
    ),
    varm_mapping = FALSE,
    obsp_mapping = FALSE,
    varp_mapping = FALSE,
    uns_mapping = FALSE
)


if (!file.exists(
    OUTPUT_H5AD
)) {

    stop(
        "The H5AD output file was not created."
    )
}


cat(
    "H5AD successfully written:\n",
    OUTPUT_H5AD,
    "\n"
)

################################################################################
# Read-back validation
################################################################################

section(
    "Validating written H5AD"
)

sce_check <- anndataR::read_h5ad(
    OUTPUT_H5AD,
    as = "SingleCellExperiment"
)


cat(
    "Read-back object:\n"
)

print(
    sce_check
)


cat(
    "\nRead-back dimensions [genes x spots]:\n"
)

print(
    dim(sce_check)
)


cat(
    "\nRead-back assays:\n"
)

print(
    assayNames(sce_check)
)


cat(
    "\nRead-back reduced dimensions:\n"
)

print(
    reducedDimNames(sce_check)
)


if (nrow(sce_check) != nrow(sce_clean)) {

    stop(
        "Read-back validation failed: gene count mismatch."
    )
}


if (ncol(sce_check) != ncol(sce_clean)) {

    stop(
        "Read-back validation failed: spot count mismatch."
    )
}


if (!(
    "spatial"
    %in%
    reducedDimNames(sce_check)
)) {

    stop(
        "Read-back validation failed: spatial coordinates are missing."
    )
}

################################################################################
# Summary CSV
################################################################################

summary_df <- data.frame(
    metric = c(
        "input_rds",
        "output_h5ad",
        "n_genes",
        "n_spots",
        "X_source",
        "counts_layer_source",
        "spatial_source",
        "has_spatial_after_readback",
        "output_size_GB",
        "dropped_colData_columns",
        "dropped_rowData_columns"
    ),
    value = c(
        INPUT_RDS,
        OUTPUT_H5AD,
        as.character(
            nrow(sce_check)
        ),
        as.character(
            ncol(sce_check)
        ),
        "logcounts",
        "counts",
        spatial_source,
        as.character(
            "spatial"
            %in%
            reducedDimNames(sce_check)
        ),
        as.character(
            round(
                file.info(
                    OUTPUT_H5AD
                )$size / 1024^3,
                4
            )
        ),
        paste(
            coldata_result$dropped,
            collapse = ";"
        ),
        paste(
            rowdata_result$dropped,
            collapse = ";"
        )
    ),
    stringsAsFactors = FALSE
)


write.csv(
    summary_df,
    SUMMARY_CSV,
    row.names = FALSE
)


cat(
    "\nConversion summary:\n"
)

print(
    summary_df
)


section(
    "DONE"
)

cat(
    "H5AD:\n",
    OUTPUT_H5AD,
    "\n"
)

cat(
    "\nSummary CSV:\n",
    SUMMARY_CSV,
    "\n"
)
