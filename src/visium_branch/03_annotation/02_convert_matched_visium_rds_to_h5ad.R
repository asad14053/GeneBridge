#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(scuttle)
    library(anndataR)
})

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

################################################################################
# Configuration
################################################################################

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "visium",
    "visium_N24_matched_layer_annotated.rds"
)

OUTPUT_H5AD <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "visium",
    "visium_N24_matched_layer_annotated.h5ad"
)

SUMMARY_CSV <- file.path(
    PROJECT_ROOT,
    "outputs",
    "visium_branch",
    "03_annotation",
    "layerwise_matched_N24",
    "tables",
    "visium_N24_h5ad_conversion_summary.csv"
)

XENIUM_H5AD <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "xenium",
    "xenium_N24_layer_celltype_annotated.h5ad"
)

################################################################################
# Helper functions
################################################################################

section <- function(title) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


clean_metadata <- function(x, name = "metadata") {

    df <- as.data.frame(
        x,
        stringsAsFactors = FALSE,
        optional = TRUE
    )

    if (ncol(df) == 0) {
        return(df)
    }

    keep <- logical(ncol(df))
    names(keep) <- colnames(df)

    dropped <- character()

    for (nm in colnames(df)) {

        v <- df[[nm]]

        if (
            is.character(v) ||
            is.numeric(v) ||
            is.integer(v) ||
            is.logical(v)
        ) {

            keep[nm] <- TRUE

        } else if (is.factor(v)) {

            df[[nm]] <- as.character(v)
            keep[nm] <- TRUE

        } else {

            # Try converting simple unusual atomic classes to character.
            converted <- tryCatch(
                as.character(v),
                error = function(e) NULL
            )

            if (
                !is.null(converted) &&
                length(converted) == nrow(df)
            ) {

                df[[nm]] <- converted
                keep[nm] <- TRUE

            } else {

                keep[nm] <- FALSE
                dropped <- c(dropped, nm)
            }
        }
    }

    if (length(dropped) > 0) {

        cat(
            "\nWARNING: Dropping unsupported ",
            name,
            " columns:\n",
            sep = ""
        )

        print(dropped)
    }

    df <- df[
        ,
        keep,
        drop = FALSE
    ]

    return(df)
}


add_summary <- function(rows, metric, value) {

    rows[[length(rows) + 1]] <- data.frame(
        metric = as.character(metric),
        value = paste(value, collapse = ";"),
        stringsAsFactors = FALSE
    )

    return(rows)
}

################################################################################
# Start
################################################################################

section("Convert matched N24 Visium RDS to H5AD")

cat("Input RDS:\n")
cat(INPUT_RDS, "\n")

cat("\nOutput H5AD:\n")
cat(OUTPUT_H5AD, "\n")

cat("\nSummary CSV:\n")
cat(SUMMARY_CSV, "\n")

cat("\nXenium H5AD:\n")
cat(XENIUM_H5AD, "\n")

################################################################################
# Check input
################################################################################

if (!file.exists(INPUT_RDS)) {

    stop(
        "Input RDS does not exist:\n",
        INPUT_RDS
    )
}

################################################################################
# Read RDS
################################################################################

section("Reading Visium RDS")

spe <- readRDS(INPUT_RDS)

cat("Class:\n")
print(class(spe))

cat("\nDimensions [genes x spots]:\n")
print(dim(spe))

cat("\nAssays:\n")
print(assayNames(spe))

cat("\nNumber of colData columns:\n")
print(ncol(colData(spe)))

cat("\nNumber of rowData columns:\n")
print(ncol(rowData(spe)))

################################################################################
# Check object type
################################################################################

if (!inherits(spe, "SingleCellExperiment")) {

    stop(
        "Input object is not a SingleCellExperiment/SpatialExperiment."
    )
}

################################################################################
# Prepare gene and spot names
################################################################################

section("Preparing gene and spot names")

if (is.null(rownames(spe))) {

    gene_names <- paste0(
        "gene_",
        seq_len(nrow(spe))
    )

} else {

    gene_names <- make.unique(
        as.character(
            rownames(spe)
        )
    )
}


if (is.null(colnames(spe))) {

    spot_names <- paste0(
        "spot_",
        seq_len(ncol(spe))
    )

} else {

    spot_names <- make.unique(
        as.character(
            colnames(spe)
        )
    )
}


cat("Genes:", length(gene_names), "\n")
cat("Spots:", length(spot_names), "\n")

################################################################################
# Check counts
################################################################################

section("Checking assays")

if (!("counts" %in% assayNames(spe))) {

    stop(
        "Required assay 'counts' is missing."
    )
}


cat("counts assay found.\n")

################################################################################
# Create logcounts if necessary
################################################################################

if (!("logcounts" %in% assayNames(spe))) {

    cat(
        "\nlogcounts not found.\n"
    )

    cat(
        "Creating logcounts using scuttle::logNormCounts()...\n"
    )

    spe <- scuttle::logNormCounts(spe)

} else {

    cat(
        "logcounts assay already exists.\n"
    )
}


cat("\nAvailable assays:\n")
print(assayNames(spe))

################################################################################
# Extract assays WITHOUT dimnames
#
# This directly addresses the previous error:
#
# "the rownames and colnames of the supplied assay(s) must be NULL
#  or identical..."
#
# We deliberately strip assay dimnames first.
################################################################################

section("Extracting expression matrices")

counts_matrix <- assay(
    spe,
    "counts",
    withDimnames = FALSE
)

logcounts_matrix <- assay(
    spe,
    "logcounts",
    withDimnames = FALSE
)


cat("counts dimensions:\n")
print(dim(counts_matrix))

cat("\nlogcounts dimensions:\n")
print(dim(logcounts_matrix))


if (!identical(
    dim(counts_matrix),
    dim(logcounts_matrix)
)) {

    stop(
        "counts and logcounts dimensions do not match."
    )
}

################################################################################
# Extract metadata
################################################################################

section("Preparing metadata")

cd <- clean_metadata(
    colData(spe),
    "colData"
)

rd <- clean_metadata(
    rowData(spe),
    "rowData"
)


cat(
    "Retained colData columns:",
    ncol(cd),
    "\n"
)

cat(
    "Retained rowData columns:",
    ncol(rd),
    "\n"
)

################################################################################
# Extract spatial coordinates
################################################################################

section("Preparing spatial coordinates")

coords <- tryCatch(
    spatialCoords(spe),
    error = function(e) NULL
)


if (is.null(coords)) {

    stop(
        "spatialCoords(spe) is missing."
    )
}


coords <- as.matrix(coords)


cat(
    "Original spatial dimensions:\n"
)

print(
    dim(coords)
)


if (nrow(coords) != ncol(spe)) {

    stop(
        "Spatial coordinate row count does not match number of spots.\n",
        "Spatial rows: ",
        nrow(coords),
        "\nSpots: ",
        ncol(spe)
    )
}


if (ncol(coords) < 2) {

    stop(
        "spatialCoords contains fewer than two coordinate columns."
    )
}


# Keep the first two spatial-coordinate columns.
coords <- coords[
    ,
    1:2,
    drop = FALSE
]


rownames(coords) <- spot_names


cat(
    "\nFinal spatial dimensions [spots x 2]:\n"
)

print(
    dim(coords)
)


cat(
    "\nFirst spatial coordinates:\n"
)

print(
    head(coords)
)

################################################################################
# Build CLEAN SingleCellExperiment
#
# Important sequence:
#
# 1. Build SCE using assays with NO dimnames
# 2. Assign gene/spot names
# 3. Attach rowData and colData
# 4. Attach spatial coordinates
#
# This avoids the previous SummarizedExperiment dimname mismatch.
################################################################################

section("Building clean SingleCellExperiment")

sce <- SingleCellExperiment(
    assays = list(
        counts = counts_matrix,
        logcounts = logcounts_matrix
    )
)


rownames(sce) <- gene_names

colnames(sce) <- spot_names


# Ensure metadata row names match the SCE exactly.
rownames(rd) <- gene_names

rownames(cd) <- spot_names


rowData(sce) <- S4Vectors::DataFrame(
    rd
)

colData(sce) <- S4Vectors::DataFrame(
    cd
)


reducedDim(
    sce,
    "spatial"
) <- coords


cat("\nClean SCE:\n")
print(sce)


cat("\nAssays:\n")
print(assayNames(sce))


cat("\nReduced dimensions:\n")
print(reducedDimNames(sce))

################################################################################
# Validate clean SCE
################################################################################

section("Validating clean SCE")

if (nrow(sce) != length(gene_names)) {

    stop(
        "Gene count mismatch in clean SCE."
    )
}


if (ncol(sce) != length(spot_names)) {

    stop(
        "Spot count mismatch in clean SCE."
    )
}


if (
    nrow(
        reducedDim(
            sce,
            "spatial"
        )
    ) != ncol(sce)
) {

    stop(
        "Spatial coordinate count does not match SCE spot count."
    )
}


cat(
    "Clean SCE validation PASSED.\n"
)

################################################################################
# Create directories
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

################################################################################
# Remove previous H5AD
################################################################################

if (file.exists(OUTPUT_H5AD)) {

    cat(
        "\nRemoving existing output:\n"
    )

    cat(
        OUTPUT_H5AD,
        "\n"
    )

    unlink(
        OUTPUT_H5AD
    )
}

################################################################################
# Write H5AD using anndataR
#
# Mapping:
#
# AnnData X
#     <- SCE assay "logcounts"
#
# AnnData layers["counts"]
#     <- SCE assay "counts"
#
# AnnData obs
#     <- SCE colData
#
# AnnData var
#     <- SCE rowData
#
# AnnData obsm["spatial"]
#     <- SCE reducedDim "spatial"
################################################################################

section("Writing H5AD with anndataR")

cat(
    "Writing:\n",
    OUTPUT_H5AD,
    "\n"
)


anndataR::write_h5ad(
    object = sce,
    path = OUTPUT_H5AD,
    compression = "lzf",
    x_mapping = "logcounts",
    layers_mapping = c(
        counts = "counts"
    ),
    obs_mapping = TRUE,
    var_mapping = TRUE,
    obsm_mapping = c(
        spatial = "spatial"
    ),
    varm_mapping = FALSE,
    obsp_mapping = FALSE,
    varp_mapping = FALSE,
    uns_mapping = FALSE
)


if (!file.exists(OUTPUT_H5AD)) {

    stop(
        "Output H5AD was not created."
    )
}


cat(
    "\nH5AD successfully written.\n"
)

cat(
    "File size GB:",
    round(
        file.info(OUTPUT_H5AD)$size / 1024^3,
        4
    ),
    "\n"
)

################################################################################
# Read H5AD back as SingleCellExperiment
#
# This uses anndataR's documented reader and avoids depending on R6 method names.
################################################################################

section("Validating written H5AD")

check <- anndataR::read_h5ad(
    OUTPUT_H5AD,
    as = "SingleCellExperiment"
)


cat(
    "Read-back object:\n"
)

print(
    check
)


cat(
    "\nRead-back dimensions [genes x spots]:\n"
)

print(
    dim(check)
)


cat(
    "\nRead-back assays:\n"
)

print(
    assayNames(check)
)


cat(
    "\nRead-back reducedDims:\n"
)

print(
    reducedDimNames(check)
)


cat(
    "\nRead-back colData columns:\n"
)

print(
    colnames(
        colData(check)
    )
)


cat(
    "\nRead-back rowData columns:\n"
)

print(
    colnames(
        rowData(check)
    )
)

################################################################################
# Validate dimensions
################################################################################

if (nrow(check) != nrow(sce)) {

    stop(
        "Validation failed: gene count changed after H5AD conversion."
    )
}


if (ncol(check) != ncol(sce)) {

    stop(
        "Validation failed: spot count changed after H5AD conversion."
    )
}

################################################################################
# Determine read-back mappings
################################################################################

has_counts <- (
    "counts"
    %in%
    assayNames(check)
)


has_spatial <- (
    "spatial"
    %in%
    reducedDimNames(check)
)


has_brnum_matched <- (
    "BrNum_matched"
    %in%
    colnames(
        colData(check)
    )
)


has_layer_annotation <- (
    "visium_layer_annotation"
    %in%
    colnames(
        colData(check)
    )
)

################################################################################
# Print annotation checks
################################################################################

cat(
    "\nValidation checks:\n"
)

cat(
    "counts assay:",
    has_counts,
    "\n"
)

cat(
    "spatial reducedDim:",
    has_spatial,
    "\n"
)

cat(
    "BrNum_matched:",
    has_brnum_matched,
    "\n"
)

cat(
    "visium_layer_annotation:",
    has_layer_annotation,
    "\n"
)

################################################################################
# Build summary CSV
################################################################################

section("Creating conversion summary CSV")

summary_rows <- list()


summary_rows <- add_summary(
    summary_rows,
    "input_rds",
    INPUT_RDS
)


summary_rows <- add_summary(
    summary_rows,
    "output_h5ad",
    OUTPUT_H5AD
)


summary_rows <- add_summary(
    summary_rows,
    "xenium_h5ad",
    XENIUM_H5AD
)


summary_rows <- add_summary(
    summary_rows,
    "xenium_h5ad_exists",
    file.exists(XENIUM_H5AD)
)


summary_rows <- add_summary(
    summary_rows,
    "visium_h5ad_exists",
    file.exists(OUTPUT_H5AD)
)


summary_rows <- add_summary(
    summary_rows,
    "n_spots",
    ncol(check)
)


summary_rows <- add_summary(
    summary_rows,
    "n_genes",
    nrow(check)
)


summary_rows <- add_summary(
    summary_rows,
    "AnnData_X_source",
    "logcounts"
)


summary_rows <- add_summary(
    summary_rows,
    "AnnData_counts_layer_source",
    "counts"
)


summary_rows <- add_summary(
    summary_rows,
    "has_counts_after_readback",
    has_counts
)


summary_rows <- add_summary(
    summary_rows,
    "has_spatial_after_readback",
    has_spatial
)


summary_rows <- add_summary(
    summary_rows,
    "has_BrNum_matched",
    has_brnum_matched
)


summary_rows <- add_summary(
    summary_rows,
    "has_visium_layer_annotation",
    has_layer_annotation
)


summary_rows <- add_summary(
    summary_rows,
    "h5ad_file_size_GB",
    round(
        file.info(OUTPUT_H5AD)$size / 1024^3,
        4
    )
)

################################################################################
# Donor summary
################################################################################

if (
    "BrNum_matched"
    %in%
    colnames(cd)
) {

    donors <- unique(
        as.character(
            cd$BrNum_matched
        )
    )

    donors <- donors[
        !is.na(donors) &
        donors != ""
    ]


    summary_rows <- add_summary(
        summary_rows,
        "n_donors",
        length(donors)
    )


    summary_rows <- add_summary(
        summary_rows,
        "donors",
        sort(donors)
    )
}

################################################################################
# Layer summary
################################################################################

if (
    "visium_layer_annotation"
    %in%
    colnames(cd)
) {

    layer_counts <- table(
        as.character(
            cd$visium_layer_annotation
        ),
        useNA = "ifany"
    )


    for (
        layer_name in names(layer_counts)
    ) {

        metric_name <- paste0(
            "layer_count_",
            layer_name
        )


        summary_rows <- add_summary(
            summary_rows,
            metric_name,
            as.integer(
                layer_counts[
                    layer_name
                ]
            )
        )
    }
}

################################################################################
# Save summary
################################################################################

summary_df <- do.call(
    rbind,
    summary_rows
)


write.csv(
    summary_df,
    SUMMARY_CSV,
    row.names = FALSE
)


cat(
    "\nSummary:\n"
)

print(
    summary_df
)


cat(
    "\nSaved summary CSV:\n"
)

cat(
    SUMMARY_CSV,
    "\n"
)

################################################################################
# Done
################################################################################

section("DONE")

cat(
    "Visium H5AD:\n"
)

cat(
    OUTPUT_H5AD,
    "\n"
)


cat(
    "\nXenium H5AD:\n"
)

cat(
    XENIUM_H5AD,
    "\n"
)


cat(
    "\nConversion summary:\n"
)

cat(
    SUMMARY_CSV,
    "\n"
)

