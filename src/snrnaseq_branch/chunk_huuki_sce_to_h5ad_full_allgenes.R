#!/usr/bin/env Rscript

############################:contentReference[oaicite:1]{index=1}##############
# Huuki-Myers snRNA-seq SCE -> H5AD chunks
#
# Strategy:
#   - Keep ALL genes
#   - Split CELLS into chunks
#   - Read HDF5-backed SCE
#   - Explicitly transpose each assay:
#
#       SCE     = genes x cells
#       AnnData = cells x genes
#
#   - counts    -> AnnData X
#   - logcounts -> AnnData layers["logcounts"]
#   - colData   -> AnnData obs
#   - rowData   -> AnnData var
#
# No zellkonverter.
# No Python installation.
#
# Usage:
#
# Test:
# Rscript src/snrnaseq_branch/chunk_huuki_sce_to_h5ad_full_allgenes.R test 50 200
#
# Full:
# Rscript src/snrnaseq_branch/chunk_huuki_sce_to_h5ad_full_allgenes.R full 200 0
################################################################################


################################################################################
# 0. Settings
################################################################################

options(stringsAsFactors = FALSE)
options(width = 120)

Sys.setenv(
    HDF5_USE_FILE_LOCKING = "FALSE"
)


################################################################################
# 1. Arguments
################################################################################

args <- commandArgs(
    trailingOnly = TRUE
)

MODE <- if (
    length(args) >= 1
) {
    args[[1]]
} else {
    "full"
}

CELLS_PER_CHUNK <- if (
    length(args) >= 2
) {
    as.integer(args[[2]])
} else {
    200
}

MAX_CELLS <- if (
    length(args) >= 3
) {
    as.integer(args[[3]])
} else {
    0
}


if (!MODE %in% c("test", "full")) {
    stop("MODE must be either 'test' or 'full'")
}

if (
    is.na(CELLS_PER_CHUNK) ||
    CELLS_PER_CHUNK <= 0
) {
    stop("CELLS_PER_CHUNK must be > 0")
}

if (
    is.na(MAX_CELLS) ||
    MAX_CELLS < 0
) {
    stop("MAX_CELLS must be >= 0")
}


################################################################################
# 2. Paths
################################################################################

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

SCE_DIR <- file.path(
    PROJECT_ROOT,
    "data",
    "processed",
    "snrnaseq",
    "sce_DLPFC_annotated"
)

SE_RDS <- file.path(
    SCE_DIR,
    "se.rds"
)

ASSAYS_H5 <- file.path(
    SCE_DIR,
    "assays.h5"
)

OUT_DIR <- file.path(
    SCE_DIR,
    paste0(
        "huuki_h5ad_chunks_allgenes_",
        MODE,
        "_maxcells",
        MAX_CELLS,
        "_chunk",
        CELLS_PER_CHUNK
    )
)

MANIFEST <- file.path(
    OUT_DIR,
    "chunk_manifest.csv"
)

LOG_FILE <- file.path(
    OUT_DIR,
    "chunk_conversion_log.txt"
)


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


################################################################################
# 3. Logging
################################################################################

sink(
    LOG_FILE,
    split = TRUE
)


cat(
    "============================================================\n"
)

cat(
    "Huuki SCE/RDS -> H5AD chunks, all genes\n"
)

cat(
    "Manual genes x cells -> cells x genes conversion\n"
)

cat(
    "============================================================\n"
)

cat(
    "PROJECT_ROOT:    ", PROJECT_ROOT, "\n"
)

cat(
    "MODE:            ", MODE, "\n"
)

cat(
    "CELLS_PER_CHUNK: ", CELLS_PER_CHUNK, "\n"
)

cat(
    "MAX_CELLS:       ", MAX_CELLS, "\n"
)

cat(
    "SE_RDS:          ", SE_RDS, "\n"
)

cat(
    "ASSAYS_H5:       ", ASSAYS_H5, "\n"
)

cat(
    "OUT_DIR:         ", OUT_DIR, "\n"
)

cat(
    "MANIFEST:        ", MANIFEST, "\n"
)

cat(
    "LOG_FILE:        ", LOG_FILE, "\n"
)


################################################################################
# 4. Check inputs
################################################################################

if (!file.exists(SE_RDS)) {
    stop(
        "Missing se.rds: ",
        SE_RDS
    )
}

if (!file.exists(ASSAYS_H5)) {
    stop(
        "Missing assays.h5: ",
        ASSAYS_H5
    )
}


cat(
    "\nInput files found.\n"
)

cat(
    "se.rds size GB: ",
    round(file.info(SE_RDS)$size / 1e9, 3),
    "\n"
)

cat(
    "assays.h5 size GB: ",
    round(file.info(ASSAYS_H5)$size / 1e9, 3),
    "\n"
)


################################################################################
# 5. Set working directory
#
# Required because se.rds references assays.h5 relatively.
################################################################################

setwd(
    SCE_DIR
)

cat(
    "\nWorking directory: ",
    getwd(),
    "\n"
)


################################################################################
# 6. Required packages
################################################################################

required_pkgs <- c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "HDF5Array",
    "DelayedArray",
    "anndataR",
    "rhdf5"
)


missing_pkgs <- required_pkgs[
    !vapply(
        required_pkgs,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (
    length(missing_pkgs) > 0
) {

    stop(
        "Missing R packages: ",
        paste(
            missing_pkgs,
            collapse = ", "
        )
    )
}


suppressPackageStartupMessages({

    library(SingleCellExperiment)

    library(SummarizedExperiment)

    library(HDF5Array)

    library(DelayedArray)

    library(anndataR)

})


cat(
    "\nRequired packages loaded.\n"
)

cat(
    "anndataR version: ",
    as.character(
        packageVersion("anndataR")
    ),
    "\n"
)


################################################################################
# 7. DelayedArray configuration
################################################################################

DelayedArray::setAutoBlockSize(
    50e6
)


################################################################################
# 8. Read SCE
################################################################################

cat(
    "\nReading Huuki SCE...\n"
)

sce <- readRDS(
    SE_RDS
)


cat(
    "\nSCE summary:\n"
)

print(
    sce
)


cat(
    "\nOriginal dimensions [genes x cells]:\n"
)

print(
    dim(sce)
)


cat(
    "\nAssays:\n"
)

print(
    assayNames(sce)
)


################################################################################
# 9. Determine assays
################################################################################

assay_names <- assayNames(
    sce
)

if (
    length(assay_names) == 0
) {
    stop(
        "No assays found."
    )
}


if (
    "counts" %in% assay_names
) {

    assay_to_use <- "counts"

} else if (
    "logcounts" %in% assay_names
) {

    assay_to_use <- "logcounts"

} else {

    assay_to_use <- assay_names[[1]]
}


other_assays <- setdiff(
    assay_names,
    assay_to_use
)


cat(
    "\nAnnData X assay: ",
    assay_to_use,
    "\n"
)


cat(
    "AnnData layers: ",
    if (
        length(other_assays) > 0
    ) {
        paste(
            other_assays,
            collapse = ", "
        )
    } else {
        "None"
    },
    "\n"
)


################################################################################
# 10. Gene names
################################################################################

rowdata <- as.data.frame(
    rowData(sce)
)


symbol_candidates <- c(
    "gene_name",
    "gene_symbol",
    "symbol",
    "gene",
    "external_gene_name"
)


symbol_cols <- intersect(
    symbol_candidates,
    colnames(rowdata)
)


if (
    length(symbol_cols) > 0
) {

    symbol_col <- symbol_cols[[1]]

    symbols <- as.character(
        rowdata[[symbol_col]]
    )

    valid <- (
        !is.na(symbols) &
        symbols != ""
    )


    cat(
        "\nGene symbol column: ",
        symbol_col,
        "\n"
    )

    cat(
        "Fraction valid: ",
        round(mean(valid), 3),
        "\n"
    )


    if (
        mean(valid) > 0.8
    ) {

        rownames(sce) <- make.unique(
            symbols
        )

        cat(
            "Using gene symbols as gene names.\n"
        )
    }
}


################################################################################
# 11. Validate names
################################################################################

if (
    is.null(rownames(sce)) ||
    any(is.na(rownames(sce))) ||
    any(rownames(sce) == "")
) {

    rownames(sce) <- paste0(
        "gene_",
        seq_len(nrow(sce))
    )
}


if (
    is.null(colnames(sce)) ||
    any(is.na(colnames(sce))) ||
    any(colnames(sce) == "")
) {

    colnames(sce) <- paste0(
        "cell_",
        seq_len(ncol(sce))
    )
}


rownames(sce) <- make.unique(
    as.character(
        rownames(sce)
    )
)

colnames(sce) <- make.unique(
    as.character(
        colnames(sce)
    )
)


cat(
    "\nFinal gene count: ",
    nrow(sce),
    "\n"
)

cat(
    "Final cell count: ",
    ncol(sce),
    "\n"
)


################################################################################
# 12. Select cells
################################################################################

all_cells <- colnames(
    sce
)


if (
    MAX_CELLS > 0
) {

    set.seed(
        42
    )

    n_keep <- min(
        MAX_CELLS,
        length(all_cells)
    )

    selected_cells <- sample(
        all_cells,
        size = n_keep,
        replace = FALSE
    )

} else {

    selected_cells <- all_cells
}


cat(
    "\nSelected cells: ",
    length(selected_cells),
    "\n"
)

cat(
    "All genes kept: ",
    nrow(sce),
    "\n"
)


################################################################################
# 13. Number of chunks
################################################################################

n_chunks <- ceiling(
    length(selected_cells) /
    CELLS_PER_CHUNK
)


cat(
    "Cells per chunk: ",
    CELLS_PER_CHUNK,
    "\n"
)

cat(
    "Number of chunks: ",
    n_chunks,
    "\n"
)


################################################################################
# 14. Manifest
################################################################################

manifest <- data.frame(

    chunk_id = integer(),

    file = character(),

    n_genes = integer(),

    n_cells = integer(),

    X_assay = character(),

    layer_assays = character(),

    file_size_mb = numeric(),

    stringsAsFactors = FALSE

)


################################################################################
# 15. Write chunks
################################################################################

for (
    i in seq_len(
        n_chunks
    )
) {


    ##########################################################################
    # Cell indices
    ##########################################################################

    start_idx <- (
        (i - 1) *
        CELLS_PER_CHUNK
    ) + 1


    end_idx <- min(
        i *
        CELLS_PER_CHUNK,
        length(selected_cells)
    )


    cells_i <- selected_cells[
        start_idx:end_idx
    ]


    ##########################################################################
    # Output
    ##########################################################################

    chunk_file <- file.path(

        OUT_DIR,

        sprintf(
            "huuki_allgenes_chunk_%04d.h5ad",
            i
        )

    )


    cat(
        "\n------------------------------------------------------------\n"
    )

    cat(
        "Writing chunk ",
        i,
        "/",
        n_chunks,
        "\n"
    )

    cat(
        "Cells: ",
        length(cells_i),
        "\n"
    )

    cat(
        "Genes: ",
        nrow(sce),
        "\n"
    )

    cat(
        "File: ",
        chunk_file,
        "\n"
    )

    cat(
        "------------------------------------------------------------\n"
    )


    ##########################################################################
    # Subset cells
    #
    # SCE orientation:
    #
    # genes x cells
    ##########################################################################

    sce_i <- sce[
        ,
        cells_i
    ]


    cat(
        "SCE dimensions [genes x cells]: ",
        nrow(sce_i),
        " x ",
        ncol(sce_i),
        "\n"
    )


    ##########################################################################
    # Cell metadata -> AnnData obs
    ##########################################################################

    obs_df <- as.data.frame(
        colData(sce_i)
    )


    rownames(
        obs_df
    ) <- colnames(
        sce_i
    )


    ##########################################################################
    # Gene metadata -> AnnData var
    ##########################################################################

    var_df <- as.data.frame(
        rowData(sce_i)
    )


    rownames(
        var_df
    ) <- rownames(
        sce_i
    )


    ##########################################################################
    # CRITICAL FIX
    #
    # SCE assay:
    #
    #     genes x cells
    #     36601 x 200
    #
    # AnnData X:
    #
    #     cells x genes
    #     200 x 36601
    #
    # Explicit transpose avoids the anndataR SCE conversion shape error.
    ##########################################################################

    cat(
        "Reading and transposing ",
        assay_to_use,
        "...\n"
    )


    X_mat <- as.matrix(

        t(

            assay(
                sce_i,
                assay_to_use,
                withDimnames = FALSE
            )

        )

    )


    rownames(
        X_mat
    ) <- colnames(
        sce_i
    )


    colnames(
        X_mat
    ) <- rownames(
        sce_i
    )


    cat(
        "AnnData X dimensions [cells x genes]: ",
        nrow(X_mat),
        " x ",
        ncol(X_mat),
        "\n"
    )


    ##########################################################################
    # Verify orientation
    ##########################################################################

    expected_shape <- c(
        length(cells_i),
        nrow(sce_i)
    )


    if (
        !identical(
            as.integer(dim(X_mat)),
            as.integer(expected_shape)
        )
    ) {

        stop(
            "Unexpected X shape. Expected ",
            paste(
                expected_shape,
                collapse = " x "
            ),
            ", got ",
            paste(
                dim(X_mat),
                collapse = " x "
            )
        )
    }


    ##########################################################################
    # Additional assays -> AnnData layers
    ##########################################################################

    layers_list <- list()


    if (
        length(other_assays) > 0
    ) {


        for (
            layer_name in other_assays
        ) {


            cat(
                "Reading and transposing layer: ",
                layer_name,
                "\n"
            )


            layer_mat <- as.matrix(

                t(

                    assay(
                        sce_i,
                        layer_name,
                        withDimnames = FALSE
                    )

                )

            )


            rownames(
                layer_mat
            ) <- colnames(
                sce_i
            )


            colnames(
                layer_mat
            ) <- rownames(
                sce_i
            )


            if (
                !identical(
                    dim(layer_mat),
                    dim(X_mat)
                )
            ) {

                stop(
                    "Layer ",
                    layer_name,
                    " has incorrect dimensions."
                )
            }


            layers_list[[
                layer_name
            ]] <- layer_mat

        }

    }


    ##########################################################################
    # Remove old output
    ##########################################################################

    if (
        file.exists(
            chunk_file
        )
    ) {

        file.remove(
            chunk_file
        )

    }


    gc()


    ##########################################################################
    # Manually construct AnnData
    #
    # This bypasses:
    #
    #   SingleCellExperiment
    #           ->
    #   as_AnnData()
    #
    # which caused the shape mismatch.
    ##########################################################################

    cat(
        "Creating AnnData object...\n"
    )


    adata <- anndataR::AnnData(

        X = X_mat,

        obs = obs_df,

        var = var_df,

        layers = layers_list

    )


    cat(
        "AnnData dimensions [cells x genes]: ",
        adata$n_obs(),
        " x ",
        adata$n_vars(),
        "\n"
    )


    ##########################################################################
    # Final validation before writing
    ##########################################################################

    if (
        adata$n_obs() != length(cells_i)
    ) {

        stop(
            "AnnData cell count mismatch."
        )

    }


    if (
        adata$n_vars() != nrow(sce_i)
    ) {

        stop(
            "AnnData gene count mismatch."
        )

    }


    ##########################################################################
    # Write H5AD
    ##########################################################################

    cat(
        "Writing H5AD...\n"
    )


    anndataR::write_h5ad(

        adata,

        chunk_file

    )


    ##########################################################################
    # Verify
    ##########################################################################

    if (
        !file.exists(
            chunk_file
        )
    ) {

        stop(
            "Chunk write failed: ",
            chunk_file
        )

    }


    chunk_size_mb <- file.info(
        chunk_file
    )$size / 1e6


    cat(
        "Chunk completed successfully.\n"
    )

    cat(
        "Chunk size MB: ",
        round(
            chunk_size_mb,
            2
        ),
        "\n"
    )


    ##########################################################################
    # Update manifest
    ##########################################################################

    manifest <- rbind(

        manifest,

        data.frame(

            chunk_id = i,

            file = chunk_file,

            n_genes = nrow(sce_i),

            n_cells = ncol(sce_i),

            X_assay = assay_to_use,

            layer_assays = paste(
                other_assays,
                collapse = ";"
            ),

            file_size_mb = round(
                chunk_size_mb,
                2
            ),

            stringsAsFactors = FALSE

        )

    )


    write.csv(

        manifest,

        MANIFEST,

        row.names = FALSE

    )


    ##########################################################################
    # Clean memory
    ##########################################################################

    rm(
        sce_i,
        obs_df,
        var_df,
        X_mat,
        layers_list,
        adata
    )


    if (
        exists(
            "layer_mat"
        )
    ) {

        rm(
            layer_mat
        )

    }


    gc()

}


################################################################################
# 16. Done
################################################################################

cat(
    "\n============================================================\n"
)

cat(
    "DONE writing all H5AD chunks.\n"
)

cat(
    "============================================================\n"
)

cat(
    "Total chunks: ",
    n_chunks,
    "\n"
)

cat(
    "Total cells: ",
    length(selected_cells),
    "\n"
)

cat(
    "Genes per chunk: ",
    nrow(sce),
    "\n"
)

cat(
    "AnnData X: ",
    assay_to_use,
    "\n"
)

cat(
    "AnnData layers: ",
    paste(
        other_assays,
        collapse = ", "
    ),
    "\n"
)

cat(
    "Manifest: ",
    MANIFEST,
    "\n"
)

cat(
    "Output directory: ",
    OUT_DIR,
    "\n"
)

cat(
    "\nConversion completed successfully.\n"
)


sink()