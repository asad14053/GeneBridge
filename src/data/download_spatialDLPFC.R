#!/usr/bin/env Rscript

####################################################################################################
# download_spatialDLPFC.R
#
# HGCC version
#
# Purpose:
#   Download Huuki-Myers / spatialDLPFC datasets from spatialLIBD.
#
# Outputs:
#
#   1. Huuki-Myers Visium:
#      data/processed/visium/spatialDLPFC_Visium_sce_light.RDS
#
#      Contains:
#        - counts
#        - logcounts
#        - rowData
#        - colData
#        - spatial coordinates
#
#      Does NOT contain SpatialExperiment image objects.
#
#   2. Huuki-Myers snRNA-seq:
#      data/processed/snrnaseq/sce_DLPFC_annotated/se.rds
#      data/processed/snrnaseq/sce_DLPFC_annotated/assays.h5
#
#   3. Summary:
#      outputs/huuki_myers/tables/huuki_download_summary.txt
#
#
# Run:
#
#   Rscript src/data/download_spatialDLPFC.R
#
# Force rerun:
#
#   Rscript src/data/download_spatialDLPFC.R --force
#
####################################################################################################


####################################################################################################
# 0. Basic HGCC settings
####################################################################################################

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

args <- commandArgs(trailingOnly = TRUE)

force_run <- "--force" %in% args


####################################################################################################
# Project paths
####################################################################################################

project_root <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

r_lib <- file.path(
    project_root,
    "software",
    "R_libs"
)

bio_cache <- file.path(
    project_root,
    "data",
    "cache",
    "bioconductor"
)

xdg_cache <- file.path(
    project_root,
    "data",
    "cache",
    "xdg"
)


####################################################################################################
# Create directories
####################################################################################################

dir.create(
    project_root,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    r_lib,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    bio_cache,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    xdg_cache,
    recursive = TRUE,
    showWarnings = FALSE
)


####################################################################################################
# R library path
####################################################################################################

.libPaths(
    c(
        r_lib,
        .libPaths()
    )
)


####################################################################################################
# Bioconductor / cache settings
####################################################################################################

Sys.setenv(
    XDG_CACHE_HOME = xdg_cache,
    BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "FALSE",
    EXPERIMENT_HUB_CACHE = bio_cache,
    ANNOTATION_HUB_CACHE = bio_cache,
    BIOCFILECACHE_CACHE = bio_cache
)

options(
    repos = c(
        CRAN = "https://cloud.r-project.org"
    ),
    timeout = 100000
)


####################################################################################################
# Output directories
####################################################################################################

visium_dir <- file.path(
    project_root,
    "data",
    "processed",
    "visium"
)

snrna_dir <- file.path(
    project_root,
    "data",
    "processed",
    "snrnaseq"
)

out_table_dir <- file.path(
    project_root,
    "outputs",
    "huuki_myers",
    "tables"
)


dir.create(
    visium_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    snrna_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    out_table_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


####################################################################################################
# Output files
####################################################################################################

visium_rds <- file.path(
    visium_dir,
    "spatialDLPFC_Visium_sce_light.RDS"
)

sce_extract_dir <- file.path(
    snrna_dir,
    "sce_DLPFC_annotated"
)

se_rds <- file.path(
    sce_extract_dir,
    "se.rds"
)

assays_h5 <- file.path(
    sce_extract_dir,
    "assays.h5"
)

summary_path <- file.path(
    out_table_dir,
    "huuki_download_summary.txt"
)


####################################################################################################
# Helper functions
####################################################################################################

section <- function(title) {

    cat(
        "\n",
        paste(
            rep("=", 100),
            collapse = ""
        ),
        "\n",
        sep = ""
    )

    cat(
        title,
        "\n"
    )

    cat(
        paste(
            rep("=", 100),
            collapse = ""
        ),
        "\n",
        sep = ""
    )
}


file_size_gb <- function(path) {

    if (!file.exists(path)) {
        return(NA_real_)
    }

    file.info(path)$size / 1e9
}


####################################################################################################
# 1. Install and load packages
####################################################################################################

section(
    "0. HGCC setup"
)

cat(
    "Project root:\n",
    project_root,
    "\n"
)

cat(
    "R library path:\n",
    r_lib,
    "\n"
)

cat(
    "Bioconductor cache:\n",
    bio_cache,
    "\n"
)

cat(
    "Force run:\n",
    force_run,
    "\n"
)

cat(
    "\n.libPaths():\n"
)

print(
    .libPaths()
)


####################################################################################################
# Install BiocManager if missing
####################################################################################################

if (!requireNamespace(
    "BiocManager",
    quietly = TRUE
)) {

    install.packages(
        "BiocManager",
        lib = r_lib,
        repos = "https://cloud.r-project.org"
    )
}


####################################################################################################
# Required packages
####################################################################################################

pkgs <- c(
    "spatialLIBD",
    "SpatialExperiment",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "HDF5Array",
    "rhdf5",
    "S4Vectors"
)


for (pkg in pkgs) {

    if (!requireNamespace(
        pkg,
        quietly = TRUE
    )) {

        cat(
            "\nInstalling missing package:",
            pkg,
            "\n"
        )

        BiocManager::install(
            pkg,
            ask = FALSE,
            update = FALSE,
            lib = r_lib
        )
    }
}


####################################################################################################
# Load packages
####################################################################################################

suppressPackageStartupMessages({

    library(spatialLIBD)

    library(SpatialExperiment)

    library(SingleCellExperiment)

    library(SummarizedExperiment)

    library(HDF5Array)

    library(rhdf5)

    library(S4Vectors)

})


####################################################################################################
# Explicit Bioconductor cache locations
####################################################################################################

if (requireNamespace(
    "ExperimentHub",
    quietly = TRUE
)) {

    try(
        ExperimentHub::setExperimentHubOption(
            "CACHE",
            bio_cache
        ),
        silent = TRUE
    )
}


if (requireNamespace(
    "AnnotationHub",
    quietly = TRUE
)) {

    try(
        AnnotationHub::setAnnotationHubOption(
            "CACHE",
            bio_cache
        ),
        silent = TRUE
    )
}


####################################################################################################
# 2. Huuki-Myers Visium
####################################################################################################

section(
    "1. Huuki-Myers Visium SpatialExperiment"
)


####################################################################################################
# Check existing output
####################################################################################################

visium_rds_ready <-

    file.exists(visium_rds) &&

    file.info(visium_rds)$size > 0 &&

    !force_run


if (visium_rds_ready) {

    cat(
        "Visium RDS already exists. Skipping Visium fetch/save:\n"
    )

    cat(
        visium_rds,
        "\n"
    )

    cat(
        "Size GB:\n"
    )

    print(
        file_size_gb(visium_rds)
    )


    cat(
        "\nLoading existing Visium RDS...\n"
    )

    sce_visium <- readRDS(
        visium_rds
    )


} else {


    ################################################################################################
    # Remove partial output
    ################################################################################################

    if (file.exists(visium_rds)) {

        cat(
            "Removing old/partial Visium RDS:\n"
        )

        cat(
            visium_rds,
            "\n"
        )

        file.remove(
            visium_rds
        )
    }


    ################################################################################################
    # Fetch Visium
    ################################################################################################

    cat(
        "Fetching Visium with spatialLIBD::fetch_data(type = 'spatialDLPFC_Visium')\n"
    )


    spe <- spatialLIBD::fetch_data(
        type = "spatialDLPFC_Visium"
    )


    ################################################################################################
    # Inspect original object
    ################################################################################################

    cat(
        "\nOriginal Visium object:\n"
    )

    print(
        spe
    )


    cat(
        "\nDimensions:\n"
    )

    print(
        dim(spe)
    )


    cat(
        "\nAssays:\n"
    )

    print(
        assayNames(spe)
    )


    cat(
        "\nrowData columns:\n"
    )

    print(
        colnames(
            rowData(spe)
        )
    )


    cat(
        "\nFirst 50 colData columns:\n"
    )

    print(
        colnames(
            colData(spe)
        )[
            1:min(
                50,
                ncol(
                    colData(spe)
                )
            )
        ]
    )


    cat(
        "\nspatialCoords columns:\n"
    )

    print(
        colnames(
            spatialCoords(spe)
        )
    )


    cat(
        "\nImage data in original SpatialExperiment:\n"
    )

    print(
        imgData(spe)
    )


    ################################################################################################
    # Build clean SingleCellExperiment manually
    #
    # IMPORTANT:
    #
    # We do NOT use:
    #
    #     as(spe, "SingleCellExperiment")
    #
    # and we do NOT directly construct a SummarizedExperiment with
    # assays that still contain conflicting dimnames.
    #
    # Instead:
    #
    #   1. Save the correct gene names.
    #   2. Save the correct spot names.
    #   3. Extract assays without top-level dimnames.
    #   4. Explicitly remove assay dimnames.
    #   5. Construct SingleCellExperiment from assays only.
    #   6. Restore gene and spot names.
    #   7. Add rowData and colData.
    #   8. Add spatial coordinates.
    #
    ################################################################################################

    section(
        "2. Building full clean Visium SingleCellExperiment without images"
    )


    ################################################################################################
    # Save names
    ################################################################################################

    gene_names <- as.character(
        rownames(spe)
    )

    spot_names <- as.character(
        colnames(spe)
    )


    ################################################################################################
    # Check names
    ################################################################################################

    if (is.null(gene_names)) {

        gene_names <- paste0(
            "gene_",
            seq_len(
                nrow(spe)
            )
        )
    }


    if (is.null(spot_names)) {

        spot_names <- paste0(
            "spot_",
            seq_len(
                ncol(spe)
            )
        )
    }


    ################################################################################################
    # Make names unique
    ################################################################################################

    gene_names <- make.unique(
        gene_names
    )

    spot_names <- make.unique(
        spot_names
    )


    ################################################################################################
    # Save metadata
    ################################################################################################

    rd <- rowData(
        spe
    )

    cd <- colData(
        spe
    )


    ################################################################################################
    # Force metadata rownames to exactly match new object
    ################################################################################################

    rownames(rd) <- gene_names

    rownames(cd) <- spot_names


    ################################################################################################
    # Save spatial coordinates
    ################################################################################################

    spatial_mat <- as.matrix(
        spatialCoords(spe)
    )


    cat(
        "Spatial matrix dimensions:\n"
    )

    print(
        dim(spatial_mat)
    )


    ################################################################################################
    # Extract counts
    ################################################################################################

    cat(
        "Adding counts assay...\n"
    )


    counts_mat <- assay(
        spe,
        "counts",
        withDimnames = FALSE
    )


    ################################################################################################
    # IMPORTANT FIX
    #
    # Remove assay-level dimnames completely.
    #
    # This prevents:
    #
    # "the rownames and colnames of the supplied assay(s)
    #  must be NULL or identical"
    #
    ################################################################################################

    dimnames(counts_mat) <- list(
        NULL,
        NULL
    )


    ################################################################################################
    # Extract logcounts
    ################################################################################################

    cat(
        "Adding logcounts assay...\n"
    )


    logcounts_mat <- assay(
        spe,
        "logcounts",
        withDimnames = FALSE
    )


    dimnames(logcounts_mat) <- list(
        NULL,
        NULL
    )


    ################################################################################################
    # Validate dimensions before constructing object
    ################################################################################################

    cat(
        "\nCounts dimensions:\n"
    )

    print(
        dim(counts_mat)
    )


    cat(
        "\nLogcounts dimensions:\n"
    )

    print(
        dim(logcounts_mat)
    )


    cat(
        "\nrowData dimensions:\n"
    )

    print(
        dim(rd)
    )


    cat(
        "\ncolData dimensions:\n"
    )

    print(
        dim(cd)
    )


    if (
        nrow(counts_mat) != nrow(spe) ||
        ncol(counts_mat) != ncol(spe)
    ) {

        stop(
            "Counts dimensions do not match original SpatialExperiment."
        )
    }


    if (
        nrow(logcounts_mat) != nrow(spe) ||
        ncol(logcounts_mat) != ncol(spe)
    ) {

        stop(
            "Logcounts dimensions do not match original SpatialExperiment."
        )
    }


    if (
        nrow(rd) != nrow(spe)
    ) {

        stop(
            "rowData dimensions do not match number of genes."
        )
    }


    if (
        nrow(cd) != ncol(spe)
    ) {

        stop(
            "colData dimensions do not match number of spots."
        )
    }


    ################################################################################################
    # Construct SCE FROM ASSAYS ONLY
    #
    # This is the key fix.
    ################################################################################################

    cat(
        "\nConstructing clean SingleCellExperiment...\n"
    )


    sce_visium <- SingleCellExperiment::SingleCellExperiment(

        assays = S4Vectors::SimpleList(

            counts = counts_mat,

            logcounts = logcounts_mat

        )

    )


    ################################################################################################
    # Restore gene and spot names
    ################################################################################################

    cat(
        "Restoring gene and spot names...\n"
    )


    rownames(
        sce_visium
    ) <- gene_names


    colnames(
        sce_visium
    ) <- spot_names


    ################################################################################################
    # Add rowData
    ################################################################################################

    cat(
        "Adding rowData...\n"
    )


    rowData(
        sce_visium
    ) <- rd


    ################################################################################################
    # Add colData
    ################################################################################################

    cat(
        "Adding colData...\n"
    )


    colData(
        sce_visium
    ) <- cd


    ################################################################################################
    # Add spatial coordinates
    ################################################################################################

    if (

        !is.null(spatial_mat) &&

        nrow(spatial_mat) == ncol(sce_visium) &&

        ncol(spatial_mat) >= 2

    ) {

        rownames(
            spatial_mat
        ) <- spot_names


        SingleCellExperiment::reducedDim(

            sce_visium,

            "spatial"

        ) <- spatial_mat


        cat(
            "Spatial coordinates added to reducedDim(sce_visium, 'spatial').\n"
        )


    } else {

        warning(
            "Spatial coordinates were not added because dimensions did not match."
        )

    }


    ################################################################################################
    # Remove inherited/extra metadata
    #
    # The new object was manually constructed, so image data are already excluded.
    ################################################################################################

    metadata(
        sce_visium
    ) <- list()


    ################################################################################################
    # Validate final object
    ################################################################################################

    cat(
        "\nClean Visium SingleCellExperiment:\n"
    )

    print(
        sce_visium
    )


    cat(
        "\nDimensions:\n"
    )

    print(
        dim(sce_visium)
    )


    cat(
        "\nAssays:\n"
    )

    print(
        assayNames(sce_visium)
    )


    cat(
        "\nrowData columns:\n"
    )

    print(
        colnames(
            rowData(sce_visium)
        )
    )


    cat(
        "\nFirst 50 colData columns:\n"
    )

    print(
        colnames(
            colData(sce_visium)
        )[
            1:min(
                50,
                ncol(
                    colData(sce_visium)
                )
            )
        ]
    )


    cat(
        "\nReduced dimensions:\n"
    )

    print(
        reducedDimNames(sce_visium)
    )


    ################################################################################################
    # Final consistency checks
    ################################################################################################

    stopifnot(

        nrow(
            sce_visium
        ) == length(
            gene_names
        ),

        ncol(
            sce_visium
        ) == length(
            spot_names
        ),

        identical(
            rownames(
                sce_visium
            ),
            rownames(
                rowData(
                    sce_visium
                )
            )
        ),

        identical(
            colnames(
                sce_visium
            ),
            rownames(
                colData(
                    sce_visium
                )
            )
        )

    )


    cat(
        "\nFinal dimname consistency checks passed.\n"
    )


    ################################################################################################
    # Save Visium RDS
    ################################################################################################

    section(
        "3. Saving Huuki-Myers Visium RDS"
    )


    cat(
        "Saving Visium RDS to:\n"
    )

    cat(
        visium_rds,
        "\n"
    )


    saveRDS(

        sce_visium,

        file = visium_rds,

        compress = FALSE

    )


    ################################################################################################
    # Verify save
    ################################################################################################

    if (!file.exists(visium_rds)) {

        stop(
            "Visium RDS was not created."
        )

    }


    cat(
        "\nSaved Visium RDS successfully.\n"
    )


    cat(
        "File:\n"
    )

    cat(
        visium_rds,
        "\n"
    )


    cat(
        "Size GB:\n"
    )

    print(
        file_size_gb(
            visium_rds
        )
    )


    ################################################################################################
    # Free original object
    ################################################################################################

    rm(
        spe,
        counts_mat,
        logcounts_mat,
        rd,
        cd,
        spatial_mat
    )

    gc()

}


####################################################################################################
# 3. Huuki-Myers snRNA-seq
####################################################################################################

section(
    "4. Huuki-Myers snRNA-seq"
)


####################################################################################################
# Check whether extracted data already exist
####################################################################################################

snrna_ready <-

    file.exists(se_rds) &&

    file.exists(assays_h5) &&

    !force_run


if (snrna_ready) {

    assays_h5_size <- file.info(
        assays_h5
    )$size


    if (

        !is.na(assays_h5_size) &&

        assays_h5_size >= 3e9

    ) {

        cat(
            "snRNA-seq extracted files already exist and look complete.\n"
        )

        cat(
            "Skipping snRNA-seq download/extraction.\n"
        )


        cat(
            "\nse.rds:\n",
            se_rds,
            "\n"
        )


        cat(
            "\nassays.h5:\n",
            assays_h5,
            "\n"
        )


        cat(
            "\nassays.h5 size GB:\n"
        )

        print(
            assays_h5_size / 1e9
        )


    } else {

        cat(
            "Existing assays.h5 looks too small.\n"
        )

        cat(
            "Will refetch and re-extract snRNA-seq.\n"
        )

        snrna_ready <- FALSE

    }

}


####################################################################################################
# Download snRNA-seq if necessary
####################################################################################################

if (!snrna_ready) {


    cat(
        "Fetching snRNA-seq with spatialLIBD::fetch_data(type = 'spatialDLPFC_snRNAseq')\n"
    )


    sce_path_zip <- spatialLIBD::fetch_data(

        type = "spatialDLPFC_snRNAseq"

    )


    cat(
        "\nsnRNA-seq zip path:\n"
    )

    cat(
        sce_path_zip,
        "\n"
    )


    if (!file.exists(sce_path_zip)) {

        stop(
            "snRNA-seq zip file does not exist: ",
            sce_path_zip
        )

    }


    cat(
        "\nZip file size GB:\n"
    )

    print(
        file.info(
            sce_path_zip
        )$size / 1e9
    )


    ################################################################################################
    # Extract with Linux unzip
    ################################################################################################

    section(
        "5. Extracting snRNA-seq using Linux unzip"
    )


    unzip_bin <- Sys.which(
        "unzip"
    )


    if (unzip_bin == "") {

        stop(
            "Linux unzip was not found on HGCC."
        )

    }


    cat(
        "Using unzip:\n"
    )

    cat(
        unzip_bin,
        "\n"
    )


    ################################################################################################
    # Remove previous extracted directories
    ################################################################################################

    macosx_dir <- file.path(
        snrna_dir,
        "__MACOSX"
    )


    if (dir.exists(sce_extract_dir)) {

        cat(
            "Removing old extracted snRNA-seq directory:\n"
        )

        cat(
            sce_extract_dir,
            "\n"
        )

        unlink(
            sce_extract_dir,
            recursive = TRUE
        )

    }


    if (dir.exists(macosx_dir)) {

        unlink(
            macosx_dir,
            recursive = TRUE
        )

    }


    ################################################################################################
    # Extract
    ################################################################################################

    unzip_log_path <- file.path(
        snrna_dir,
        "unzip_snrna.log"
    )


    unzip_log <- system2(

        command = unzip_bin,

        args = c(
            "-o",
            sce_path_zip,
            "-d",
            snrna_dir
        ),

        stdout = TRUE,

        stderr = TRUE

    )


    writeLines(
        unzip_log,
        unzip_log_path
    )


    unzip_status <- attr(
        unzip_log,
        "status"
    )


    if (is.null(unzip_status)) {

        unzip_status <- 0

    }


    if (unzip_status != 0) {

        stop(
            "Linux unzip failed. Check log: ",
            unzip_log_path
        )

    }


    cat(
        "Unzip completed successfully.\n"
    )


    cat(
        "Unzip log:\n"
    )

    cat(
        unzip_log_path,
        "\n"
    )

}


####################################################################################################
# 4. Validate snRNA-seq files
####################################################################################################

section(
    "6. Validating extracted snRNA-seq files"
)


cat(
    "Expected snRNA-seq files:\n"
)


cat(
    se_rds,
    "\n"
)


cat(
    assays_h5,
    "\n"
)


####################################################################################################
# Validate se.rds
####################################################################################################

if (!file.exists(se_rds)) {

    stop(
        "Missing se.rds: ",
        se_rds
    )

}


####################################################################################################
# Validate assays.h5
####################################################################################################

if (!file.exists(assays_h5)) {

    stop(
        "Missing assays.h5: ",
        assays_h5
    )

}


assays_h5_size <- file.info(
    assays_h5
)$size


cat(
    "\nassays.h5 size GB:\n"
)


print(
    assays_h5_size / 1e9
)


if (

    is.na(assays_h5_size) ||

    assays_h5_size < 3e9

) {

    stop(

        "assays.h5 looks truncated. Expected about 3.9 GB, but got ",

        assays_h5_size / 1e9,

        " GB."

    )

}


####################################################################################################
# Inspect HDF5
####################################################################################################

cat(
    "\nTesting assays.h5 with rhdf5::h5ls():\n"
)


print(

    rhdf5::h5ls(

        assays_h5,

        recursive = FALSE

    )

)


####################################################################################################
# 5. Load snRNA-seq object
####################################################################################################

section(
    "7. Loading snRNA-seq object using readRDS from inside extracted directory"
)


####################################################################################################
# Important:
#
# The HDF5-backed RDS references assays.h5 relatively.
# Therefore read se.rds while working inside sce_DLPFC_annotated.
####################################################################################################

old_wd <- getwd()


setwd(
    sce_extract_dir
)


sce_snrna <- readRDS(
    "se.rds"
)


cat(
    "\nsnRNA-seq object:\n"
)


print(
    sce_snrna
)


cat(
    "\nAssays:\n"
)


print(
    assayNames(
        sce_snrna
    )
)


cat(
    "\nDimensions:\n"
)


print(
    dim(
        sce_snrna
    )
)


cat(
    "\nFirst 50 colData columns:\n"
)


print(

    colnames(

        colData(
            sce_snrna
        )

    )[
        1:min(
            50,
            ncol(
                colData(
                    sce_snrna
                )
            )
        )
    ]

)


cat(
    "\nrowData columns:\n"
)


print(

    colnames(

        rowData(
            sce_snrna
        )

    )

)


####################################################################################################
# Test counts access
####################################################################################################

cat(
    "\nTesting counts assay access:\n"
)


print(

    assay(

        sce_snrna,

        "counts",

        withDimnames = FALSE

    )[1:5, 1:5]

)


####################################################################################################
# Test logcounts access
####################################################################################################

cat(
    "\nTesting logcounts assay access:\n"
)


print(

    assay(

        sce_snrna,

        "logcounts",

        withDimnames = FALSE

    )[1:5, 1:5]

)


####################################################################################################
# Restore original working directory
####################################################################################################

setwd(
    old_wd
)


####################################################################################################
# 6. Save summary
####################################################################################################

section(
    "8. Saving download summary"
)


summary_lines <- c(

    "Huuki-Myers / spatialDLPFC download summary",

    "",

    "Cluster:",

    "HGCC",

    "",

    "Project root:",

    project_root,

    "",

    "R library path:",

    r_lib,

    "",

    "Bioconductor cache:",

    bio_cache,

    "",

    "Huuki Visium RDS:",

    visium_rds,

    "",

    paste0(
        "Huuki Visium RDS exists: ",
        file.exists(
            visium_rds
        )
    ),

    paste0(
        "Huuki Visium RDS size GB: ",
        round(
            file_size_gb(
                visium_rds
            ),
            3
        )
    ),

    "",

    "Huuki snRNA-seq outputs:",

    se_rds,

    assays_h5,

    "",

    paste0(
        "snRNA-seq assays.h5 size GB: ",
        round(
            assays_h5_size / 1e9,
            3
        )
    ),

    "",

    "Important notes:",

    "1. HGCC version using BeEGFS project paths.",

    "2. R libraries are stored under project_root/software/R_libs.",

    "3. Bioconductor cache is stored under project_root/data/cache/bioconductor.",

    "4. Original Visium SpatialExperiment images are excluded from the saved RDS.",

    "5. Visium counts and logcounts assays are retained.",

    "6. Visium rowData and colData are retained.",

    "7. Visium spatial coordinates are stored in reducedDim(sce_visium, 'spatial').",

    "8. Visium assay dimnames are removed during SCE construction and restored afterward to avoid SummarizedExperiment dimname mismatch.",

    "9. snRNA-seq zip is extracted using Linux unzip.",

    "10. snRNA-seq se.rds is loaded from inside sce_DLPFC_annotated so assays.h5 relative links resolve correctly.",

    "",

    "Visium SingleCellExperiment:",

    capture.output(
        print(
            sce_visium
        )
    ),

    "",

    "snRNA-seq SingleCellExperiment:",

    capture.output(
        print(
            sce_snrna
        )
    )

)


dir.create(

    dirname(
        summary_path
    ),

    recursive = TRUE,

    showWarnings = FALSE

)


writeLines(

    summary_lines,

    summary_path

)


cat(
    "Saved summary to:\n"
)


cat(
    summary_path,
    "\n"
)


####################################################################################################
# DONE
####################################################################################################

section(
    "DONE"
)


cat(
    "\nMain outputs:\n"
)


cat(
    "\nHuuki Visium RDS:\n"
)


cat(
    visium_rds,
    "\n"
)


cat(
    "\nHuuki snRNA-seq se.rds:\n"
)


cat(
    se_rds,
    "\n"
)


cat(
    "\nHuuki snRNA-seq assays.h5:\n"
)


cat(
    assays_h5,
    "\n"
)


cat(
    "\nDownload summary:\n"
)


cat(
    summary_path,
    "\n"
)


cat(
    "\nAll tasks completed successfully.\n"
)