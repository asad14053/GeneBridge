#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SummarizedExperiment)
})

# ================================================================================================
# 01b_compare_original_pseudobulk_to_paper.R
#
# PURPOSE
# -------
# Validate the GeneBridge donor × spatial-domain pseudobulk against the
# published Xenium pseudobulk and summarize ENVI OOF-imputed 300-gene
# pseudobulk performance.
#
# Comparisons:
#
#   A. Published paper Xenium 300
#          vs
#      GeneBridge original Xenium 300
#
#   B. GeneBridge original Xenium 300
#          vs
#      ENVI OOF-imputed Xenium 300
#
#
# OUTPUT LEVELS
# -------------
#   1. Per donor × domain : 161 rows
#   2. Per donor          : 23 rows
#   3. Per domain         : 7 rows
#   4. Overall            : 1 row
#
# IMPORTANT
# ---------
# Paper object contains N24.
# GeneBridge analysis contains N23.
# Br6432 is therefore removed from the paper object before comparison.
# ================================================================================================


# ================================================================================================
# CONFIGURATION
# ================================================================================================

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"


PAPER_RDS <- file.path(
    ROOT,
    "data/reference/paper_xenium",
    "spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"
)


OUR_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "original_xenium_300gene_donor_layer_pseudobulk_counts.csv.gz"
)


OUR_META <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)


# NOTE:
# Legacy filename says "measured300".
# In this comparison this file is expected to represent the ENVI OOF-imputed
# predictions for the 300 benchmark genes.
#
# The script explicitly checks whether this matrix is numerically identical
# to original Xenium and raises a warning if so.

ENVI_OOF300_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
)


ENVI_OOF300_META <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_donor_layer_pseudobulk_metadata.csv"
)


OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "paper_raw_comparison"
)


dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


# ================================================================================================
# HELPERS
# ================================================================================================

section <- function(title) {

    cat(
        "\n",
        paste(rep("=", 100), collapse = ""),
        "\n",
        title,
        "\n",
        paste(rep("=", 100), collapse = ""),
        "\n",
        sep = ""
    )
}


safe_cor <- function(
    x,
    y,
    method = "pearson"
) {

    if (
        length(x) < 2 ||
        length(unique(x)) < 2 ||
        length(unique(y)) < 2
    ) {
        return(NA_real_)
    }

    suppressWarnings(
        cor(
            x,
            y,
            method = method
        )
    )
}


read_pb_matrix <- function(path) {

    x <- read.csv(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )


    first_col <- colnames(x)[1]
    first_values <- as.character(x[[1]])


    looks_like_id <- (
        first_col %in% c(
            "pseudobulk_id",
            "X",
            "Unnamed..0",
            "Unnamed: 0",
            "row.names"
        )
        ||
        any(grepl("__spd", first_values))
        ||
        any(grepl("::spd", first_values))
    )


    if (looks_like_id) {

        rownames(x) <- first_values
        x[[1]] <- NULL
    }


    m <- as.matrix(x)

    storage.mode(m) <- "double"

    return(m)
}


write_csv_gz <- function(
    x,
    path
) {

    con <- gzfile(
        path,
        open = "wt"
    )

    write.csv(
        x,
        con,
        row.names = FALSE
    )

    close(con)
}


check_file <- function(path) {

    if (!file.exists(path)) {
        stop(
            "Missing required file:\n",
            path
        )
    }
}


# ================================================================================================
# CHECK INPUT FILES
# ================================================================================================

section("CHECK INPUT FILES")


for (
    f in c(
        PAPER_RDS,
        OUR_COUNTS,
        OUR_META,
        ENVI_OOF300_COUNTS,
        ENVI_OOF300_META
    )
) {

    check_file(f)

    cat(
        "FOUND:",
        f,
        "\n"
    )
}


# ================================================================================================
# STEP 1 — LOAD PAPER PSEUDOBULK
# ================================================================================================

section("1. LOAD PAPER PSEUDOBULK")


paper <- readRDS(
    PAPER_RDS
)


cat(
    "Paper object class:\n"
)

print(
    class(paper)
)


cat(
    "\nPaper dimensions:\n"
)

print(
    dim(paper)
)


paper_counts_gene_by_pb <- assay(
    paper,
    "counts"
)


cat(
    "\nPaper raw count matrix:",
    nrow(paper_counts_gene_by_pb),
    "genes x",
    ncol(paper_counts_gene_by_pb),
    "pseudobulks\n"
)


paper_meta <- as.data.frame(
    colData(paper)
)


paper_rowdata <- as.data.frame(
    rowData(paper)
)


cat(
    "\nPaper colData columns:\n"
)

print(
    colnames(paper_meta)
)


cat(
    "\nPaper rowData columns:\n"
)

print(
    colnames(paper_rowdata)
)


# ================================================================================================
# STEP 2 — DETERMINE PAPER GENE SYMBOLS
# ================================================================================================

section("2. DETERMINE PAPER GENE SYMBOLS")


gene_candidates <- c(
    "Symbol",
    "symbol",
    "gene",
    "gene_name",
    "gene_symbol"
)


gene_col <- gene_candidates[
    gene_candidates %in%
    colnames(paper_rowdata)
][1]


if (is.na(gene_col)) {

    stop(
        "Could not identify gene-symbol column in paper rowData."
    )
}


paper_genes <- as.character(
    paper_rowdata[[gene_col]]
)


if (anyDuplicated(paper_genes)) {

    stop(
        "Duplicated paper gene symbols detected."
    )
}


cat(
    "Using paper rowData gene column:",
    gene_col,
    "\n"
)


# Convert paper matrix:
#
# genes × pseudobulks
#        ↓
# pseudobulks × genes

paper_matrix_all <- t(
    as.matrix(
        paper_counts_gene_by_pb
    )
)


storage.mode(
    paper_matrix_all
) <- "double"


colnames(
    paper_matrix_all
) <- paper_genes


# ================================================================================================
# STEP 3 — LOAD GENEBRIDGE ORIGINAL PSEUDOBULK
# ================================================================================================

section("3. LOAD GENEBRIDGE ORIGINAL PSEUDOBULK")


our_matrix <- read_pb_matrix(
    OUR_COUNTS
)


our_meta <- read.csv(
    OUR_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat(
    "Our raw count matrix:",
    nrow(our_matrix),
    "pseudobulks x",
    ncol(our_matrix),
    "genes\n"
)


cat(
    "Our metadata:",
    nrow(our_meta),
    "rows\n"
)


required_meta <- c(
    "pseudobulk_id",
    "BrNum",
    "predictions_smooth",
    "layer_annotation",
    "n_cells"
)


missing_meta <- setdiff(
    required_meta,
    colnames(our_meta)
)


if (length(missing_meta) > 0) {

    stop(
        "Missing GeneBridge metadata columns: ",
        paste(
            missing_meta,
            collapse = ", "
        )
    )
}


# Align original count matrix rows to metadata.

if (
    !is.null(rownames(our_matrix)) &&
    all(
        our_meta$pseudobulk_id
        %in%
        rownames(our_matrix)
    )
) {

    our_matrix <- our_matrix[
        match(
            our_meta$pseudobulk_id,
            rownames(our_matrix)
        ),
        ,
        drop = FALSE
    ]

} else {

    if (
        nrow(our_matrix)
        !=
        nrow(our_meta)
    ) {

        stop(
            "Original matrix and metadata row numbers differ."
        )
    }


    rownames(
        our_matrix
    ) <- our_meta$pseudobulk_id
}


our_genes <- colnames(
    our_matrix
)


# ================================================================================================
# STEP 4 — ALIGN 300 GENES
# ================================================================================================

section("4. ALIGN 300 GENES")


shared_genes <- intersect(
    paper_genes,
    our_genes
)


paper_only <- setdiff(
    paper_genes,
    our_genes
)


our_only <- setdiff(
    our_genes,
    paper_genes
)


cat(
    "Paper genes     :",
    length(paper_genes),
    "\n"
)

cat(
    "GeneBridge genes:",
    length(our_genes),
    "\n"
)

cat(
    "Shared genes    :",
    length(shared_genes),
    "\n"
)

cat(
    "Paper-only genes:",
    length(paper_only),
    "\n"
)

cat(
    "Our-only genes  :",
    length(our_only),
    "\n"
)


if (
    length(shared_genes) != 300 ||
    length(paper_only) != 0 ||
    length(our_only) != 0
) {

    stop(
        "The 300-gene universes are not identical."
    )
}


# Use GeneBridge gene order everywhere.

paper_matrix_all <- paper_matrix_all[
    ,
    our_genes,
    drop = FALSE
]


# ================================================================================================
# STEP 5 — SUBSET PAPER N24 TO GENEBRIDGE N23
# ================================================================================================

section("5. SUBSET PAPER N24 → GENEBRIDGE N23")


paper_donors_total <- sort(
    unique(
        as.character(
            paper_meta$BrNum
        )
    )
)


our_donors <- sort(
    unique(
        as.character(
            our_meta$BrNum
        )
    )
)


paper_excluded_donors <- setdiff(
    paper_donors_total,
    our_donors
)


keep <- (
    as.character(
        paper_meta$BrNum
    )
    %in%
    our_donors
)


paper_meta_n23 <- paper_meta[
    keep,
    ,
    drop = FALSE
]


paper_matrix <- paper_matrix_all[
    keep,
    ,
    drop = FALSE
]


cat(
    "Paper donors total     :",
    length(paper_donors_total),
    "\n"
)

cat(
    "GeneBridge donors      :",
    length(our_donors),
    "\n"
)

cat(
    "Paper donors excluded from N23 comparison:",
    paste(
        paper_excluded_donors,
        collapse = ", "
    ),
    "\n"
)

cat(
    "Paper pseudobulks after N23 subset:",
    nrow(paper_matrix),
    "\n"
)

cat(
    "Our pseudobulks:",
    nrow(our_matrix),
    "\n"
)


# ================================================================================================
# STEP 6 — ALIGN DONOR × SPATIAL-DOMAIN KEYS
# ================================================================================================

section("6. ALIGN DONOR × SPATIAL DOMAIN")


paper_key <- paste(
    as.character(
        paper_meta_n23$BrNum
    ),
    as.character(
        paper_meta_n23$predictions_smooth
    ),
    sep = "::"
)


our_key <- paste(
    as.character(
        our_meta$BrNum
    ),
    as.character(
        our_meta$predictions_smooth
    ),
    sep = "::"
)


if (anyDuplicated(paper_key)) {

    stop(
        "Duplicated paper donor-domain keys."
    )
}


if (anyDuplicated(our_key)) {

    stop(
        "Duplicated GeneBridge donor-domain keys."
    )
}


paper_only_keys <- setdiff(
    paper_key,
    our_key
)


our_only_keys <- setdiff(
    our_key,
    paper_key
)


cat(
    "Paper N23 donor-domain keys:",
    length(paper_key),
    "\n"
)

cat(
    "Our donor-domain keys      :",
    length(our_key),
    "\n"
)

cat(
    "Paper-only keys:",
    length(paper_only_keys),
    "\n"
)

cat(
    "Our-only keys  :",
    length(our_only_keys),
    "\n"
)


if (
    length(paper_only_keys) > 0 ||
    length(our_only_keys) > 0
) {

    stop(
        "Donor-domain key mismatch."
    )
}


# Reorder paper rows to exactly match GeneBridge.

paper_idx <- match(
    our_key,
    paper_key
)


paper_matrix <- paper_matrix[
    paper_idx,
    ,
    drop = FALSE
]


paper_meta_n23 <- paper_meta_n23[
    paper_idx,
    ,
    drop = FALSE
]


# ================================================================================================
# STEP 7 — COMPARE SPATIAL LABELS
# ================================================================================================

section("7. COMPARE SPATIAL LABELS")


spd_identical <- all(
    as.character(
        paper_meta_n23$predictions_smooth
    )
    ==
    as.character(
        our_meta$predictions_smooth
    )
)


layer_identical <- TRUE


if (
    "domain_annotations"
    %in%
    colnames(paper_meta_n23)
) {

    layer_identical <- all(
        as.character(
            paper_meta_n23$domain_annotations
        )
        ==
        as.character(
            our_meta$layer_annotation
        )
    )
}


spatial_labels_identical <- (
    spd_identical &&
    layer_identical
)


cat(
    "predictions_smooth identical:",
    spd_identical,
    "\n"
)

cat(
    "Layer/domain annotation identical:",
    layer_identical,
    "\n"
)

cat(
    "Spatial-domain labels identical:",
    spatial_labels_identical,
    "\n"
)


# ================================================================================================
# STEP 8 — COMPARE NUMBER OF CELLS
# ================================================================================================

section("8. COMPARE NUMBER OF CELLS PER PSEUDOBULK")


if (
    "ncells"
    %in%
    colnames(paper_meta_n23)
) {

    paper_ncells <- as.numeric(
        paper_meta_n23$ncells
    )

} else if (
    "n_cells"
    %in%
    colnames(paper_meta_n23)
) {

    paper_ncells <- as.numeric(
        paper_meta_n23$n_cells
    )

} else {

    paper_ncells <- rep(
        NA_real_,
        nrow(paper_meta_n23)
    )
}


our_ncells <- as.numeric(
    our_meta$n_cells
)


ncells_table <- data.frame(

    pseudobulk_id =
        our_meta$pseudobulk_id,

    BrNum =
        our_meta$BrNum,

    predictions_smooth =
        our_meta$predictions_smooth,

    layer_annotation =
        our_meta$layer_annotation,

    paper_ncells =
        paper_ncells,

    our_ncells =
        our_ncells,

    difference =
        our_ncells -
        paper_ncells,

    exact_match =
        our_ncells ==
        paper_ncells,

    stringsAsFactors = FALSE
)


ncells_exact_n <- sum(
    ncells_table$exact_match,
    na.rm = TRUE
)


ncells_identical <- all(
    ncells_table$exact_match,
    na.rm = TRUE
)


cat(
    "Exact n_cells matches:",
    ncells_exact_n,
    "/",
    nrow(ncells_table),
    "\n"
)

cat(
    "All donor-domain n_cells identical:",
    ncells_identical,
    "\n"
)


write.csv(
    ncells_table,
    file.path(
        OUT,
        "paper_vs_genebridge_ncells.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# STEP 9 — PAPER vs ORIGINAL RAW PSEUDOBULK COUNTS
# ================================================================================================

section("9. PAPER vs GENEBRIDGE ORIGINAL RAW PSEUDOBULK COUNTS")


if (
    !identical(
        dim(paper_matrix),
        dim(our_matrix)
    )
) {

    stop(
        "Paper and original matrices have different dimensions."
    )
}


count_difference <- (
    our_matrix -
    paper_matrix
)


exact_matrix <- (
    our_matrix ==
    paper_matrix
)


compared_entries <- length(
    exact_matrix
)


exact_count_entries <- sum(
    exact_matrix
)


mismatched_count_entries <- (
    compared_entries -
    exact_count_entries
)


exact_count_fraction <- mean(
    exact_matrix
)


max_abs_count_difference <- max(
    abs(
        count_difference
    )
)


total_abs_count_difference <- sum(
    abs(
        count_difference
    )
)


pearson_raw_counts <- safe_cor(
    as.vector(
        paper_matrix
    ),
    as.vector(
        our_matrix
    ),
    method = "pearson"
)


spearman_raw_counts <- safe_cor(
    as.vector(
        paper_matrix
    ),
    as.vector(
        our_matrix
    ),
    method = "spearman"
)


raw_counts_identical <- (
    mismatched_count_entries == 0
)


cat(
    "Entries compared:",
    compared_entries,
    "\n"
)

cat(
    "Exact entries:",
    exact_count_entries,
    "\n"
)

cat(
    "Mismatched entries:",
    mismatched_count_entries,
    "\n"
)

cat(
    "Exact fraction:",
    sprintf(
        "%.8f",
        exact_count_fraction
    ),
    "\n"
)

cat(
    "Maximum absolute count difference:",
    max_abs_count_difference,
    "\n"
)

cat(
    "Total absolute count difference:",
    total_abs_count_difference,
    "\n"
)

cat(
    "Pearson:",
    sprintf(
        "%.12f",
        pearson_raw_counts
    ),
    "\n"
)

cat(
    "Spearman:",
    sprintf(
        "%.12f",
        spearman_raw_counts
    ),
    "\n"
)

cat(
    "\nRAW COUNTS IDENTICAL:",
    raw_counts_identical,
    "\n"
)


# ================================================================================================
# STEP 10 — PER-PSEUDOBULK PAPER vs ORIGINAL SUMMARY
# ================================================================================================

section("10. PAPER vs ORIGINAL — PER DONOR × DOMAIN")


paper_original_pb <- data.frame(

    pseudobulk_id =
        our_meta$pseudobulk_id,

    BrNum =
        our_meta$BrNum,

    Dx =
        if (
            "Dx"
            %in%
            colnames(our_meta)
        ) {
            our_meta$Dx
        } else {
            NA
        },

    predictions_smooth =
        our_meta$predictions_smooth,

    layer_annotation =
        our_meta$layer_annotation,

    paper_ncells =
        paper_ncells,

    our_ncells =
        our_ncells,

    paper_total_counts =
        rowSums(
            paper_matrix
        ),

    our_total_counts =
        rowSums(
            our_matrix
        ),

    total_difference =
        rowSums(
            count_difference
        ),

    n_gene_mismatches =
        rowSums(
            !exact_matrix
        ),

    max_abs_gene_difference =
        apply(
            abs(
                count_difference
            ),
            1,
            max
        ),

    exact_match =
        rowSums(
            !exact_matrix
        )
        == 0,

    stringsAsFactors = FALSE
)


write.csv(
    paper_original_pb,
    file.path(
        OUT,
        "paper_vs_genebridge_per_pseudobulk_summary.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# STEP 11 — PER-GENE PAPER vs ORIGINAL SUMMARY
# ================================================================================================

section("11. PAPER vs ORIGINAL — PER-GENE SUMMARY")


paper_original_gene <- data.frame(

    gene =
        colnames(
            our_matrix
        ),

    n_pseudobulk_mismatches =
        colSums(
            !exact_matrix
        ),

    total_abs_difference =
        colSums(
            abs(
                count_difference
            )
        ),

    max_abs_difference =
        apply(
            abs(
                count_difference
            ),
            2,
            max
        ),

    stringsAsFactors = FALSE
)


write.csv(
    paper_original_gene,
    file.path(
        OUT,
        "paper_vs_genebridge_per_gene_summary.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# STEP 12 — SAVE RAW MISMATCHES
# ================================================================================================

section("12. SAVE PAPER vs ORIGINAL RAW MISMATCHES")


mismatch_idx <- which(
    !exact_matrix,
    arr.ind = TRUE
)


if (
    nrow(mismatch_idx) > 0
) {

    mismatch_table <- data.frame(

        pseudobulk_id =
            our_meta$pseudobulk_id[
                mismatch_idx[, 1]
            ],

        BrNum =
            our_meta$BrNum[
                mismatch_idx[, 1]
            ],

        predictions_smooth =
            our_meta$predictions_smooth[
                mismatch_idx[, 1]
            ],

        layer_annotation =
            our_meta$layer_annotation[
                mismatch_idx[, 1]
            ],

        gene =
            colnames(
                our_matrix
            )[
                mismatch_idx[, 2]
            ],

        paper_count =
            paper_matrix[
                mismatch_idx
            ],

        our_count =
            our_matrix[
                mismatch_idx
            ],

        difference =
            count_difference[
                mismatch_idx
            ],

        stringsAsFactors = FALSE
    )


    write_csv_gz(
        mismatch_table,
        file.path(
            OUT,
            "paper_vs_genebridge_raw_count_mismatches.csv.gz"
        )
    )


    cat(
        "Saved",
        nrow(mismatch_table),
        "mismatched entries.\n"
    )

} else {

    cat(
        "No raw-count mismatches.\n"
    )
}


# ================================================================================================
# STEP 13 — LOAD ENVI OOF-IMPUTED 300 PSEUDOBULK
# ================================================================================================

section("13. LOAD ENVI OOF-IMPUTED 300 PSEUDOBULK")


envi_matrix <- read_pb_matrix(
    ENVI_OOF300_COUNTS
)


envi_meta <- read.csv(
    ENVI_OOF300_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat(
    "ENVI OOF-300 matrix:",
    nrow(envi_matrix),
    "pseudobulks x",
    ncol(envi_matrix),
    "genes\n"
)

cat(
    "ENVI OOF-300 metadata:",
    nrow(envi_meta),
    "rows\n"
)


# ================================================================================================
# STEP 14 — ALIGN ENVI GENES
# ================================================================================================

section("14. ALIGN ENVI OOF-300 GENES")


if (
    !setequal(
        colnames(envi_matrix),
        our_genes
    )
) {

    stop(
        "ENVI OOF-300 gene universe does not match original 300 genes."
    )
}


envi_matrix <- envi_matrix[
    ,
    our_genes,
    drop = FALSE
]


cat(
    "ENVI genes aligned:",
    ncol(envi_matrix),
    "/ 300\n"
)


# ================================================================================================
# STEP 15 — ALIGN ENVI DONOR × DOMAIN ROWS
# ================================================================================================

section("15. ALIGN ENVI OOF-300 DONOR × DOMAIN")


if (
    "pseudobulk_id"
    %in%
    colnames(envi_meta)
) {

    envi_pb_ids <- as.character(
        envi_meta$pseudobulk_id
    )

} else {

    envi_pb_ids <- paste(
        as.character(
            envi_meta$BrNum
        ),
        as.character(
            envi_meta$predictions_smooth
        ),
        sep = "::"
    )
}


if (
    !is.null(
        rownames(envi_matrix)
    )
    &&
    all(
        our_meta$pseudobulk_id
        %in%
        rownames(envi_matrix)
    )
) {

    envi_idx <- match(
        our_meta$pseudobulk_id,
        rownames(envi_matrix)
    )

} else {

    envi_key <- paste(
        as.character(
            envi_meta$BrNum
        ),
        as.character(
            envi_meta$predictions_smooth
        ),
        sep = "::"
    )


    envi_idx <- match(
        our_key,
        envi_key
    )
}


if (
    any(
        is.na(
            envi_idx
        )
    )
) {

    stop(
        "Could not align all ENVI donor-domain pseudobulks."
    )
}


envi_matrix <- envi_matrix[
    envi_idx,
    ,
    drop = FALSE
]


if (
    !identical(
        dim(envi_matrix),
        dim(our_matrix)
    )
) {

    stop(
        "ENVI and original matrices have different dimensions."
    )
}


cat(
    "Aligned ENVI donor-domain pseudobulks:",
    nrow(envi_matrix),
    "\n"
)


# ================================================================================================
# STEP 16 — ORIGINAL vs ENVI OOF-300
# ================================================================================================

section("16. ORIGINAL vs ENVI OOF-IMPUTED 300")


envi_difference <- (
    envi_matrix -
    our_matrix
)


envi_exact <- (
    envi_matrix ==
    our_matrix
)


envi_exact_entries <- sum(
    envi_exact
)


envi_total_entries <- length(
    envi_exact
)


envi_exact_fraction <- mean(
    envi_exact
)


envi_pearson <- safe_cor(
    as.vector(
        our_matrix
    ),
    as.vector(
        envi_matrix
    ),
    method = "pearson"
)


envi_spearman <- safe_cor(
    as.vector(
        our_matrix
    ),
    as.vector(
        envi_matrix
    ),
    method = "spearman"
)


envi_mae <- mean(
    abs(
        envi_difference
    )
)


envi_rmse <- sqrt(
    mean(
        envi_difference^2
    )
)


cat(
    "Entries compared:",
    envi_total_entries,
    "\n"
)

cat(
    "Exact entries:",
    envi_exact_entries,
    "\n"
)

cat(
    "Exact fraction:",
    sprintf(
        "%.12f",
        envi_exact_fraction
    ),
    "\n"
)

cat(
    "Pearson:",
    sprintf(
        "%.12f",
        envi_pearson
    ),
    "\n"
)

cat(
    "Spearman:",
    sprintf(
        "%.12f",
        envi_spearman
    ),
    "\n"
)

cat(
    "MAE:",
    sprintf(
        "%.12f",
        envi_mae
    ),
    "\n"
)

cat(
    "RMSE:",
    sprintf(
        "%.12f",
        envi_rmse
    ),
    "\n"
)


# Critical sanity check.

envi_identical_to_original <- (
    envi_exact_fraction == 1
)


if (
    envi_identical_to_original
) {

    cat(
        "\n",
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
        "WARNING: ENVI OOF-300 matrix is EXACTLY IDENTICAL to original Xenium.\n",
        "All ", envi_total_entries, " donor-domain × gene values are identical.\n",
        "If this file is intended to contain OOF predictions, verify its upstream source.\n",
        "Legacy filename currently used:\n",
        ENVI_OOF300_COUNTS,
        "\n",
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
        sep = ""
    )
}


# ================================================================================================
# STEP 17 — 161-ROW DONOR × DOMAIN THREE-WAY SUMMARY
# ================================================================================================

section("17. DONOR × DOMAIN SUMMARY — 161 ROWS")


donor_domain_summary <- data.frame(

    pseudobulk_id =
        our_meta$pseudobulk_id,

    BrNum =
        our_meta$BrNum,

    Dx =
        if (
            "Dx"
            %in%
            colnames(our_meta)
        ) {
            our_meta$Dx
        } else {
            NA
        },

    predictions_smooth =
        our_meta$predictions_smooth,

    layer_annotation =
        our_meta$layer_annotation,

    paper_ncells =
        paper_ncells,

    original_ncells =
        our_ncells,

    paper_total_counts =
        rowSums(
            paper_matrix
        ),

    original_total_counts =
        rowSums(
            our_matrix
        ),

    envi_oof300_total =
        rowSums(
            envi_matrix
        ),

    paper_vs_original_n_mismatch_genes =
        rowSums(
            paper_matrix
            !=
            our_matrix
        ),

    paper_vs_original_max_abs_diff =
        apply(
            abs(
                paper_matrix -
                our_matrix
            ),
            1,
            max
        ),

    original_vs_envi_oof300_pearson =
        vapply(
            seq_len(
                nrow(our_matrix)
            ),
            function(i) {

                safe_cor(
                    our_matrix[i, ],
                    envi_matrix[i, ],
                    method = "pearson"
                )
            },
            numeric(1)
        ),

    original_vs_envi_oof300_spearman =
        vapply(
            seq_len(
                nrow(our_matrix)
            ),
            function(i) {

                safe_cor(
                    our_matrix[i, ],
                    envi_matrix[i, ],
                    method = "spearman"
                )
            },
            numeric(1)
        ),

    original_vs_envi_oof300_MAE =
        rowMeans(
            abs(
                envi_difference
            )
        ),

    original_vs_envi_oof300_RMSE =
        sqrt(
            rowMeans(
                envi_difference^2
            )
        ),

    stringsAsFactors = FALSE
)


write.csv(
    donor_domain_summary,
    file.path(
        OUT,
        "paper_vs_original_vs_ENVI_OOF300_donor_domain_summary.csv"
    ),
    row.names = FALSE
)


cat(
    "Saved 161-row donor × domain summary.\n"
)


# ================================================================================================
# STEP 17B — PAPER / ORIGINAL / ENVI OOF-300 SIDE-BY-SIDE TABLE
# ================================================================================================

section(
    "17B. PAPER vs ORIGINAL vs ENVI OOF-300 — SIDE BY SIDE"
)


side_by_side <- data.frame(

    Donor =
        as.character(
            donor_domain_summary$BrNum
        ),

    Dx =
        as.character(
            donor_domain_summary$Dx
        ),

    Domain =
        as.character(
            donor_domain_summary$layer_annotation
        ),

    Paper =
        donor_domain_summary$paper_total_counts,

    Original =
        donor_domain_summary$original_total_counts,

    ENVI_OOF300 =
        donor_domain_summary$envi_oof300_total,

    Paper_minus_Original =
        donor_domain_summary$paper_total_counts -
        donor_domain_summary$original_total_counts,

    ENVI_minus_Original =
        donor_domain_summary$envi_oof300_total -
        donor_domain_summary$original_total_counts,

    Paper_mismatch_genes =
        donor_domain_summary$paper_vs_original_n_mismatch_genes,

    ENVI_Pearson =
        donor_domain_summary$original_vs_envi_oof300_pearson,

    ENVI_MAE =
        donor_domain_summary$original_vs_envi_oof300_MAE,

    ENVI_RMSE =
        donor_domain_summary$original_vs_envi_oof300_RMSE,

    stringsAsFactors = FALSE
)


# Save full side-by-side table.

side_by_side_file <- file.path(
    OUT,
    "paper_original_ENVI_OOF300_side_by_side.csv"
)


write.csv(
    side_by_side,
    side_by_side_file,
    row.names = FALSE
)


cat(
    "Saved side-by-side table:\n",
    side_by_side_file,
    "\n\n",
    sep = ""
)


# --------------------------------------------------------------------------------
# Terminal-friendly fixed-width table.
#
# This avoids R wrapping a wide data.frame into several separate blocks.
# --------------------------------------------------------------------------------

cat(
    sprintf(
        "%-8s %-4s %-7s %14s %14s %14s %12s %12s %9s\n",
        "Donor",
        "Dx",
        "Domain",
        "Paper",
        "Original",
        "ENVI_OOF300",
        "P-O",
        "E-O",
        "MisGenes"
    )
)


cat(
    paste(
        rep("-", 108),
        collapse = ""
    ),
    "\n",
    sep = ""
)


for (
    i in seq_len(
        nrow(side_by_side)
    )
) {

    cat(
        sprintf(
            "%-8s %-4s %-7s %14.0f %14.0f %14.0f %12.0f %12.0f %9d\n",

            side_by_side$Donor[i],

            side_by_side$Dx[i],

            side_by_side$Domain[i],

            side_by_side$Paper[i],

            side_by_side$Original[i],

            side_by_side$ENVI_OOF300[i],

            side_by_side$Paper_minus_Original[i],

            side_by_side$ENVI_minus_Original[i],

            side_by_side$Paper_mismatch_genes[i]
        )
    )
}


# ================================================================================================
# STEP 17C — ONLY NON-IDENTICAL ROWS
# ================================================================================================

section(
    "17C. NON-IDENTICAL DONOR × DOMAIN ROWS"
)


non_identical <- side_by_side[
    (
        side_by_side$Paper_minus_Original != 0
        |
        side_by_side$ENVI_minus_Original != 0
    ),
    ,
    drop = FALSE
]


if (
    nrow(non_identical) == 0
) {

    cat(
        "All three matrices have identical donor-domain total counts.\n"
    )

} else {

    cat(
        "Number of donor-domain rows with any total-count difference:",
        nrow(non_identical),
        "\n\n"
    )


    cat(
        sprintf(
            "%-8s %-4s %-7s %14s %14s %14s %12s %12s %9s\n",
            "Donor",
            "Dx",
            "Domain",
            "Paper",
            "Original",
            "ENVI_OOF300",
            "P-O",
            "E-O",
            "MisGenes"
        )
    )


    cat(
        paste(
            rep("-", 108),
            collapse = ""
        ),
        "\n",
        sep = ""
    )


    for (
        i in seq_len(
            nrow(non_identical)
        )
    ) {

        cat(
            sprintf(
                "%-8s %-4s %-7s %14.0f %14.0f %14.0f %12.0f %12.0f %9d\n",

                non_identical$Donor[i],

                non_identical$Dx[i],

                non_identical$Domain[i],

                non_identical$Paper[i],

                non_identical$Original[i],

                non_identical$ENVI_OOF300[i],

                non_identical$Paper_minus_Original[i],

                non_identical$ENVI_minus_Original[i],

                non_identical$Paper_mismatch_genes[i]
            )
        )
    }
}


write.csv(
    non_identical,
    file.path(
        OUT,
        "paper_original_ENVI_OOF300_nonidentical_rows.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# STEP 17D — THREE-WAY CORRELATION / ERROR SUMMARY
# ================================================================================================

section(
    "17D. THREE-WAY COMPARISON SUMMARY"
)


cat(
    sprintf(
        "%-28s %18s\n",
        "Metric",
        "Value"
    )
)


cat(
    paste(
        rep("-", 48),
        collapse = ""
    ),
    "\n",
    sep = ""
)


cat(
    sprintf(
        "%-28s %18d\n",
        "Total pseudobulks",
        nrow(side_by_side)
    )
)


cat(
    sprintf(
        "%-28s %18d\n",
        "Genes",
        ncol(our_matrix)
    )
)


cat(
    sprintf(
        "%-28s %18d\n",
        "Gene x PB entries",
        length(our_matrix)
    )
)


cat(
    sprintf(
        "%-28s %18d\n",
        "Paper=Original PBs",
        sum(
            donor_domain_summary$
                paper_vs_original_n_mismatch_genes
            == 0
        )
    )
)


cat(
    sprintf(
        "%-28s %17.8f%%\n",
        "Paper=Original entries",
        100 *
        mean(
            paper_matrix
            ==
            our_matrix
        )
    )
)


cat(
    sprintf(
        "%-28s %18.12f\n",
        "Paper vs Original Pearson",
        safe_cor(
            as.vector(paper_matrix),
            as.vector(our_matrix),
            "pearson"
        )
    )
)


cat(
    sprintf(
        "%-28s %18.12f\n",
        "Original vs ENVI Pearson",
        safe_cor(
            as.vector(our_matrix),
            as.vector(envi_matrix),
            "pearson"
        )
    )
)


cat(
    sprintf(
        "%-28s %18.12f\n",
        "Original vs ENVI Spearman",
        safe_cor(
            as.vector(our_matrix),
            as.vector(envi_matrix),
            "spearman"
        )
    )
)


cat(
    sprintf(
        "%-28s %18.12f\n",
        "Original vs ENVI MAE",
        mean(
            abs(
                envi_matrix -
                our_matrix
            )
        )
    )
)


cat(
    sprintf(
        "%-28s %18.12f\n",
        "Original vs ENVI RMSE",
        sqrt(
            mean(
                (
                    envi_matrix -
                    our_matrix
                )^2
            )
        )
    )
)



# ================================================================================================
# STEP 18 — 23-ROW DONOR SUMMARY
# ================================================================================================

section("18. DONOR SUMMARY — 23 ROWS")


donor_groups <- split(
    seq_len(
        nrow(our_matrix)
    ),
    as.character(
        our_meta$BrNum
    )
)


donor_summary <- do.call(
    rbind,
    lapply(
        names(
            donor_groups
        ),
        function(donor) {

            idx <- donor_groups[[donor]]


            original_vector <- as.vector(
                our_matrix[
                    idx,
                    ,
                    drop = FALSE
                ]
            )


            envi_vector <- as.vector(
                envi_matrix[
                    idx,
                    ,
                    drop = FALSE
                ]
            )


            paper_vector <- as.vector(
                paper_matrix[
                    idx,
                    ,
                    drop = FALSE
                ]
            )


            data.frame(

                BrNum =
                    donor,

                Dx =
                    if (
                        "Dx"
                        %in%
                        colnames(our_meta)
                    ) {
                        unique(
                            our_meta$Dx[idx]
                        )[1]
                    } else {
                        NA
                    },

                n_domains =
                    length(idx),

                gene_domain_entries =
                    length(
                        original_vector
                    ),

                paper_vs_original_exact_fraction =
                    mean(
                        paper_vector
                        ==
                        original_vector
                    ),

                original_vs_envi_oof300_pearson =
                    safe_cor(
                        original_vector,
                        envi_vector,
                        "pearson"
                    ),

                original_vs_envi_oof300_spearman =
                    safe_cor(
                        original_vector,
                        envi_vector,
                        "spearman"
                    ),

                original_vs_envi_oof300_MAE =
                    mean(
                        abs(
                            envi_vector -
                            original_vector
                        )
                    ),

                original_vs_envi_oof300_RMSE =
                    sqrt(
                        mean(
                            (
                                envi_vector -
                                original_vector
                            )^2
                        )
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)


rownames(
    donor_summary
) <- NULL


write.csv(
    donor_summary,
    file.path(
        OUT,
        "paper_vs_original_vs_ENVI_OOF300_donor_summary.csv"
    ),
    row.names = FALSE
)


print(
    donor_summary,
    row.names = FALSE
)


# ================================================================================================
# STEP 19 — 7-ROW DOMAIN SUMMARY
# ================================================================================================

section("19. DOMAIN SUMMARY — 7 ROWS")


domain_groups <- split(
    seq_len(
        nrow(our_matrix)
    ),
    as.character(
        our_meta$layer_annotation
    )
)


domain_summary <- do.call(
    rbind,
    lapply(
        names(
            domain_groups
        ),
        function(domain) {

            idx <- domain_groups[[domain]]


            original_vector <- as.vector(
                our_matrix[
                    idx,
                    ,
                    drop = FALSE
                ]
            )


            envi_vector <- as.vector(
                envi_matrix[
                    idx,
                    ,
                    drop = FALSE
                ]
            )


            paper_vector <- as.vector(
                paper_matrix[
                    idx,
                    ,
                    drop = FALSE
                ]
            )


            data.frame(

                layer =
                    domain,

                donors =
                    length(idx),

                paper_vs_original_exact_pseudobulks =
                    sum(
                        rowSums(
                            paper_matrix[
                                idx,
                                ,
                                drop = FALSE
                            ]
                            !=
                            our_matrix[
                                idx,
                                ,
                                drop = FALSE
                            ]
                        )
                        == 0
                    ),

                paper_vs_original_exact_fraction =
                    mean(
                        paper_vector
                        ==
                        original_vector
                    ),

                original_vs_envi_oof300_pearson =
                    safe_cor(
                        original_vector,
                        envi_vector,
                        "pearson"
                    ),

                original_vs_envi_oof300_spearman =
                    safe_cor(
                        original_vector,
                        envi_vector,
                        "spearman"
                    ),

                original_vs_envi_oof300_MAE =
                    mean(
                        abs(
                            envi_vector -
                            original_vector
                        )
                    ),

                original_vs_envi_oof300_RMSE =
                    sqrt(
                        mean(
                            (
                                envi_vector -
                                original_vector
                            )^2
                        )
                    ),

                stringsAsFactors = FALSE
            )
        }
    )
)


rownames(
    domain_summary
) <- NULL


write.csv(
    domain_summary,
    file.path(
        OUT,
        "paper_vs_original_vs_ENVI_OOF300_domain_summary.csv"
    ),
    row.names = FALSE
)


print(
    domain_summary,
    row.names = FALSE
)


# ================================================================================================
# STEP 20 — ONE-ROW OVERALL SUMMARY
# ================================================================================================

section("20. OVERALL SUMMARY")


overall_summary <- data.frame(

    paper_donors_total =
        length(
            paper_donors_total
        ),

    compared_donors =
        length(
            our_donors
        ),

    paper_excluded_donors =
        paste(
            paper_excluded_donors,
            collapse = ";"
        ),

    compared_domains =
        length(
            unique(
                our_meta$layer_annotation
            )
        ),

    compared_pseudobulks =
        nrow(
            our_matrix
        ),

    compared_genes =
        ncol(
            our_matrix
        ),

    compared_entries =
        length(
            our_matrix
        ),

    paper_vs_original_exact_entries =
        exact_count_entries,

    paper_vs_original_exact_fraction =
        exact_count_fraction,

    paper_vs_original_exact_pseudobulks =
        sum(
            paper_original_pb$exact_match
        ),

    paper_vs_original_pearson =
        pearson_raw_counts,

    paper_vs_original_spearman =
        spearman_raw_counts,

    paper_vs_original_max_abs_diff =
        max_abs_count_difference,

    paper_vs_original_total_abs_diff =
        total_abs_count_difference,

    original_vs_envi_oof300_exact_entries =
        envi_exact_entries,

    original_vs_envi_oof300_exact_fraction =
        envi_exact_fraction,

    original_vs_envi_oof300_pearson =
        envi_pearson,

    original_vs_envi_oof300_spearman =
        envi_spearman,

    original_vs_envi_oof300_MAE =
        envi_mae,

    original_vs_envi_oof300_RMSE =
        envi_rmse,

    spatial_labels_identical =
        spatial_labels_identical,

    ncells_identical =
        ncells_identical,

    raw_counts_identical =
        raw_counts_identical,

    envi_exactly_identical_to_original =
        envi_identical_to_original,

    stringsAsFactors = FALSE
)


write.csv(
    overall_summary,
    file.path(
        OUT,
        "paper_vs_original_vs_ENVI_OOF300_overall_summary.csv"
    ),
    row.names = FALSE
)


print(
    overall_summary,
    row.names = FALSE
)


# ================================================================================================
# STEP 21 — FINAL INTERPRETATION
# ================================================================================================

section("21. FINAL INTERPRETATION")


cat(
    "\nPAPER → GENEBRIDGE ORIGINAL\n"
)

cat(
    "----------------------------------------\n"
)

cat(
    "Genes aligned              :",
    length(shared_genes),
    "/ 300\n"
)

cat(
    "Donor-domain keys aligned  :",
    length(our_key),
    "/ 161\n"
)

cat(
    "Spatial labels identical   :",
    spatial_labels_identical,
    "\n"
)

cat(
    "n_cells exact pseudobulks  :",
    ncells_exact_n,
    "/ 161\n"
)

cat(
    "Raw-count exact entries    :",
    exact_count_entries,
    "/",
    compared_entries,
    "\n"
)

cat(
    "Raw-count exact fraction   :",
    sprintf(
        "%.6f%%",
        100 * exact_count_fraction
    ),
    "\n"
)

cat(
    "Exact pseudobulks          :",
    sum(
        paper_original_pb$exact_match
    ),
    "/ 161\n"
)


cat(
    "\nORIGINAL → ENVI OOF-300\n"
)

cat(
    "----------------------------------------\n"
)

cat(
    "Exact entries              :",
    envi_exact_entries,
    "/",
    envi_total_entries,
    "\n"
)

cat(
    "Pearson                    :",
    sprintf(
        "%.12f",
        envi_pearson
    ),
    "\n"
)

cat(
    "Spearman                   :",
    sprintf(
        "%.12f",
        envi_spearman
    ),
    "\n"
)

cat(
    "MAE                        :",
    sprintf(
        "%.12f",
        envi_mae
    ),
    "\n"
)

cat(
    "RMSE                       :",
    sprintf(
        "%.12f",
        envi_rmse
    ),
    "\n"
)


cat(
    "\nSTATUS\n"
)

cat(
    "----------------------------------------\n"
)


if (
    raw_counts_identical
) {

    cat(
        "PAPER RECONSTRUCTION: PASS — exact raw pseudobulk reproduction.\n"
    )

} else {

    cat(
        "PAPER RECONSTRUCTION: NEAR-MATCH — inspect localized cell-universe differences.\n"
    )
}


if (
    envi_identical_to_original
) {

    cat(
        "ENVI OOF-300: WARNING — matrix is exactly identical to original Xenium.\n"
    )

    cat(
        "Verify that ENVI_OOF300_COUNTS really contains OOF predictions.\n"
    )

} else {

    cat(
        "ENVI OOF-300: predictions differ from original ground truth as expected.\n"
    )
}


cat(
    "\nOutput directory:\n",
    OUT,
    "\n",
    sep = ""
)


cat(
    "\nFINAL STATUS: COMPLETE\n"
)
