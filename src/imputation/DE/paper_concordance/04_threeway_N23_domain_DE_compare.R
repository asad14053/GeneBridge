#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(spatialLIBD)
    library(edgeR)
    library(limma)
})


# ================================================================================================
# PURPOSE
# ================================================================================================
#
# Canonical N23 three-way donor x spatial-domain DE comparison:
#
#   1. Published Paper Xenium pseudobulk
#   2. GeneBridge original/before-imputation Xenium pseudobulk
#   3. ENVI measured-300 preserved Xenium pseudobulk
#
# All analyses:
#
#   23 donors only
#   7 spatial domains
#   161 donor-domain pseudobulks
#   300 Xenium genes
#
# Paper-style model:
#
#   expression ~ Dx + predictions_smooth + Age + Sex + slide_id
#
# with:
#
#   block = BrNum
#
# IMPORTANT:
#
#   run_date is NOT included in the DE model because the paper's
#   02_pseudobulk_de.R does not include run_date in its DE covariates.
#
#   ENVI measured300 is a technical preservation control.
#   It is NOT OOF-imputed expression.
#
# ================================================================================================


ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"


PAPER_RDS <- file.path(
    ROOT,
    "data/reference/paper_xenium",
    "spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"
)


ORIGINAL_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "original_xenium_300gene_donor_layer_pseudobulk_counts.csv.gz"
)


ORIGINAL_META <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)


ENVI_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
)


ENVI_META <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_donor_layer_pseudobulk_metadata.csv"
)


OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "threeway_N23_domain_DE"
)


dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


# ================================================================================================
# HELPERS
# ================================================================================================


section <- function(x) {

    cat(
        "\n",
        paste(rep("=", 110), collapse = ""),
        "\n",
        x,
        "\n",
        paste(rep("=", 110), collapse = ""),
        "\n",
        sep = ""
    )
}


check_file <- function(x) {

    if (!file.exists(x)) {

        stop(
            "Missing input:\n",
            x
        )
    }

    cat(
        "FOUND: ",
        x,
        "\n",
        sep = ""
    )
}


make_key <- function(meta) {

    paste(
        as.character(meta$BrNum),
        as.character(meta$predictions_smooth),
        sep = "::"
    )
}


safe_cor <- function(
    x,
    y,
    method = "pearson"
) {

    suppressWarnings(
        cor(
            x,
            y,
            method = method,
            use = "complete.obs"
        )
    )
}


# ================================================================================================
# INPUT CHECK
# ================================================================================================


section("INPUT CHECK")


for (
    f in c(
        PAPER_RDS,
        ORIGINAL_COUNTS,
        ORIGINAL_META,
        ENVI_COUNTS,
        ENVI_META
    )
) {

    check_file(f)
}


# ================================================================================================
# LOAD GENEBRIDGE ORIGINAL
# ================================================================================================


section("LOAD GENEBRIDGE ORIGINAL N23")


original_pb <- read.csv(
    ORIGINAL_COUNTS,
    row.names = 1,
    check.names = FALSE
)


original_pb <- as.matrix(
    original_pb
)


storage.mode(
    original_pb
) <- "double"


original_meta <- read.csv(
    ORIGINAL_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (
    "pseudobulk_id"
    %in%
    colnames(original_meta)
) {

    idx <- match(
        rownames(original_pb),
        original_meta$pseudobulk_id
    )

    if (any(is.na(idx))) {

        stop(
            "Could not align original metadata to pseudobulk matrix."
        )
    }

    original_meta <- original_meta[
        idx,
        ,
        drop = FALSE
    ]
}


original_key <- make_key(
    original_meta
)


if (anyDuplicated(original_key)) {

    stop(
        "Duplicate donor-domain keys in GeneBridge original."
    )
}


rownames(original_pb) <- original_key
rownames(original_meta) <- original_key


cat(
    "Original matrix:",
    nrow(original_pb),
    "PB x",
    ncol(original_pb),
    "genes\n"
)


cat(
    "Original donors:",
    length(unique(original_meta$BrNum)),
    "\n"
)


cat(
    "Original domains:",
    length(unique(original_meta$predictions_smooth)),
    "\n"
)


# ================================================================================================
# DEFINE CANONICAL N23 UNIVERSE
# ================================================================================================


section("CANONICAL N23 UNIVERSE")


canonical_keys <- original_key


canonical_donors <- sort(
    unique(
        as.character(
            original_meta$BrNum
        )
    )
)


canonical_genes <- colnames(
    original_pb
)


cat(
    "Canonical donors:",
    length(canonical_donors),
    "\n"
)


cat(
    "Canonical PBs:",
    length(canonical_keys),
    "\n"
)


cat(
    "Canonical genes:",
    length(canonical_genes),
    "\n"
)


if (
    length(canonical_donors)
    != 23
) {

    stop(
        "Expected exactly 23 donors."
    )
}


if (
    length(canonical_keys)
    != 161
) {

    stop(
        "Expected exactly 161 donor-domain pseudobulks."
    )
}


if (
    length(canonical_genes)
    != 300
) {

    stop(
        "Expected exactly 300 Xenium genes."
    )
}


if (
    "Br6432"
    %in%
    canonical_donors
) {

    stop(
        "Br6432 is present. Canonical analysis must be N23."
    )
}


cat(
    "Br6432 excluded: TRUE\n"
)


# ================================================================================================
# LOAD ENVI MEASURED-300
# ================================================================================================


section("LOAD ENVI MEASURED-300 N23")


envi_pb <- read.csv(
    ENVI_COUNTS,
    row.names = 1,
    check.names = FALSE
)


envi_pb <- as.matrix(
    envi_pb
)


storage.mode(
    envi_pb
) <- "double"


envi_meta <- read.csv(
    ENVI_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (
    "pseudobulk_id"
    %in%
    colnames(envi_meta)
) {

    idx <- match(
        rownames(envi_pb),
        envi_meta$pseudobulk_id
    )

    if (any(is.na(idx))) {

        stop(
            "Could not align ENVI metadata."
        )
    }

    envi_meta <- envi_meta[
        idx,
        ,
        drop = FALSE
    ]
}


envi_key <- make_key(
    envi_meta
)


if (anyDuplicated(envi_key)) {

    stop(
        "Duplicate ENVI donor-domain keys."
    )
}


rownames(envi_pb) <- envi_key
rownames(envi_meta) <- envi_key


if (
    !setequal(
        canonical_keys,
        rownames(envi_pb)
    )
) {

    stop(
        "ENVI donor-domain universe differs from canonical N23."
    )
}


envi_pb <- envi_pb[
    canonical_keys,
    ,
    drop = FALSE
]


envi_meta <- envi_meta[
    canonical_keys,
    ,
    drop = FALSE
]


if (
    !setequal(
        canonical_genes,
        colnames(envi_pb)
    )
) {

    stop(
        "ENVI measured300 gene universe differs."
    )
}


envi_pb <- envi_pb[
    ,
    canonical_genes,
    drop = FALSE
]


cat(
    "ENVI matrix:",
    nrow(envi_pb),
    "PB x",
    ncol(envi_pb),
    "genes\n"
)


cat(
    "Original vs ENVI measured300 exact:",
    identical(
        original_pb,
        envi_pb
    ),
    "\n"
)


cat(
    "Original vs ENVI measured300 all values equal:",
    all(
        original_pb
        ==
        envi_pb
    ),
    "\n"
)


# ================================================================================================
# LOAD PAPER
# ================================================================================================


section("LOAD PAPER PSEUDOBULK")


paper <- readRDS(
    PAPER_RDS
)


paper_meta_all <- as.data.frame(
    colData(paper)
)


paper_gene_info <- as.data.frame(
    rowData(paper)
)


required_gene_cols <- c(
    "ID",
    "Symbol"
)


if (
    !all(
        required_gene_cols
        %in%
        colnames(paper_gene_info)
    )
) {

    stop(
        "Paper rowData is missing ID and/or Symbol."
    )
}


paper_pb_all <- t(
    as.matrix(
        counts(paper)
    )
)


storage.mode(
    paper_pb_all
) <- "double"


colnames(
    paper_pb_all
) <- as.character(
    paper_gene_info$Symbol
)


paper_key_all <- make_key(
    paper_meta_all
)


rownames(
    paper_pb_all
) <- paper_key_all


rownames(
    paper_meta_all
) <- paper_key_all


cat(
    "Paper full matrix:",
    nrow(paper_pb_all),
    "PB x",
    ncol(paper_pb_all),
    "genes\n"
)


# ================================================================================================
# PAPER N24 -> CANONICAL N23
# ================================================================================================


section("SUBSET PAPER TO EXACT SAME 23 DONORS")


if (
    !all(
        canonical_keys
        %in%
        rownames(paper_pb_all)
    )
) {

    missing_keys <- setdiff(
        canonical_keys,
        rownames(paper_pb_all)
    )

    cat(
        "Missing keys:\n"
    )

    print(
        missing_keys
    )

    stop(
        "Paper does not contain every canonical donor-domain key."
    )
}


paper_pb <- paper_pb_all[
    canonical_keys,
    ,
    drop = FALSE
]


paper_meta <- paper_meta_all[
    canonical_keys,
    ,
    drop = FALSE
]


if (
    !setequal(
        canonical_genes,
        colnames(paper_pb)
    )
) {

    stop(
        "Paper 300-gene universe differs from GeneBridge."
    )
}


paper_pb <- paper_pb[
    ,
    canonical_genes,
    drop = FALSE
]


cat(
    "Paper N23:",
    nrow(paper_pb),
    "PB x",
    ncol(paper_pb),
    "genes\n"
)


cat(
    "Paper N23 donors:",
    length(unique(paper_meta$BrNum)),
    "\n"
)


# ================================================================================================
# BUILD GENE INFO IN CANONICAL ORDER
# ================================================================================================


section("GENE INFORMATION")


gene_idx <- match(
    canonical_genes,
    as.character(
        paper_gene_info$Symbol
    )
)


if (
    any(
        is.na(gene_idx)
    )
) {

    stop(
        "Could not map all 300 genes to paper rowData."
    )
}


gene_info <- data.frame(

    ID = as.character(
        paper_gene_info$ID[
            gene_idx
        ]
    ),

    Symbol = canonical_genes,

    stringsAsFactors = FALSE
)


cat(
    "Gene information rows:",
    nrow(gene_info),
    "\n"
)


# ================================================================================================
# VERIFY METADATA
# ================================================================================================


section("VERIFY THREE-WAY METADATA")


required_meta <- c(
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "predictions_smooth"
)


for (
    nm in c(
        "Paper",
        "GeneBridge",
        "ENVI"
    )
) {

    m <- switch(
        nm,

        "Paper" = paper_meta,

        "GeneBridge" = original_meta,

        "ENVI" = envi_meta
    )


    missing_cols <- setdiff(
        required_meta,
        colnames(m)
    )


    if (
        length(
            missing_cols
        )
        > 0
    ) {

        stop(
            nm,
            " metadata missing: ",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }
}


cat(
    "\nDiagnosis counts by donor:\n"
)


print(
    unique(
        original_meta[
            ,
            c(
                "BrNum",
                "Dx"
            )
        ]
    )$Dx
    |>
        table()
)


cat(
    "\nSpatial-domain counts:\n"
)


print(
    table(
        original_meta$predictions_smooth
    )
)


# ================================================================================================
# FUNCTION: PREPARE PAPER-STYLE PSEUDOBULK
# ================================================================================================


prepare_pseudobulk <- function(
    pb,
    meta,
    gene_info,
    label
) {

    section(
        paste0(
            "PREPARE: ",
            label
        )
    )


    if (
        nrow(pb)
        !=
        nrow(meta)
    ) {

        stop(
            label,
            ": matrix / metadata rows differ."
        )
    }


    sce <- SingleCellExperiment(

        assays = list(

            counts = t(
                pb
            )
        ),

        colData = DataFrame(
            meta
        ),

        rowData = DataFrame(
            gene_info
        )
    )


    rownames(sce) <- gene_info$Symbol

    colnames(sce) <- rownames(meta)


    # ----------------------------------------------------------------
    # Match paper factor handling.
    # ----------------------------------------------------------------

    sce$Dx <- factor(
        as.character(
            sce$Dx
        ),
        levels = c(
            "NTC",
            "SCZ"
        )
    )


    sce$predictions_smooth <- factor(
        as.character(
            sce$predictions_smooth
        )
    )


    sce$Age <- as.numeric(
        as.character(
            sce$Age
        )
    )


    sce$Sex <- factor(
        as.character(
            sce$Sex
        )
    )


    sce$slide_id <- factor(
        as.character(
            sce$slide_id
        )
    )


    sce$BrNum <- factor(
        as.character(
            sce$BrNum
        )
    )


    # ----------------------------------------------------------------
    # Recompute filterByExpr AFTER restricting to the same N23.
    #
    # Paper registration_pseudobulk() uses the domain as the group.
    # predictions_smooth is one-to-one with the seven domains.
    # ----------------------------------------------------------------

    keep_expr <- edgeR::filterByExpr(
        sce,
        group = sce$predictions_smooth
    )


    cat(
        label,
        " filterByExpr retained:",
        sum(keep_expr),
        "/",
        length(keep_expr),
        "genes\n"
    )


    sce <- sce[
        keep_expr,
        ,
        drop = FALSE
    ]


    # ----------------------------------------------------------------
    # Exact spatialLIBD normalization:
    #
    # calcNormFactors -> cpm(log=TRUE, prior.count=1)
    # ----------------------------------------------------------------

    logcounts(sce) <- edgeR::cpm(

        edgeR::calcNormFactors(
            sce
        ),

        log = TRUE,

        prior.count = 1
    )


    cat(
        label,
        " normalized shape:",
        nrow(sce),
        "genes x",
        ncol(sce),
        "PBs\n"
    )


    return(
        sce
    )
}


# ================================================================================================
# PREPARE THREE OBJECTS
# ================================================================================================


paper_sce <- prepare_pseudobulk(
    paper_pb,
    paper_meta,
    gene_info,
    "PAPER N23"
)


original_sce <- prepare_pseudobulk(
    original_pb,
    original_meta,
    gene_info,
    "GENEBRIDGE BEFORE N23"
)


envi_sce <- prepare_pseudobulk(
    envi_pb,
    envi_meta,
    gene_info,
    "ENVI MEASURED300 N23"
)


# ================================================================================================
# ENSURE SAME FILTERED GENE UNIVERSE
# ================================================================================================


section("CHECK FILTERED GENE UNIVERSE")


cat(
    "Paper genes:",
    nrow(paper_sce),
    "\n"
)


cat(
    "GeneBridge genes:",
    nrow(original_sce),
    "\n"
)


cat(
    "ENVI measured genes:",
    nrow(envi_sce),
    "\n"
)


if (
    !setequal(
        rownames(paper_sce),
        rownames(original_sce)
    )
    ||
    !setequal(
        rownames(paper_sce),
        rownames(envi_sce)
    )
) {

    stop(
        "filterByExpr produced different three-way gene universes."
    )
}


common_genes <- rownames(
    paper_sce
)


original_sce <- original_sce[
    common_genes,
    ,
    drop = FALSE
]


envi_sce <- envi_sce[
    common_genes,
    ,
    drop = FALSE
]


# ================================================================================================
# FUNCTION: RUN EXACT PAPER-STYLE DE
# ================================================================================================


run_domain_de <- function(
    sce,
    label
) {

    section(
        paste0(
            "DOMAIN-ADJUSTED DE: ",
            label
        )
    )


    covars <- c(
        "predictions_smooth",
        "Age",
        "Sex",
        "slide_id"
    )


    # ----------------------------------------------------------------
    # Same registration_model() call as paper.
    # ----------------------------------------------------------------

    dx_mod <- spatialLIBD::registration_model(

        sce,

        covars = covars,

        var_registration = "Dx"
    )


    cat(
        "Design:",
        nrow(dx_mod),
        "rows x",
        ncol(dx_mod),
        "columns\n"
    )


    cat(
        "Design rank:",
        qr(dx_mod)$rank,
        "/",
        ncol(dx_mod),
        "\n"
    )


    cat(
        "Design columns:\n"
    )


    print(
        colnames(dx_mod)
    )


    # ----------------------------------------------------------------
    # Same donor blocking as paper.
    # ----------------------------------------------------------------

    block_cor <- spatialLIBD::registration_block_cor(

        sce,

        registration_model = dx_mod,

        var_sample_id = "BrNum"
    )


    cat(
        "Block correlation:",
        sprintf(
            "%.10f",
            block_cor
        ),
        "\n"
    )


    # ----------------------------------------------------------------
    # Same enrichment statistics as paper.
    # ----------------------------------------------------------------

    dx_res <- spatialLIBD::registration_stats_enrichment(

        sce,

        block_cor = block_cor,

        covars = covars,

        var_registration = "Dx",

        var_sample_id = "BrNum",

        gene_ensembl = "ID",

        gene_name = "Symbol"
    )


    required_res <- c(
        "gene",
        "ensembl",
        "logFC_SCZ",
        "t_stat_SCZ",
        "p_value_SCZ",
        "fdr_SCZ"
    )


    missing_res <- setdiff(
        required_res,
        colnames(dx_res)
    )


    if (
        length(missing_res)
        > 0
    ) {

        stop(
            label,
            " result missing: ",
            paste(
                missing_res,
                collapse = ", "
            )
        )
    }


    cat(
        "Genes tested:",
        nrow(dx_res),
        "\n"
    )


    cat(
        "Nominal P < 0.05:",
        sum(
            dx_res$p_value_SCZ
            <
            0.05,
            na.rm = TRUE
        ),
        "\n"
    )


    cat(
        "FDR < 0.10:",
        sum(
            dx_res$fdr_SCZ
            <
            0.10,
            na.rm = TRUE
        ),
        "\n"
    )


    cat(
        "FDR < 0.05:",
        sum(
            dx_res$fdr_SCZ
            <
            0.05,
            na.rm = TRUE
        ),
        "\n"
    )


    return(
        list(
            result = dx_res,
            block_cor = block_cor,
            design = dx_mod
        )
    )
}


# ================================================================================================
# RUN THREE DE MODELS
# ================================================================================================


paper_fit <- run_domain_de(
    paper_sce,
    "PAPER N23"
)


original_fit <- run_domain_de(
    original_sce,
    "GENEBRIDGE BEFORE N23"
)


envi_fit <- run_domain_de(
    envi_sce,
    "ENVI MEASURED300 N23"
)


# ================================================================================================
# SAVE INDIVIDUAL RESULTS
# ================================================================================================


section("SAVE INDIVIDUAL DE RESULTS")


write.csv(
    paper_fit$result,
    file.path(
        OUT,
        "paper_N23_domain_adjusted_DE.csv"
    ),
    row.names = FALSE
)


write.csv(
    original_fit$result,
    file.path(
        OUT,
        "genebridge_before_N23_domain_adjusted_DE.csv"
    ),
    row.names = FALSE
)


write.csv(
    envi_fit$result,
    file.path(
        OUT,
        "envi_measured300_N23_domain_adjusted_DE.csv"
    ),
    row.names = FALSE
)


saveRDS(
    paper_sce,
    file.path(
        OUT,
        "paper_N23_reprocessed_pseudobulk.rds"
    )
)


saveRDS(
    original_sce,
    file.path(
        OUT,
        "genebridge_before_N23_reprocessed_pseudobulk.rds"
    )
)


saveRDS(
    envi_sce,
    file.path(
        OUT,
        "envi_measured300_N23_reprocessed_pseudobulk.rds"
    )
)


# ================================================================================================
# MASTER THREE-WAY TABLE
# ================================================================================================


section("BUILD THREE-WAY MASTER TABLE")


extract_result <- function(
    x,
    prefix
) {

    y <- x[
        ,
        c(
            "gene",
            "ensembl",
            "logFC_SCZ",
            "t_stat_SCZ",
            "p_value_SCZ",
            "fdr_SCZ"
        )
    ]


    colnames(y) <- c(
        "gene",
        paste0(prefix, "_ensembl"),
        paste0(prefix, "_logFC"),
        paste0(prefix, "_t"),
        paste0(prefix, "_P"),
        paste0(prefix, "_FDR")
    )


    return(y)
}


paper_res <- extract_result(
    paper_fit$result,
    "paper"
)


original_res <- extract_result(
    original_fit$result,
    "genebridge"
)


envi_res <- extract_result(
    envi_fit$result,
    "envi_measured"
)


master <- merge(
    paper_res,
    original_res,
    by = "gene",
    all = TRUE
)


master <- merge(
    master,
    envi_res,
    by = "gene",
    all = TRUE
)


master$paper_vs_genebridge_direction_match <- (
    sign(master$paper_logFC)
    ==
    sign(master$genebridge_logFC)
)


master$paper_vs_envi_direction_match <- (
    sign(master$paper_logFC)
    ==
    sign(master$envi_measured_logFC)
)


master$genebridge_vs_envi_direction_match <- (
    sign(master$genebridge_logFC)
    ==
    sign(master$envi_measured_logFC)
)


master$paper_FDR10 <- (
    master$paper_FDR
    <
    0.10
)


master$genebridge_FDR10 <- (
    master$genebridge_FDR
    <
    0.10
)


master$envi_measured_FDR10 <- (
    master$envi_measured_FDR
    <
    0.10
)


master$all_three_FDR10 <- (
    master$paper_FDR10
    &
    master$genebridge_FDR10
    &
    master$envi_measured_FDR10
)


master$paper_genebridge_abs_logFC_diff <- abs(
    master$paper_logFC
    -
    master$genebridge_logFC
)


master$paper_envi_abs_logFC_diff <- abs(
    master$paper_logFC
    -
    master$envi_measured_logFC
)


master <- master[
    order(
        master$paper_FDR,
        master$genebridge_FDR
    ),
    ,
    drop = FALSE
]


MASTER_FILE <- file.path(
    OUT,
    "threeway_N23_domain_DE_master_300genes.csv"
)


write.csv(
    master,
    MASTER_FILE,
    row.names = FALSE
)


cat(
    "Master table:",
    MASTER_FILE,
    "\n"
)


# ================================================================================================
# SUMMARY
# ================================================================================================


section("THREE-WAY DE SUMMARY")


summary_df <- data.frame(

    dataset = c(
        "Paper_N23",
        "GeneBridge_before_N23",
        "ENVI_measured300_N23"
    ),

    donors = c(
        length(unique(paper_sce$BrNum)),
        length(unique(original_sce$BrNum)),
        length(unique(envi_sce$BrNum))
    ),

    pseudobulks = c(
        ncol(paper_sce),
        ncol(original_sce),
        ncol(envi_sce)
    ),

    genes_tested = c(
        nrow(paper_fit$result),
        nrow(original_fit$result),
        nrow(envi_fit$result)
    ),

    block_correlation = c(
        paper_fit$block_cor,
        original_fit$block_cor,
        envi_fit$block_cor
    ),

    nominal_P05 = c(
        sum(paper_fit$result$p_value_SCZ < 0.05),
        sum(original_fit$result$p_value_SCZ < 0.05),
        sum(envi_fit$result$p_value_SCZ < 0.05)
    ),

    FDR10 = c(
        sum(paper_fit$result$fdr_SCZ < 0.10),
        sum(original_fit$result$fdr_SCZ < 0.10),
        sum(envi_fit$result$fdr_SCZ < 0.10)
    ),

    FDR05 = c(
        sum(paper_fit$result$fdr_SCZ < 0.05),
        sum(original_fit$result$fdr_SCZ < 0.05),
        sum(envi_fit$result$fdr_SCZ < 0.05)
    ),

    stringsAsFactors = FALSE
)


print(
    summary_df,
    row.names = FALSE
)


write.csv(
    summary_df,
    file.path(
        OUT,
        "threeway_N23_domain_DE_summary.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# CONCORDANCE METRICS
# ================================================================================================


section("CONCORDANCE METRICS")


metrics <- data.frame(

    comparison = c(

        "Paper_vs_GeneBridge_logFC",
        "Paper_vs_ENVI_measured_logFC",

        "Paper_vs_GeneBridge_t",
        "Paper_vs_ENVI_measured_t",

        "GeneBridge_vs_ENVI_measured_logFC",
        "GeneBridge_vs_ENVI_measured_t"
    ),

    Pearson = c(

        safe_cor(
            master$paper_logFC,
            master$genebridge_logFC
        ),

        safe_cor(
            master$paper_logFC,
            master$envi_measured_logFC
        ),

        safe_cor(
            master$paper_t,
            master$genebridge_t
        ),

        safe_cor(
            master$paper_t,
            master$envi_measured_t
        ),

        safe_cor(
            master$genebridge_logFC,
            master$envi_measured_logFC
        ),

        safe_cor(
            master$genebridge_t,
            master$envi_measured_t
        )
    ),

    Spearman = c(

        safe_cor(
            master$paper_logFC,
            master$genebridge_logFC,
            "spearman"
        ),

        safe_cor(
            master$paper_logFC,
            master$envi_measured_logFC,
            "spearman"
        ),

        safe_cor(
            master$paper_t,
            master$genebridge_t,
            "spearman"
        ),

        safe_cor(
            master$paper_t,
            master$envi_measured_t,
            "spearman"
        ),

        safe_cor(
            master$genebridge_logFC,
            master$envi_measured_logFC,
            "spearman"
        ),

        safe_cor(
            master$genebridge_t,
            master$envi_measured_t,
            "spearman"
        )
    ),

    stringsAsFactors = FALSE
)


print(
    metrics,
    row.names = FALSE
)


write.csv(
    metrics,
    file.path(
        OUT,
        "threeway_N23_domain_DE_concordance_metrics.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# PLOT HELPERS
# ================================================================================================


scatter_panel <- function(
    x,
    y,
    xlab,
    ylab,
    title
) {

    lim <- range(
        c(
            x,
            y
        ),
        finite = TRUE
    )


    plot(
        x,
        y,
        xlim = lim,
        ylim = lim,
        pch = 16,
        cex = 0.75,
        xlab = xlab,
        ylab = ylab,
        main = title
    )


    abline(
        a = 0,
        b = 1,
        lty = 2,
        lwd = 2
    )


    abline(
        h = 0,
        v = 0,
        lty = 3
    )


    legend(
        "topleft",

        legend = c(

            paste0(
                "Pearson r = ",
                sprintf(
                    "%.4f",
                    safe_cor(
                        x,
                        y
                    )
                )
            ),

            paste0(
                "Spearman rho = ",
                sprintf(
                    "%.4f",
                    safe_cor(
                        x,
                        y,
                        "spearman"
                    )
                )
            )
        ),

        bty = "n",
        cex = 0.82
    )
}


make_concordance_figure <- function() {

    par(
        mfrow = c(
            2,
            2
        ),

        mar = c(
            5.5,
            5.8,
            4.2,
            1.3
        ),

        oma = c(
            0,
            0,
            3,
            0
        ),

        mgp = c(
            3.6,
            1,
            0
        ),

        las = 1
    )


    scatter_panel(

        master$paper_logFC,
        master$genebridge_logFC,

        "Paper SCZ logFC",
        "GeneBridge before SCZ logFC",

        "A. Paper vs GeneBridge"
    )


    scatter_panel(

        master$paper_logFC,
        master$envi_measured_logFC,

        "Paper SCZ logFC",
        "ENVI measured-300 SCZ logFC",

        "B. Paper vs ENVI measured"
    )


    scatter_panel(

        master$paper_t,
        master$genebridge_t,

        "Paper SCZ t-statistic",
        "GeneBridge before SCZ t-statistic",

        "C. Paper vs GeneBridge"
    )


    scatter_panel(

        master$paper_t,
        master$envi_measured_t,

        "Paper SCZ t-statistic",
        "ENVI measured-300 SCZ t-statistic",

        "D. Paper vs ENVI measured"
    )


    mtext(
        "N23 domain-adjusted SCZ differential-expression concordance",
        outer = TRUE,
        font = 2,
        cex = 1.25
    )
}


# ================================================================================================
# SAVE CONCORDANCE FIGURE
# ================================================================================================


section("PLOT 1: DE CONCORDANCE")


pdf(
    file.path(
        OUT,
        "04A_threeway_N23_logFC_t_concordance.pdf"
    ),
    width = 12,
    height = 10
)


make_concordance_figure()


dev.off()


png(
    file.path(
        OUT,
        "04A_threeway_N23_logFC_t_concordance.png"
    ),
    width = 3600,
    height = 3000,
    res = 300
)


make_concordance_figure()


dev.off()


# ================================================================================================
# VOLCANO PLOTS
# ================================================================================================


volcano_panel <- function(
    logfc,
    fdr,
    title
) {

    y <- -log10(
        pmax(
            fdr,
            .Machine$double.xmin
        )
    )


    sig <- (
        fdr
        <
        0.10
    )


    plot(
        logfc,
        y,
        pch = ifelse(
            sig,
            16,
            1
        ),
        cex = 0.8,
        xlab = "SCZ logFC",
        ylab = "-log10(FDR)",
        main = paste0(
            title,
            "\nFDR < 0.10: ",
            sum(
                sig,
                na.rm = TRUE
            )
        )
    )


    abline(
        h = -log10(0.10),
        lty = 2
    )


    abline(
        v = 0,
        lty = 3
    )
}


make_volcano_figure <- function() {

    par(
        mfrow = c(
            1,
            3
        ),

        mar = c(
            5.5,
            5.5,
            4.5,
            1
        ),

        oma = c(
            0,
            0,
            3,
            0
        ),

        mgp = c(
            3.5,
            1,
            0
        ),

        las = 1
    )


    volcano_panel(

        master$paper_logFC,
        master$paper_FDR,

        "Paper N23"
    )


    volcano_panel(

        master$genebridge_logFC,
        master$genebridge_FDR,

        "GeneBridge before N23"
    )


    volcano_panel(

        master$envi_measured_logFC,
        master$envi_measured_FDR,

        "ENVI measured-300 N23"
    )


    mtext(
        "Three-way domain-adjusted SCZ differential expression",
        outer = TRUE,
        font = 2,
        cex = 1.25
    )
}


section("PLOT 2: THREE-WAY VOLCANO")


pdf(
    file.path(
        OUT,
        "04B_threeway_N23_volcano.pdf"
    ),
    width = 15,
    height = 5.5
)


make_volcano_figure()


dev.off()


png(
    file.path(
        OUT,
        "04B_threeway_N23_volcano.png"
    ),
    width = 4500,
    height = 1650,
    res = 300
)


make_volcano_figure()


dev.off()


# ================================================================================================
# P-VALUE HISTOGRAM + QQ PLOTS
# ================================================================================================


qq_panel <- function(
    p,
    title
) {

    p <- p[
        is.finite(p)
        &
        !is.na(p)
    ]


    p <- sort(
        p
    )


    n <- length(
        p
    )


    expected <- -log10(
        ppoints(
            n
        )
    )


    observed <- -log10(
        pmax(
            p,
            .Machine$double.xmin
        )
    )


    lim <- range(
        c(
            expected,
            observed
        ),
        finite = TRUE
    )


    plot(
        expected,
        observed,
        xlim = lim,
        ylim = lim,
        pch = 16,
        cex = 0.7,
        xlab = "Expected -log10(P)",
        ylab = "Observed -log10(P)",
        main = title
    )


    abline(
        a = 0,
        b = 1,
        lty = 2,
        lwd = 2
    )
}


make_pvalue_figure <- function() {

    par(
        mfrow = c(
            2,
            3
        ),

        mar = c(
            5,
            5,
            4,
            1
        ),

        oma = c(
            0,
            0,
            3,
            0
        ),

        mgp = c(
            3.2,
            1,
            0
        ),

        las = 1
    )


    hist(
        master$paper_P,
        breaks = 20,
        main = "Paper N23",
        xlab = "P-value"
    )


    hist(
        master$genebridge_P,
        breaks = 20,
        main = "GeneBridge before",
        xlab = "P-value"
    )


    hist(
        master$envi_measured_P,
        breaks = 20,
        main = "ENVI measured-300",
        xlab = "P-value"
    )


    qq_panel(
        master$paper_P,
        "Paper N23 QQ"
    )


    qq_panel(
        master$genebridge_P,
        "GeneBridge before QQ"
    )


    qq_panel(
        master$envi_measured_P,
        "ENVI measured-300 QQ"
    )


    mtext(
        "N23 domain-adjusted P-value diagnostics",
        outer = TRUE,
        font = 2,
        cex = 1.25
    )
}


section("PLOT 3: P-VALUE + QQ DIAGNOSTICS")


pdf(
    file.path(
        OUT,
        "04C_threeway_N23_pvalue_QQ.pdf"
    ),
    width = 15,
    height = 10
)


make_pvalue_figure()


dev.off()


png(
    file.path(
        OUT,
        "04C_threeway_N23_pvalue_QQ.png"
    ),
    width = 4500,
    height = 3000,
    res = 300
)


make_pvalue_figure()


dev.off()


# ================================================================================================
# SIGNIFICANT-GENE COMPARISON
# ================================================================================================


section("FDR < 0.10 GENE COMPARISON")


sig_union <- master[
    master$paper_FDR10
    |
    master$genebridge_FDR10
    |
    master$envi_measured_FDR10,
    ,
    drop = FALSE
]


sig_union <- sig_union[
    order(
        sig_union$paper_FDR,
        sig_union$genebridge_FDR
    ),
    ,
    drop = FALSE
]


print_cols <- c(
    "gene",

    "paper_logFC",
    "paper_FDR",

    "genebridge_logFC",
    "genebridge_FDR",

    "envi_measured_logFC",
    "envi_measured_FDR",

    "paper_vs_genebridge_direction_match",
    "paper_vs_envi_direction_match"
)


print(
    sig_union[
        ,
        print_cols,
        drop = FALSE
    ],
    row.names = FALSE
)


write.csv(
    sig_union,
    file.path(
        OUT,
        "threeway_N23_FDR10_union_genes.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# FINAL STATUS
# ================================================================================================


section("FINAL STATUS")


cat(
    "Analysis universe:\n"
)


cat(
    "  Donors      : 23\n"
)


cat(
    "  Br6432      : excluded\n"
)


cat(
    "  Domains     : 7\n"
)


cat(
    "  Pseudobulks : 161\n"
)


cat(
    "  Model       : Dx + predictions_smooth + Age + Sex + slide_id\n"
)


cat(
    "  Block       : BrNum\n"
)


cat(
    "  run_date    : metadata only; NOT DE covariate\n"
)


cat(
    "\nThree-way comparison:\n"
)


cat(
    "  Paper N23\n"
)


cat(
    "  GeneBridge before N23\n"
)


cat(
    "  ENVI measured300 N23 [technical control, NOT OOF]\n"
)


cat(
    "\nMain outputs:\n"
)


cat(
    "  threeway_N23_domain_DE_master_300genes.csv\n"
)


cat(
    "  threeway_N23_domain_DE_summary.csv\n"
)


cat(
    "  threeway_N23_domain_DE_concordance_metrics.csv\n"
)


cat(
    "  threeway_N23_FDR10_union_genes.csv\n"
)


cat(
    "  04A_threeway_N23_logFC_t_concordance.png/pdf\n"
)


cat(
    "  04B_threeway_N23_volcano.png/pdf\n"
)


cat(
    "  04C_threeway_N23_pvalue_QQ.png/pdf\n"
)


cat(
    "\nFINAL STATUS: THREE-WAY N23 DOMAIN DE COMPLETE\n"
)
