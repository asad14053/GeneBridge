#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(edgeR)
    library(limma)
    library(spatialLIBD)
})


# ================================================================================================
# PURPOSE
# ================================================================================================
#
# N23 ENVI full-transcriptome domain-adjusted SCZ DE.
#
# Input:
#   161 donor x SpD pseudobulks
#   34,987 total genes
#
#       300 measured_xenium
#    34,687 envi_imputed
#
# Model:
#
#   expression ~ Dx + predictions_smooth + Age + Sex + slide_id
#
# Repeated donor-domain observations:
#
#   block = BrNum
#
# run_date is metadata only and is NOT included in the DE model,
# matching the Lieber Xenium domain-level DE.
#
# ================================================================================================


ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"


PB_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_full34987_donor_layer_pseudobulk_countscale.csv.gz"
)


META_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_donor_layer_pseudobulk_metadata.csv"
)


GENE_INFO_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_full34987_gene_info.csv"
)


OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_full_N23_domain_DE"
)


dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


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


# ================================================================================================
# INPUT CHECK
# ================================================================================================


section("INPUT CHECK")


for (
    f in c(
        PB_FILE,
        META_FILE,
        GENE_INFO_FILE
    )
) {

    if (!file.exists(f)) {

        stop(
            "Missing input:\n",
            f
        )
    }

    cat(
        "FOUND: ",
        f,
        "\n",
        sep = ""
    )
}


# ================================================================================================
# LOAD PSEUDOBULK
# ================================================================================================


section("LOAD ENVI FULL PSEUDOBULK")


pb <- read.csv(
    PB_FILE,
    row.names = 1,
    check.names = FALSE
)


pb <- as.matrix(
    pb
)


storage.mode(
    pb
) <- "double"


cat(
    "Pseudobulk matrix:",
    nrow(pb),
    "samples x",
    ncol(pb),
    "genes\n"
)


if (
    nrow(pb)
    != 161
) {

    stop(
        "Expected exactly 161 donor-domain pseudobulks."
    )
}


if (
    ncol(pb)
    != 34987
) {

    stop(
        "Expected exactly 34,987 genes."
    )
}


if (
    any(!is.finite(pb))
) {

    stop(
        "Non-finite expression values found."
    )
}


if (
    any(pb < 0)
) {

    stop(
        "Negative count-scale expression values found."
    )
}


# ================================================================================================
# LOAD METADATA
# ================================================================================================


section("LOAD METADATA")


meta <- read.csv(
    META_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (
    !"pseudobulk_id"
    %in%
    colnames(meta)
) {

    stop(
        "pseudobulk_id missing from metadata."
    )
}


idx <- match(
    rownames(pb),
    meta$pseudobulk_id
)


if (
    any(is.na(idx))
) {

    stop(
        "Could not align pseudobulk metadata."
    )
}


meta <- meta[
    idx,
    ,
    drop = FALSE
]


rownames(meta) <- rownames(pb)


required_meta <- c(
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "predictions_smooth"
)


missing_meta <- setdiff(
    required_meta,
    colnames(meta)
)


if (
    length(missing_meta)
    > 0
) {

    stop(
        "Missing metadata: ",
        paste(
            missing_meta,
            collapse = ", "
        )
    )
}


cat(
    "Donors:",
    length(unique(meta$BrNum)),
    "\n"
)


cat(
    "Pseudobulks:",
    nrow(meta),
    "\n"
)


cat(
    "Domains:",
    length(unique(meta$predictions_smooth)),
    "\n"
)


if (
    length(unique(meta$BrNum))
    != 23
) {

    stop(
        "Expected exactly 23 donors."
    )
}


if (
    "Br6432"
    %in%
    as.character(meta$BrNum)
) {

    stop(
        "Br6432 must be excluded."
    )
}


cat(
    "Br6432 excluded: TRUE\n"
)


cat(
    "\nDiagnosis by donor:\n"
)


donor_dx <- unique(
    meta[
        ,
        c(
            "BrNum",
            "Dx"
        )
    ]
)


print(
    table(
        donor_dx$Dx
    )
)


cat(
    "\nDomain pseudobulks:\n"
)


print(
    table(
        meta$predictions_smooth
    )
)


# ================================================================================================
# LOAD GENE INFORMATION
# ================================================================================================


section("LOAD GENE INFORMATION")


gene_info <- read.csv(
    GENE_INFO_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_gene <- c(
    "gene",
    "expression_source"
)


missing_gene <- setdiff(
    required_gene,
    colnames(gene_info)
)


if (
    length(missing_gene)
    > 0
) {

    stop(
        "Gene info missing: ",
        paste(
            missing_gene,
            collapse = ", "
        )
    )
}


gene_idx <- match(
    colnames(pb),
    gene_info$gene
)


if (
    any(is.na(gene_idx))
) {

    stop(
        "Could not align all genes to gene_info."
    )
}


gene_info <- gene_info[
    gene_idx,
    ,
    drop = FALSE
]


if (
    !identical(
        as.character(gene_info$gene),
        colnames(pb)
    )
) {

    stop(
        "Gene order mismatch."
    )
}


cat(
    "\nExpression source:\n"
)


print(
    table(
        gene_info$expression_source
    )
)


n_measured <- sum(
    gene_info$expression_source
    ==
    "measured_xenium"
)


n_imputed <- sum(
    gene_info$expression_source
    ==
    "envi_imputed"
)


if (
    n_measured
    != 300
) {

    stop(
        "Expected 300 measured genes."
    )
}


if (
    n_imputed
    != 34687
) {

    stop(
        "Expected 34,687 imputed genes."
    )
}


# ================================================================================================
# BUILD SINGLECELLEXPERIMENT
# ================================================================================================


section("BUILD PSEUDOBULK OBJECT")


row_info <- DataFrame(

    Symbol = as.character(
        gene_info$gene
    ),

    expression_source = as.character(
        gene_info$expression_source
    )
)


if (
    "ensembl"
    %in%
    colnames(gene_info)
) {

    row_info$ID <- as.character(
        gene_info$ensembl
    )

} else if (
    "ID"
    %in%
    colnames(gene_info)
) {

    row_info$ID <- as.character(
        gene_info$ID
    )

} else {

    row_info$ID <- as.character(
        gene_info$gene
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

    rowData = row_info
)


rownames(sce) <- gene_info$gene
colnames(sce) <- rownames(meta)


# ================================================================================================
# FACTORS
# ================================================================================================


section("PREPARE MODEL VARIABLES")


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


# ================================================================================================
# FILTER BY EXPRESSION
# ================================================================================================


section("FILTER BY EXPRESSION")


keep_expr <- edgeR::filterByExpr(

    sce,

    group = sce$predictions_smooth
)


cat(
    "Retained:",
    sum(keep_expr),
    "/",
    length(keep_expr),
    "genes\n"
)


filter_table <- data.frame(

    expression_source = rowData(sce)$expression_source,

    retained = keep_expr,

    stringsAsFactors = FALSE
)


cat(
    "\nFilter result by expression source:\n"
)


print(
    with(
        filter_table,
        table(
            expression_source,
            retained
        )
    )
)


write.csv(
    data.frame(
        gene = rownames(sce),
        expression_source = rowData(sce)$expression_source,
        filterByExpr_retained = keep_expr
    ),

    file.path(
        OUT,
        "ENVI_full34987_filterByExpr_status.csv"
    ),

    row.names = FALSE
)


sce <- sce[
    keep_expr,
    ,
    drop = FALSE
]


# ================================================================================================
# NORMALIZE
# ================================================================================================


section("TMM + LOG2 CPM NORMALIZATION")


logcounts(sce) <- edgeR::cpm(

    edgeR::calcNormFactors(
        sce
    ),

    log = TRUE,

    prior.count = 1
)


cat(
    "Normalized object:",
    nrow(sce),
    "genes x",
    ncol(sce),
    "pseudobulks\n"
)


# ================================================================================================
# MODEL
# ================================================================================================


section("PAPER-STYLE DOMAIN-ADJUSTED MODEL")


covars <- c(
    "predictions_smooth",
    "Age",
    "Sex",
    "slide_id"
)


dx_mod <- spatialLIBD::registration_model(

    sce,

    covars = covars,

    var_registration = "Dx"
)


cat(
    "Design:",
    nrow(dx_mod),
    "x",
    ncol(dx_mod),
    "\n"
)


cat(
    "Design rank:",
    qr(dx_mod)$rank,
    "/",
    ncol(dx_mod),
    "\n"
)


cat(
    "\nDesign columns:\n"
)


print(
    colnames(dx_mod)
)


# ================================================================================================
# DONOR BLOCKING
# ================================================================================================


section("DONOR BLOCK CORRELATION")


dx_block_cor <- spatialLIBD::registration_block_cor(

    sce,

    registration_model = dx_mod,

    var_sample_id = "BrNum"
)


cat(
    "Block correlation:",
    sprintf(
        "%.10f",
        dx_block_cor
    ),
    "\n"
)


# ================================================================================================
# DE
# ================================================================================================


section("RUN DOMAIN-ADJUSTED SCZ DE")


dx_res <- spatialLIBD::registration_stats_enrichment(

    sce,

    block_cor = dx_block_cor,

    covars = covars,

    var_registration = "Dx",

    var_sample_id = "BrNum",

    gene_ensembl = "ID",

    gene_name = "Symbol"
)


source_map <- data.frame(

    gene = rownames(sce),

    expression_source = as.character(
        rowData(sce)$expression_source
    ),

    stringsAsFactors = FALSE
)


dx_res <- merge(

    dx_res,

    source_map,

    by = "gene",

    all.x = TRUE,

    sort = FALSE
)


# ================================================================================================
# SPLIT MEASURED / IMPUTED
# ================================================================================================


section("SPLIT RESULTS")


measured_res <- dx_res[
    dx_res$expression_source
    ==
    "measured_xenium",
    ,
    drop = FALSE
]


imputed_res <- dx_res[
    dx_res$expression_source
    ==
    "envi_imputed",
    ,
    drop = FALSE
]


cat(
    "All genes tested:",
    nrow(dx_res),
    "\n"
)


cat(
    "Measured genes tested:",
    nrow(measured_res),
    "\n"
)


cat(
    "Imputed genes tested:",
    nrow(imputed_res),
    "\n"
)


# ================================================================================================
# SUMMARY FUNCTION
# ================================================================================================


summarize_de <- function(
    x,
    label
) {

    data.frame(

        gene_set = label,

        genes_tested = nrow(x),

        nominal_P05 = sum(
            x$p_value_SCZ
            <
            0.05,
            na.rm = TRUE
        ),

        FDR10 = sum(
            x$fdr_SCZ
            <
            0.10,
            na.rm = TRUE
        ),

        FDR05 = sum(
            x$fdr_SCZ
            <
            0.05,
            na.rm = TRUE
        ),

        stringsAsFactors = FALSE
    )
}


summary_df <- rbind(

    summarize_de(
        dx_res,
        "ENVI_full"
    ),

    summarize_de(
        measured_res,
        "measured_xenium"
    ),

    summarize_de(
        imputed_res,
        "envi_imputed"
    )
)


section("DE SUMMARY")


print(
    summary_df,
    row.names = FALSE
)


# ================================================================================================
# SAVE TABLES
# ================================================================================================


section("SAVE RESULTS")


write.csv(
    dx_res,
    file.path(
        OUT,
        "ENVI_full_domain_adjusted_DE_SCZ_vs_NTC.csv"
    ),
    row.names = FALSE
)


write.csv(
    measured_res,
    file.path(
        OUT,
        "ENVI_measured_domain_adjusted_DE_SCZ_vs_NTC.csv"
    ),
    row.names = FALSE
)


write.csv(
    imputed_res,
    file.path(
        OUT,
        "ENVI_imputed_domain_adjusted_DE_SCZ_vs_NTC.csv"
    ),
    row.names = FALSE
)


write.csv(
    imputed_res[
        imputed_res$p_value_SCZ
        <
        0.05,
        ,
        drop = FALSE
    ],

    file.path(
        OUT,
        "ENVI_imputed_domain_adjusted_DE_nominalP05.csv"
    ),

    row.names = FALSE
)


write.csv(
    imputed_res[
        imputed_res$fdr_SCZ
        <
        0.10,
        ,
        drop = FALSE
    ],

    file.path(
        OUT,
        "ENVI_imputed_domain_adjusted_DE_FDR10.csv"
    ),

    row.names = FALSE
)


write.csv(
    imputed_res[
        imputed_res$fdr_SCZ
        <
        0.05,
        ,
        drop = FALSE
    ],

    file.path(
        OUT,
        "ENVI_imputed_domain_adjusted_DE_FDR05.csv"
    ),

    row.names = FALSE
)


write.csv(
    summary_df,
    file.path(
        OUT,
        "ENVI_full_domain_adjusted_DE_summary.csv"
    ),
    row.names = FALSE
)


# ================================================================================================
# P-VALUE HISTOGRAM
# ================================================================================================


section("PLOT P-VALUE HISTOGRAM")


pdf(
    file.path(
        OUT,
        "05A_ENVI_imputed34687_domain_DE_pvalue_histogram.pdf"
    ),
    width = 7,
    height = 6
)


hist(
    imputed_res$p_value_SCZ,

    breaks = 50,

    main = paste0(
        "ENVI-imputed genes: P-value distribution\n",
        "Domain-adjusted N23 SCZ vs NTC"
    ),

    xlab = "P-value",

    ylab = "Number of genes"
)


dev.off()


png(
    file.path(
        OUT,
        "05A_ENVI_imputed34687_domain_DE_pvalue_histogram.png"
    ),
    width = 2100,
    height = 1800,
    res = 300
)


hist(
    imputed_res$p_value_SCZ,

    breaks = 50,

    main = paste0(
        "ENVI-imputed genes: P-value distribution\n",
        "Domain-adjusted N23 SCZ vs NTC"
    ),

    xlab = "P-value",

    ylab = "Number of genes"
)


dev.off()


# ================================================================================================
# QQ PLOT
# ================================================================================================


section("PLOT QQ")


p <- imputed_res$p_value_SCZ


p <- p[
    is.finite(p)
    &
    !is.na(p)
]


p <- sort(
    p
)


expected <- -log10(
    ppoints(
        length(p)
    )
)


observed <- -log10(
    pmax(
        p,
        .Machine$double.xmin
    )
)


qq_lim <- range(
    c(
        expected,
        observed
    ),
    finite = TRUE
)


make_qq <- function() {

    plot(
        expected,
        observed,

        xlim = qq_lim,
        ylim = qq_lim,

        pch = 16,
        cex = 0.45,

        xlab = "Expected -log10(P)",

        ylab = "Observed -log10(P)",

        main = paste0(
            "QQ plot: ENVI-imputed genes\n",
            "Domain-adjusted N23 SCZ vs NTC"
        )
    )


    abline(
        a = 0,
        b = 1,
        lty = 2,
        lwd = 2
    )
}


pdf(
    file.path(
        OUT,
        "05B_ENVI_imputed34687_domain_DE_QQ.pdf"
    ),
    width = 7,
    height = 7
)


make_qq()


dev.off()


png(
    file.path(
        OUT,
        "05B_ENVI_imputed34687_domain_DE_QQ.png"
    ),
    width = 2100,
    height = 2100,
    res = 300
)


make_qq()


dev.off()


# ================================================================================================
# VOLCANO
# ================================================================================================


section("PLOT VOLCANO")


volcano_y <- -log10(
    pmax(
        imputed_res$fdr_SCZ,
        .Machine$double.xmin
    )
)


sig10 <- (
    imputed_res$fdr_SCZ
    <
    0.10
)


make_volcano <- function() {

    plot(
        imputed_res$logFC_SCZ,
        volcano_y,

        pch = ifelse(
            sig10,
            16,
            1
        ),

        cex = 0.45,

        xlab = "SCZ logFC",

        ylab = "-log10(FDR)",

        main = paste0(
            "ENVI-imputed domain-adjusted DE\n",
            "FDR < 0.10: ",
            sum(
                sig10,
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


pdf(
    file.path(
        OUT,
        "05C_ENVI_imputed34687_domain_DE_volcano.pdf"
    ),
    width = 8,
    height = 7
)


make_volcano()


dev.off()


png(
    file.path(
        OUT,
        "05C_ENVI_imputed34687_domain_DE_volcano.png"
    ),
    width = 2400,
    height = 2100,
    res = 300
)


make_volcano()


dev.off()


# ================================================================================================
# FINAL
# ================================================================================================


section("FINAL STATUS")


cat(
    "Donors                : 23\n"
)


cat(
    "Pseudobulks           : 161\n"
)


cat(
    "Input genes           : 34,987\n"
)


cat(
    "Measured genes input  : 300\n"
)


cat(
    "Imputed genes input   : 34,687\n"
)


cat(
    "Genes after filtering :",
    nrow(sce),
    "\n"
)


cat(
    "Model                 : Dx + predictions_smooth + Age + Sex + slide_id\n"
)


cat(
    "Block                 : BrNum\n"
)


cat(
    "Block correlation     :",
    sprintf(
        "%.10f",
        dx_block_cor
    ),
    "\n"
)


cat(
    "\n"
)


print(
    summary_df,
    row.names = FALSE
)


cat(
    "\nFINAL STATUS: ENVI FULL N23 DOMAIN DE COMPLETE\n"
)
