#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(S4Vectors)
    library(edgeR)
    library(limma)
    library(spatialLIBD)
})

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PAPER_RDS <- file.path(
    ROOT,
    "data/reference/spatialDLPFC_SCZ_XENIUM/07_cell_type_de",
    "spe_pseudo_donor_celltype_clust_M0_lam0.1_k50_res0.7.rds"
)

PAPER_OFFICIAL_CSV <- file.path(
    ROOT,
    "data/reference/spatialDLPFC_SCZ_XENIUM/07_cell_type_de",
    "donor_cell_type_level_pseudobulk_Dx_DEGs_paper_N24.csv"
)

PB_DIR <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "celltype_N23/pseudobulk"
)

META_FILE <- file.path(
    PB_DIR,
    "N23_donor_celltype_pseudobulk_metadata.csv"
)

META23_FILE <- file.path(
    ROOT,
    "data/metadata/xenium_DE_metadata_23.csv"
)

GENE_INFO_FILE <- file.path(
    PB_DIR,
    "ENVI_full34987_gene_info.csv"
)

OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "celltype_N23/global_celltype_adjusted_DE"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


section <- function(x) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(x, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


# =============================================================================
# CORE PAPER DE FUNCTION
# =============================================================================

run_paper_de <- function(
    sce,
    label,
    renormalize = TRUE
) {

    section(paste("RUN:", label))

    sce$Dx <- factor(
        as.character(sce$Dx),
        levels = c("NTC", "SCZ")
    )

    sce$annots <- droplevels(
        factor(as.character(sce$annots))
    )

    sce$Age <- as.numeric(sce$Age)
    sce$Sex <- droplevels(factor(sce$Sex))
    sce$slide_id <- droplevels(factor(sce$slide_id))
    sce$BrNum <- droplevels(factor(sce$BrNum))

    sce$registration_sample_id <- sce$BrNum

    # Paper script creates this although it is not used in the model.
    sce$predictions_smooth <- factor(sce$annots)


    if (renormalize) {

        section(paste(label, "- FILTER + NORMALIZE"))

        n_before <- nrow(sce)

        keep <- edgeR::filterByExpr(
            sce,
            group = sce$annots
        )

        sce <- sce[
            keep,
            ,
            drop = FALSE
        ]

        cat("Genes before:", n_before, "\n")
        cat("Genes retained:", nrow(sce), "\n")
        cat("Genes filtered:", n_before - nrow(sce), "\n")

        logcounts(sce) <- edgeR::cpm(
            edgeR::calcNormFactors(sce),
            log = TRUE,
            prior.count = 1
        )

    } else {

        cat(
            "Using stored Paper N24 logcounts without renormalization.\n"
        )
    }


    # Exact paper Step-02 model
    dx_mod <- registration_model(
        sce,
        covars = c(
            "annots",
            "Age",
            "Sex",
            "slide_id"
        ),
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


    dx_block_cor <- registration_block_cor(
        sce,
        registration_model = dx_mod,
        var_sample_id = "BrNum"
    )


    cat(
        "Block correlation:",
        dx_block_cor,
        "\n"
    )


    dx_res <- registration_stats_enrichment(
        sce,
        block_cor = dx_block_cor,
        covars = c(
            "annots",
            "Age",
            "Sex",
            "slide_id"
        ),
        var_registration = "Dx",
        var_sample_id = "BrNum",
        gene_ensembl = "ID",
        gene_name = "Symbol"
    )


    cat(
        "Genes tested:",
        nrow(dx_res),
        "\n"
    )

    cat(
        "Nominal P<0.05:",
        sum(
            dx_res$p_value_SCZ < 0.05,
            na.rm = TRUE
        ),
        "\n"
    )

    cat(
        "FDR<0.10:",
        sum(
            dx_res$fdr_SCZ < 0.10,
            na.rm = TRUE
        ),
        "\n"
    )

    cat(
        "FDR<0.05:",
        sum(
            dx_res$fdr_SCZ < 0.05,
            na.rm = TRUE
        ),
        "\n"
    )


    list(
        sce = sce,
        result = dx_res,
        block_cor = dx_block_cor
    )
}


# =============================================================================
# LOAD EXACT PAPER N24
# =============================================================================

section("A. EXACT PAPER N24")

paper24 <- readRDS(
    PAPER_RDS
)

cat(
    "Dimensions:",
    dim(paper24),
    "\n"
)

cat(
    "Donors:",
    length(unique(paper24$BrNum)),
    "\n"
)

cat(
    "Cell types:",
    length(unique(paper24$annots)),
    "\n"
)


# Run exact stored N24 object using its stored logcounts.
paper24_run <- run_paper_de(
    paper24,
    "Paper_official_source_N24",
    renormalize = FALSE
)


# =============================================================================
# VALIDATE PAPER N24 AGAINST OFFICIAL GITHUB RESULT
# =============================================================================

section("A2. PAPER N24 REPRODUCTION CHECK")

official <- read.csv(
    PAPER_OFFICIAL_CSV,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

rerun <- paper24_run$result


official <- official[
    match(
        rerun$gene,
        official$gene
    ),
    ,
    drop = FALSE
]


if (any(is.na(official$gene))) {
    stop(
        "Official Paper N24 gene alignment failed."
    )
}


check_cols <- c(
    "t_stat_SCZ",
    "p_value_SCZ",
    "fdr_SCZ",
    "logFC_SCZ"
)


official_matrix <- as.matrix(
    official[, check_cols]
)

rerun_matrix <- as.matrix(
    rerun[, check_cols]
)


paper24_maxdiff <- max(
    abs(
        official_matrix -
        rerun_matrix
    ),
    na.rm = TRUE
)


paper24_logfc_cor <- cor(
    official$logFC_SCZ,
    rerun$logFC_SCZ
)

paper24_t_cor <- cor(
    official$t_stat_SCZ,
    rerun$t_stat_SCZ
)


cat(
    "Official N24 vs rerun max abs difference:",
    paper24_maxdiff,
    "\n"
)

cat(
    "Official N24 vs rerun logFC Pearson:",
    paper24_logfc_cor,
    "\n"
)

cat(
    "Official N24 vs rerun t Pearson:",
    paper24_t_cor,
    "\n"
)


write.csv(
    rerun,
    file.path(
        OUT,
        "A_paper_rerun_N24_celltype_adjusted_DE.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# B. PAPER-SOURCE N23
# =============================================================================

section("B. PAPER-SOURCE N23")

meta23 <- read.csv(
    META23_FILE,
    stringsAsFactors = FALSE
)

donors23 <- as.character(
    meta23$BrNum
)


paper23 <- paper24[
    ,
    as.character(paper24$BrNum) %in% donors23
]


if (ncol(paper23) != 276) {
    stop(
        "Paper-source N23 expected 276 pseudobulks."
    )
}


# IMPORTANT:
# We subset counts to N23 and recompute filtering/TMM/logCPM.
paper23_run <- run_paper_de(
    paper23,
    "Paper_source_N23",
    renormalize = TRUE
)


write.csv(
    paper23_run$result,
    file.path(
        OUT,
        "B_paper_source_N23_celltype_adjusted_DE.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# HELPERS FOR GENEBRIDGE MATRICES
# =============================================================================

gb_meta <- read.csv(
    META_FILE,
    stringsAsFactors = FALSE
)


paper_gene_map <- data.frame(
    Symbol = as.character(
        rowData(paper24)$Symbol
    ),
    ID = as.character(
        rowData(paper24)$ID
    ),
    stringsAsFactors = FALSE
)


build_sce_from_csv <- function(
    matrix_file,
    expected_genes,
    label
) {

    section(
        paste("LOAD MATRIX:", label)
    )

    pb <- read.csv(
        matrix_file,
        row.names = 1,
        check.names = FALSE
    )


    if (
        !all(
            dim(pb) ==
            c(276, expected_genes)
        )
    ) {
        stop(
            label,
            ": unexpected matrix dimensions: ",
            paste(dim(pb), collapse = " x ")
        )
    }


    idx <- match(
        rownames(pb),
        gb_meta$pseudobulk_id
    )


    if (any(is.na(idx))) {
        stop(
            label,
            ": metadata alignment failed."
        )
    }


    md <- gb_meta[
        idx,
        ,
        drop = FALSE
    ]


    if (
        !identical(
            rownames(pb),
            as.character(md$pseudobulk_id)
        )
    ) {
        stop(
            label,
            ": pseudobulk ordering mismatch."
        )
    }


    count_mat <- t(
        as.matrix(pb)
    )

    storage.mode(
        count_mat
    ) <- "double"


    genes <- rownames(
        count_mat
    )


    ens <- paper_gene_map$ID[
        match(
            genes,
            paper_gene_map$Symbol
        )
    ]


    # Full ENVI contains genes outside the 300-gene panel.
    # Preserve gene symbol as ID when Ensembl mapping is not in paper panel.
    ens[
        is.na(ens)
    ] <- genes[
        is.na(ens)
    ]


    sce <- SingleCellExperiment(
        assays = list(
            counts = count_mat
        ),
        rowData = DataFrame(
            ID = ens,
            Symbol = genes,
            row.names = genes
        ),
        colData = DataFrame(
            md
        )
    )


    sce
}


# =============================================================================
# C. GENEBRIDGE-BEFORE N23
# =============================================================================

section("C. GENEBRIDGE-BEFORE N23")

before_sce <- build_sce_from_csv(
    file.path(
        PB_DIR,
        "original_xenium_300gene_donor_celltype_pseudobulk.csv.gz"
    ),
    300,
    "GeneBridge_before_N23"
)


before_run <- run_paper_de(
    before_sce,
    "GeneBridge_before_N23",
    renormalize = TRUE
)


write.csv(
    before_run$result,
    file.path(
        OUT,
        "C_GeneBridge_before_N23_celltype_adjusted_DE.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# D. ENVI-MEASURED N23
# =============================================================================

section("D. ENVI-MEASURED N23")

measured_sce <- build_sce_from_csv(
    file.path(
        PB_DIR,
        "ENVI_measured300_donor_celltype_pseudobulk.csv.gz"
    ),
    300,
    "ENVI_measured_N23"
)


measured_run <- run_paper_de(
    measured_sce,
    "ENVI_measured_N23",
    renormalize = TRUE
)


write.csv(
    measured_run$result,
    file.path(
        OUT,
        "D_ENVI_measured_N23_celltype_adjusted_DE.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# BEFORE vs ENVI-MEASURED MUST BE EXACT
# =============================================================================

section("C vs D TECHNICAL CONTROL")

b <- before_run$result[
    order(before_run$result$gene),
    ,
    drop = FALSE
]

m <- measured_run$result[
    order(measured_run$result$gene),
    ,
    drop = FALSE
]


if (!identical(b$gene, m$gene)) {
    stop(
        "Before vs ENVI-measured gene alignment failed."
    )
}


tech_maxdiff <- max(
    abs(
        as.matrix(
            b[, check_cols]
        ) -
        as.matrix(
            m[, check_cols]
        )
    ),
    na.rm = TRUE
)


cat(
    "Before vs ENVI-measured max abs DE difference:",
    tech_maxdiff,
    "\n"
)


if (tech_maxdiff != 0) {
    stop(
        "Before and ENVI-measured DE are not exact."
    )
}


cat(
    "GeneBridge-before == ENVI-measured DE: EXACT\n"
)


# =============================================================================
# B vs C PAPER-SOURCE N23 CONCORDANCE
# =============================================================================

section("B vs C PAPER-SOURCE vs BEFORE N23")

p23 <- paper23_run$result[
    order(paper23_run$result$gene),
    ,
    drop = FALSE
]

gb <- before_run$result[
    order(before_run$result$gene),
    ,
    drop = FALSE
]


common <- intersect(
    p23$gene,
    gb$gene
)


p23c <- p23[
    match(common, p23$gene),
    ,
    drop = FALSE
]

gbc <- gb[
    match(common, gb$gene),
    ,
    drop = FALSE
]


paper_before_logfc_cor <- cor(
    p23c$logFC_SCZ,
    gbc$logFC_SCZ
)

paper_before_t_cor <- cor(
    p23c$t_stat_SCZ,
    gbc$t_stat_SCZ
)


cat(
    "Common tested genes:",
    length(common),
    "\n"
)

cat(
    "Paper-source N23 vs Before logFC Pearson:",
    paper_before_logfc_cor,
    "\n"
)

cat(
    "Paper-source N23 vs Before t Pearson:",
    paper_before_t_cor,
    "\n"
)


# =============================================================================
# E. ENVI-AFTER FULL N23
# =============================================================================

section("E. ENVI-AFTER FULL N23")

full_sce <- build_sce_from_csv(
    file.path(
        PB_DIR,
        "ENVI_full34987_donor_celltype_pseudobulk_countscale.csv.gz"
    ),
    34987,
    "ENVI_after_full_N23"
)


full_run <- run_paper_de(
    full_sce,
    "ENVI_after_full_N23",
    renormalize = TRUE
)


gene_info <- read.csv(
    GENE_INFO_FILE,
    stringsAsFactors = FALSE
)


full_res <- full_run$result


full_res$expression_source <- gene_info$expression_source[
    match(
        full_res$gene,
        gene_info$gene
    )
]


if (any(is.na(full_res$expression_source))) {
    stop(
        "Could not map ENVI expression_source."
    )
}


full_measured <- full_res[
    full_res$expression_source ==
    "measured_xenium",
    ,
    drop = FALSE
]

full_imputed <- full_res[
    full_res$expression_source ==
    "envi_imputed",
    ,
    drop = FALSE
]


write.csv(
    full_res,
    file.path(
        OUT,
        "E_ENVI_after_full_N23_celltype_adjusted_DE.csv"
    ),
    row.names = FALSE
)

write.csv(
    full_measured,
    file.path(
        OUT,
        "E1_ENVI_after_full_measured_subset_N23_DE.csv"
    ),
    row.names = FALSE
)

write.csv(
    full_imputed,
    file.path(
        OUT,
        "E2_ENVI_after_full_imputed_subset_N23_DE.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# SUMMARY
# =============================================================================

section("FINAL PAPER vs BEFORE vs AFTER SUMMARY")


make_summary <- function(
    arm,
    cohort,
    result,
    block_cor
) {

    data.frame(
        arm = arm,
        cohort = cohort,
        genes_tested = nrow(result),
        block_correlation = block_cor,
        nominal_p_lt_0.05 =
            sum(
                result$p_value_SCZ < 0.05,
                na.rm = TRUE
            ),
        fdr_lt_0.10 =
            sum(
                result$fdr_SCZ < 0.10,
                na.rm = TRUE
            ),
        fdr_lt_0.05 =
            sum(
                result$fdr_SCZ < 0.05,
                na.rm = TRUE
            )
    )
}


summary_df <- rbind(

    make_summary(
        "Paper_official_N24_rerun",
        "N24",
        paper24_run$result,
        paper24_run$block_cor
    ),

    make_summary(
        "Paper_source_N23",
        "N23",
        paper23_run$result,
        paper23_run$block_cor
    ),

    make_summary(
        "GeneBridge_before_N23",
        "N23",
        before_run$result,
        before_run$block_cor
    ),

    make_summary(
        "ENVI_measured_N23",
        "N23",
        measured_run$result,
        measured_run$block_cor
    ),

    make_summary(
        "ENVI_after_full_N23",
        "N23",
        full_res,
        full_run$block_cor
    ),

    make_summary(
        "ENVI_after_imputed_only",
        "N23",
        full_imputed,
        full_run$block_cor
    )
)


write.csv(
    summary_df,
    file.path(
        OUT,
        "paper_before_after_global_celltype_DE_summary.csv"
    ),
    row.names = FALSE
)


print(
    summary_df,
    row.names = FALSE
)


cat("\nOfficial N24 validation:\n")
cat(
    "  max abs difference:",
    paper24_maxdiff,
    "\n"
)
cat(
    "  logFC Pearson:",
    paper24_logfc_cor,
    "\n"
)
cat(
    "  t Pearson:",
    paper24_t_cor,
    "\n"
)


cat("\nPaper-source N23 vs Before:\n")
cat(
    "  logFC Pearson:",
    paper_before_logfc_cor,
    "\n"
)
cat(
    "  t Pearson:",
    paper_before_t_cor,
    "\n"
)


cat("\nBefore vs ENVI measured:\n")
cat(
    "  max difference:",
    tech_maxdiff,
    "\n"
)


cat(
    "\nFull ENVI measured genes retained:",
    nrow(full_measured),
    "\n"
)

cat(
    "Full ENVI imputed genes retained:",
    nrow(full_imputed),
    "\n"
)

cat(
    "Full ENVI imputed FDR<0.10:",
    sum(
        full_imputed$fdr_SCZ < 0.10,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Full ENVI imputed FDR<0.05:",
    sum(
        full_imputed$fdr_SCZ < 0.05,
        na.rm = TRUE
    ),
    "\n"
)


cat(
    "\nFINAL STATUS: GLOBAL CELL-TYPE PAPER-vs-BEFORE-vs-AFTER DE PASS\n"
)
