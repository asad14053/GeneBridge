#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(S4Vectors)
    library(edgeR)
    library(spatialLIBD)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: 03_layer_adjusted_de.R <original|envi_measured>")
}

dataset <- args[1]

if (!dataset %in% c("original", "envi_measured")) {
    stop("dataset must be original or envi_measured")
}


ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance/layer_adjusted_de",
    dataset
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


if (dataset == "original") {

    counts_file <- file.path(
        ROOT,
        "outputs/imputation_full/DE/paper_concordance",
        "original_layer_pseudobulk",
        "original_xenium_300gene_donor_layer_pseudobulk_counts.csv.gz"
    )

    meta_file <- file.path(
        ROOT,
        "outputs/imputation_full/DE/paper_concordance",
        "original_layer_pseudobulk",
        "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
    )

} else {

    counts_file <- file.path(
        ROOT,
        "outputs/imputation_full/DE/paper_concordance",
        "envi_layer_pseudobulk",
        "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
    )

    meta_file <- file.path(
        ROOT,
        "outputs/imputation_full/DE/paper_concordance",
        "envi_layer_pseudobulk",
        "ENVI_donor_layer_pseudobulk_metadata.csv"
    )
}


cat(strrep("=", 100), "\n", sep = "")
cat("PAPER-MATCHED XENIUM LAYER-ADJUSTED SCZ DE\n")
cat(strrep("=", 100), "\n", sep = "")

cat("Dataset:", dataset, "\n")
cat("Counts :", counts_file, "\n")
cat("Meta   :", meta_file, "\n")


# =============================================================================
# Read pseudobulk data
# =============================================================================

pb <- read.csv(
    counts_file,
    row.names = 1,
    check.names = FALSE
)

meta <- read.csv(
    meta_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat("\nRaw pseudobulk shape:\n")
cat("samples =", nrow(pb), "\n")
cat("genes   =", ncol(pb), "\n")


stopifnot(nrow(pb) == 161)
stopifnot(ncol(pb) == 300)
stopifnot(nrow(meta) == 161)


if (!identical(
    rownames(pb),
    meta$pseudobulk_id
)) {
    stop("Pseudobulk count rows and metadata are not in identical order.")
}


# =============================================================================
# Metadata validation
# =============================================================================

required_cols <- c(
    "pseudobulk_id",
    "BrNum",
    "Dx",
    "Age",
    "Sex",
    "slide_id",
    "predictions_smooth",
    "layer_annotation",
    "n_cells"
)

missing_cols <- setdiff(
    required_cols,
    colnames(meta)
)

if (length(missing_cols) > 0) {
    stop(
        "Missing metadata columns: ",
        paste(missing_cols, collapse = ", ")
    )
}


stopifnot(length(unique(meta$BrNum)) == 23)
stopifnot(length(unique(meta$predictions_smooth)) == 7)

stopifnot(
    all(
        table(meta$BrNum) == 7
    )
)

stopifnot(
    all(meta$n_cells >= 10)
)


# Important: NTC reference / SCZ effect.
meta$Dx <- factor(
    meta$Dx,
    levels = c("NTC", "SCZ")
)

meta$predictions_smooth <- factor(
    meta$predictions_smooth,
    levels = c(
        "spd01",
        "spd02",
        "spd03",
        "spd04",
        "spd05",
        "spd06",
        "spd07"
    )
)

meta$Age <- as.numeric(
    meta$Age
)

meta$Sex <- factor(
    meta$Sex
)

meta$slide_id <- factor(
    meta$slide_id
)

meta$BrNum <- factor(
    meta$BrNum
)


cat("\nDiagnosis by donor:\n")

donor_dx <- unique(
    meta[, c("BrNum", "Dx")]
)

print(
    table(donor_dx$Dx)
)


cat("\nSpatial-domain samples:\n")

print(
    table(meta$predictions_smooth)
)


# =============================================================================
# Create pseudobulk SingleCellExperiment
#
# pb is samples × genes.
# SCE expects genes × samples.
# =============================================================================

count_matrix <- t(
    as.matrix(pb)
)

storage.mode(
    count_matrix
) <- "double"


if (any(!is.finite(count_matrix))) {
    stop("Non-finite pseudobulk values.")
}

if (any(count_matrix < 0)) {
    stop("Negative pseudobulk values.")
}


sce <- SingleCellExperiment(
    assays = list(
        counts = count_matrix
    ),
    colData = DataFrame(
        meta,
        row.names = meta$pseudobulk_id
    )
)


rownames(sce) <- colnames(pb)

rowData(sce)$ID <- rownames(sce)
rowData(sce)$Symbol <- rownames(sce)


# =============================================================================
# Match registration_pseudobulk() post-aggregation behavior
#
# Published spatialLIBD workflow:
#   edgeR::filterByExpr(group = spatial-domain variable)
#   edgeR::calcNormFactors()
#   edgeR::cpm(log=TRUE, prior.count=1)
# =============================================================================

cat("\n", strrep("=", 100), "\n", sep = "")
cat("FILTER LOWLY EXPRESSED GENES\n")
cat(strrep("=", 100), "\n", sep = "")


keep_expr <- edgeR::filterByExpr(
    sce,
    group = sce$predictions_smooth
)


cat(
    "Genes before filterByExpr:",
    nrow(sce),
    "\n"
)

cat(
    "Genes retained:",
    sum(keep_expr),
    "\n"
)

cat(
    "Genes removed:",
    sum(!keep_expr),
    "\n"
)


filter_table <- data.frame(
    gene = rownames(sce),
    retained = keep_expr
)

write.csv(
    filter_table,
    file.path(
        OUT,
        paste0(
            dataset,
            "_filterByExpr_300genes.csv"
        )
    ),
    row.names = FALSE
)


sce <- sce[
    keep_expr,
]


# =============================================================================
# Exact spatialLIBD-style normalization
# =============================================================================

cat("\n", strrep("=", 100), "\n", sep = "")
cat("TMM + LOG2-CPM NORMALIZATION\n")
cat(strrep("=", 100), "\n", sep = "")


dge <- edgeR::calcNormFactors(
    sce
)

logcounts(sce) <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 1
)


if (any(!is.finite(logcounts(sce)))) {
    stop("Non-finite log2-CPM values.")
}


cat(
    "Normalized matrix:",
    nrow(sce),
    "genes ×",
    ncol(sce),
    "pseudobulks\n"
)


# =============================================================================
# Paper model
#
# Published:
#
# registration_model(
#   spe_pseudo,
#   covars = c(
#       "predictions_smooth",
#       "Age",
#       "Sex",
#       "slide_id"
#   ),
#   var_registration = "Dx"
# )
#
# block by BrNum.
# =============================================================================

cat("\n", strrep("=", 100), "\n", sep = "")
cat("BUILD PAPER-MATCHED MODEL\n")
cat(strrep("=", 100), "\n", sep = "")


covars <- c(
    "predictions_smooth",
    "Age",
    "Sex",
    "slide_id"
)


dx_mod <- registration_model(
    sce,
    covars = covars,
    var_registration = "Dx"
)


cat("\nModel matrix dimensions:\n")
print(
    dim(dx_mod)
)

cat("\nModel columns:\n")
print(
    colnames(dx_mod)
)

cat(
    "\nModel rank:",
    qr(dx_mod)$rank,
    "/",
    ncol(dx_mod),
    "\n"
)


if (qr(dx_mod)$rank != ncol(dx_mod)) {
    stop("DE model is rank deficient.")
}


# =============================================================================
# Repeated donor / duplicate correlation
# =============================================================================

cat("\n", strrep("=", 100), "\n", sep = "")
cat("ESTIMATE WITHIN-DONOR BLOCK CORRELATION\n")
cat(strrep("=", 100), "\n", sep = "")


dx_block_cor <- registration_block_cor(
    sce,
    registration_model = dx_mod,
    var_sample_id = "BrNum"
)


cat(
    "\nEstimated donor block correlation:",
    dx_block_cor,
    "\n"
)


if (!is.finite(dx_block_cor)) {
    stop("Block correlation is not finite. Do not continue with paper-matched analysis.")
}


# =============================================================================
# SCZ enrichment / differential expression
# =============================================================================

cat("\n", strrep("=", 100), "\n", sep = "")
cat("RUN SCZ DIFFERENTIAL EXPRESSION\n")
cat(strrep("=", 100), "\n", sep = "")


dx_res <- registration_stats_enrichment(
    sce,
    block_cor = dx_block_cor,
    covars = covars,
    var_registration = "Dx",
    var_sample_id = "BrNum",
    gene_ensembl = "ID",
    gene_name = "Symbol"
)


required_result_cols <- c(
    "gene",
    "logFC_SCZ",
    "t_stat_SCZ",
    "p_value_SCZ",
    "fdr_SCZ"
)

missing_result_cols <- setdiff(
    required_result_cols,
    colnames(dx_res)
)

if (length(missing_result_cols) > 0) {

    cat("\nResult columns returned:\n")
    print(colnames(dx_res))

    stop(
        "Expected result columns missing: ",
        paste(
            missing_result_cols,
            collapse = ", "
        )
    )
}


dx_res <- dx_res[
    order(dx_res$fdr_SCZ),
]


write.csv(
    dx_res,
    file.path(
        OUT,
        paste0(
            dataset,
            "_layer_adjusted_SCZ_DE_all.csv"
        )
    ),
    row.names = FALSE
)


# =============================================================================
# Significance tables
# =============================================================================

sig10 <- dx_res[
    dx_res$fdr_SCZ < 0.10,
]

sig05 <- dx_res[
    dx_res$fdr_SCZ < 0.05,
]

nominal <- dx_res[
    dx_res$p_value_SCZ < 0.05,
]


write.csv(
    sig10,
    file.path(
        OUT,
        paste0(
            dataset,
            "_layer_adjusted_SCZ_DE_FDR10.csv"
        )
    ),
    row.names = FALSE
)

write.csv(
    sig05,
    file.path(
        OUT,
        paste0(
            dataset,
            "_layer_adjusted_SCZ_DE_FDR05.csv"
        )
    ),
    row.names = FALSE
)


# =============================================================================
# Summary
# =============================================================================

summary_df <- data.frame(
    dataset = dataset,
    donors = length(unique(meta$BrNum)),
    donor_layer_samples = ncol(sce),
    genes_input = 300,
    genes_after_filterByExpr = nrow(sce),
    block_correlation = dx_block_cor,
    nominal_p05 = nrow(nominal),
    fdr_lt_010 = nrow(sig10),
    fdr_lt_005 = nrow(sig05),
    up_fdr010 = sum(
        sig10$logFC_SCZ > 0
    ),
    down_fdr010 = sum(
        sig10$logFC_SCZ < 0
    )
)


write.csv(
    summary_df,
    file.path(
        OUT,
        paste0(
            dataset,
            "_layer_adjusted_SCZ_DE_summary.csv"
        )
    ),
    row.names = FALSE
)


cat("\n", strrep("=", 100), "\n", sep = "")
cat("FINAL SUMMARY\n")
cat(strrep("=", 100), "\n", sep = "")

print(
    summary_df,
    row.names = FALSE
)


cat("\nTop 20 genes by FDR:\n")

print(
    head(
        dx_res[
            ,
            c(
                "gene",
                "logFC_SCZ",
                "t_stat_SCZ",
                "p_value_SCZ",
                "fdr_SCZ"
            )
        ],
        20
    ),
    row.names = FALSE
)


cat("\nFDR < 0.10 genes:\n")

if (nrow(sig10) == 0) {

    cat("NONE\n")

} else {

    print(
        sig10[
            ,
            c(
                "gene",
                "logFC_SCZ",
                "p_value_SCZ",
                "fdr_SCZ"
            )
        ],
        row.names = FALSE
    )
}


cat(
    "\nFINAL STATUS: PASS — ",
    dataset,
    " paper-matched layer-adjusted SCZ DE completed.\n",
    sep = ""
)
