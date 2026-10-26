#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(edgeR)
    library(SingleCellExperiment)
    library(S4Vectors)
    library(spatialLIBD)
})

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PB_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_preimputation_baseline/domain_pseudobulk",
    "original_xenium_N24_donor_domain_pseudobulk_300genes.csv.gz"
)

META_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_preimputation_baseline/domain_pseudobulk",
    "original_xenium_N24_donor_domain_pseudobulk_metadata.csv"
)

OUT_ROOT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_preimputation_baseline/domain_de"
)

dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)


section <- function(x) {
    cat("\n")
    cat(paste(rep("=", 110), collapse=""), "\n")
    cat(x, "\n")
    cat(paste(rep("=", 110), collapse=""), "\n")
}


# ==========================================================================================
# LOAD ORIGINAL MEASURED PSEUDOBULK
# ==========================================================================================

section("ORIGINAL MEASURED XENIUM — DONOR × DOMAIN DE")

pb <- read.csv(
    PB_FILE,
    check.names=FALSE
)

meta <- read.csv(
    META_FILE,
    stringsAsFactors=FALSE,
    check.names=FALSE
)


rownames(pb) <- pb$pseudobulk_id
pb$pseudobulk_id <- NULL

rownames(meta) <- meta$pseudobulk_id


# Ensure identical order
meta <- meta[rownames(pb), , drop=FALSE]

stopifnot(
    identical(
        rownames(pb),
        rownames(meta)
    )
)


cat("PB matrix :", nrow(pb), "samples x", ncol(pb), "genes\n")
cat("Metadata  :", nrow(meta), "rows\n")


# ==========================================================================================
# VALIDATE DESIGN
# ==========================================================================================

stopifnot(nrow(pb) == 168)
stopifnot(ncol(pb) == 300)
stopifnot(length(unique(meta$BrNum)) == 24)
stopifnot(length(unique(meta$predictions_smooth)) == 7)

ntc_n <- length(unique(meta$BrNum[meta$Dx == "NTC"]))
scz_n <- length(unique(meta$BrNum[meta$Dx == "SCZ"]))

stopifnot(ntc_n == 12)
stopifnot(scz_n == 12)


cat("Donors  :", length(unique(meta$BrNum)), "\n")
cat("NTC     :", ntc_n, "\n")
cat("SCZ     :", scz_n, "\n")
cat("Domains :", length(unique(meta$predictions_smooth)), "\n")
cat("Slides  :", length(unique(meta$slide_id)), "\n")


# ==========================================================================================
# GENE × PSEUDOBULK COUNT MATRIX
# ==========================================================================================

counts <- t(
    as.matrix(pb)
)

storage.mode(counts) <- "numeric"

if (any(counts < 0, na.rm=TRUE)) {
    stop("Negative values found in measured pseudobulk counts.")
}


# ==========================================================================================
# GITHUB / spatialLIBD PREPROCESSING
#
# registration_pseudobulk() performs:
#   filterByExpr
#   TMM
#   log2 CPM
#
# registration variable for donor-domain PB = spatial domain.
# ==========================================================================================

section("GITHUB-COMPATIBLE FILTER + TMM + LOG2 CPM")

dge <- edgeR::DGEList(
    counts=counts
)

keep <- edgeR::filterByExpr(
    dge,
    group=factor(meta$predictions_smooth)
)

genes_before <- nrow(dge)
genes_after  <- sum(keep)
genes_drop   <- genes_before - genes_after

cat("Genes before filter :", genes_before, "\n")
cat("Genes retained      :", genes_after, "\n")
cat("Genes dropped       :", genes_drop, "\n")


dge <- dge[
    keep,
    ,
    keep.lib.sizes=FALSE
]

dge <- edgeR::calcNormFactors(
    dge,
    method="TMM"
)

logcpm <- edgeR::cpm(
    dge,
    log=TRUE,
    prior.count=1
)


# ==========================================================================================
# METADATA TYPES
# ==========================================================================================

meta$Dx <- factor(
    meta$Dx,
    levels=c("NTC", "SCZ")
)

meta$predictions_smooth <- factor(
    meta$predictions_smooth
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

# Explicit repeated-measure donor field for provenance.
meta$registration_sample_id <- meta$BrNum


# ==========================================================================================
# BUILD SCE
# ==========================================================================================

rd <- DataFrame(
    ID=rownames(logcpm),
    Symbol=rownames(logcpm)
)

rownames(rd) <- rownames(logcpm)

sce <- SingleCellExperiment(
    assays=list(
        counts=as.matrix(dge$counts),
        logcounts=logcpm
    ),
    rowData=rd,
    colData=DataFrame(meta)
)


# ==========================================================================================
# FUNCTION: RUN ONE DE MODEL
# ==========================================================================================

run_de <- function(
    label,
    covars,
    out_dir
) {

    section(
        paste0(
            "ORIGINAL XENIUM DOMAIN DE — ",
            label
        )
    )

    dir.create(
        out_dir,
        recursive=TRUE,
        showWarnings=FALSE
    )


    cat(
        "Covariates:",
        paste(covars, collapse=" + "),
        "\n"
    )

    cat(
        "Repeated-measure block: BrNum\n"
    )


    # ----------------------------------------------------------------------
    # GitHub:
    #
    # registration_model(
    #   covars = c("predictions_smooth","Age","Sex","slide_id"),
    #   var_registration = "Dx"
    # )
    # ----------------------------------------------------------------------

    dx_mod <- spatialLIBD::registration_model(
        sce,
        covars=covars,
        var_registration="Dx"
    )


    # Explicit BrNum because this SCE was reconstructed manually,
    # rather than created directly by registration_pseudobulk().
    dx_block_cor <- spatialLIBD::registration_block_cor(
        sce,
        registration_model=dx_mod,
        var_sample_id="BrNum"
    )


    dx_res <- spatialLIBD::registration_stats_enrichment(
        sce,
        block_cor=dx_block_cor,
        covars=covars,
        var_registration="Dx",
        var_sample_id="BrNum",
        gene_ensembl="ID",
        gene_name="Symbol"
    )


    result_file <- file.path(
        out_dir,
        "original_xenium_domain_SCZ_DE_dx_results.csv"
    )

    write.csv(
        dx_res,
        result_file,
        row.names=FALSE
    )


    sig <- dx_res[
        !is.na(dx_res$fdr_SCZ) &
        dx_res$fdr_SCZ < 0.10,
        ,
        drop=FALSE
    ]


    n_up <- sum(
        sig$logFC_SCZ > 0,
        na.rm=TRUE
    )

    n_down <- sum(
        sig$logFC_SCZ < 0,
        na.rm=TRUE
    )


    cat("\nBlock correlation :", dx_block_cor, "\n")
    cat("Genes tested       :", nrow(dx_res), "\n")
    cat("FDR < 0.10         :", nrow(sig), "\n")
    cat("  UP               :", n_up, "\n")
    cat("  DOWN             :", n_down, "\n")

    cat("\nSignificant genes:\n")

    if (nrow(sig) == 0) {

        cat("NONE\n")

    } else {

        print(
            sig[
                order(sig$fdr_SCZ),
                c(
                    "gene",
                    "logFC_SCZ",
                    "p_value_SCZ",
                    "fdr_SCZ"
                ),
                drop=FALSE
            ],
            row.names=FALSE
        )
    }


    summary <- data.frame(
        analysis=label,
        genes_before_filter=genes_before,
        genes_tested=nrow(dx_res),
        genes_dropped=genes_drop,
        block_correlation=dx_block_cor,
        FDR_lt_010=nrow(sig),
        up_FDR010=n_up,
        down_FDR010=n_down,
        result_file=result_file
    )


    write.csv(
        summary,
        file.path(
            out_dir,
            "original_xenium_domain_SCZ_DE_summary.csv"
        ),
        row.names=FALSE
    )


    cat("\nOutput:\n")
    cat(result_file, "\n")


    return(summary)
}


# ==========================================================================================
# A. PAPER-COMPATIBLE WITH SLIDE
# ==========================================================================================

with_slide <- run_de(
    label="WITH_SLIDE",
    covars=c(
        "predictions_smooth",
        "Age",
        "Sex",
        "slide_id"
    ),
    out_dir=file.path(
        OUT_ROOT,
        "with_slide"
    )
)


# ==========================================================================================
# B. SENSITIVITY — NO SLIDE
# ==========================================================================================

no_slide <- run_de(
    label="NO_SLIDE",
    covars=c(
        "predictions_smooth",
        "Age",
        "Sex"
    ),
    out_dir=file.path(
        OUT_ROOT,
        "no_slide"
    )
)


# ==========================================================================================
# COMBINED SUMMARY
# ==========================================================================================

combined <- rbind(
    with_slide,
    no_slide
)

section("ORIGINAL XENIUM DOMAIN DE — FINAL SUMMARY")

print(
    combined,
    row.names=FALSE
)

write.csv(
    combined,
    file.path(
        OUT_ROOT,
        "original_xenium_domain_DE_both_models_summary.csv"
    ),
    row.names=FALSE
)


section("FINAL STATUS")

cat(
    "PASS — ORIGINAL MEASURED XENIUM DOMAIN DE COMPLETE\n"
)
