#!/usr/bin/env Rscript

# ================================================================================================
# 05_celltype_adjusted_de_nofilter_all_models.R
#
# TASK 2.3 — DONOR × CELL-TYPE ADJUSTED SCZ DE
#
# PURPOSE
# -------
# Perform SCZ vs NTC differential-expression analysis using donor × cell-type
# pseudobulk expression generated from NTC-trained OOF-imputed count_scale values.
#
# REFERENCE WORKFLOW
# ------------------
# LieberInstitute/spatialDLPFC_SCZ_XENIUM
#
# code/analysis/07_cell_type_de/
#   01_create_cell_type_pseudobulk_spe.R
#   02_cell_type_pseudobulk_de.R
#
# GITHUB MODEL
# ------------
# var_registration = Dx
#
# covariates:
#   annots
#   Age
#   Sex
#   slide_id
#
# repeated-measures donor:
#   BrNum
#
# Preprocessing equivalent to registration_pseudobulk():
#   filterByExpr
#   TMM normalization
#   log2 CPM, prior.count = 1
#
# INPUT
# -----
# 24 donors
# 12 NTC + 12 SCZ
# 12 cell types
# 288 donor × cell-type pseudobulks/model
# 300 OOF genes before expression filtering
#
# MODELS
# ------
# VISTA
# gimVI
# ENVI
# TransImp
# Tangram
# SpaGE
#
# NOTE
# ----
# run_date is carried in pseudobulk metadata but is NOT included in the
# GitHub cell-type-adjusted DE model.
#
# ================================================================================================


suppressPackageStartupMessages({

    library(edgeR)
    library(limma)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(spatialLIBD)
})


# ================================================================================================
# CONFIGURATION
# ================================================================================================

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PB_ROOT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_ntc_dysregulation/celltype_pseudobulk"
)

OUT_ROOT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_ntc_dysregulation/celltype_adjusted_de_nofilter_noslide"
)

dir.create(
    OUT_ROOT,
    recursive = TRUE,
    showWarnings = FALSE
)


MODELS <- c(
    "VISTA",
    "gimVI",
    "ENVI",
    "TransImp",
    "Tangram",
    "SpaGE"
)


EXPECTED_DONORS <- 24
EXPECTED_NTC <- 12
EXPECTED_SCZ <- 12
EXPECTED_CELLTYPES <- 12
EXPECTED_PB <- 288
EXPECTED_GENES <- 300


# ================================================================================================
# HELPERS
# ================================================================================================

section <- function(title) {

    cat(
        "\n",
        paste(rep("=", 110), collapse = ""),
        "\n",
        title,
        "\n",
        paste(rep("=", 110), collapse = ""),
        "\n",
        sep = ""
    )
}


get_unique_donor_table <- function(meta) {

    unique(
        meta[
            ,
            c(
                "BrNum",
                "Dx",
                "Age",
                "Sex",
                "slide_id",
                "run_date"
            ),
            drop = FALSE
        ]
    )
}


# ================================================================================================
# GLOBAL SUMMARY
# ================================================================================================

summary_list <- list()


section(
    "TASK 2.3 — DONOR × CELL-TYPE ADJUSTED DE / ALL SIX MODELS"
)


# ================================================================================================
# MODEL LOOP
# ================================================================================================

for (model in MODELS) {


    section(
        paste0(
            "MODEL: ",
            model
        )
    )


    # ============================================================================================
    # INPUT PATHS
    # ============================================================================================

    pb_file <- file.path(
        PB_ROOT,
        model,
        paste0(
            model,
            "_NTCtrained_OOF_donor_celltype_pseudobulk_300genes.csv.gz"
        )
    )


    meta_file <- file.path(
        PB_ROOT,
        model,
        paste0(
            model,
            "_NTCtrained_OOF_donor_celltype_pseudobulk_metadata.csv"
        )
    )


    if (!file.exists(pb_file)) {

        stop(
            model,
            ": missing pseudobulk matrix:\n",
            pb_file
        )
    }


    if (!file.exists(meta_file)) {

        stop(
            model,
            ": missing pseudobulk metadata:\n",
            meta_file
        )
    }


    cat(
        "\nCounts:\n",
        pb_file,
        "\n",
        sep = ""
    )


    cat(
        "\nMetadata:\n",
        meta_file,
        "\n",
        sep = ""
    )


    # ============================================================================================
    # LOAD
    # ============================================================================================

    pb <- read.csv(
        pb_file,
        row.names = 1,
        check.names = FALSE
    )


    meta <- read.csv(
        meta_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )


    cat(
        "\nPB matrix:",
        nrow(pb),
        "samples x",
        ncol(pb),
        "genes\n"
    )


    cat(
        "Metadata:",
        nrow(meta),
        "rows\n"
    )


    # ============================================================================================
    # ALIGN MATRIX AND METADATA
    # ============================================================================================

    if (!"pseudobulk_id" %in% colnames(meta)) {

        stop(
            model,
            ": metadata missing pseudobulk_id"
        )
    }


    if (anyDuplicated(meta$pseudobulk_id)) {

        stop(
            model,
            ": duplicated pseudobulk_id in metadata"
        )
    }


    if (anyDuplicated(rownames(pb))) {

        stop(
            model,
            ": duplicated pseudobulk IDs in matrix"
        )
    }


    if (!setequal(
        rownames(pb),
        meta$pseudobulk_id
    )) {

        stop(
            model,
            ": pseudobulk IDs differ between matrix and metadata"
        )
    }


    meta <- meta[
        match(
            rownames(pb),
            meta$pseudobulk_id
        ),
        ,
        drop = FALSE
    ]


    stopifnot(
        identical(
            rownames(pb),
            meta$pseudobulk_id
        )
    )


    # ============================================================================================
    # COHORT VALIDATION
    # ============================================================================================

    required_meta <- c(
        "BrNum",
        "Dx",
        "Age",
        "Sex",
        "slide_id",
        "run_date",
        "annots",
        "n_cells"
    )


    missing_meta <- setdiff(
        required_meta,
        colnames(meta)
    )


    if (length(missing_meta) > 0) {

        stop(
            model,
            ": missing metadata columns: ",
            paste(
                missing_meta,
                collapse = ", "
            )
        )
    }


    donor_meta <- get_unique_donor_table(
        meta
    )


    n_donors <- length(
        unique(meta$BrNum)
    )


    n_ntc <- sum(
        donor_meta$Dx == "NTC"
    )


    n_scz <- sum(
        donor_meta$Dx == "SCZ"
    )


    n_celltypes <- length(
        unique(meta$annots)
    )


    cat(
        "\nDonors :",
        n_donors,
        "\n"
    )

    cat(
        "NTC    :",
        n_ntc,
        "\n"
    )

    cat(
        "SCZ    :",
        n_scz,
        "\n"
    )

    cat(
        "Cell types:",
        n_celltypes,
        "\n"
    )

    cat(
        "Pseudobulks:",
        nrow(meta),
        "\n"
    )

    cat(
        "Slides:",
        length(
            unique(meta$slide_id)
        ),
        "\n"
    )


    if (n_donors != EXPECTED_DONORS) {

        stop(
            model,
            ": expected 24 donors"
        )
    }


    if (n_ntc != EXPECTED_NTC) {

        stop(
            model,
            ": expected 12 NTC donors"
        )
    }


    if (n_scz != EXPECTED_SCZ) {

        stop(
            model,
            ": expected 12 SCZ donors"
        )
    }


    if (n_celltypes != EXPECTED_CELLTYPES) {

        stop(
            model,
            ": expected 12 cell types"
        )
    }


    if (nrow(meta) != EXPECTED_PB) {

        stop(
            model,
            ": expected 288 donor-celltype pseudobulks"
        )
    }


    if (ncol(pb) != EXPECTED_GENES) {

        stop(
            model,
            ": expected 300 genes"
        )
    }


    # ============================================================================================
    # PREPARE COUNTS
    #
    # Python output:
    # rows    = pseudobulks
    # columns = genes
    #
    # edgeR/spatialLIBD:
    # rows    = genes
    # columns = pseudobulks
    # ============================================================================================

    counts <- t(
        as.matrix(pb)
    )


    storage.mode(
        counts
    ) <- "double"


    if (any(!is.finite(counts))) {

        stop(
            model,
            ": non-finite pseudobulk values detected"
        )
    }


    if (any(counts < 0)) {

        stop(
            model,
            ": negative pseudobulk values detected"
        )
    }


    # ============================================================================================
    # PAPER/GITHUB EXPRESSION FILTERING
    #
    # registration_pseudobulk() uses filterByExpr with registration_variable.
    #
    # Here:
    # registration_variable = annots
    # ============================================================================================

    section(
        paste0(
            model,
            " FILTER + TMM + LOG2 CPM"
        )
    )


    dge <- edgeR::DGEList(
        counts = counts
    )


    # ============================================================================================
    # FIXED 300-GENE BENCHMARK — NO EXPRESSION FILTERING
    #
    # For cross-model comparison, retain the complete Xenium 300-gene universe.
    # The original filterByExpr() analysis is preserved separately in:
    #   05_celltype_adjusted_de_nofilter_all_models.R
    # ============================================================================================

    keep <- rep(
        TRUE,
        nrow(dge)
    )

    names(keep) <- rownames(dge)


    genes_before <- nrow(
        dge
    )


    genes_after <- sum(
        keep
    )


    dropped_genes <- rownames(
        dge
    )[
        !keep
    ]


    cat(
        "Genes before filter:",
        genes_before,
        "\n"
    )


    cat(
        "Genes retained     :",
        genes_after,
        "\n"
    )


    cat(
        "Genes dropped      :",
        length(
            dropped_genes
        ),
        "\n"
    )


    if (length(dropped_genes) > 0) {

        cat(
            "\nDropped genes:\n"
        )

        print(
            dropped_genes
        )
    }


    dge <- dge[
        keep,
        ,
        keep.lib.sizes = FALSE
    ]


    # ============================================================================================
    # TMM NORMALIZATION
    # ============================================================================================

    dge <- edgeR::calcNormFactors(
        dge,
        method = "TMM"
    )


    # ============================================================================================
    # LOG2 CPM
    #
    # Matches registration_pseudobulk:
    # edgeR::cpm(..., log=TRUE, prior.count=1)
    # ============================================================================================

    logcpm <- edgeR::cpm(
        dge,
        log = TRUE,
        prior.count = 1
    )


    # ============================================================================================
    # METADATA TYPES
    # ============================================================================================

    meta$Dx <- factor(
        meta$Dx,
        levels = c(
            "NTC",
            "SCZ"
        )
    )


    meta$annots <- factor(
        meta$annots
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


    # This reproduces the internal field created by
    # spatialLIBD::registration_pseudobulk().
    #
    # We also explicitly pass BrNum below so donor blocking
    # is unambiguous.

    meta$registration_sample_id <- meta$BrNum


    # ============================================================================================
    # CONSTRUCT SCE
    # ============================================================================================

    sce <- SingleCellExperiment(
        assays = list(
            counts =
                dge$counts,

            logcounts =
                logcpm
        ),

        colData =
            S4Vectors::DataFrame(
                meta
            )
    )


    rowData(
        sce
    )$ID <- rownames(
        sce
    )


    rowData(
        sce
    )$Symbol <- rownames(
        sce
    )


    # ============================================================================================
    # GITHUB CELL-TYPE ADJUSTED DE MODEL
    #
    # covars = c("annots", "Age", "Sex")  # NO-SLIDE sensitivity analysis
    #
    # run_date intentionally NOT included.
    # ============================================================================================

    section(
        paste0(
            model,
            " CELL-TYPE-ADJUSTED SCZ vs NTC DE"
        )
    )


    covars <- c(
        "annots",
        "Age",
        "Sex"
    )


    dx_mod <- spatialLIBD::registration_model(
        sce,
        covars = covars,
        var_registration = "Dx"
    )


    cat(
        "\nDesign matrix:\n"
    )


    print(
        colnames(
            dx_mod
        )
    )


    if (
        qr(
            dx_mod
        )$rank != ncol(
            dx_mod
        )
    ) {

        stop(
            model,
            ": DE design matrix is not full rank."
        )
    }


    # ============================================================================================
    # DONOR CORRELATION
    # ============================================================================================

    dx_block_cor <- spatialLIBD::registration_block_cor(
        sce,
        registration_model = dx_mod,
        var_sample_id = "BrNum"
    )


    cat(
        "\nBlock correlation:",
        dx_block_cor,
        "\n"
    )


    # ============================================================================================
    # DIFFERENTIAL EXPRESSION
    #
    # IMPORTANT:
    # Explicitly use BrNum here.
    #
    # The GitHub object created by registration_pseudobulk()
    # already contains registration_sample_id.
    #
    # Our SCE was reconstructed manually from CSV.
    # ============================================================================================

    dx_res <- spatialLIBD::registration_stats_enrichment(
        sce,
        block_cor = dx_block_cor,
        covars = covars,
        var_registration = "Dx",
        var_sample_id = "BrNum",
        gene_ensembl = "ID",
        gene_name = "Symbol"
    )


    cat(
        "\nDE result columns:\n"
    )


    print(
        colnames(
            dx_res
        )
    )


    # ============================================================================================
    # DE SUMMARY
    # ============================================================================================

    n_nominal <- sum(
        dx_res$p_value_SCZ < 0.05,
        na.rm = TRUE
    )


    n_fdr010 <- sum(
        dx_res$fdr_SCZ < 0.10,
        na.rm = TRUE
    )


    n_fdr005 <- sum(
        dx_res$fdr_SCZ < 0.05,
        na.rm = TRUE
    )


    n_up010 <- sum(
        dx_res$fdr_SCZ < 0.10 &
        dx_res$logFC_SCZ > 0,
        na.rm = TRUE
    )


    n_down010 <- sum(
        dx_res$fdr_SCZ < 0.10 &
        dx_res$logFC_SCZ < 0,
        na.rm = TRUE
    )


    cat(
        "\nGenes tested      :",
        nrow(
            dx_res
        ),
        "\n"
    )


    cat(
        "Nominal p < 0.05  :",
        n_nominal,
        "\n"
    )


    cat(
        "FDR < 0.10        :",
        n_fdr010,
        "\n"
    )


    cat(
        "FDR < 0.05        :",
        n_fdr005,
        "\n"
    )


    cat(
        "Up FDR < 0.10     :",
        n_up010,
        "\n"
    )


    cat(
        "Down FDR < 0.10   :",
        n_down010,
        "\n"
    )


    # ============================================================================================
    # SAVE MODEL OUTPUT
    # ============================================================================================

    model_out <- file.path(
        OUT_ROOT,
        model
    )


    dir.create(
        model_out,
        recursive = TRUE,
        showWarnings = FALSE
    )


    de_file <- file.path(
        model_out,
        paste0(
            model,
            "_celltype_adjusted_SCZ_DE_dx_results.csv"
        )
    )


    filter_file <- file.path(
        model_out,
        paste0(
            model,
            "_celltype_adjusted_filter_summary.csv"
        )
    )


    dropped_file <- file.path(
        model_out,
        paste0(
            model,
            "_celltype_adjusted_dropped_genes.csv"
        )
    )


    write.csv(
        dx_res,
        de_file,
        row.names = FALSE
    )


    write.csv(
        data.frame(
            model =
                model,

            genes_before_filter =
                genes_before,

            genes_after_filter =
                genes_after,

            genes_dropped =
                length(
                    dropped_genes
                ),

            block_correlation =
                dx_block_cor,

            nominal_p05 =
                n_nominal,

            FDR_lt_010 =
                n_fdr010,

            FDR_lt_005 =
                n_fdr005,

            up_FDR010 =
                n_up010,

            down_FDR010 =
                n_down010
        ),
        filter_file,
        row.names = FALSE
    )


    write.csv(
        data.frame(
            gene =
                dropped_genes
        ),
        dropped_file,
        row.names = FALSE
    )


    # ============================================================================================
    # STORE COMBINED SUMMARY
    # ============================================================================================

    summary_list[[model]] <- data.frame(

        model =
            model,

        input_donors =
            n_donors,

        NTC_donors =
            n_ntc,

        SCZ_donors =
            n_scz,

        cell_types =
            n_celltypes,

        pseudobulks =
            nrow(
                meta
            ),

        input_genes =
            genes_before,

        genes_tested =
            genes_after,

        genes_dropped =
            genes_before -
            genes_after,

        block_correlation =
            dx_block_cor,

        nominal_p05 =
            n_nominal,

        FDR_lt_010 =
            n_fdr010,

        FDR_lt_005 =
            n_fdr005,

        up_FDR010 =
            n_up010,

        down_FDR010 =
            n_down010
    )


    cat(
        "\nOutput:\n",
        de_file,
        "\n",
        sep = ""
    )


    cat(
        "\nSTATUS: ",
        model,
        " PASS\n",
        sep = ""
    )
}


# ================================================================================================
# COMBINED SIX-MODEL SUMMARY
# ================================================================================================

section(
    "ALL SIX MODELS — CELL-TYPE-ADJUSTED DE SUMMARY"
)


summary_df <- do.call(
    rbind,
    summary_list
)


rownames(
    summary_df
) <- NULL


print(
    summary_df
)


summary_file <- file.path(
    OUT_ROOT,
    "all_models_celltype_adjusted_DE_summary.csv"
)


write.csv(
    summary_df,
    summary_file,
    row.names = FALSE
)


# ================================================================================================
# FINAL STATUS
# ================================================================================================

section(
    "FINAL STATUS"
)


cat(
    "Models completed :",
    nrow(
        summary_df
    ),
    "/ 6\n"
)


cat(
    "Donors/model     : 24 (12 NTC + 12 SCZ)\n"
)


cat(
    "Cell types/model : 12\n"
)


cat(
    "Pseudobulks      : 288 / model\n"
)


cat(
    "OOF genes        : 300 / model before GitHub-style expression filtering\n"
)


cat(
    "\nCombined summary:\n",
    summary_file,
    "\n",
    sep = ""
)


if (
    nrow(
        summary_df
    ) == 6
) {

    cat(
        "\nFINAL STATUS: PASS — ALL SIX CELL-TYPE-ADJUSTED DE MODELS COMPLETE\n"
    )

} else {

    stop(
        "Not all six models completed."
    )
}
