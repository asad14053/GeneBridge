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
# TASK 2.3
#
# Can NTC-trained imputation models recover SCZ-associated dysregulation?
#
# INPUT
# -----
# Six models:
#   VISTA, gimVI, ENVI, TransImp, Tangram, SpaGE
#
# Each model:
#   24 donors = 12 NTC + 12 SCZ
#   7 spatial domains / donor
#   168 donor × domain pseudobulks
#   300 OOF-imputed genes
#
# MODEL
# -----
# SCZ vs NTC
# Adjust for:
#   predictions_smooth (SpD)
#   Age
#   Sex
#   slide_id
#
# Repeated donor structure:
#   block = BrNum
#
# Paper-compatible workflow:
#   filterByExpr
#   TMM
#   log2 CPM
#   spatialLIBD::registration_model
#   spatialLIBD::registration_block_cor
#   spatialLIBD::registration_stats_enrichment
#
# OUTPUT
# ------
# 1. Full 300-gene DE result per model
# 2. Reported-19-gene recovery table per model
# 3. Six-model recovery summary
# ================================================================================================


# ================================================================================================
# CONFIGURATION
# ================================================================================================

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PB_ROOT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_ntc_dysregulation/pseudobulk"
)

OUT_ROOT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/task2_3_ntc_dysregulation/layer_adjusted_de_nofilter_noslide"
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

PAPER_DE <- file.path(
    ROOT,
    "data/reference/paper_xenium",
    "paper_xenium_layer_adjusted_24donor.csv"
)

FDR_THRESHOLD <- 0.10


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


check_file <- function(path) {

    if (!file.exists(path)) {
        stop(
            "Missing required file:\n",
            path
        )
    }
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

    m
}


pick_col <- function(
    x,
    candidates,
    label
) {

    found <- candidates[
        candidates %in% colnames(x)
    ]

    if (length(found) == 0) {

        stop(
            "Could not identify ",
            label,
            " column.\nAvailable columns:\n",
            paste(colnames(x), collapse = ", ")
        )
    }

    found[1]
}


safe_cor <- function(
    x,
    y,
    method = "pearson"
) {

    ok <- (
        is.finite(x) &
        is.finite(y)
    )

    x <- x[ok]
    y <- y[ok]

    if (
        length(x) < 3 ||
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


# ================================================================================================
# PAPER REFERENCE
# ================================================================================================

section("LOAD OFFICIAL N24 PAPER DE REFERENCE")


if (!file.exists(PAPER_DE)) {

    candidates <- list.files(
        file.path(
            ROOT,
            "data/reference/paper_xenium"
        ),
        pattern = "layer.*adjusted.*\\.csv$",
        full.names = TRUE,
        ignore.case = TRUE
    )

    cat(
        "Expected file not found:\n",
        PAPER_DE,
        "\n\n",
        sep = ""
    )

    cat(
        "Candidate files:\n"
    )

    print(candidates)

    stop(
        "Please verify PAPER_DE path."
    )
}


paper <- read.csv(
    PAPER_DE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


paper_gene_col <- pick_col(
    paper,
    c(
        "gene",
        "Symbol",
        "symbol",
        "gene_name",
        "gene_symbol"
    ),
    "paper gene"
)


paper_fdr_col <- pick_col(
    paper,
    c(
        "fdr_SCZ",
        "FDR_SCZ",
        "fdr",
        "FDR",
        "adj.P.Val"
    ),
    "paper SCZ FDR"
)


paper_logfc_col <- pick_col(
    paper,
    c(
        "logFC_SCZ",
        "LogFC_SCZ",
        "logFC",
        "log2FC",
        "estimate"
    ),
    "paper SCZ logFC"
)


paper[[paper_gene_col]] <- as.character(
    paper[[paper_gene_col]]
)


paper19 <- paper[
    !is.na(paper[[paper_fdr_col]]) &
    paper[[paper_fdr_col]] < FDR_THRESHOLD,
    ,
    drop = FALSE
]


cat(
    "Paper file     :",
    PAPER_DE,
    "\n"
)

cat(
    "Gene column    :",
    paper_gene_col,
    "\n"
)

cat(
    "FDR column     :",
    paper_fdr_col,
    "\n"
)

cat(
    "logFC column   :",
    paper_logfc_col,
    "\n"
)

cat(
    "FDR < 0.10 genes:",
    nrow(paper19),
    "\n"
)


if (nrow(paper19) != 19) {

    stop(
        "Expected exactly 19 reported genes at FDR < 0.10, found ",
        nrow(paper19),
        ". Verify the official N24 paper DE file."
    )
}


paper19_std <- data.frame(

    gene =
        as.character(
            paper19[[paper_gene_col]]
        ),

    paper_logFC_SCZ =
        as.numeric(
            paper19[[paper_logfc_col]]
        ),

    paper_fdr_SCZ =
        as.numeric(
            paper19[[paper_fdr_col]]
        ),

    stringsAsFactors = FALSE
)


if (anyDuplicated(paper19_std$gene)) {

    stop(
        "Duplicated gene symbols in reported 19."
    )
}


cat(
    "\nReported 19 genes:\n"
)

print(
    paper19_std,
    row.names = FALSE
)


# ================================================================================================
# STORAGE
# ================================================================================================

all_summary <- list()
all_recovery <- list()


# ================================================================================================
# PROCESS EACH MODEL
# ================================================================================================

for (MODEL in MODELS) {

    section(
        paste(
            "MODEL:",
            MODEL
        )
    )


    MODEL_OUT <- file.path(
        OUT_ROOT,
        MODEL
    )

    dir.create(
        MODEL_OUT,
        recursive = TRUE,
        showWarnings = FALSE
    )


    COUNTS_FILE <- file.path(
        PB_ROOT,
        MODEL,
        paste0(
            MODEL,
            "_NTCtrained_OOF_donor_domain_pseudobulk_300genes.csv.gz"
        )
    )


    META_FILE <- file.path(
        PB_ROOT,
        MODEL,
        paste0(
            MODEL,
            "_NTCtrained_OOF_donor_domain_pseudobulk_metadata.csv"
        )
    )


    check_file(
        COUNTS_FILE
    )

    check_file(
        META_FILE
    )


    # ============================================================================================
    # LOAD
    # ============================================================================================

    cat(
        "\nCounts:\n",
        COUNTS_FILE,
        "\n",
        sep = ""
    )

    cat(
        "\nMetadata:\n",
        META_FILE,
        "\n",
        sep = ""
    )


    pb <- read_pb_matrix(
        COUNTS_FILE
    )


    meta <- read.csv(
        META_FILE,
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
    # STRUCTURAL VALIDATION
    # ============================================================================================

    required_meta <- c(
        "pseudobulk_id",
        "BrNum",
        "Dx",
        "Age",
        "Sex",
        "slide_id",
        "predictions_smooth",
        "domain_annotations",
        "n_cells"
    )


    missing_meta <- setdiff(
        required_meta,
        colnames(meta)
    )


    if (length(missing_meta) > 0) {

        stop(
            MODEL,
            ": missing metadata columns: ",
            paste(
                missing_meta,
                collapse = ", "
            )
        )
    }


    if (
        nrow(pb) != 168 ||
        nrow(meta) != 168
    ) {

        stop(
            MODEL,
            ": expected 168 pseudobulks."
        )
    }


    if (ncol(pb) != 300) {

        stop(
            MODEL,
            ": expected 300 genes, found ",
            ncol(pb)
        )
    }


    if (anyDuplicated(colnames(pb))) {

        stop(
            MODEL,
            ": duplicated genes."
        )
    }


    # Align matrix rows to metadata.

    if (
        !is.null(rownames(pb)) &&
        all(
            meta$pseudobulk_id
            %in%
            rownames(pb)
        )
    ) {

        pb <- pb[
            match(
                meta$pseudobulk_id,
                rownames(pb)
            ),
            ,
            drop = FALSE
        ]

    } else {

        if (nrow(pb) != nrow(meta)) {

            stop(
                MODEL,
                ": matrix / metadata mismatch."
            )
        }

        rownames(pb) <- meta$pseudobulk_id
    }


    if (
        any(!is.finite(pb)) ||
        any(pb < 0)
    ) {

        stop(
            MODEL,
            ": invalid count-scale values."
        )
    }


    meta$BrNum <- as.character(
        meta$BrNum
    )

    meta$Dx <- factor(
        meta$Dx,
        levels = c(
            "NTC",
            "SCZ"
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

    meta$predictions_smooth <- factor(
        meta$predictions_smooth,
        levels = c(
            "spd07",
            "spd06",
            "spd02",
            "spd05",
            "spd03",
            "spd01",
            "spd04"
        )
    )

    meta$BrNum <- factor(
        meta$BrNum
    )


    if (anyNA(meta$Dx)) {

        stop(
            MODEL,
            ": invalid Dx labels."
        )
    }


    donor_dx <- unique(
        data.frame(
            BrNum = meta$BrNum,
            Dx = meta$Dx
        )
    )


    dx_counts <- table(
        donor_dx$Dx
    )


    if (
        as.integer(dx_counts["NTC"]) != 12 ||
        as.integer(dx_counts["SCZ"]) != 12
    ) {

        stop(
            MODEL,
            ": expected 12 NTC + 12 SCZ donors."
        )
    }


    if (
        length(unique(meta$BrNum)) != 24 ||
        length(unique(meta$predictions_smooth)) != 7
    ) {

        stop(
            MODEL,
            ": expected 24 donors and 7 SpDs."
        )
    }


    if (any(meta$n_cells < 10)) {

        stop(
            MODEL,
            ": pseudobulk with fewer than 10 cells detected."
        )
    }


    cat(
        "Donors :",
        length(unique(meta$BrNum)),
        "\n"
    )

    cat(
        "NTC    :",
        dx_counts["NTC"],
        "\n"
    )

    cat(
        "SCZ    :",
        dx_counts["SCZ"],
        "\n"
    )

    cat(
        "SpDs   :",
        length(unique(meta$predictions_smooth)),
        "\n"
    )


    # ============================================================================================
    # PAPER-COMPATIBLE NORMALIZATION
    #
    # registration_pseudobulk() normally performs:
    #   filterByExpr
    #   calcNormFactors (TMM)
    #   edgeR::cpm(log=TRUE, prior.count=1)
    #
    # We are already pseudobulked, so reproduce those steps directly.
    # ============================================================================================

    section(
        paste(
            MODEL,
            "NO FILTER + TMM + LOG2 CPM"
        )
    )


    counts_gene_by_pb <- t(
        pb
    )


    dge <- edgeR::DGEList(
        counts = counts_gene_by_pb
    )


    # ========================================================================================
    # NO EXPRESSION FILTERING
    #
    # Task 2.3 sensitivity analysis:
    # retain all 300 OOF-imputed genes for DE.
    # ========================================================================================

    cat(
        "Genes before DE    :",
        nrow(dge),
        "\n"
    )

    cat(
        "Genes retained     :",
        nrow(dge),
        "/ 300 (filterByExpr DISABLED)\n"
    )


    dge <- edgeR::calcNormFactors(
        dge,
        method = "TMM"
    )


    logcpm <- edgeR::cpm(
        dge,
        log = TRUE,
        prior.count = 1
    )


    if (any(!is.finite(logcpm))) {

        stop(
            MODEL,
            ": non-finite logCPM."
        )
    }


    # ============================================================================================
    # BUILD PSEUDOBULK SCE
    # ============================================================================================

    genes_tested <- rownames(
        dge$counts
    )


    rd <- S4Vectors::DataFrame(

        ID =
            genes_tested,

        Symbol =
            genes_tested,

        row.names =
            genes_tested
    )


    cd <- S4Vectors::DataFrame(
        meta,
        row.names = meta$pseudobulk_id
    )


    spe_pseudo <- SingleCellExperiment::SingleCellExperiment(

        assays = list(

            counts =
                dge$counts,

            logcounts =
                logcpm
        ),

        rowData =
            rd,

        colData =
            cd
    )


    # ============================================================================================
    # EXACT PAPER-STYLE LAYER-ADJUSTED MODEL
    # ============================================================================================

    section(
        paste(
            MODEL,
            "LAYER-ADJUSTED SCZ vs NTC DE"
        )
    )


    covars <- c(
        "predictions_smooth",
        "Age",
        "Sex"
    )


    dx_mod <- spatialLIBD::registration_model(

        spe_pseudo,

        covars =
            covars,

        var_registration =
            "Dx"
    )


    cat(
        "\nDesign matrix:\n"
    )

    print(
        colnames(dx_mod)
    )


    if (
        qr(dx_mod)$rank
        !=
        ncol(dx_mod)
    ) {

        stop(
            MODEL,
            ": DE design matrix is not full rank."
        )
    }


    dx_block_cor <- spatialLIBD::registration_block_cor(

        spe_pseudo,

        registration_model =
            dx_mod,

        var_sample_id =
            "BrNum"
    )


    cat(
        "\nBlock correlation:",
        dx_block_cor,
        "\n"
    )


    dx_res <- spatialLIBD::registration_stats_enrichment(

        spe_pseudo,

        block_cor =
            dx_block_cor,

        covars =
            covars,

        var_registration =
            "Dx",

        var_sample_id =
            "BrNum",

        gene_ensembl =
            "ID",

        gene_name =
            "Symbol"
    )


    cat(
        "\nDE result columns:\n"
    )

    print(
        colnames(dx_res)
    )


    # ============================================================================================
    # IDENTIFY SCZ OUTPUT COLUMNS
    # ============================================================================================

    result_gene_col <- pick_col(
        dx_res,
        c(
            "gene",
            "Symbol",
            "symbol"
        ),
        paste(
            MODEL,
            "DE gene"
        )
    )


    result_logfc_col <- pick_col(
        dx_res,
        c(
            "logFC_SCZ",
            "logFC.SCZ",
            "SCZ_logFC"
        ),
        paste(
            MODEL,
            "DE SCZ logFC"
        )
    )


    result_p_col <- pick_col(
        dx_res,
        c(
            "p_value_SCZ",
            "p_SCZ",
            "p.SCZ",
            "PValue_SCZ",
            "pvalue_SCZ"
        ),
        paste(
            MODEL,
            "DE SCZ p"
        )
    )


    result_fdr_col <- pick_col(
        dx_res,
        c(
            "fdr_SCZ",
            "fdr.SCZ",
            "FDR_SCZ"
        ),
        paste(
            MODEL,
            "DE SCZ FDR"
        )
    )


    result_t_col_candidates <- c(
        "t_stat_SCZ",
        "t_SCZ",
        "t.SCZ",
        "SCZ_t"
    )


    result_t_found <- result_t_col_candidates[
        result_t_col_candidates
        %in%
        colnames(dx_res)
    ]


    if (length(result_t_found) > 0) {

        result_t_col <- result_t_found[1]

    } else {

        result_t_col <- NA_character_
    }


    # ============================================================================================
    # SAVE FULL DE RESULTS
    # ============================================================================================

    dx_res$model <- MODEL

    dx_res$block_correlation <- dx_block_cor


    FULL_DE_FILE <- file.path(
        MODEL_OUT,
        paste0(
            MODEL,
            "_layer_adjusted_SCZ_DE_dx_results.csv"
        )
    )


    write.csv(
        dx_res,
        FULL_DE_FILE,
        row.names = FALSE
    )


    # ============================================================================================
    # STANDARDIZED RESULT
    # ============================================================================================

    model_de <- data.frame(

        gene =
            as.character(
                dx_res[[result_gene_col]]
            ),

        model_logFC_SCZ =
            as.numeric(
                dx_res[[result_logfc_col]]
            ),

        model_p_SCZ =
            as.numeric(
                dx_res[[result_p_col]]
            ),

        model_fdr_SCZ =
            as.numeric(
                dx_res[[result_fdr_col]]
            ),

        stringsAsFactors = FALSE
    )


    if (!is.na(result_t_col)) {

        model_de$model_t_SCZ <- as.numeric(
            dx_res[[result_t_col]]
        )

    } else {

        model_de$model_t_SCZ <- NA_real_
    }


    # NO-FILTER analysis must retain all 300 genes.
    if (nrow(model_de) != 300) {

        stop(
            MODEL,
            ": NO-FILTER DE expected 300 genes but returned ",
            nrow(model_de)
        )
    }



    # ============================================================================================
    # REPORTED 19-GENE RECOVERY
    # ============================================================================================

    section(
        paste(
            MODEL,
            "REPORTED 19-GENE RECOVERY"
        )
    )


    missing19 <- setdiff(
        paper19_std$gene,
        model_de$gene
    )


    if (length(missing19) > 0) {

        cat(
            "\nNOTE:",
            MODEL,
            "reported gene(s) removed by filterByExpr():",
            paste(
                missing19,
                collapse = ", "
            ),
            "\n"
        )

        cat(
            "These genes will be counted as NOT RECOVERED, not treated as an error.\n"
        )
    }


    recovery <- merge(
        paper19_std,
        model_de,
        by = "gene",
        all.x = TRUE,
        sort = FALSE
    )


    recovery <- recovery[
        match(
            paper19_std$gene,
            recovery$gene
        ),
        ,
        drop = FALSE
    ]


    recovery$model <- MODEL


    recovery$recovered_FDR010 <- (
        !is.na(recovery$model_fdr_SCZ) &
        recovery$model_fdr_SCZ < FDR_THRESHOLD
    )


    recovery$nominal_p05 <- (
        !is.na(recovery$model_p_SCZ) &
        recovery$model_p_SCZ < 0.05
    )


    recovery$direction_match <- (
        is.finite(recovery$paper_logFC_SCZ) &
        is.finite(recovery$model_logFC_SCZ) &
        sign(recovery$paper_logFC_SCZ)
        ==
        sign(recovery$model_logFC_SCZ)
    )


    recovery$recovered_same_direction <- (
        recovery$recovered_FDR010 &
        recovery$direction_match
    )


    RECOVERY_FILE <- file.path(
        MODEL_OUT,
        paste0(
            MODEL,
            "_reported19_gene_recovery.csv"
        )
    )


    write.csv(
        recovery,
        RECOVERY_FILE,
        row.names = FALSE
    )


    # ============================================================================================
    # MODEL SUMMARY
    # ============================================================================================

    pvals <- model_de$model_p_SCZ
    fdrs <- model_de$model_fdr_SCZ
    logfc <- model_de$model_logFC_SCZ


    model_summary <- data.frame(

        model =
            MODEL,

        input_donors =
            24,

        NTC_donors =
            12,

        SCZ_donors =
            12,

        spatial_domains =
            7,

        pseudobulks =
            168,

        input_genes =
            300,

        genes_tested =
            nrow(model_de),

        block_correlation =
            dx_block_cor,

        nominal_p05 =
            sum(
                pvals < 0.05,
                na.rm = TRUE
            ),

        FDR_lt_010 =
            sum(
                fdrs < 0.10,
                na.rm = TRUE
            ),

        FDR_lt_005 =
            sum(
                fdrs < 0.05,
                na.rm = TRUE
            ),

        up_FDR010 =
            sum(
                fdrs < 0.10 &
                logfc > 0,
                na.rm = TRUE
            ),

        down_FDR010 =
            sum(
                fdrs < 0.10 &
                logfc < 0,
                na.rm = TRUE
            ),

        reported_genes =
            19,

        reported19_nominal_p05 =
            sum(
                recovery$nominal_p05,
                na.rm = TRUE
            ),

        reported19_recovered_FDR010 =
            sum(
                recovery$recovered_FDR010,
                na.rm = TRUE
            ),

        reported19_direction_match =
            sum(
                recovery$direction_match,
                na.rm = TRUE
            ),

        reported19_recovered_same_direction =
            sum(
                recovery$recovered_same_direction,
                na.rm = TRUE
            ),

        reported19_recovery_rate =
            mean(
                recovery$recovered_FDR010,
                na.rm = TRUE
            ),

        paper_vs_model_logFC_Pearson_19 =
            safe_cor(
                recovery$paper_logFC_SCZ,
                recovery$model_logFC_SCZ,
                "pearson"
            ),

        paper_vs_model_logFC_Spearman_19 =
            safe_cor(
                recovery$paper_logFC_SCZ,
                recovery$model_logFC_SCZ,
                "spearman"
            ),

        stringsAsFactors = FALSE
    )


    SUMMARY_FILE <- file.path(
        MODEL_OUT,
        paste0(
            MODEL,
            "_layer_adjusted_SCZ_DE_summary.csv"
        )
    )


    write.csv(
        model_summary,
        SUMMARY_FILE,
        row.names = FALSE
    )


    cat(
        "\n",
        MODEL,
        " summary:\n",
        sep = ""
    )

    print(
        model_summary,
        row.names = FALSE
    )


    cat(
        "\nReported 19 recovery:\n"
    )

    print(
        recovery[
            ,
            c(
                "gene",
                "paper_logFC_SCZ",
                "paper_fdr_SCZ",
                "model_logFC_SCZ",
                "model_p_SCZ",
                "model_fdr_SCZ",
                "recovered_FDR010",
                "direction_match"
            )
        ],
        row.names = FALSE
    )


    all_summary[[MODEL]] <- model_summary
    all_recovery[[MODEL]] <- recovery


    cat(
        "\nSTATUS:",
        MODEL,
        "PASS\n"
    )
}


# ================================================================================================
# SIX-MODEL COMBINED SUMMARY
# ================================================================================================

section(
    "ALL SIX MODELS — TASK 2.3 SUMMARY"
)


summary_all <- do.call(
    rbind,
    all_summary
)

rownames(
    summary_all
) <- NULL


recovery_all <- do.call(
    rbind,
    all_recovery
)

rownames(
    recovery_all
) <- NULL


SUMMARY_ALL_FILE <- file.path(
    OUT_ROOT,
    "all_models_reported19_recovery_summary.csv"
)


RECOVERY_ALL_FILE <- file.path(
    OUT_ROOT,
    "all_models_reported19_gene_recovery.csv"
)


write.csv(
    summary_all,
    SUMMARY_ALL_FILE,
    row.names = FALSE
)


write.csv(
    recovery_all,
    RECOVERY_ALL_FILE,
    row.names = FALSE
)


print(
    summary_all[
        ,
        c(
            "model",
            "genes_tested",
            "block_correlation",
            "FDR_lt_010",
            "reported19_nominal_p05",
            "reported19_recovered_FDR010",
            "reported19_direction_match",
            "reported19_recovered_same_direction",
            "reported19_recovery_rate",
            "paper_vs_model_logFC_Pearson_19",
            "paper_vs_model_logFC_Spearman_19"
        )
    ],
    row.names = FALSE
)


section(
    "FINAL STATUS"
)


cat(
    "Models completed: 6 / 6\n"
)

cat(
    "Donors/model   : 24 (12 NTC + 12 SCZ)\n"
)

cat(
    "Domains/model  : 7\n"
)

cat(
    "Pseudobulks    : 168 / model\n"
)

cat(
    "OOF genes      : 300 / model; filterByExpr disabled\n"
)

cat(
    "Reported genes : 19\n"
)

cat(
    "\nCombined summary:\n",
    SUMMARY_ALL_FILE,
    "\n",
    sep = ""
)

cat(
    "\nCombined 19-gene table:\n",
    RECOVERY_ALL_FILE,
    "\n",
    sep = ""
)

cat(
    "\nFINAL STATUS: PASS — ALL SIX NTC-TRAINED OOF MODELS COMPLETE\n"
)
