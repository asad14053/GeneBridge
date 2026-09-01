#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SummarizedExperiment)
})

# ================================================================================================
# 01g_plot_pseudobulk_concordance_qc.R
#
# PURPOSE
# -------
# Visual QC of donor × spatial-domain pseudobulk concordance:
#
#   A. Published paper Xenium vs GeneBridge original Xenium
#      48,300 gene × pseudobulk values
#
#   B. Paper vs GeneBridge donor-domain TOTAL counts
#      161 donor × domain pseudobulks
#
#   C. GeneBridge original vs ENVI measured-300 preserved values
#      161 donor × domain pseudobulks
#
# IMPORTANT
# ---------
# The current ENVI measured-300 matrix is NOT OOF prediction.
# It contains the 300 measured Xenium genes preserved in the
# full-transcriptome production ENVI object.
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


SIDE_BY_SIDE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "paper_raw_comparison",
    "paper_original_ENVI_OOF300_side_by_side.csv"
)


OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "paper_raw_comparison",
    "plots"
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
        paste(rep("=", 100), collapse = ""),
        "\n",
        x,
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

    suppressWarnings(
        cor(
            x,
            y,
            method = method
        )
    )
}


# ================================================================================================
# CHECK INPUTS
# ================================================================================================

section("CHECK INPUTS")


for (
    f in c(
        PAPER_RDS,
        OUR_COUNTS,
        OUR_META,
        SIDE_BY_SIDE
    )
) {

    if (!file.exists(f)) {

        stop(
            "Missing input:\n",
            f
        )
    }

    cat(
        "FOUND:",
        f,
        "\n"
    )
}


# ================================================================================================
# LOAD PAPER
# ================================================================================================

section("LOAD PAPER")


paper <- readRDS(
    PAPER_RDS
)


paper_counts <- assay(
    paper,
    "counts"
)


paper_meta <- as.data.frame(
    colData(paper)
)


paper_rowdata <- as.data.frame(
    rowData(paper)
)


paper_genes <- as.character(
    paper_rowdata$Symbol
)


paper_matrix_all <- t(
    as.matrix(
        paper_counts
    )
)


storage.mode(
    paper_matrix_all
) <- "double"


colnames(
    paper_matrix_all
) <- paper_genes


cat(
    "Paper:",
    nrow(paper_matrix_all),
    "pseudobulks x",
    ncol(paper_matrix_all),
    "genes\n"
)


# ================================================================================================
# LOAD GENEBRIDGE ORIGINAL
# ================================================================================================

section("LOAD GENEBRIDGE ORIGINAL")


our_matrix <- read.csv(
    OUR_COUNTS,
    row.names = 1,
    check.names = FALSE
)


our_matrix <- as.matrix(
    our_matrix
)


storage.mode(
    our_matrix
) <- "double"


our_meta <- read.csv(
    OUR_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat(
    "Original:",
    nrow(our_matrix),
    "pseudobulks x",
    ncol(our_matrix),
    "genes\n"
)


# ================================================================================================
# ALIGN GENES
# ================================================================================================

section("ALIGN GENES")


if (
    !setequal(
        colnames(our_matrix),
        paper_genes
    )
) {

    stop(
        "Paper and GeneBridge gene universes differ."
    )
}


paper_matrix_all <- paper_matrix_all[
    ,
    colnames(our_matrix),
    drop = FALSE
]


cat(
    "Shared genes:",
    ncol(our_matrix),
    "\n"
)


# ================================================================================================
# SUBSET PAPER N24 -> N23
# ================================================================================================

section("SUBSET PAPER N24 -> N23")


our_donors <- unique(
    as.character(
        our_meta$BrNum
    )
)


keep <- (
    as.character(
        paper_meta$BrNum
    )
    %in%
    our_donors
)


paper_meta <- paper_meta[
    keep,
    ,
    drop = FALSE
]


paper_matrix <- paper_matrix_all[
    keep,
    ,
    drop = FALSE
]


paper_key <- paste(
    as.character(
        paper_meta$BrNum
    ),
    as.character(
        paper_meta$predictions_smooth
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


idx <- match(
    our_key,
    paper_key
)


if (any(is.na(idx))) {

    stop(
        "Could not align paper donor-domain pseudobulks."
    )
}


paper_matrix <- paper_matrix[
    idx,
    ,
    drop = FALSE
]


paper_meta <- paper_meta[
    idx,
    ,
    drop = FALSE
]


cat(
    "Aligned pseudobulks:",
    nrow(paper_matrix),
    "\n"
)


# ================================================================================================
# RAW COUNT QC METRICS
# ================================================================================================

section("RAW COUNT QC")


paper_vec <- as.vector(
    paper_matrix
)


our_vec <- as.vector(
    our_matrix
)


exact <- (
    paper_vec
    ==
    our_vec
)


pearson <- safe_cor(
    paper_vec,
    our_vec,
    "pearson"
)


spearman <- safe_cor(
    paper_vec,
    our_vec,
    "spearman"
)


exact_fraction <- mean(
    exact
)


cat(
    "Entries:",
    length(paper_vec),
    "\n"
)


cat(
    "Exact:",
    sum(exact),
    "\n"
)


cat(
    "Exact fraction:",
    sprintf(
        "%.8f",
        exact_fraction
    ),
    "\n"
)


cat(
    "Pearson:",
    sprintf(
        "%.12f",
        pearson
    ),
    "\n"
)


cat(
    "Spearman:",
    sprintf(
        "%.12f",
        spearman
    ),
    "\n"
)


# ================================================================================================
# LOAD SIDE-BY-SIDE DONOR-DOMAIN TOTALS
# ================================================================================================

section("LOAD DONOR-DOMAIN TOTALS")


side <- read.csv(
    SIDE_BY_SIDE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


# Current CSV header may still say ENVI_OOF300.
# Scientifically this is measured/preserved ENVI 300.

if (
    !"ENVI_OOF300"
    %in%
    colnames(side)
) {

    stop(
        "Expected ENVI_OOF300 legacy column not found."
    )
}


side$ENVI_measured300 <- side$ENVI_OOF300


cat(
    "Donor-domain rows:",
    nrow(side),
    "\n"
)


# ================================================================================================
# IDENTIFY MISMATCHED DONOR-DOMAIN PSEUDOBULKS
# ================================================================================================

paper_original_bad <- (
    side$Paper
    !=
    side$Original
)


cat(
    "Paper-vs-original exact pseudobulks:",
    sum(!paper_original_bad),
    "/",
    nrow(side),
    "\n"
)


cat(
    "\nMismatched donor-domain rows:\n"
)


print(
    side[
        paper_original_bad,
        c(
            "Donor",
            "Dx",
            "Domain",
            "Paper",
            "Original",
            "ENVI_measured300"
        )
    ],
    row.names = FALSE
)


# ================================================================================================
# PLOTTING FUNCTION
# ================================================================================================

make_plot <- function() {

    oldpar <- par(
        no.readonly = TRUE
    )

    on.exit(
        par(oldpar)
    )


    par(
        mfrow = c(1, 3),
        mar = c(5.3, 5.3, 4.3, 1.5),
        oma = c(0, 0, 3, 0),
        las = 1
    )


    # ============================================================================================
    # PANEL A
    # ============================================================================================

    x <- log1p(
        paper_vec
    )

    y <- log1p(
        our_vec
    )


    plot(
        x,
        y,
        pch = 16,
        cex = 0.23,
        xlab = "Published Xenium pseudobulk\nlog1p(raw count)",
        ylab = "GeneBridge original pseudobulk\nlog1p(raw count)",
        main = "A  Gene-level concordance"
    )


    abline(
        a = 0,
        b = 1,
        lwd = 2,
        lty = 2
    )


    legend(
        "topleft",
        legend = c(
            paste0(
                "n = ",
                format(
                    length(paper_vec),
                    big.mark = ","
                )
            ),

            paste0(
                "Pearson r = ",
                sprintf(
                    "%.9f",
                    pearson
                )
            ),

            paste0(
                "Spearman rho = ",
                sprintf(
                    "%.9f",
                    spearman
                )
            ),

            paste0(
                "Exact = ",
                sprintf(
                    "%.3f%%",
                    100 * exact_fraction
                )
            )
        ),
        bty = "n",
        cex = 0.78
    )


    # ============================================================================================
    # PANEL B
    # ============================================================================================

    x2 <- side$Paper
    y2 <- side$Original


    plot(
        x2,
        y2,
        pch = 16,
        cex = 0.8,
        xlab = "Published Xenium\npseudobulk total count",
        ylab = "GeneBridge original\npseudobulk total count",
        main = "B  Donor-domain totals"
    )


    abline(
        a = 0,
        b = 1,
        lwd = 2,
        lty = 2
    )


    # Highlight mismatches.

    points(
        x2[
            paper_original_bad
        ],
        y2[
            paper_original_bad
        ],
        pch = 1,
        cex = 2,
        lwd = 2
    )


    if (
        sum(
            paper_original_bad
        ) > 0
    ) {

        labels_bad <- paste0(
            side$Donor[
                paper_original_bad
            ],
            " ",
            side$Domain[
                paper_original_bad
            ]
        )


        text(
            x2[
                paper_original_bad
            ],
            y2[
                paper_original_bad
            ],
            labels = labels_bad,
            pos = 4,
            cex = 0.72
        )
    }


    donor_total_r <- safe_cor(
        x2,
        y2,
        "pearson"
    )


    legend(
        "topleft",
        legend = c(
            "161 donor-domain pseudobulks",

            paste0(
                "Exact = ",
                sum(
                    !paper_original_bad
                ),
                "/161"
            ),

            paste0(
                "Pearson r = ",
                sprintf(
                    "%.9f",
                    donor_total_r
                )
            )
        ),
        bty = "n",
        cex = 0.78
    )


    # ============================================================================================
    # PANEL C
    # ============================================================================================

    x3 <- side$Original
    y3 <- side$ENVI_measured300


    plot(
        x3,
        y3,
        pch = 16,
        cex = 0.8,
        xlab = "GeneBridge original Xenium\npseudobulk total count",
        ylab = "ENVI preserved measured-300\npseudobulk total count",
        main = "C  ENVI measured-gene preservation"
    )


    abline(
        a = 0,
        b = 1,
        lwd = 2,
        lty = 2
    )


    envi_exact <- (
        x3
        ==
        y3
    )


    envi_r <- safe_cor(
        x3,
        y3,
        "pearson"
    )


    legend(
        "topleft",
        legend = c(
            "Not OOF prediction",

            paste0(
                "Exact = ",
                sum(envi_exact),
                "/161"
            ),

            paste0(
                "Pearson r = ",
                sprintf(
                    "%.9f",
                    envi_r
                )
            ),

            paste0(
                "MAE = ",
                sprintf(
                    "%.3f",
                    mean(
                        abs(
                            y3 - x3
                        )
                    )
                )
            )
        ),
        bty = "n",
        cex = 0.78
    )


    mtext(
        "Xenium donor × spatial-domain pseudobulk concordance QC",
        outer = TRUE,
        cex = 1.3,
        font = 2
    )
}


# ================================================================================================
# SAVE PDF
# ================================================================================================

section("SAVE PDF")


PDF_FILE <- file.path(
    OUT,
    "01g_pseudobulk_concordance_qc.pdf"
)


pdf(
    PDF_FILE,
    width = 16,
    height = 5.7
)


make_plot()


dev.off()


cat(
    "Saved:",
    PDF_FILE,
    "\n"
)


# ================================================================================================
# SAVE PNG
# ================================================================================================

section("SAVE PNG")


PNG_FILE <- file.path(
    OUT,
    "01g_pseudobulk_concordance_qc.png"
)


png(
    PNG_FILE,
    width = 4800,
    height = 1700,
    res = 300
)


make_plot()


dev.off()


cat(
    "Saved:",
    PNG_FILE,
    "\n"
)


# ================================================================================================
# SAVE PLOT SUMMARY
# ================================================================================================

summary <- data.frame(

    metric = c(
        "gene_pseudobulk_entries",
        "paper_original_exact_entries",
        "paper_original_exact_fraction",
        "paper_original_pearson",
        "paper_original_spearman",
        "paper_original_exact_pseudobulks",
        "paper_original_total_pseudobulks",
        "envi_measured300_exact_pseudobulks",
        "envi_measured300_total_pseudobulks"
    ),

    value = c(
        length(paper_vec),
        sum(exact),
        exact_fraction,
        pearson,
        spearman,
        sum(!paper_original_bad),
        nrow(side),
        sum(
            side$Original
            ==
            side$ENVI_measured300
        ),
        nrow(side)
    ),

    stringsAsFactors = FALSE
)


SUMMARY_FILE <- file.path(
    OUT,
    "01g_pseudobulk_concordance_qc_summary.csv"
)


write.csv(
    summary,
    SUMMARY_FILE,
    row.names = FALSE
)


section("FINAL SUMMARY")


print(
    summary,
    row.names = FALSE
)


cat(
    "\nFigure interpretation:\n"
)


cat(
    "A = paper vs GeneBridge original at 48,300 gene-by-pseudobulk values.\n"
)


cat(
    "B = paper vs GeneBridge original at 161 donor-domain totals.\n"
)


cat(
    "C = original vs ENVI preserved measured-300 control; this is NOT OOF imputation.\n"
)


cat(
    "\nFINAL STATUS: PLOT COMPLETE\n"
)
