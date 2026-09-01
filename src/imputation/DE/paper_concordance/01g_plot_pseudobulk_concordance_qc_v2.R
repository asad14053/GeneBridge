#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SummarizedExperiment)
})

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


safe_cor <- function(x, y, method = "pearson") {
    suppressWarnings(
        cor(
            x,
            y,
            method = method
        )
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


# ================================================================================================
# LOAD ORIGINAL
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


# ================================================================================================
# ALIGN GENES
# ================================================================================================

if (
    !setequal(
        colnames(our_matrix),
        paper_genes
    )
) {
    stop("Gene sets differ.")
}

paper_matrix_all <- paper_matrix_all[
    ,
    colnames(our_matrix),
    drop = FALSE
]


# ================================================================================================
# PAPER N24 -> N23
# ================================================================================================

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
    paper_meta$BrNum,
    paper_meta$predictions_smooth,
    sep = "::"
)

our_key <- paste(
    our_meta$BrNum,
    our_meta$predictions_smooth,
    sep = "::"
)

idx <- match(
    our_key,
    paper_key
)

if (any(is.na(idx))) {
    stop("Could not align donor-domain keys.")
}

paper_matrix <- paper_matrix[
    idx,
    ,
    drop = FALSE
]


# ================================================================================================
# METRICS
# ================================================================================================

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


# ================================================================================================
# SIDE-BY-SIDE TOTALS
# ================================================================================================

side <- read.csv(
    SIDE_BY_SIDE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (
    !"ENVI_OOF300"
    %in%
    colnames(side)
) {
    stop(
        "Legacy ENVI_OOF300 column not found."
    )
}

# Correct scientific interpretation.
side$ENVI_measured300 <- side$ENVI_OOF300

bad <- (
    side$Paper
    !=
    side$Original
)


# ================================================================================================
# PLOT FUNCTION
# ================================================================================================

make_plot <- function() {

    oldpar <- par(
        no.readonly = TRUE
    )

    on.exit(
        par(oldpar)
    )


    # ----------------------------------------------------------------
    # IMPORTANT FIXES:
    #
    # mar:
    #   much larger LEFT margin
    #
    # mgp:
    #   move axis title farther from tick labels
    #
    # cex.axis:
    #   slightly smaller tick labels
    #
    # las=1:
    #   keep tick numbers horizontal
    # ----------------------------------------------------------------

    par(
        mfrow = c(1, 3),

        mar = c(
            5.5,   # bottom
            7.0,   # LEFT -- increased substantially
            4.5,   # top
            1.5    # right
        ),

        oma = c(
            0,
            0,
            3,
            0
        ),

        mgp = c(
            4.6,   # axis title distance
            1.4,   # tick-label distance
            0
        ),

        las = 1,

        cex.axis = 0.82,

        cex.lab = 0.92,

        cex.main = 1.05
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
        cex = 0.22,

        xlab = "Published Xenium\nlog1p(raw count)",

        ylab = "GeneBridge original\nlog1p(raw count)",

        main = "A. Gene-level concordance"
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
        cex = 0.76
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
        cex = 0.75,

        xlab = "Published Xenium\npseudobulk total",

        ylab = "GeneBridge original\npseudobulk total",

        main = "B. Donor-domain totals"
    )


    abline(
        a = 0,
        b = 1,
        lwd = 2,
        lty = 2
    )


    points(
        x2[bad],
        y2[bad],

        pch = 1,
        cex = 2.0,
        lwd = 2
    )


    if (
        sum(bad) > 0
    ) {

        labels_bad <- paste0(
            side$Donor[bad],
            " ",
            side$Domain[bad]
        )


        text(
            x2[bad],
            y2[bad],

            labels = labels_bad,

            pos = 4,

            offset = 0.5,

            cex = 0.72
        )
    }


    r2 <- safe_cor(
        x2,
        y2,
        "pearson"
    )


    legend(
        "topleft",

        legend = c(
            "161 donor-domain PBs",

            paste0(
                "Exact = ",
                sum(!bad),
                "/161"
            ),

            paste0(
                "Pearson r = ",
                sprintf(
                    "%.9f",
                    r2
                )
            )
        ),

        bty = "n",
        cex = 0.76
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
        cex = 0.75,

        xlab = "GeneBridge original\npseudobulk total",

        ylab = "ENVI measured-300\npseudobulk total",

        main = "C. Measured-gene preservation"
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
            "Measured genes; not OOF",

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
        cex = 0.76
    )


    mtext(
        "Xenium donor × spatial-domain pseudobulk concordance",

        outer = TRUE,

        cex = 1.35,

        font = 2
    )
}


# ================================================================================================
# SAVE PDF
# ================================================================================================

section("SAVE PDF")

PDF_FILE <- file.path(
    OUT,
    "01g_pseudobulk_concordance_qc_v2.pdf"
)

pdf(
    PDF_FILE,

    width = 19,

    height = 6.5,

    pointsize = 11
)

make_plot()

dev.off()

cat(
    "Saved PDF:\n",
    PDF_FILE,
    "\n",
    sep = ""
)


# ================================================================================================
# SAVE PNG
# ================================================================================================

section("SAVE PNG")

PNG_FILE <- file.path(
    OUT,
    "01g_pseudobulk_concordance_qc_v2.png"
)

png(
    PNG_FILE,

    width = 5700,

    height = 1950,

    res = 300,

    pointsize = 11
)

make_plot()

dev.off()

cat(
    "Saved PNG:\n",
    PNG_FILE,
    "\n",
    sep = ""
)


# ================================================================================================
# FINAL
# ================================================================================================

section("FINAL SUMMARY")

cat(
    "Panel A entries              :",
    length(paper_vec),
    "\n"
)

cat(
    "Paper-original exact entries:",
    sum(exact),
    "/",
    length(exact),
    "\n"
)

cat(
    "Paper-original exact PBs    :",
    sum(!bad),
    "/",
    nrow(side),
    "\n"
)

cat(
    "ENVI measured exact PBs     :",
    sum(
        side$Original
        ==
        side$ENVI_measured300
    ),
    "/",
    nrow(side),
    "\n"
)

cat(
    "\nFINAL STATUS: PLOT COMPLETE\n"
)
