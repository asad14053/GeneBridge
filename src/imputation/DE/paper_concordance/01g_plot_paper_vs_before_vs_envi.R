#!/usr/bin/env Rscript

# ================================================================================================
# Plot:
#
# X axis = published-paper donor × domain pseudobulk total
#
# Y axis:
#   1. GeneBridge original/before-imputation pseudobulk total
#   2. ENVI 300-gene pseudobulk total
#
# NOTE:
# Current ENVI matrix is measured/preserved 300 genes, NOT true OOF predictions.
# Therefore ENVI and GeneBridge original currently overlap exactly.
# ================================================================================================


ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"


INPUT <- file.path(
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


if (!file.exists(INPUT)) {
    stop(
        "Missing input:\n",
        INPUT
    )
}


# ================================================================================================
# LOAD
# ================================================================================================

x <- read.csv(
    INPUT,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat(
    "Rows:",
    nrow(x),
    "\n"
)


cat(
    "Columns:\n"
)

print(
    colnames(x)
)


required <- c(
    "Donor",
    "Domain",
    "Paper",
    "Original",
    "ENVI_OOF300"
)


missing <- setdiff(
    required,
    colnames(x)
)


if (length(missing) > 0) {
    stop(
        "Missing columns: ",
        paste(
            missing,
            collapse = ", "
        )
    )
}


# ================================================================================================
# METRICS
# ================================================================================================

paper <- as.numeric(
    x$Paper
)


before <- as.numeric(
    x$Original
)


# Legacy column name only.
# Current data are actually preserved/measured 300 genes.
envi300 <- as.numeric(
    x$ENVI_OOF300
)


r_before <- cor(
    paper,
    before,
    method = "pearson"
)


rho_before <- cor(
    paper,
    before,
    method = "spearman"
)


r_envi <- cor(
    paper,
    envi300,
    method = "pearson"
)


rho_envi <- cor(
    paper,
    envi300,
    method = "spearman"
)


mae_before <- mean(
    abs(
        before - paper
    )
)


mae_envi <- mean(
    abs(
        envi300 - paper
    )
)


rmse_before <- sqrt(
    mean(
        (
            before - paper
        )^2
    )
)


rmse_envi <- sqrt(
    mean(
        (
            envi300 - paper
        )^2
    )
)


cat(
    "\nPaper vs GeneBridge-before Pearson:",
    sprintf("%.12f", r_before),
    "\n"
)

cat(
    "Paper vs current ENVI-300 Pearson:",
    sprintf("%.12f", r_envi),
    "\n"
)


cat(
    "\nOriginal vs current ENVI identical:",
    identical(
        before,
        envi300
    ),
    "\n"
)


# ================================================================================================
# PLOT FUNCTION
# ================================================================================================

make_plot <- function() {

    par(
        mar = c(
            6,
            7,
            4,
            2
        ),

        mgp = c(
            4.5,
            1.3,
            0
        ),

        las = 1
    )


    lim <- range(
        c(
            paper,
            before,
            envi300
        ),
        finite = TRUE
    )


    # ----------------------------------------------------------------
    # First series:
    # Paper -> GeneBridge BEFORE
    # ----------------------------------------------------------------

    plot(
        paper,
        before,

        xlim = lim,
        ylim = lim,

        pch = 16,
        cex = 0.85,

        xlab = "Published paper donor-domain pseudobulk total",

        ylab = "GeneBridge / ENVI donor-domain pseudobulk total",

        main = "Paper-level pseudobulk concordance"
    )


    # ----------------------------------------------------------------
    # Identity reference
    # ----------------------------------------------------------------

    abline(
        a = 0,
        b = 1,
        lty = 2,
        lwd = 2
    )


    # ----------------------------------------------------------------
    # Second series:
    # Paper -> ENVI 300
    #
    # Open circles are drawn slightly larger so that if they are
    # exactly on top of GeneBridge points, both can still be recognized.
    # ----------------------------------------------------------------

    points(
        paper,
        envi300,

        pch = 1,
        cex = 1.35,
        lwd = 1.4
    )


    # ----------------------------------------------------------------
    # Highlight paper-vs-original mismatched pseudobulks
    # ----------------------------------------------------------------

    bad <- (
        paper != before
    )


    if (any(bad)) {

        points(
            paper[bad],
            before[bad],

            pch = 4,
            cex = 1.7,
            lwd = 2
        )


        labels_bad <- paste0(
            x$Donor[bad],
            " ",
            x$Domain[bad]
        )


        text(
            paper[bad],
            before[bad],

            labels = labels_bad,

            pos = 4,
            offset = 0.6,

            cex = 0.8
        )
    }


    # ----------------------------------------------------------------
    # Legend
    # ----------------------------------------------------------------

    legend(
        "topleft",

        legend = c(

            paste0(
                "GeneBridge before: r = ",
                sprintf(
                    "%.9f",
                    r_before
                )
            ),

            paste0(
                "ENVI 300 current: r = ",
                sprintf(
                    "%.9f",
                    r_envi
                )
            ),

            "Dashed line: y = x",

            "ENVI current = measured/preserved, not OOF"
        ),

        pch = c(
            16,
            1,
            NA,
            NA
        ),

        lty = c(
            NA,
            NA,
            2,
            NA
        ),

        bty = "n",

        cex = 0.82
    )
}


# ================================================================================================
# SAVE PDF
# ================================================================================================

PDF <- file.path(
    OUT,
    "01g_paper_vs_genebridge_before_vs_ENVI300.pdf"
)


pdf(
    PDF,
    width = 8.5,
    height = 7.5
)


make_plot()


dev.off()


# ================================================================================================
# SAVE PNG
# ================================================================================================

PNG <- file.path(
    OUT,
    "01g_paper_vs_genebridge_before_vs_ENVI300.png"
)


png(
    PNG,
    width = 2550,
    height = 2250,
    res = 300
)


make_plot()


dev.off()


# ================================================================================================
# SUMMARY TABLE
# ================================================================================================

summary <- data.frame(

    comparison = c(
        "Paper_vs_GeneBridge_before",
        "Paper_vs_ENVI300_current"
    ),

    Pearson = c(
        r_before,
        r_envi
    ),

    Spearman = c(
        rho_before,
        rho_envi
    ),

    MAE = c(
        mae_before,
        mae_envi
    ),

    RMSE = c(
        rmse_before,
        rmse_envi
    )
)


SUMMARY <- file.path(
    OUT,
    "01g_paper_vs_before_vs_ENVI300_metrics.csv"
)


write.csv(
    summary,
    SUMMARY,
    row.names = FALSE
)


cat(
    "\n========================================\n"
)

cat(
    "COMPARISON SUMMARY\n"
)

cat(
    "========================================\n"
)


print(
    summary,
    row.names = FALSE
)


cat(
    "\nOriginal and ENVI current identical:",
    all(
        before == envi300
    ),
    "\n"
)


cat(
    "\nSaved:\n",
    PDF,
    "\n",
    PNG,
    "\n",
    SUMMARY,
    "\n",
    sep = ""
)


cat(
    "\nFINAL STATUS: COMPLETE\n"
)
