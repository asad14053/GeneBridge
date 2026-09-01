#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(SummarizedExperiment)
})

# ================================================================================================
# 01h_plot_hist_qq_pca_qc.R
#
# PURPOSE
# -------
# Additional QC for pseudobulk concordance:
#
# A. Histogram:
#       log1p(GeneBridge original) - log1p(Paper)
#
# B. Q-Q plot:
#       non-zero log1p residuals only
#
# C. PCA:
#       Paper
#       GeneBridge original
#       ENVI measured-300 preserved control
#
# PCA uses log2(CPM + 1) before PCA to reduce library-size effects.
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


ENVI_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
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


read_pb <- function(path) {

    x <- read.csv(
        path,
        row.names = 1,
        check.names = FALSE
    )

    x <- as.matrix(x)

    storage.mode(x) <- "double"

    x
}


logcpm <- function(x) {

    lib <- rowSums(x)

    if (any(lib <= 0)) {
        stop("Found pseudobulk with zero library size.")
    }

    cpm <- sweep(
        x,
        1,
        lib,
        "/"
    ) * 1e6

    log2(
        cpm + 1
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
        ENVI_COUNTS
    )
) {

    if (!file.exists(f)) {
        stop(
            "Missing file:\n",
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


paper_meta <- as.data.frame(
    colData(paper)
)


paper_rd <- as.data.frame(
    rowData(paper)
)


paper_counts <- assay(
    paper,
    "counts"
)


paper_genes <- as.character(
    paper_rd$Symbol
)


paper_all <- t(
    as.matrix(
        paper_counts
    )
)


storage.mode(
    paper_all
) <- "double"


colnames(
    paper_all
) <- paper_genes


cat(
    "Paper:",
    nrow(paper_all),
    "pseudobulks x",
    ncol(paper_all),
    "genes\n"
)


# ================================================================================================
# LOAD GENEBRIDGE + ENVI
# ================================================================================================

section("LOAD GENEBRIDGE AND ENVI")


original <- read_pb(
    OUR_COUNTS
)


envi <- read_pb(
    ENVI_COUNTS
)


meta <- read.csv(
    OUR_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


cat(
    "Original:",
    nrow(original),
    "x",
    ncol(original),
    "\n"
)


cat(
    "ENVI measured-300:",
    nrow(envi),
    "x",
    ncol(envi),
    "\n"
)


# ================================================================================================
# ALIGN ORIGINAL ROWS
# ================================================================================================

section("ALIGN PSEUDOBULKS")


if (
    all(
        meta$pseudobulk_id
        %in%
        rownames(original)
    )
) {

    original <- original[
        meta$pseudobulk_id,
        ,
        drop = FALSE
    ]
}


if (
    all(
        meta$pseudobulk_id
        %in%
        rownames(envi)
    )
) {

    envi <- envi[
        meta$pseudobulk_id,
        ,
        drop = FALSE
    ]
}


# ================================================================================================
# ALIGN GENES
# ================================================================================================

section("ALIGN GENES")


if (
    !setequal(
        colnames(original),
        paper_genes
    )
) {
    stop("Paper and Original 300-gene sets differ.")
}


if (
    !setequal(
        colnames(original),
        colnames(envi)
    )
) {
    stop("Original and ENVI 300-gene sets differ.")
}


genes <- colnames(
    original
)


paper_all <- paper_all[
    ,
    genes,
    drop = FALSE
]


envi <- envi[
    ,
    genes,
    drop = FALSE
]


cat(
    "Aligned genes:",
    length(genes),
    "\n"
)


# ================================================================================================
# SUBSET PAPER N24 -> N23
# ================================================================================================

section("SUBSET PAPER N24 -> N23")


our_donors <- unique(
    as.character(
        meta$BrNum
    )
)


keep <- (
    as.character(
        paper_meta$BrNum
    )
    %in%
    our_donors
)


paper_meta23 <- paper_meta[
    keep,
    ,
    drop = FALSE
]


paper23 <- paper_all[
    keep,
    ,
    drop = FALSE
]


paper_key <- paste(
    paper_meta23$BrNum,
    paper_meta23$predictions_smooth,
    sep = "::"
)


our_key <- paste(
    meta$BrNum,
    meta$predictions_smooth,
    sep = "::"
)


idx <- match(
    our_key,
    paper_key
)


if (any(is.na(idx))) {

    stop(
        "Could not align all 161 paper pseudobulks."
    )
}


paper23 <- paper23[
    idx,
    ,
    drop = FALSE
]


paper_meta23 <- paper_meta23[
    idx,
    ,
    drop = FALSE
]


cat(
    "Paper N23 aligned:",
    nrow(paper23),
    "x",
    ncol(paper23),
    "\n"
)


# ================================================================================================
# HISTOGRAM / RESIDUAL DATA
# ================================================================================================

section("CALCULATE RESIDUALS")


# Use log1p so huge-expression genes do not dominate the difference distribution.

residual_matrix <- (
    log1p(original)
    -
    log1p(paper23)
)


residuals <- as.vector(
    residual_matrix
)


nonzero <- residuals[
    residuals != 0
]


cat(
    "Total residuals:",
    length(residuals),
    "\n"
)


cat(
    "Exactly zero:",
    sum(residuals == 0),
    "\n"
)


cat(
    "Non-zero:",
    length(nonzero),
    "\n"
)


cat(
    "Median absolute residual:",
    median(
        abs(residuals)
    ),
    "\n"
)


cat(
    "Maximum absolute residual:",
    max(
        abs(residuals)
    ),
    "\n"
)


# ================================================================================================
# PCA INPUT
# ================================================================================================

section("PREPARE PCA")


paper_norm <- logcpm(
    paper23
)


original_norm <- logcpm(
    original
)


envi_norm <- logcpm(
    envi
)


# Stack the three sources:
#
# 161 Paper
# 161 Original
# 161 ENVI measured preserved
#
# total = 483 profiles.

pca_input <- rbind(
    paper_norm,
    original_norm,
    envi_norm
)


source <- factor(
    rep(
        c(
            "Paper",
            "GeneBridge original",
            "ENVI measured-300"
        ),
        each = nrow(original)
    ),
    levels = c(
        "Paper",
        "GeneBridge original",
        "ENVI measured-300"
    )
)


pb_id <- rep(
    meta$pseudobulk_id,
    times = 3
)


donor <- rep(
    meta$BrNum,
    times = 3
)


domain <- rep(
    meta$layer_annotation,
    times = 3
)


cat(
    "PCA profiles:",
    nrow(pca_input),
    "\n"
)


cat(
    "Genes:",
    ncol(pca_input),
    "\n"
)


# ================================================================================================
# PCA
# ================================================================================================

section("RUN PCA")


pca <- prcomp(
    pca_input,
    center = TRUE,
    scale. = TRUE
)


variance <- (
    pca$sdev^2
    /
    sum(
        pca$sdev^2
    )
)


pc <- data.frame(

    PC1 =
        pca$x[, 1],

    PC2 =
        pca$x[, 2],

    PC3 =
        pca$x[, 3],

    source =
        source,

    pseudobulk_id =
        pb_id,

    BrNum =
        donor,

    layer =
        domain,

    stringsAsFactors = FALSE
)


PCA_CSV <- file.path(
    OUT,
    "01h_pca_coordinates.csv"
)


write.csv(
    pc,
    PCA_CSV,
    row.names = FALSE
)


cat(
    "PC1 variance:",
    sprintf(
        "%.2f%%",
        100 * variance[1]
    ),
    "\n"
)


cat(
    "PC2 variance:",
    sprintf(
        "%.2f%%",
        100 * variance[2]
    ),
    "\n"
)


# ================================================================================================
# PCA MATCHED-DISTANCE QC
# ================================================================================================

section("PCA MATCHED PROFILE DISTANCES")


n <- nrow(
    original
)


paper_pc <- pca$x[
    1:n,
    1:10,
    drop = FALSE
]


original_pc <- pca$x[
    (n + 1):(2 * n),
    1:10,
    drop = FALSE
]


envi_pc <- pca$x[
    (2 * n + 1):(3 * n),
    1:10,
    drop = FALSE
]


paper_original_distance <- sqrt(
    rowSums(
        (
            paper_pc -
            original_pc
        )^2
    )
)


original_envi_distance <- sqrt(
    rowSums(
        (
            original_pc -
            envi_pc
        )^2
    )
)


cat(
    "Median Paper -> Original PCA distance:",
    median(
        paper_original_distance
    ),
    "\n"
)


cat(
    "Max Paper -> Original PCA distance:",
    max(
        paper_original_distance
    ),
    "\n"
)


cat(
    "Median Original -> ENVI PCA distance:",
    median(
        original_envi_distance
    ),
    "\n"
)


cat(
    "Max Original -> ENVI PCA distance:",
    max(
        original_envi_distance
    ),
    "\n"
)


# ================================================================================================
# DRAW FIGURE
# ================================================================================================

make_figure <- function() {

    oldpar <- par(
        no.readonly = TRUE
    )

    on.exit(
        par(oldpar)
    )


    par(
        mfrow = c(1, 3),
        mar = c(5.2, 5.5, 4.2, 1.4),
        oma = c(0, 0, 2.5, 0),
        mgp = c(3.4, 1.0, 0),
        las = 1
    )


    # ============================================================================================
    # A. HISTOGRAM
    # ============================================================================================

    hist(
        residuals,

        breaks = 100,

        main = "A. Difference distribution",

        xlab = expression(
            log(1 + GeneBridge) -
            log(1 + Paper)
        ),

        ylab = "Frequency"
    )


    abline(
        v = 0,
        lwd = 2,
        lty = 2
    )


    legend(
        "topright",

        legend = c(
            paste0(
                "Zero = ",
                format(
                    sum(residuals == 0),
                    big.mark = ","
                )
            ),

            paste0(
                "Non-zero = ",
                format(
                    length(nonzero),
                    big.mark = ","
                )
            )
        ),

        bty = "n",
        cex = 0.8
    )


    # ============================================================================================
    # B. QQ PLOT OF NON-ZERO RESIDUALS
    # ============================================================================================

    if (
        length(nonzero) >= 3
    ) {

        qqnorm(
            nonzero,

            main = "B. Q-Q of non-zero differences",

            xlab = "Theoretical normal quantiles",

            ylab = "Observed log-scale differences",

            pch = 16,

            cex = 0.55
        )


        qqline(
            nonzero,
            lwd = 2,
            lty = 2
        )


        legend(
            "topleft",

            legend = paste0(
                "n = ",
                length(nonzero),
                " non-zero values"
            ),

            bty = "n",
            cex = 0.8
        )

    } else {

        plot.new()

        title(
            "B. Q-Q plot"
        )

        text(
            0.5,
            0.5,
            "Too few non-zero differences"
        )
    }


    # ============================================================================================
    # C. PCA
    # ============================================================================================

    pch_map <- c(
        1,
        16,
        4
    )


    plot(
        pc$PC1,
        pc$PC2,

        type = "n",

        xlab = paste0(
            "PC1 (",
            sprintf(
                "%.1f",
                100 * variance[1]
            ),
            "%)"
        ),

        ylab = paste0(
            "PC2 (",
            sprintf(
                "%.1f",
                100 * variance[2]
            ),
            "%)"
        ),

        main = "C. 300-gene pseudobulk PCA"
    )


    for (
        k in seq_along(
            levels(source)
        )
    ) {

        use <- (
            pc$source
            ==
            levels(source)[k]
        )


        points(
            pc$PC1[use],
            pc$PC2[use],

            pch = pch_map[k],

            cex = 0.75
        )
    }


    legend(
        "topright",

        legend = levels(
            source
        ),

        pch = pch_map,

        bty = "n",

        cex = 0.78
    )


    mtext(
        "Additional pseudobulk concordance diagnostics",

        outer = TRUE,

        font = 2,

        cex = 1.25
    )
}


# ================================================================================================
# SAVE PNG
# ================================================================================================

section("SAVE PNG")


PNG_FILE <- file.path(
    OUT,
    "01h_histogram_qq_pca_qc.png"
)


png(
    PNG_FILE,
    width = 5700,
    height = 1900,
    res = 300
)


make_figure()


dev.off()


cat(
    "Saved:",
    PNG_FILE,
    "\n"
)


# ================================================================================================
# SAVE PDF
# ================================================================================================

section("SAVE PDF")


PDF_FILE <- file.path(
    OUT,
    "01h_histogram_qq_pca_qc.pdf"
)


pdf(
    PDF_FILE,
    width = 19,
    height = 6.3
)


make_figure()


dev.off()


cat(
    "Saved:",
    PDF_FILE,
    "\n"
)


# ================================================================================================
# SAVE SUMMARY
# ================================================================================================

SUMMARY_FILE <- file.path(
    OUT,
    "01h_histogram_qq_pca_qc_summary.csv"
)


summary <- data.frame(

    metric = c(
        "total_gene_pseudobulk_values",
        "exact_zero_log_residuals",
        "nonzero_log_residuals",
        "pc1_variance_percent",
        "pc2_variance_percent",
        "median_paper_original_pc_distance",
        "max_paper_original_pc_distance",
        "median_original_envi_pc_distance",
        "max_original_envi_pc_distance"
    ),

    value = c(
        length(residuals),
        sum(residuals == 0),
        length(nonzero),
        100 * variance[1],
        100 * variance[2],
        median(paper_original_distance),
        max(paper_original_distance),
        median(original_envi_distance),
        max(original_envi_distance)
    ),

    stringsAsFactors = FALSE
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
    "\nInterpretation:\n"
)


cat(
    "Histogram: should be sharply concentrated at zero if reconstruction is concordant.\n"
)


cat(
    "Q-Q: diagnostic only; evaluates the shape of the non-zero residual distribution.\n"
)


cat(
    "PCA: matched Paper / Original profiles should overlap closely.\n"
)


cat(
    "Current ENVI measured-300 is a preservation control, so Original and ENVI should overlap exactly.\n"
)


cat(
    "\nFINAL STATUS: HISTOGRAM / QQ / PCA COMPLETE\n"
)
