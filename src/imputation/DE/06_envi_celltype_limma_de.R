#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(limma)
})


ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PB_DIR <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_celltype/pseudobulk"
)

OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_celltype/results"
)

META_FILE <- file.path(
    ROOT,
    "data/metadata/patient_xenium_visium_24_common_with_dx.csv"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


cat(strrep("=", 110), "\n")
cat("CELL-TYPE-SPECIFIC ENVI DE: SCZ vs NTC\n")
cat(strrep("=", 110), "\n")


# =============================================================================
# Files
# =============================================================================

file_manifest <- read.csv(
    file.path(
        PB_DIR,
        "celltype_file_manifest.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)


gene_info <- read.csv(
    file.path(
        PB_DIR,
        "ENVI_full34987_gene_info.csv"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
)


stopifnot(
    nrow(gene_info) == 34987
)

stopifnot(
    sum(
        gene_info$expression_source
        == "measured_xenium"
    ) == 300
)

stopifnot(
    sum(
        gene_info$expression_source
        == "envi_imputed"
    ) == 34687
)


# =============================================================================
# Metadata
# =============================================================================

meta <- read.csv(
    META_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (!"patient_id" %in% colnames(meta)) {
    stop("patient_id missing")
}


age_col <- if (
    "AGE" %in% colnames(meta)
) {
    "AGE"
} else if (
    "Age" %in% colnames(meta)
) {
    "Age"
} else {
    stop("Age column missing")
}


sex_col <- if (
    "SEX" %in% colnames(meta)
) {
    "SEX"
} else if (
    "Sex" %in% colnames(meta)
) {
    "Sex"
} else {
    stop("Sex column missing")
}


meta$patient_id <- trimws(
    as.character(
        meta$patient_id
    )
)


meta$Dx <- toupper(
    trimws(
        as.character(
            meta$Dx
        )
    )
)


meta$Age_DE <- as.numeric(
    meta[[age_col]]
)


meta$Sex_DE <- toupper(
    trimws(
        as.character(
            meta[[sex_col]]
        )
    )
)


# =============================================================================
# Result containers
# =============================================================================

result_list <- list()

full_result_list <- list()


# =============================================================================
# Cell type loop
# =============================================================================

for (
    i in seq_len(
        nrow(
            file_manifest
        )
    )
) {

    cell_type <- file_manifest$cell_type[i]

    slug <- file_manifest$slug[i]

    pb_file <- file_manifest$pseudobulk_file[i]


    cat(
        "\n",
        strrep("=", 110),
        "\n",
        sep = ""
    )

    cat(
        sprintf(
            "[%02d/%02d] %s\n",
            i,
            nrow(file_manifest),
            cell_type
        )
    )

    cat(
        strrep("=", 110),
        "\n"
    )


    # -------------------------------------------------------------------------
    # Pseudobulk
    # -------------------------------------------------------------------------

    pb <- read.csv(
        pb_file,
        row.names = 1,
        check.names = FALSE
    )


    stopifnot(
        nrow(pb) == 23
    )

    stopifnot(
        ncol(pb) == 34987
    )

    stopifnot(
        identical(
            colnames(pb),
            gene_info$gene
        )
    )


    donor_meta <- meta[
        match(
            rownames(pb),
            meta$patient_id
        ),
        ,
        drop = FALSE
    ]


    if (
        any(
            is.na(
                donor_meta$patient_id
            )
        )
    ) {

        stop(
            paste(
                cell_type,
                ": donor metadata match failed"
            )
        )
    }


    stopifnot(
        identical(
            donor_meta$patient_id,
            rownames(pb)
        )
    )


    donor_meta$Dx <- factor(
        donor_meta$Dx,
        levels = c(
            "NTC",
            "SCZ"
        )
    )


    donor_meta$Sex_DE <- factor(
        donor_meta$Sex_DE
    )


    if (
        any(
            is.na(
                donor_meta$Age_DE
            )
        )
    ) {

        stop(
            paste(
                cell_type,
                ": missing age"
            )
        )
    }


    stopifnot(
        sum(
            donor_meta$Dx
            == "SCZ"
        ) == 12
    )

    stopifnot(
        sum(
            donor_meta$Dx
            == "NTC"
        ) == 11
    )


    # -------------------------------------------------------------------------
    # Normalize within this cell type
    #
    # These are ENVI expected/count-scale values, not raw sequencing counts.
    # -------------------------------------------------------------------------

    counts <- as.matrix(
        pb
    )

    storage.mode(
        counts
    ) <- "double"


    if (
        any(
            !is.finite(
                counts
            )
        )
    ) {

        stop(
            paste(
                cell_type,
                ": non-finite values"
            )
        )
    }


    if (
        any(
            counts < 0
        )
    ) {

        stop(
            paste(
                cell_type,
                ": negative values"
            )
        )
    }


    library_size <- rowSums(
        counts
    )


    if (
        any(
            library_size <= 0
        )
    ) {

        stop(
            paste(
                cell_type,
                ": zero library size"
            )
        )
    }


    cpm <- sweep(
        counts,
        1,
        library_size,
        "/"
    ) * 1e6


    log2cpm <- log2(
        cpm + 1
    )


    # -------------------------------------------------------------------------
    # Model
    #
    # DxSCZ = SCZ - NTC
    # -------------------------------------------------------------------------

    Age_centered <- (
        donor_meta$Age_DE
        - mean(
            donor_meta$Age_DE
        )
    )


    design <- model.matrix(
        ~ Dx + Age_centered + Sex_DE,
        data = donor_meta
    )


    if (
        !"DxSCZ"
        %in%
        colnames(
            design
        )
    ) {

        stop(
            paste(
                cell_type,
                ": DxSCZ coefficient missing"
            )
        )
    }


    if (
        qr(
            design
        )$rank
        !=
        ncol(
            design
        )
    ) {

        stop(
            paste(
                cell_type,
                ": design is not full rank"
            )
        )
    }


    # -------------------------------------------------------------------------
    # limma
    # -------------------------------------------------------------------------

    Y <- t(
        log2cpm
    )


    fit <- lmFit(
        Y,
        design
    )


    fit <- eBayes(
        fit,
        trend = TRUE,
        robust = TRUE
    )


    tt <- topTable(
        fit,
        coef = "DxSCZ",
        number = Inf,
        sort.by = "none",
        adjust.method = "none"
    )


    tt <- tt[
        colnames(pb),
        ,
        drop = FALSE
    ]


    stopifnot(
        identical(
            rownames(tt),
            colnames(pb)
        )
    )


    # -------------------------------------------------------------------------
    # Full 34,987
    # -------------------------------------------------------------------------

    full <- data.frame(
        cell_type =
            cell_type,

        gene =
            rownames(tt),

        expression_source =
            gene_info$expression_source,

        mean_log2CPM_NTC =
            colMeans(
                log2cpm[
                    donor_meta$Dx
                    == "NTC",
                    ,
                    drop = FALSE
                ]
            ),

        mean_log2CPM_SCZ =
            colMeans(
                log2cpm[
                    donor_meta$Dx
                    == "SCZ",
                    ,
                    drop = FALSE
                ]
            ),

        log2FC =
            tt$logFC,

        P.value =
            tt$P.Value,

        stringsAsFactors = FALSE
    )


    full$FDR_within_celltype <- p.adjust(
        full$P.value,
        method = "BH"
    )


    full$direction <- ifelse(
        full$log2FC > 0,
        "Higher_in_SCZ",
        ifelse(
            full$log2FC < 0,
            "Higher_in_NTC",
            "No_change"
        )
    )


    full_result_list[[cell_type]] <- full


    # -------------------------------------------------------------------------
    # Primary: 34,687 imputed-only genes
    #
    # Recalculate BH only among imputed genes.
    # -------------------------------------------------------------------------

    imp <- full[
        full$expression_source
        == "envi_imputed",
        ,
        drop = FALSE
    ]


    stopifnot(
        nrow(imp)
        == 34687
    )


    imp$FDR_within_celltype <- p.adjust(
        imp$P.value,
        method = "BH"
    )


    result_list[[cell_type]] <- imp


    cat(
        "Genes tested       :",
        nrow(imp),
        "\n"
    )

    cat(
        "Higher in SCZ      :",
        sum(
            imp$log2FC > 0
        ),
        "\n"
    )

    cat(
        "Higher in NTC      :",
        sum(
            imp$log2FC < 0
        ),
        "\n"
    )

    nominal_now <- imp[
        imp$P.value < 0.05,
        ,
        drop = FALSE
    ]

    fdr_now <- imp[
        imp$FDR_within_celltype < 0.05,
        ,
        drop = FALSE
    ]

    cat(
        "Nominal P < 0.05   :",
        nrow(nominal_now),
        "\n"
    )

    cat(
        "  Higher in SCZ    :",
        sum(nominal_now$log2FC > 0),
        "\n"
    )

    cat(
        "  Higher in NTC    :",
        sum(nominal_now$log2FC < 0),
        "\n"
    )

    cat(
        "Within-cell FDR<.05:",
        nrow(fdr_now),
        "\n"
    )

    cat(
        "  Higher in SCZ    :",
        sum(fdr_now$log2FC > 0),
        "\n"
    )

    cat(
        "  Higher in NTC    :",
        sum(fdr_now$log2FC < 0),
        "\n"
    )

}


# =============================================================================
# Combine all 10 cell types
# =============================================================================

combined <- do.call(
    rbind,
    result_list
)


rownames(
    combined
) <- NULL


expected_tests <- (
    10
    * 34687
)


stopifnot(
    nrow(combined)
    == expected_tests
)


# Strict correction across ALL cell-type × gene tests.
combined$FDR_global_10celltypes <- p.adjust(
    combined$P.value,
    method = "BH"
)


# =============================================================================
# Save per-cell-type results AFTER global FDR exists
# =============================================================================

summary_rows <- list()


for (
    i in seq_len(
        nrow(
            file_manifest
        )
    )
) {

    cell_type <- file_manifest$cell_type[i]

    slug <- file_manifest$slug[i]


    x <- combined[
        combined$cell_type
        == cell_type,
        ,
        drop = FALSE
    ]


    x <- x[
        order(
            x$FDR_within_celltype,
            x$P.value
        ),
        ,
        drop = FALSE
    ]


    # Internal rich table.
    write.csv(
        x,
        file.path(
            OUT,
            paste0(
                "ENVI_",
                slug,
                "_imputed34687_DE_SCZ_vs_NTC_full.csv"
            )
        ),
        row.names = FALSE
    )


    # Advisor table: requested four columns.
    advisor <- x[
        ,
        c(
            "gene",
            "P.value",
            "FDR_within_celltype",
            "log2FC"
        )
    ]


    colnames(
        advisor
    ) <- c(
        "gene",
        "P-value",
        "FDR",
        "log2FC"
    )


    write.csv(
        advisor,
        file.path(
            OUT,
            paste0(
                "ENVI_",
                slug,
                "_imputed34687_DE_advisor_table.csv"
            )
        ),
        row.names = FALSE
    )


    nominal <- x[
        x$P.value < 0.05,
        ,
        drop = FALSE
    ]


    sig_within <- x[
        x$FDR_within_celltype < 0.05,
        ,
        drop = FALSE
    ]


    sig_global <- x[
        x$FDR_global_10celltypes < 0.05,
        ,
        drop = FALSE
    ]


    summary_rows[[i]] <- data.frame(
        cell_type =
            cell_type,

        genes_tested =
            nrow(x),

        higher_SCZ =
            sum(
                x$log2FC > 0
            ),

        higher_NTC =
            sum(
                x$log2FC < 0
            ),

        nominal_P05 =
            nrow(
                nominal
            ),

        nominal_SCZ =
            sum(
                nominal$log2FC > 0
            ),

        nominal_NTC =
            sum(
                nominal$log2FC < 0
            ),

        FDR05_within_celltype =
            nrow(
                sig_within
            ),

        FDR05_global_10celltypes =
            nrow(
                sig_global
            ),

        smallest_P =
            min(
                x$P.value
            ),

        top_gene =
            x$gene[
                which.min(
                    x$P.value
                )
            ],

        stringsAsFactors = FALSE
    )

}


summary <- do.call(
    rbind,
    summary_rows
)


# =============================================================================
# Combined outputs
# =============================================================================

combined <- combined[
    order(
        combined$FDR_global_10celltypes,
        combined$P.value
    ),
    ,
    drop = FALSE
]


write.csv(
    combined,
    file.path(
        OUT,
        "ENVI_10celltypes_imputed34687_DE_combined.csv.gz"
    ),
    row.names = FALSE
)


write.csv(
    summary,
    file.path(
        OUT,
        "ENVI_10celltypes_DE_summary.csv"
    ),
    row.names = FALSE
)


# Significant across all cell-type × gene tests.
global_sig <- combined[
    combined$FDR_global_10celltypes < 0.05,
    ,
    drop = FALSE
]


write.csv(
    global_sig,
    file.path(
        OUT,
        "ENVI_10celltypes_DE_globalFDR05.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# Print summary
# =============================================================================

cat(
    "\n",
    strrep("=", 130),
    "\n",
    sep = ""
)

cat(
    "CELL-TYPE DE SUMMARY\n"
)

cat(
    strrep("=", 130),
    "\n"
)


print(
    summary,
    row.names = FALSE
)


cat(
    "\nTotal cell-type × gene tests:",
    nrow(combined),
    "\n"
)


cat(
    "Global FDR < 0.05:",
    nrow(global_sig),
    "\n"
)


if (
    nrow(global_sig) > 0
) {

    cat(
        "\nTop globally significant results:\n"
    )

    print(
        head(
            global_sig[
                ,
                c(
                    "cell_type",
                    "gene",
                    "P.value",
                    "FDR_within_celltype",
                    "FDR_global_10celltypes",
                    "log2FC",
                    "direction"
                )
            ],
            30
        ),
        row.names = FALSE
    )
}


cat(
    "\n",
    strrep("=", 130),
    "\n",
    sep = ""
)

cat(
    "SUCCESS\n"
)

cat(
    strrep("=", 130),
    "\n"
)

cat(
    "Summary: ",
    file.path(
        OUT,
        "ENVI_10celltypes_DE_summary.csv"
    ),
    "\n",
    sep = ""
)
