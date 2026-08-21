#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(limma)
})

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PB_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_transcriptome",
    "ENVI_full34987_pseudobulk_countscale_23donors.csv.gz"
)

GENE_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_transcriptome",
    "ENVI_full34987_gene_info.csv"
)

META_FILE <- file.path(
    ROOT,
    "data/metadata/patient_xenium_visium_24_common_with_dx.csv"
)

OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_transcriptome"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)

cat(strrep("=", 100), "\n")
cat("ENVI TRANSCRIPTOME-WIDE SCZ vs NTC DIFFERENTIAL EXPRESSION\n")
cat(strrep("=", 100), "\n")


# =============================================================================
# Load donor pseudobulk
# =============================================================================

pb <- read.csv(
    PB_FILE,
    row.names = 1,
    check.names = FALSE
)

cat("\nPseudobulk shape:", nrow(pb), "donors x", ncol(pb), "genes\n")

stopifnot(
    nrow(pb) == 23,
    ncol(pb) == 34987
)


# =============================================================================
# Gene source
# =============================================================================

gene_info <- read.csv(
    GENE_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

stopifnot(
    nrow(gene_info) == 34987
)

stopifnot(
    identical(
        colnames(pb),
        gene_info$gene
    )
)

cat("\nGene source:\n")
print(
    table(
        gene_info$expression_source
    )
)


# =============================================================================
# Donor metadata
# =============================================================================

meta <- read.csv(
    META_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!"patient_id" %in% colnames(meta)) {
    stop("patient_id column not found")
}

if (!"Dx" %in% colnames(meta)) {
    stop("Dx column not found")
}


age_col <- if ("AGE" %in% colnames(meta)) {
    "AGE"
} else if ("Age" %in% colnames(meta)) {
    "Age"
} else {
    stop("Age column not found")
}


sex_col <- if ("SEX" %in% colnames(meta)) {
    "SEX"
} else if ("Sex" %in% colnames(meta)) {
    "Sex"
} else {
    stop("Sex column not found")
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

meta <- meta[
    meta$patient_id %in% rownames(pb),
    ,
    drop = FALSE
]

meta <- meta[
    match(
        rownames(pb),
        meta$patient_id
    ),
    ,
    drop = FALSE
]


stopifnot(
    identical(
        meta$patient_id,
        rownames(pb)
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


if (any(is.na(meta$Age_DE))) {
    stop(
        "Missing Age detected."
    )
}


meta$Dx <- factor(
    meta$Dx,
    levels = c(
        "NTC",
        "SCZ"
    )
)


meta$Sex_DE <- factor(
    meta$Sex_DE
)


cat("\nDiagnosis:\n")
print(
    table(
        meta$Dx
    )
)

cat("\nSex:\n")
print(
    table(
        meta$Sex_DE
    )
)

cat("\nAge summary:\n")
print(
    summary(
        meta$Age_DE
    )
)


stopifnot(
    sum(meta$Dx == "NTC") == 11,
    sum(meta$Dx == "SCZ") == 12
)


# =============================================================================
# Normalize donor pseudobulk
#
# Important:
# ENVI count_scale values are continuous expected expression.
# We use full 34,987-gene library normalization, then log2(CPM + 1).
# =============================================================================

counts <- as.matrix(pb)

storage.mode(counts) <- "double"


if (any(!is.finite(counts))) {
    stop(
        "Non-finite pseudobulk values detected."
    )
}


if (any(counts < 0)) {
    stop(
        "Negative pseudobulk values detected."
    )
}


library_size <- rowSums(
    counts
)


if (any(library_size <= 0)) {
    stop(
        "Zero donor library detected."
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


write.csv(
    data.frame(
        donor = rownames(pb),
        Dx = meta$Dx,
        Age = meta$Age_DE,
        Sex = meta$Sex_DE,
        full_library = library_size
    ),
    file.path(
        OUT,
        "ENVI_full34987_DE_donor_metadata.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# Design
#
# DxSCZ coefficient = adjusted SCZ - NTC difference
# =============================================================================

Age_centered <- (
    meta$Age_DE
    - mean(
        meta$Age_DE
    )
)


design <- model.matrix(
    ~ Dx + Age_centered + Sex_DE,
    data = meta
)


cat("\nDesign matrix:\n")
print(design)

cat(
    "\nDesign rank:",
    qr(design)$rank,
    "/",
    ncol(design),
    "\n"
)


write.csv(
    design,
    file.path(
        OUT,
        "ENVI_full34987_DE_design_matrix.csv"
    )
)


if (!"DxSCZ" %in% colnames(design)) {
    stop(
        paste(
            "DxSCZ coefficient missing.",
            "Design columns:",
            paste(
                colnames(design),
                collapse = ", "
            )
        )
    )
}


# =============================================================================
# Limma
#
# Matrix for limma = genes x donors
# =============================================================================

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


coef_name <- "DxSCZ"


raw <- topTable(
    fit,
    coef = coef_name,
    number = Inf,
    sort.by = "none",
    adjust.method = "none"
)


# Restore original canonical gene order.
raw <- raw[
    colnames(pb),
    ,
    drop = FALSE
]


stopifnot(
    identical(
        rownames(raw),
        colnames(pb)
    )
)


# =============================================================================
# Full 34,987-gene results
# =============================================================================

full <- data.frame(
    gene = rownames(raw),
    expression_source = gene_info$expression_source,
    mean_log2CPM_NTC = colMeans(
        log2cpm[
            meta$Dx == "NTC",
            ,
            drop = FALSE
        ]
    ),
    mean_log2CPM_SCZ = colMeans(
        log2cpm[
            meta$Dx == "SCZ",
            ,
            drop = FALSE
        ]
    ),
    log2FC = raw$logFC,
    P.value = raw$P.Value,
    stringsAsFactors = FALSE
)


full$FDR <- p.adjust(
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


full <- full[
    order(
        full$FDR,
        full$P.value
    ),
]


write.csv(
    full,
    file.path(
        OUT,
        "ENVI_full34987_DE_SCZ_vs_NTC_full.csv"
    ),
    row.names = FALSE
)


full_advisor <- full[
    ,
    c(
        "gene",
        "P.value",
        "FDR",
        "log2FC"
    )
]


colnames(full_advisor)[2] <- "P-value"


write.csv(
    full_advisor,
    file.path(
        OUT,
        "ENVI_full34987_DE_SCZ_vs_NTC_advisor_table.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# PRIMARY analysis:
# 34,687 ENVI-imputed genes only
#
# Important:
# FDR is recalculated across the 34,687 imputed genes.
# =============================================================================

imputed <- full[
    full$expression_source == "envi_imputed",
    ,
    drop = FALSE
]


stopifnot(
    nrow(imputed) == 34687
)


imputed$FDR <- p.adjust(
    imputed$P.value,
    method = "BH"
)


imputed <- imputed[
    order(
        imputed$FDR,
        imputed$P.value
    ),
]


write.csv(
    imputed,
    file.path(
        OUT,
        "ENVI_imputed34687_DE_SCZ_vs_NTC_full.csv"
    ),
    row.names = FALSE
)


imputed_advisor <- imputed[
    ,
    c(
        "gene",
        "P.value",
        "FDR",
        "log2FC"
    )
]


colnames(imputed_advisor)[2] <- "P-value"


write.csv(
    imputed_advisor,
    file.path(
        OUT,
        "ENVI_imputed34687_DE_SCZ_vs_NTC_advisor_table.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# Significant / nominal tables
# =============================================================================

imputed_fdr05 <- imputed[
    imputed$FDR < 0.05,
    ,
    drop = FALSE
]


imputed_nominal <- imputed[
    imputed$P.value < 0.05,
    ,
    drop = FALSE
]


write.csv(
    imputed_fdr05,
    file.path(
        OUT,
        "ENVI_imputed34687_DE_SCZ_vs_NTC_FDR05.csv"
    ),
    row.names = FALSE
)


write.csv(
    imputed_nominal,
    file.path(
        OUT,
        "ENVI_imputed34687_DE_SCZ_vs_NTC_nominalP05.csv"
    ),
    row.names = FALSE
)


# =============================================================================
# Summary function
# =============================================================================

report_summary <- function(df, title) {

    nominal <- df[
        df$P.value < 0.05,
        ,
        drop = FALSE
    ]

    sig <- df[
        df$FDR < 0.05,
        ,
        drop = FALSE
    ]


    cat("\n", strrep("=", 100), "\n", sep = "")
    cat(title, "\n")
    cat(strrep("=", 100), "\n")


    cat("\nALL GENES\n")

    cat(
        "Total         :",
        nrow(df),
        "\n"
    )

    cat(
        "Higher in SCZ :",
        sum(
            df$log2FC > 0
        ),
        "\n"
    )

    cat(
        "Higher in NTC :",
        sum(
            df$log2FC < 0
        ),
        "\n"
    )


    cat("\nNOMINAL P < 0.05\n")

    cat(
        "Total         :",
        nrow(nominal),
        "\n"
    )

    cat(
        "Higher in SCZ :",
        sum(
            nominal$log2FC > 0
        ),
        "\n"
    )

    cat(
        "Higher in NTC :",
        sum(
            nominal$log2FC < 0
        ),
        "\n"
    )


    cat("\nFDR < 0.05\n")

    cat(
        "Total         :",
        nrow(sig),
        "\n"
    )

    cat(
        "Higher in SCZ :",
        sum(
            sig$log2FC > 0
        ),
        "\n"
    )

    cat(
        "Higher in NTC :",
        sum(
            sig$log2FC < 0
        ),
        "\n"
    )


    cat("\nTop 20:\n")

    print(
        head(
            df[
                ,
                c(
                    "gene",
                    "P.value",
                    "FDR",
                    "log2FC",
                    "direction"
                )
            ],
            20
        ),
        row.names = FALSE
    )
}


report_summary(
    imputed,
    "PRIMARY: 34,687 ENVI-IMPUTED GENES — SCZ vs NTC"
)


report_summary(
    full,
    "FULL TRANSCRIPTOME: 34,987 GENES — SCZ vs NTC"
)


cat("\n", strrep("=", 100), "\n", sep = "")
cat("SUCCESS\n")
cat(strrep("=", 100), "\n")

cat(
    "\nPrimary advisor table:\n",
    file.path(
        OUT,
        "ENVI_imputed34687_DE_SCZ_vs_NTC_advisor_table.csv"
    ),
    "\n",
    sep = ""
)
