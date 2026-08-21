#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("clusterProfiler is not installed.")
    }

    if (!requireNamespace("msigdbr", quietly = TRUE)) {
        stop("msigdbr is not installed.")
    }

    library(clusterProfiler)
    library(msigdbr)
})


ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_DIR <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_celltype/results/",
    "significant_candidates/L6_validation"
)

GENE_INFO_FILE <- file.path(
    ROOT,
    "outputs/imputation_full/DE/envi_celltype/pseudobulk/",
    "ENVI_full34987_gene_info.csv"
)

OUT <- file.path(
    INPUT_DIR,
    "pathway_enrichment"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


cat(strrep("=", 110), "\n")
cat("L6 EX PATHWAY ENRICHMENT\n")
cat(strrep("=", 110), "\n")


# =============================================================================
# Input
# =============================================================================

SCZ_FILE <- file.path(
    INPUT_DIR,
    "L6_Ex_77_higher_in_SCZ_genes.txt"
)

NTC_FILE <- file.path(
    INPUT_DIR,
    "L6_Ex_64_higher_in_NTC_genes.txt"
)


for (x in c(
    SCZ_FILE,
    NTC_FILE,
    GENE_INFO_FILE
)) {

    if (!file.exists(x)) {
        stop(
            paste(
                "Missing input:",
                x
            )
        )
    }
}


scz_genes <- unique(
    trimws(
        readLines(
            SCZ_FILE
        )
    )
)

ntc_genes <- unique(
    trimws(
        readLines(
            NTC_FILE
        )
    )
)

scz_genes <- scz_genes[
    nzchar(
        scz_genes
    )
]

ntc_genes <- ntc_genes[
    nzchar(
        ntc_genes
    )
]


stopifnot(
    length(scz_genes) == 77
)

stopifnot(
    length(ntc_genes) == 64
)


# =============================================================================
# Background universe
#
# IMPORTANT:
# Only genes that could have been discovered in the ENVI-imputed DE analysis
# belong in the enrichment universe.
# =============================================================================

gene_info <- read.csv(
    GENE_INFO_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

background <- unique(
    gene_info$gene[
        gene_info$expression_source
        == "envi_imputed"
    ]
)

stopifnot(
    length(background) == 34687
)


cat("\nGene sets:\n")
cat(
    "Higher in SCZ :",
    length(scz_genes),
    "\n"
)

cat(
    "Higher in NTC :",
    length(ntc_genes),
    "\n"
)

cat(
    "Background    :",
    length(background),
    "\n"
)


# =============================================================================
# Robust msigdbr helper
#
# Handles both newer:
#     collection / subcollection
#
# and older:
#     category / subcategory
# APIs.
# =============================================================================

get_msig <- function(
    collection,
    subcollection = NULL
) {

    f <- names(
        formals(
            msigdbr::msigdbr
        )
    )


    if (
        "collection"
        %in%
        f
    ) {

        args <- list(
            species = "Homo sapiens",
            collection = collection
        )

        if (
            !is.null(
                subcollection
            )
        ) {
            args$subcollection <- subcollection
        }

    } else {

        args <- list(
            species = "Homo sapiens",
            category = collection
        )

        if (
            !is.null(
                subcollection
            )
        ) {
            args$subcategory <- subcollection
        }
    }


    x <- do.call(
        msigdbr::msigdbr,
        args
    )


    if (
        nrow(x) == 0
    ) {

        stop(
            paste(
                "No MSigDB terms returned for",
                collection,
                ifelse(
                    is.null(subcollection),
                    "",
                    subcollection
                )
            )
        )
    }


    required <- c(
        "gs_name",
        "gene_symbol"
    )


    if (
        !all(
            required
            %in%
            colnames(x)
        )
    ) {

        stop(
            paste(
                "Unexpected msigdbr columns:",
                paste(
                    colnames(x),
                    collapse = ", "
                )
            )
        )
    }


    x <- unique(
        x[
            ,
            c(
                "gs_name",
                "gene_symbol"
            )
        ]
    )


    colnames(x) <- c(
        "term",
        "gene"
    )


    return(x)
}


# =============================================================================
# Collections
# =============================================================================

cat("\nLoading MSigDB collections...\n")


hallmark <- get_msig(
    "H"
)


go_bp <- get_msig(
    "C5",
    "GO:BP"
)


reactome <- get_msig(
    "C2",
    "CP:REACTOME"
)


cat(
    "Hallmark gene sets:",
    length(
        unique(
            hallmark$term
        )
    ),
    "\n"
)

cat(
    "GO BP gene sets  :",
    length(
        unique(
            go_bp$term
        )
    ),
    "\n"
)

cat(
    "Reactome sets    :",
    length(
        unique(
            reactome$term
        )
    ),
    "\n"
)


# =============================================================================
# ORA helper
# =============================================================================

run_ora <- function(
    genes,
    direction,
    database,
    term2gene
) {

    # Keep only genes belonging to our tested background.
    genes_use <- intersect(
        genes,
        background
    )


    # Restrict annotation database itself to the tested universe.
    term2gene_use <- term2gene[
        term2gene$gene
        %in%
        background,
        ,
        drop = FALSE
    ]


    cat(
        "\n",
        strrep("-", 90),
        "\n",
        sep = ""
    )

    cat(
        direction,
        "|",
        database,
        "\n"
    )

    cat(
        "Input genes:",
        length(
            genes_use
        ),
        "\n"
    )


    fit <- clusterProfiler::enricher(
        gene = genes_use,
        universe = background,
        TERM2GENE = term2gene_use,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 5,
        maxGSSize = 500
    )


    if (
        is.null(
            fit
        )
    ) {

        cat(
            "No enrichment results.\n"
        )

        return(
            data.frame()
        )
    }


    result <- as.data.frame(
        fit
    )


    if (
        nrow(result) == 0
    ) {

        cat(
            "No enrichment results.\n"
        )

        return(
            result
        )
    }


    result$direction <- direction

    result$database <- database


    # Make names explicit for advisor/report use.
    result <- result[
        ,
        c(
            "direction",
            "database",
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "Count",
            "pvalue",
            "p.adjust",
            "qvalue",
            "geneID"
        )
    ]


    colnames(result) <- c(
        "direction",
        "database",
        "pathway_id",
        "pathway",
        "GeneRatio",
        "BackgroundRatio",
        "gene_count",
        "P-value",
        "FDR",
        "qvalue",
        "genes"
    )


    result <- result[
        order(
            result$FDR,
            result$`P-value`
        ),
        ,
        drop = FALSE
    ]


    outfile <- file.path(
        OUT,
        paste0(
            "L6_Ex_",
            direction,
            "_",
            database,
            "_enrichment.csv"
        )
    )


    write.csv(
        result,
        outfile,
        row.names = FALSE
    )


    significant <- result[
        result$FDR < 0.05,
        ,
        drop = FALSE
    ]


    cat(
        "Pathways tested :",
        nrow(result),
        "\n"
    )

    cat(
        "FDR < 0.05      :",
        nrow(significant),
        "\n"
    )


    if (
        nrow(result) > 0
    ) {

        cat(
            "\nTop pathways:\n"
        )

        print(
            head(
                result[
                    ,
                    c(
                        "pathway",
                        "gene_count",
                        "P-value",
                        "FDR",
                        "genes"
                    )
                ],
                10
            ),
            row.names = FALSE
        )
    }


    return(
        result
    )
}


# =============================================================================
# Run six analyses
# =============================================================================

results <- list()


results[["SCZ_Hallmark"]] <- run_ora(
    scz_genes,
    "Higher_in_SCZ",
    "Hallmark",
    hallmark
)


results[["SCZ_GO_BP"]] <- run_ora(
    scz_genes,
    "Higher_in_SCZ",
    "GO_BP",
    go_bp
)


results[["SCZ_Reactome"]] <- run_ora(
    scz_genes,
    "Higher_in_SCZ",
    "Reactome",
    reactome
)


results[["NTC_Hallmark"]] <- run_ora(
    ntc_genes,
    "Higher_in_NTC",
    "Hallmark",
    hallmark
)


results[["NTC_GO_BP"]] <- run_ora(
    ntc_genes,
    "Higher_in_NTC",
    "GO_BP",
    go_bp
)


results[["NTC_Reactome"]] <- run_ora(
    ntc_genes,
    "Higher_in_NTC",
    "Reactome",
    reactome
)


# =============================================================================
# Combine
# =============================================================================

nonempty <- results[
    vapply(
        results,
        nrow,
        integer(1)
    ) > 0
]


if (
    length(nonempty) > 0
) {

    combined <- do.call(
        rbind,
        nonempty
    )

    rownames(
        combined
    ) <- NULL


    combined <- combined[
        order(
            combined$direction,
            combined$database,
            combined$FDR,
            combined$`P-value`
        ),
        ,
        drop = FALSE
    ]


    write.csv(
        combined,
        file.path(
            OUT,
            "L6_Ex_pathway_enrichment_all.csv"
        ),
        row.names = FALSE
    )


    significant <- combined[
        combined$FDR < 0.05,
        ,
        drop = FALSE
    ]


    write.csv(
        significant,
        file.path(
            OUT,
            "L6_Ex_pathway_enrichment_FDR05.csv"
        ),
        row.names = FALSE
    )


} else {

    combined <- data.frame()

    significant <- data.frame()
}


# =============================================================================
# Summary
# =============================================================================

summary_rows <- list()

j <- 1


for (
    direction in c(
        "Higher_in_SCZ",
        "Higher_in_NTC"
    )
) {

    for (
        database in c(
            "Hallmark",
            "GO_BP",
            "Reactome"
        )
    ) {

        key <- if (
            direction == "Higher_in_SCZ"
        ) {
            paste0(
                "SCZ_",
                database
            )
        } else {
            paste0(
                "NTC_",
                database
            )
        }


        x <- results[
            [
                key
            ]
        ]


        if (
            nrow(x) > 0
        ) {

            summary_rows[
                [
                    j
                ]
            ] <- data.frame(
                direction = direction,
                database = database,

                pathways_tested =
                    nrow(x),

                nominal_P05 =
                    sum(
                        x$`P-value`
                        < 0.05
                    ),

                FDR05 =
                    sum(
                        x$FDR
                        < 0.05
                    ),

                smallest_P =
                    min(
                        x$`P-value`
                    ),

                smallest_FDR =
                    min(
                        x$FDR
                    ),

                top_pathway =
                    x$pathway[
                        which.min(
                            x$`P-value`
                        )
                    ],

                stringsAsFactors = FALSE
            )

        } else {

            summary_rows[
                [
                    j
                ]
            ] <- data.frame(
                direction = direction,
                database = database,
                pathways_tested = 0,
                nominal_P05 = 0,
                FDR05 = 0,
                smallest_P = NA,
                smallest_FDR = NA,
                top_pathway = NA,
                stringsAsFactors = FALSE
            )
        }


        j <- j + 1
    }
}


summary <- do.call(
    rbind,
    summary_rows
)


write.csv(
    summary,
    file.path(
        OUT,
        "L6_Ex_pathway_enrichment_summary.csv"
    ),
    row.names = FALSE
)


cat(
    "\n",
    strrep("=", 110),
    "\n",
    sep = ""
)

cat(
    "PATHWAY ENRICHMENT SUMMARY\n"
)

cat(
    strrep("=", 110),
    "\n"
)


print(
    summary,
    row.names = FALSE
)


cat(
    "\nTotal FDR < 0.05 pathways:",
    ifelse(
        nrow(significant) > 0,
        nrow(significant),
        0
    ),
    "\n"
)


cat(
    "\nIMPORTANT INTERPRETATION:\n"
)

cat(
    "These are exploratory pathway results from genes significant ",
    "at within-L6 FDR < 0.05.\n",
    sep = ""
)

cat(
    "None of the underlying gene-level results survived the stricter ",
    "global correction across all 10 cell types.\n",
    sep = ""
)


cat(
    "\nOutput:\n",
    OUT,
    "\n",
    sep = ""
)


cat(
    "\n",
    strrep("=", 110),
    "\n",
    sep = ""
)

cat(
    "SUCCESS\n"
)

cat(
    strrep("=", 110),
    "\n"
)
