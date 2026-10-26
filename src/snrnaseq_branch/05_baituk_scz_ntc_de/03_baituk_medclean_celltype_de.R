#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(edgeR)
})

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

IN_DIR <- file.path(
    ROOT,
    "outputs/snrnaseq_branch/05_baituk_scz_ntc_de"
)

COUNTS_FILE <- file.path(
    IN_DIR,
    "baituk_medclean_donor_celltype_pseudobulk_counts.csv.gz"
)

META_FILE <- file.path(
    IN_DIR,
    "baituk_medclean_donor_celltype_pseudobulk_metadata.csv"
)

OUT_DIR <- file.path(
    IN_DIR,
    "celltype_DE"
)

dir.create(
    OUT_DIR,
    recursive=TRUE,
    showWarnings=FALSE
)

MIN_CELLS <- 10


section <- function(x) {
    cat("\n", paste(rep("=", 110), collapse=""), "\n", sep="")
    cat(x, "\n")
    cat(paste(rep("=", 110), collapse=""), "\n")
}


safe_name <- function(x) {
    gsub(
        "[^A-Za-z0-9]+",
        "_",
        x
    )
}


section("LOAD DONOR × MED.CLEAN CELL-TYPE PSEUDOBULK")

counts <- read.csv(
    gzfile(COUNTS_FILE),
    row.names=1,
    check.names=FALSE
)

meta <- read.csv(
    META_FILE,
    stringsAsFactors=FALSE
)

rownames(meta) <- meta$pseudobulk_id

meta <- meta[
    colnames(counts),
    ,
    drop=FALSE
]

stopifnot(
    identical(
        rownames(meta),
        colnames(counts)
    )
)

cat("Genes:", nrow(counts), "\n")
cat("Pseudobulks:", ncol(counts), "\n")
cat("Cell types:", length(unique(meta$cell_type)), "\n")

cat("\nCell types:\n")
print(
    sort(
        unique(meta$cell_type)
    )
)


# ==============================================================================================
# LOOP CELL TYPES
# ==============================================================================================

summary_list <- list()
all_results <- list()

cell_types <- sort(
    unique(meta$cell_type)
)


for (ct in cell_types) {

    section(
        paste(
            "CELL TYPE:",
            ct
        )
    )

    m <- meta[
        meta$cell_type == ct &
        meta$n_cells >= MIN_CELLS,
        ,
        drop=FALSE
    ]


    n_ntc <- length(
        unique(
            m$donor_id[
                m$diagnosis == "NTC"
            ]
        )
    )

    n_scz <- length(
        unique(
            m$donor_id[
                m$diagnosis == "SCZ"
            ]
        )
    )


    cat("Minimum cells/PB:", MIN_CELLS, "\n")
    cat("NTC donors:", n_ntc, "\n")
    cat("SCZ donors:", n_scz, "\n")


    # Require at least 3 biological donors from each group.
    if (
        n_ntc < 3 ||
        n_scz < 3
    ) {

        cat(
            "SKIP: insufficient donor replication.\n"
        )

        summary_list[[ct]] <- data.frame(
            cell_type=ct,
            NTC_donors=n_ntc,
            SCZ_donors=n_scz,
            genes_tested=0,
            FDR010=NA,
            SCZ_up_FDR010=NA,
            SCZ_down_FDR010=NA,
            status="SKIPPED"
        )

        next
    }


    cts <- counts[
        ,
        m$pseudobulk_id,
        drop=FALSE
    ]


    m$diagnosis <- factor(
        m$diagnosis,
        levels=c("NTC", "SCZ")
    )


    y <- DGEList(
        counts=as.matrix(cts),
        group=m$diagnosis
    )


    design <- model.matrix(
        ~ diagnosis,
        data=m
    )


    keep <- filterByExpr(
        y,
        design=design
    )


    cat(
        "Genes retained:",
        sum(keep),
        "/",
        nrow(y),
        "\n"
    )


    y <- y[
        keep,
        ,
        keep.lib.sizes=FALSE
    ]

    y <- calcNormFactors(
        y,
        method="TMM"
    )

    y <- estimateDisp(
        y,
        design
    )

    fit <- glmQLFit(
        y,
        design,
        robust=TRUE
    )

    qlf <- glmQLFTest(
        fit,
        coef="diagnosisSCZ"
    )

    res <- topTags(
        qlf,
        n=Inf,
        sort.by="PValue"
    )$table


    res$gene <- rownames(res)
    res$cell_type <- ct

    res$direction <- ifelse(
        res$logFC > 0,
        "SCZ_up",
        "SCZ_down"
    )


    res <- res[
        ,
        c(
            "cell_type",
            "gene",
            "logFC",
            "logCPM",
            "F",
            "PValue",
            "FDR",
            "direction"
        )
    ]


    sig <- res[
        !is.na(res$FDR) &
        res$FDR < 0.10,
        ,
        drop=FALSE
    ]


    outfile <- file.path(
        OUT_DIR,
        paste0(
            safe_name(ct),
            "_SCZ_vs_NTC_DE.csv"
        )
    )

    write.csv(
        res,
        outfile,
        row.names=FALSE
    )


    cat("FDR < 0.10:", nrow(sig), "\n")
    cat(
        "  SCZ up:",
        sum(sig$logFC > 0),
        "\n"
    )
    cat(
        "  SCZ down:",
        sum(sig$logFC < 0),
        "\n"
    )


    if (nrow(sig) > 0) {

        cat("\nTop significant genes:\n")

        print(
            head(sig, 20),
            row.names=FALSE
        )
    }


    summary_list[[ct]] <- data.frame(
        cell_type=ct,
        NTC_donors=n_ntc,
        SCZ_donors=n_scz,
        genes_tested=nrow(res),
        FDR010=nrow(sig),
        SCZ_up_FDR010=sum(sig$logFC > 0),
        SCZ_down_FDR010=sum(sig$logFC < 0),
        status="PASS"
    )


    all_results[[ct]] <- res
}


# ==============================================================================================
# COMBINED OUTPUT
# ==============================================================================================

summary_df <- do.call(
    rbind,
    summary_list
)

rownames(summary_df) <- NULL


combined <- do.call(
    rbind,
    all_results
)

rownames(combined) <- NULL


# FDR across every gene × cell-type test in addition
# to the within-cell-type edgeR FDR.
combined$FDR_global <- p.adjust(
    combined$PValue,
    method="BH"
)


SUMMARY_FILE <- file.path(
    IN_DIR,
    "baituk_medclean_celltype_DE_summary.csv"
)

COMBINED_FILE <- file.path(
    IN_DIR,
    "baituk_medclean_celltype_SCZ_vs_NTC_DE_all.csv.gz"
)


write.csv(
    summary_df,
    SUMMARY_FILE,
    row.names=FALSE
)

con <- gzfile(
    COMBINED_FILE,
    "wt"
)

write.csv(
    combined,
    con,
    row.names=FALSE
)

close(con)


section("FINAL CELL-TYPE SUMMARY")

print(
    summary_df,
    row.names=FALSE
)

cat("\nSaved summary:\n")
cat(SUMMARY_FILE, "\n")

cat("\nSaved combined results:\n")
cat(COMBINED_FILE, "\n")

cat(
    "\nFINAL STATUS: PASS — BAITUK MED.CLEAN CELL-TYPE SCZ vs NTC DE COMPLETE\n"
)
