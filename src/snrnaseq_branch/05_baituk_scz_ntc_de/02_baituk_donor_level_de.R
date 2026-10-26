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
    "baituk_donor_pseudobulk_counts.csv.gz"
)

META_FILE <- file.path(
    IN_DIR,
    "baituk_donor_pseudobulk_metadata.csv"
)

OUT_ALL <- file.path(
    IN_DIR,
    "baituk_donor_level_SCZ_vs_NTC_DE_all.csv"
)

OUT_SIG <- file.path(
    IN_DIR,
    "baituk_donor_level_SCZ_vs_NTC_DE_FDR010.csv"
)


section <- function(x) {
    cat("\n", paste(rep("=", 110), collapse=""), "\n", sep="")
    cat(x, "\n")
    cat(paste(rep("=", 110), collapse=""), "\n")
}


section("BAITUK DONOR-LEVEL SCZ vs NTC DE")

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
cat("Donors:", length(unique(meta$donor_id)), "\n")

cat("\nDiagnosis:\n")
print(
    table(meta$diagnosis)
)


if (length(unique(meta$donor_id)) != 23) {
    stop("Expected 23 biological donors.")
}

if (length(unique(meta$donor_id[meta$diagnosis == "SCZ"])) != 9) {
    stop("Expected 9 SCZ donors.")
}

if (length(unique(meta$donor_id[meta$diagnosis == "NTC"])) != 14) {
    stop("Expected 14 NTC donors.")
}


meta$diagnosis <- factor(
    meta$diagnosis,
    levels=c("NTC", "SCZ")
)


# ==============================================================================================
# edgeR
# ==============================================================================================

section("EDGER FILTERING")

y <- DGEList(
    counts=as.matrix(counts),
    group=meta$diagnosis
)

design <- model.matrix(
    ~ diagnosis,
    data=meta
)

keep <- filterByExpr(
    y,
    design=design
)

cat("Genes before filter:", nrow(y), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")

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

res <- res[
    ,
    c(
        "gene",
        "logFC",
        "logCPM",
        "F",
        "PValue",
        "FDR"
    )
]


# Positive logFC = higher in SCZ
res$direction <- ifelse(
    res$logFC > 0,
    "SCZ_up",
    "SCZ_down"
)


write.csv(
    res,
    OUT_ALL,
    row.names=FALSE
)


sig <- res[
    !is.na(res$FDR) &
    res$FDR < 0.10,
    ,
    drop=FALSE
]

write.csv(
    sig,
    OUT_SIG,
    row.names=FALSE
)


section("RESULT")

cat("Genes tested :", nrow(res), "\n")
cat("FDR < 0.10 :", nrow(sig), "\n")
cat(
    "  SCZ up    :",
    sum(sig$logFC > 0),
    "\n"
)
cat(
    "  SCZ down  :",
    sum(sig$logFC < 0),
    "\n"
)

cat("\nTop 30 DE genes:\n")

print(
    head(res, 30),
    row.names=FALSE
)

cat("\nSaved all DE:\n")
cat(OUT_ALL, "\n")

cat("\nSaved FDR<0.10:\n")
cat(OUT_SIG, "\n")

cat(
    "\nFINAL STATUS: PASS — BAITUK DONOR-LEVEL SCZ vs NTC DE COMPLETE\n"
)
