#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SpatialExperiment)
})

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PAPER_RDS <- file.path(
    ROOT,
    "data/reference/spatialDLPFC_SCZ_XENIUM/07_cell_type_de",
    "spe_pseudo_donor_celltype_clust_M0_lam0.1_k50_res0.7.rds"
)

BEFORE_CSV <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance/celltype_N23/pseudobulk",
    "original_xenium_300gene_donor_celltype_pseudobulk.csv.gz"
)

META23 <- file.path(
    ROOT,
    "data/metadata/xenium_DE_metadata_23.csv"
)

OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance/celltype_N23/paper_source_N23"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)


section <- function(x) {
    cat("\n")
    cat(paste(rep("=", 100), collapse=""), "\n")
    cat(x, "\n")
    cat(paste(rep("=", 100), collapse=""), "\n")
}


# =============================================================================
# LOAD PAPER RDS
# =============================================================================

section("LOAD EXACT PAPER PSEUDOBULK RDS")

paper <- readRDS(PAPER_RDS)

cat("Paper N24 dimensions:", dim(paper), "\n")
cat("Paper donors:", length(unique(paper$BrNum)), "\n")
cat("Paper cell types:", length(unique(paper$annots)), "\n")

if (nrow(paper) != 300) {
    stop("Paper RDS does not contain 300 genes.")
}

if (ncol(paper) != 288) {
    stop("Paper RDS does not contain 288 pseudobulks.")
}


# =============================================================================
# CANONICAL N23 DONORS
# =============================================================================

section("SUBSET EXACT PAPER OBJECT TO CANONICAL N23")

meta23 <- read.csv(
    META23,
    stringsAsFactors = FALSE
)

donors23 <- as.character(
    meta23$BrNum
)

if (length(unique(donors23)) != 23) {
    stop("Canonical metadata does not contain 23 unique donors.")
}

if ("Br6432" %in% donors23) {
    stop("Br6432 unexpectedly present in N23 metadata.")
}


keep <- as.character(paper$BrNum) %in% donors23

paper23 <- paper[, keep]


cat("N23 dimensions:", dim(paper23), "\n")
cat("N23 donors:", length(unique(paper23$BrNum)), "\n")
cat("N23 cell types:", length(unique(paper23$annots)), "\n")
cat("Br6432 excluded:", !"Br6432" %in% as.character(paper23$BrNum), "\n")

print(table(paper23$Dx))


if (ncol(paper23) != 276) {
    stop(
        paste0(
            "Expected 276 Paper-source N23 PBs; found ",
            ncol(paper23)
        )
    )
}


pb_per_donor <- table(
    table(
        as.character(paper23$BrNum)
    )
)

cat("\nPseudobulks per donor:\n")
print(pb_per_donor)


if (
    length(pb_per_donor) != 1 ||
    names(pb_per_donor)[1] != "12" ||
    as.integer(pb_per_donor[1]) != 23
) {
    stop("Expected exactly 12 pseudobulks for each of 23 donors.")
}


# =============================================================================
# PAPER GENE IDENTIFIERS
# =============================================================================

section("PAPER GENE IDENTIFIERS")

cat("rowData columns:\n")
print(colnames(rowData(paper23)))


if (
    "Symbol" %in% colnames(rowData(paper23)) &&
    length(unique(as.character(rowData(paper23)$Symbol))) == 300
) {

    paper_genes <- as.character(
        rowData(paper23)$Symbol
    )

    cat("Using rowData$Symbol as paper gene names.\n")

} else {

    paper_genes <- rownames(paper23)

    cat("Using rownames(paper23) as paper gene names.\n")
}


cat("Paper genes:", length(paper_genes), "\n")
cat("Unique paper genes:", length(unique(paper_genes)), "\n")


if (length(unique(paper_genes)) != 300) {
    stop("Paper gene identifiers are not unique.")
}


# =============================================================================
# PAPER N23 RAW PSEUDOBULK COUNTS
# =============================================================================

section("EXTRACT PAPER-SOURCE N23 COUNTS")

paper_counts <- as.matrix(
    counts(paper23)
)

rownames(paper_counts) <- paper_genes


# build canonical donor::celltype IDs
paper_pb_ids <- paste0(
    as.character(paper23$BrNum),
    "::",
    as.character(paper23$annots)
)

if (anyDuplicated(paper_pb_ids)) {
    stop("Duplicate Paper-source donor × cell-type IDs.")
}

colnames(paper_counts) <- paper_pb_ids


# transpose to PB × gene
paper_df <- as.data.frame(
    t(paper_counts),
    check.names = FALSE
)


cat("Paper-source N23 raw counts:", dim(paper_df), "\n")


# =============================================================================
# LOAD GENEBRIDGE-BEFORE
# =============================================================================

section("LOAD GENEBRIDGE-BEFORE N23")

before <- read.csv(
    BEFORE_CSV,
    row.names = 1,
    check.names = FALSE
)


cat("GeneBridge-before dimensions:", dim(before), "\n")


if (!all(dim(before) == c(276, 300))) {
    stop("GeneBridge-before expected 276 × 300.")
}


# =============================================================================
# ALIGN PSEUDOBULKS
# =============================================================================

section("ALIGN DONOR x CELL-TYPE PSEUDOBULKS")


paper_only_pb <- setdiff(
    rownames(paper_df),
    rownames(before)
)

before_only_pb <- setdiff(
    rownames(before),
    rownames(paper_df)
)


cat("Paper-only pseudobulks:", length(paper_only_pb), "\n")
cat("Before-only pseudobulks:", length(before_only_pb), "\n")


if (
    length(paper_only_pb) > 0 ||
    length(before_only_pb) > 0
) {

    cat("\nPaper-only examples:\n")
    print(head(paper_only_pb, 20))

    cat("\nBefore-only examples:\n")
    print(head(before_only_pb, 20))

    stop("Pseudobulk identities do not match.")
}


paper_df <- paper_df[
    rownames(before),
    ,
    drop = FALSE
]


# =============================================================================
# ALIGN GENES
# =============================================================================

section("ALIGN GENES")


paper_only_gene <- setdiff(
    colnames(paper_df),
    colnames(before)
)

before_only_gene <- setdiff(
    colnames(before),
    colnames(paper_df)
)


cat("Paper-only genes:", length(paper_only_gene), "\n")
cat("Before-only genes:", length(before_only_gene), "\n")


if (
    length(paper_only_gene) > 0 ||
    length(before_only_gene) > 0
) {

    cat("\nPaper-only genes:\n")
    print(paper_only_gene)

    cat("\nBefore-only genes:\n")
    print(before_only_gene)

    stop("Paper and GeneBridge gene sets differ.")
}


paper_df <- paper_df[
    ,
    colnames(before),
    drop = FALSE
]


# =============================================================================
# COUNT-LEVEL COMPARISON
# =============================================================================

section("PAPER-SOURCE N23 vs GENEBRIDGE-BEFORE N23 COUNTS")


P <- as.matrix(paper_df)
B <- as.matrix(before)

storage.mode(P) <- "double"
storage.mode(B) <- "double"


D <- B - P


exact <- identical(
    unname(P),
    unname(B)
)

allclose <- isTRUE(
    all.equal(
        unname(P),
        unname(B),
        tolerance = 0,
        check.attributes = FALSE
    )
)


n_diff <- sum(
    D != 0,
    na.rm = TRUE
)

max_abs_diff <- max(
    abs(D),
    na.rm = TRUE
)

total_abs_diff <- sum(
    abs(D),
    na.rm = TRUE
)


cat("Exact equality:", exact, "\n")
cat("Allclose tolerance=0:", allclose, "\n")
cat("Differing PB × gene entries:", n_diff, "\n")
cat("Maximum absolute count difference:", max_abs_diff, "\n")
cat("Total absolute count difference:", total_abs_diff, "\n")


# =============================================================================
# LOCALIZE ANY DIFFERENCES
# =============================================================================

section("LOCALIZE DIFFERENCES")


pb_abs_diff <- rowSums(
    abs(D)
)

gene_abs_diff <- colSums(
    abs(D)
)


pb_diff_df <- data.frame(
    pseudobulk_id = names(pb_abs_diff),
    total_abs_diff = as.numeric(pb_abs_diff),
    stringsAsFactors = FALSE
)

pb_diff_df <- pb_diff_df[
    order(
        -pb_diff_df$total_abs_diff
    ),
    ,
    drop = FALSE
]


gene_diff_df <- data.frame(
    gene = names(gene_abs_diff),
    total_abs_diff = as.numeric(gene_abs_diff),
    stringsAsFactors = FALSE
)

gene_diff_df <- gene_diff_df[
    order(
        -gene_diff_df$total_abs_diff
    ),
    ,
    drop = FALSE
]


cat("Pseudobulks with any count difference:",
    sum(pb_diff_df$total_abs_diff > 0),
    "\n"
)

cat("Genes with any count difference:",
    sum(gene_diff_df$total_abs_diff > 0),
    "\n"
)


cat("\nTop differing pseudobulks:\n")

print(
    head(
        pb_diff_df[
            pb_diff_df$total_abs_diff > 0,
            ,
            drop = FALSE
        ],
        20
    ),
    row.names = FALSE
)


cat("\nTop differing genes:\n")

print(
    head(
        gene_diff_df[
            gene_diff_df$total_abs_diff > 0,
            ,
            drop = FALSE
        ],
        20
    ),
    row.names = FALSE
)


# =============================================================================
# SAVE PAPER-SOURCE N23
# =============================================================================

section("SAVE PAPER-SOURCE N23")


paper_counts_out <- file.path(
    OUT,
    "paper_source_N23_300gene_donor_celltype_pseudobulk_counts.csv.gz"
)

paper_meta_out <- file.path(
    OUT,
    "paper_source_N23_donor_celltype_metadata.csv"
)

pb_diff_out <- file.path(
    OUT,
    "paper_source_vs_before_pseudobulk_count_differences.csv"
)

gene_diff_out <- file.path(
    OUT,
    "paper_source_vs_before_gene_count_differences.csv"
)


write.csv(
    paper_df,
    gzfile(paper_counts_out),
    quote = TRUE
)


paper_meta <- as.data.frame(
    colData(paper23)
)

paper_meta$pseudobulk_id <- paper_pb_ids

paper_meta <- paper_meta[
    match(
        rownames(before),
        paper_meta$pseudobulk_id
    ),
    ,
    drop = FALSE
]


write.csv(
    paper_meta,
    paper_meta_out,
    row.names = FALSE
)


write.csv(
    pb_diff_df,
    pb_diff_out,
    row.names = FALSE
)


write.csv(
    gene_diff_df,
    gene_diff_out,
    row.names = FALSE
)


# Save untouched N23 SCE subset as well.
saveRDS(
    paper23,
    file.path(
        OUT,
        "paper_source_N23_spe_pseudo_donor_celltype.rds"
    )
)


cat("Saved:\n")
cat(paper_counts_out, "\n")
cat(paper_meta_out, "\n")


# =============================================================================
# FINAL
# =============================================================================

section("FINAL STATUS")

cat("Paper official source: N24\n")
cat("Paper-source comparison cohort: N23\n")
cat("Donors: 23\n")
cat("Cell types: 12\n")
cat("Pseudobulks: 276\n")
cat("Genes: 300\n")

cat(
    "Paper-source vs GeneBridge-before exact:",
    exact,
    "\n"
)

cat(
    "Differing matrix entries:",
    n_diff,
    "\n"
)

cat(
    "Maximum absolute difference:",
    max_abs_diff,
    "\n"
)

cat(
    "\nFINAL STATUS: PAPER-SOURCE N23 CELL-TYPE AUDIT PASS\n"
)
