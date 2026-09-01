#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(edgeR)
    library(limma)
    library(spatialLIBD)
})

ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

# ==========================================================================================
# INPUTS
# ==========================================================================================

PAPER_RDS <- file.path(
    ROOT,
    "data/reference/paper_xenium",
    "spe_pseudo_donor_domain_spaTransfer_k50_smoothed_predictions.rds"
)

ORIG_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "original_xenium_300gene_donor_layer_pseudobulk_counts.csv.gz"
)

ORIG_META <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "original_layer_pseudobulk",
    "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)

ENVI_MEASURED_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_measured300_donor_layer_pseudobulk.csv.gz"
)

ENVI_FULL_COUNTS <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_full34987_donor_layer_pseudobulk_countscale.csv.gz"
)

ENVI_META <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_donor_layer_pseudobulk_metadata.csv"
)

GENE_INFO <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "envi_layer_pseudobulk",
    "ENVI_full34987_gene_info.csv"
)

OUT <- file.path(
    ROOT,
    "outputs/imputation_full/DE/paper_concordance",
    "domain_stratified_N23_DE"
)

dir.create(
    OUT,
    recursive = TRUE,
    showWarnings = FALSE
)

section <- function(x) {
    cat(
        "\n",
        paste(rep("=", 110), collapse = ""),
        "\n",
        x,
        "\n",
        paste(rep("=", 110), collapse = ""),
        "\n",
        sep = ""
    )
}

# ==========================================================================================
# CHECK INPUTS
# ==========================================================================================

section("INPUT CHECK")

inputs <- c(
    PAPER_RDS,
    ORIG_COUNTS,
    ORIG_META,
    ENVI_MEASURED_COUNTS,
    ENVI_FULL_COUNTS,
    ENVI_META,
    GENE_INFO
)

for (f in inputs) {
    if (!file.exists(f)) {
        stop("Missing input:\n", f)
    }

    cat("FOUND: ", f, "\n", sep = "")
}

# ==========================================================================================
# HELPERS
# ==========================================================================================

make_key <- function(meta) {

    paste(
        as.character(meta$BrNum),
        as.character(meta$predictions_smooth),
        sep = "::"
    )
}


align_meta <- function(meta, pb) {

    if (
        "pseudobulk_id" %in% colnames(meta) &&
        all(rownames(pb) %in% meta$pseudobulk_id)
    ) {

        idx <- match(
            rownames(pb),
            meta$pseudobulk_id
        )

    } else {

        keys <- make_key(meta)

        if (!all(rownames(pb) %in% keys)) {
            stop("Unable to align metadata to pseudobulk matrix.")
        }

        idx <- match(
            rownames(pb),
            keys
        )
    }

    meta <- meta[
        idx,
        ,
        drop = FALSE
    ]

    rownames(meta) <- rownames(pb)

    meta
}


prepare_variables <- function(sce) {

    sce$Dx <- factor(
        as.character(sce$Dx),
        levels = c("NTC", "SCZ")
    )

    sce$predictions_smooth <- factor(
        as.character(sce$predictions_smooth)
    )

    sce$Age <- as.numeric(
        as.character(sce$Age)
    )

    sce$Sex <- factor(
        as.character(sce$Sex)
    )

    sce$slide_id <- factor(
        as.character(sce$slide_id)
    )

    sce$BrNum <- factor(
        as.character(sce$BrNum)
    )

    sce
}


prepare_global_expression <- function(
    sce,
    label
) {

    section(
        paste0(
            "GLOBAL FILTER/NORMALIZATION: ",
            label
        )
    )

    sce <- prepare_variables(sce)

    keep <- edgeR::filterByExpr(
        sce,
        group = sce$predictions_smooth
    )

    cat(
        "filterByExpr retained:",
        sum(keep),
        "/",
        length(keep),
        "\n"
    )

    sce <- sce[
        keep,
        ,
        drop = FALSE
    ]

    logcounts(sce) <- edgeR::cpm(
        edgeR::calcNormFactors(sce),
        log = TRUE,
        prior.count = 1
    )

    cat(
        "Normalized:",
        nrow(sce),
        "genes x",
        ncol(sce),
        "pseudobulks\n"
    )

    sce
}


run_one_domain <- function(
    sce,
    domain,
    dataset_name
) {

    keep_cols <- (
        as.character(sce$predictions_smooth)
        ==
        domain
    )

    x <- sce[
        ,
        keep_cols,
        drop = FALSE
    ]

    # Drop unused factor levels after taking one spatial domain.
    x$Dx <- droplevels(x$Dx)
    x$Sex <- droplevels(x$Sex)
    x$slide_id <- droplevels(x$slide_id)
    x$BrNum <- droplevels(x$BrNum)
    x$predictions_smooth <- droplevels(
        x$predictions_smooth
    )

    if (ncol(x) != 23) {
        stop(
            dataset_name,
            " ",
            domain,
            ": expected 23 donor pseudobulks, found ",
            ncol(x)
        )
    }

    if (
        length(unique(as.character(x$BrNum)))
        != 23
    ) {
        stop(
            dataset_name,
            " ",
            domain,
            ": expected 23 unique donors."
        )
    }

    covars <- c(
        "Age",
        "Sex",
        "slide_id"
    )

    mod <- spatialLIBD::registration_model(
        x,
        covars = covars,
        var_registration = "Dx"
    )

    design_rank <- qr(mod)$rank
    design_cols <- ncol(mod)

    cat(
        dataset_name,
        " | ",
        domain,
        " | samples=",
        ncol(x),
        " | genes=",
        nrow(x),
        " | design rank=",
        design_rank,
        "/",
        design_cols,
        "\n",
        sep = ""
    )

    if (design_rank != design_cols) {

        stop(
            dataset_name,
            " ",
            domain,
            ": design is rank deficient."
        )
    }

    # Exact paper-style layer/domain-stratified analysis:
    # one PB per donor within this domain, therefore no duplicateCorrelation block.
    res <- spatialLIBD::registration_stats_enrichment(
        x,
        block_cor = NaN,
        covars = covars,
        var_registration = "Dx",
        var_sample_id = "BrNum",
        gene_ensembl = "ID",
        gene_name = "Symbol"
    )

    res$domain <- domain
    res$dataset <- dataset_name

    list(
        result = res,
        design_rank = design_rank,
        design_cols = design_cols
    )
}


summarize_result <- function(
    res,
    dataset,
    domain
) {

    data.frame(
        dataset = dataset,
        domain = domain,
        genes_tested = nrow(res),
        nominal_P05 = sum(
            res$p_value_SCZ < 0.05,
            na.rm = TRUE
        ),
        FDR10 = sum(
            res$fdr_SCZ < 0.10,
            na.rm = TRUE
        ),
        FDR05 = sum(
            res$fdr_SCZ < 0.05,
            na.rm = TRUE
        ),
        stringsAsFactors = FALSE
    )
}

# ==========================================================================================
# LOAD ORIGINAL 300
# ==========================================================================================

section("LOAD GENEBRIDGE BEFORE 300")

orig_pb <- read.csv(
    ORIG_COUNTS,
    row.names = 1,
    check.names = FALSE
)

orig_pb <- as.matrix(orig_pb)
storage.mode(orig_pb) <- "double"

orig_meta <- read.csv(
    ORIG_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

orig_meta <- align_meta(
    orig_meta,
    orig_pb
)

cat(
    "Original:",
    nrow(orig_pb),
    "PB x",
    ncol(orig_pb),
    "genes\n"
)

if (
    nrow(orig_pb) != 161 ||
    ncol(orig_pb) != 300
) {
    stop("Original matrix must be 161 x 300.")
}

canonical_keys <- make_key(orig_meta)

if (length(unique(orig_meta$BrNum)) != 23) {
    stop("Original must contain 23 donors.")
}

if ("Br6432" %in% orig_meta$BrNum) {
    stop("Br6432 must be excluded.")
}

# ==========================================================================================
# LOAD PAPER SOURCE
# ==========================================================================================

section("LOAD PAPER-SOURCE PSEUDOBULK")

paper <- readRDS(PAPER_RDS)

paper_meta <- as.data.frame(
    colData(paper)
)

paper_keys <- make_key(paper_meta)

idx <- match(
    canonical_keys,
    paper_keys
)

if (any(is.na(idx))) {
    stop("Paper-source pseudobulks do not contain all canonical N23 keys.")
}

paper <- paper[
    ,
    idx,
    drop = FALSE
]

paper_meta <- as.data.frame(
    colData(paper)
)

cat(
    "Paper-source N23:",
    nrow(paper),
    "genes x",
    ncol(paper),
    "PB\n"
)

if (
    nrow(paper) != 300 ||
    ncol(paper) != 161
) {
    stop("Paper-source N23 must be 300 x 161.")
}

# ==========================================================================================
# ALIGN 300-GENE UNIVERSE
# ==========================================================================================

section("ALIGN 300 GENES")

paper_symbols <- if (
    "Symbol" %in% colnames(rowData(paper))
) {
    as.character(rowData(paper)$Symbol)
} else {
    rownames(paper)
}

gene_match <- match(
    colnames(orig_pb),
    paper_symbols
)

if (any(is.na(gene_match))) {

    gene_match <- match(
        colnames(orig_pb),
        rownames(paper)
    )
}

if (any(is.na(gene_match))) {
    stop("Could not align original 300 genes to paper genes.")
}

paper <- paper[
    gene_match,
    ,
    drop = FALSE
]

gene_symbols <- colnames(orig_pb)

gene_ids <- if (
    "ID" %in% colnames(rowData(paper))
) {
    as.character(rowData(paper)$ID)
} else {
    gene_symbols
}

# ==========================================================================================
# BUILD ORIGINAL SCE
# ==========================================================================================

orig_sce <- SingleCellExperiment(
    assays = list(
        counts = t(orig_pb)
    ),
    colData = DataFrame(orig_meta),
    rowData = DataFrame(
        ID = gene_ids,
        Symbol = gene_symbols
    )
)

rownames(orig_sce) <- gene_symbols

# ==========================================================================================
# LOAD ENVI MEASURED 300
# ==========================================================================================

section("LOAD ENVI MEASURED 300")

envi300_pb <- read.csv(
    ENVI_MEASURED_COUNTS,
    row.names = 1,
    check.names = FALSE
)

envi300_pb <- as.matrix(envi300_pb)
storage.mode(envi300_pb) <- "double"

envi_meta <- read.csv(
    ENVI_META,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

envi_meta_300 <- align_meta(
    envi_meta,
    envi300_pb
)

envi300_keys <- make_key(
    envi_meta_300
)

idx <- match(
    canonical_keys,
    envi300_keys
)

if (any(is.na(idx))) {
    stop("Could not align ENVI measured300 to canonical N23.")
}

envi300_pb <- envi300_pb[
    idx,
    ,
    drop = FALSE
]

envi_meta_300 <- envi_meta_300[
    idx,
    ,
    drop = FALSE
]

gene_idx <- match(
    gene_symbols,
    colnames(envi300_pb)
)

if (any(is.na(gene_idx))) {
    stop("ENVI measured300 gene mismatch.")
}

envi300_pb <- envi300_pb[
    ,
    gene_idx,
    drop = FALSE
]

cat(
    "Original vs ENVI measured300 all values equal:",
    isTRUE(
        all.equal(
            orig_pb,
            envi300_pb,
            check.attributes = FALSE
        )
    ),
    "\n"
)

envi300_sce <- SingleCellExperiment(
    assays = list(
        counts = t(envi300_pb)
    ),
    colData = DataFrame(envi_meta_300),
    rowData = DataFrame(
        ID = gene_ids,
        Symbol = gene_symbols
    )
)

rownames(envi300_sce) <- gene_symbols

# ==========================================================================================
# PREPARE PAPER N23 300
# ==========================================================================================

section("PREPARE PAPER-SOURCE N23")

paper_sce <- SingleCellExperiment(
    assays = list(
        counts = counts(paper)
    ),
    colData = colData(paper),
    rowData = DataFrame(
        ID = gene_ids,
        Symbol = gene_symbols
    )
)

rownames(paper_sce) <- gene_symbols

# ==========================================================================================
# LOAD FULL ENVI
# ==========================================================================================

section("LOAD ENVI FULL 34,987")

full_pb <- read.csv(
    ENVI_FULL_COUNTS,
    row.names = 1,
    check.names = FALSE
)

full_pb <- as.matrix(full_pb)
storage.mode(full_pb) <- "double"

full_meta <- align_meta(
    envi_meta,
    full_pb
)

full_keys <- make_key(full_meta)

idx <- match(
    canonical_keys,
    full_keys
)

if (any(is.na(idx))) {
    stop("Could not align full ENVI to canonical N23.")
}

full_pb <- full_pb[
    idx,
    ,
    drop = FALSE
]

full_meta <- full_meta[
    idx,
    ,
    drop = FALSE
]

gene_info <- read.csv(
    GENE_INFO,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

gidx <- match(
    colnames(full_pb),
    gene_info$gene
)

if (any(is.na(gidx))) {
    stop("Could not align full ENVI gene metadata.")
}

gene_info <- gene_info[
    gidx,
    ,
    drop = FALSE
]

full_id <- if (
    "ensembl" %in% colnames(gene_info)
) {
    as.character(gene_info$ensembl)
} else if (
    "ID" %in% colnames(gene_info)
) {
    as.character(gene_info$ID)
} else {
    as.character(gene_info$gene)
}

full_sce <- SingleCellExperiment(
    assays = list(
        counts = t(full_pb)
    ),
    colData = DataFrame(full_meta),
    rowData = DataFrame(
        ID = full_id,
        Symbol = as.character(gene_info$gene),
        expression_source = as.character(
            gene_info$expression_source
        )
    )
)

rownames(full_sce) <- gene_info$gene

cat(
    "Full ENVI:",
    nrow(full_sce),
    "genes x",
    ncol(full_sce),
    "PB\n"
)

print(
    table(
        rowData(full_sce)$expression_source
    )
)

# ==========================================================================================
# GLOBAL FILTERING/NORMALIZATION
# ==========================================================================================

paper_sce <- prepare_global_expression(
    paper_sce,
    "Paper-source N23 300"
)

orig_sce <- prepare_global_expression(
    orig_sce,
    "GeneBridge-before N23 300"
)

envi300_sce <- prepare_global_expression(
    envi300_sce,
    "ENVI measured300 control N23"
)

full_sce <- prepare_global_expression(
    full_sce,
    "ENVI full N23"
)

cat(
    "\nENVI full genes after filtering:\n"
)

print(
    table(
        rowData(full_sce)$expression_source
    )
)

# ==========================================================================================
# DOMAIN-STRATIFIED DE
# ==========================================================================================

section("DOMAIN-STRATIFIED DE")

domains <- sprintf(
    "spd%02d",
    1:7
)

domain_names <- c(
    spd01 = "WMtz",
    spd02 = "L3/4",
    spd03 = "L6",
    spd04 = "WM",
    spd05 = "L5",
    spd06 = "L2/3",
    spd07 = "L1/M"
)

summary_list <- list()
concordance_list <- list()

for (domain in domains) {

    section(
        paste0(
            "DOMAIN: ",
            domain,
            " / ",
            domain_names[[domain]]
        )
    )

    domain_dir <- file.path(
        OUT,
        domain
    )

    dir.create(
        domain_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    # --------------------------------------------------------------------------------------
    # PAPER-SOURCE N23
    # --------------------------------------------------------------------------------------

    paper_run <- run_one_domain(
        paper_sce,
        domain,
        "Paper_source_N23"
    )

    paper_res <- paper_run$result

    write.csv(
        paper_res,
        file.path(
            domain_dir,
            "Paper_source_N23_300gene_DE.csv"
        ),
        row.names = FALSE
    )

    summary_list[[length(summary_list) + 1]] <-
        summarize_result(
            paper_res,
            "Paper_source_N23_300",
            domain
        )

    # --------------------------------------------------------------------------------------
    # GENEBRIDGE BEFORE
    # --------------------------------------------------------------------------------------

    orig_run <- run_one_domain(
        orig_sce,
        domain,
        "GeneBridge_before_N23"
    )

    orig_res <- orig_run$result

    write.csv(
        orig_res,
        file.path(
            domain_dir,
            "GeneBridge_before_N23_300gene_DE.csv"
        ),
        row.names = FALSE
    )

    summary_list[[length(summary_list) + 1]] <-
        summarize_result(
            orig_res,
            "GeneBridge_before_N23_300",
            domain
        )

    # --------------------------------------------------------------------------------------
    # ENVI MEASURED 300 CONTROL
    # --------------------------------------------------------------------------------------

    envi300_run <- run_one_domain(
        envi300_sce,
        domain,
        "ENVI_measured300_control_N23"
    )

    envi300_res <- envi300_run$result

    write.csv(
        envi300_res,
        file.path(
            domain_dir,
            "ENVI_measured300_control_N23_DE.csv"
        ),
        row.names = FALSE
    )

    summary_list[[length(summary_list) + 1]] <-
        summarize_result(
            envi300_res,
            "ENVI_measured300_control_N23",
            domain
        )

    # --------------------------------------------------------------------------------------
    # FULL ENVI
    # --------------------------------------------------------------------------------------

    full_run <- run_one_domain(
        full_sce,
        domain,
        "ENVI_full_N23"
    )

    full_res <- full_run$result

    source_map <- data.frame(
        gene = rownames(full_sce),
        expression_source = as.character(
            rowData(full_sce)$expression_source
        ),
        stringsAsFactors = FALSE
    )

    full_res <- merge(
        full_res,
        source_map,
        by = "gene",
        all.x = TRUE,
        sort = FALSE
    )

    imputed_res <- full_res[
        full_res$expression_source == "envi_imputed",
        ,
        drop = FALSE
    ]

    measured_full_res <- full_res[
        full_res$expression_source == "measured_xenium",
        ,
        drop = FALSE
    ]

    write.csv(
        full_res,
        file.path(
            domain_dir,
            "ENVI_full_N23_domain_stratified_DE.csv"
        ),
        row.names = FALSE
    )

    write.csv(
        imputed_res,
        file.path(
            domain_dir,
            "ENVI_imputed_N23_domain_stratified_DE.csv"
        ),
        row.names = FALSE
    )

    write.csv(
        imputed_res[
            imputed_res$fdr_SCZ < 0.10,
            ,
            drop = FALSE
        ],
        file.path(
            domain_dir,
            "ENVI_imputed_N23_FDR10.csv"
        ),
        row.names = FALSE
    )

    write.csv(
        imputed_res[
            imputed_res$fdr_SCZ < 0.05,
            ,
            drop = FALSE
        ],
        file.path(
            domain_dir,
            "ENVI_imputed_N23_FDR05.csv"
        ),
        row.names = FALSE
    )

    summary_list[[length(summary_list) + 1]] <-
        summarize_result(
            full_res,
            "ENVI_full_N23",
            domain
        )

    summary_list[[length(summary_list) + 1]] <-
        summarize_result(
            measured_full_res,
            "ENVI_full_measured_subset",
            domain
        )

    summary_list[[length(summary_list) + 1]] <-
        summarize_result(
            imputed_res,
            "ENVI_imputed_N23",
            domain
        )

    # --------------------------------------------------------------------------------------
    # 300-GENE CONCORDANCE
    # --------------------------------------------------------------------------------------

    p <- paper_res[
        ,
        c(
            "gene",
            "logFC_SCZ",
            "t_stat_SCZ",
            "fdr_SCZ"
        )
    ]

    names(p)[-1] <- paste0(
        names(p)[-1],
        "_paper"
    )

    g <- orig_res[
        ,
        c(
            "gene",
            "logFC_SCZ",
            "t_stat_SCZ",
            "fdr_SCZ"
        )
    ]

    names(g)[-1] <- paste0(
        names(g)[-1],
        "_genebridge"
    )

    e <- envi300_res[
        ,
        c(
            "gene",
            "logFC_SCZ",
            "t_stat_SCZ",
            "fdr_SCZ"
        )
    ]

    names(e)[-1] <- paste0(
        names(e)[-1],
        "_envi_measured"
    )

    merged <- Reduce(
        function(x, y) merge(
            x,
            y,
            by = "gene",
            all = FALSE
        ),
        list(p, g, e)
    )

    metrics <- data.frame(
        domain = domain,
        layer = domain_names[[domain]],

        paper_vs_genebridge_logFC_Pearson =
            cor(
                merged$logFC_SCZ_paper,
                merged$logFC_SCZ_genebridge,
                method = "pearson"
            ),

        paper_vs_genebridge_logFC_Spearman =
            cor(
                merged$logFC_SCZ_paper,
                merged$logFC_SCZ_genebridge,
                method = "spearman"
            ),

        paper_vs_genebridge_t_Pearson =
            cor(
                merged$t_stat_SCZ_paper,
                merged$t_stat_SCZ_genebridge,
                method = "pearson"
            ),

        paper_vs_envi_logFC_Pearson =
            cor(
                merged$logFC_SCZ_paper,
                merged$logFC_SCZ_envi_measured,
                method = "pearson"
            ),

        paper_vs_envi_t_Pearson =
            cor(
                merged$t_stat_SCZ_paper,
                merged$t_stat_SCZ_envi_measured,
                method = "pearson"
            ),

        genebridge_vs_envi_logFC_Pearson =
            cor(
                merged$logFC_SCZ_genebridge,
                merged$logFC_SCZ_envi_measured,
                method = "pearson"
            ),

        stringsAsFactors = FALSE
    )

    concordance_list[[length(concordance_list) + 1]] <-
        metrics
}

# ==========================================================================================
# COMBINE SUMMARY
# ==========================================================================================

section("COMBINE RESULTS")

summary_df <- do.call(
    rbind,
    summary_list
)

summary_df$layer <- domain_names[
    summary_df$domain
]

summary_df <- summary_df[
    ,
    c(
        "dataset",
        "domain",
        "layer",
        "genes_tested",
        "nominal_P05",
        "FDR10",
        "FDR05"
    )
]

concordance_df <- do.call(
    rbind,
    concordance_list
)

write.csv(
    summary_df,
    file.path(
        OUT,
        "domain_stratified_DE_summary.csv"
    ),
    row.names = FALSE
)

write.csv(
    concordance_df,
    file.path(
        OUT,
        "domain_stratified_300gene_concordance.csv"
    ),
    row.names = FALSE
)

# ==========================================================================================
# PRINT ENVI IMPUTED SUMMARY
# ==========================================================================================

section("ENVI IMPUTED DOMAIN SUMMARY")

imputed_summary <- summary_df[
    summary_df$dataset == "ENVI_imputed_N23",
    ,
    drop = FALSE
]

print(
    imputed_summary,
    row.names = FALSE
)

# ==========================================================================================
# BAR PLOT
# ==========================================================================================

make_barplot <- function() {

    x <- imputed_summary

    vals <- rbind(
        x$FDR10,
        x$FDR05
    )

    colnames(vals) <- x$layer

    barplot(
        vals,
        beside = TRUE,
        names.arg = x$layer,
        xlab = "Spatial domain",
        ylab = "Number of DE genes",
        main = paste0(
            "ENVI-imputed SCZ DE by spatial domain\n",
            "N23 domain-stratified analysis"
        ),
        legend.text = c(
            "FDR < 0.10",
            "FDR < 0.05"
        )
    )
}

pdf(
    file.path(
        OUT,
        "07A_ENVI_imputed_domain_stratified_DEG_counts.pdf"
    ),
    width = 9,
    height = 6
)

make_barplot()

dev.off()


png(
    file.path(
        OUT,
        "07A_ENVI_imputed_domain_stratified_DEG_counts.png"
    ),
    width = 2700,
    height = 1800,
    res = 300
)

make_barplot()

dev.off()

# ==========================================================================================
# FINAL
# ==========================================================================================

section("FINAL STATUS")

cat(
    "Donors                  : 23\n",
    "Spatial domains         : 7\n",
    "Samples per domain      : 23\n",
    "Model within each domain: Dx + Age + Sex + slide_id\n",
    "Donor block correlation : NaN (one PB per donor/domain)\n",
    sep = ""
)

cat(
    "\nENVI imputed domain-stratified summary:\n"
)

print(
    imputed_summary,
    row.names = FALSE
)

cat(
    "\n300-gene concordance:\n"
)

print(
    concordance_df,
    row.names = FALSE
)

cat(
    "\nFINAL STATUS: N23 DOMAIN-STRATIFIED DE COMPLETE\n"
)
