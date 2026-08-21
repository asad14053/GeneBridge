#!/usr/bin/env Rscript

Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")

suppressPackageStartupMessages({
    library(conos)
})

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
    PROJECT_ROOT,
    "data",
    "raw",
    "Baituk",
    "snRNA-seq_Conos_object.RDS"
)

OUTPUT_DIR <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "04_baituk_celltype_annotation"
)

INVENTORY_TXT <- file.path(
    OUTPUT_DIR,
    "conos_object_inventory.txt"
)

CANDIDATE_SUMMARY_CSV <- file.path(
    OUTPUT_DIR,
    "conos_annotation_candidates.csv"
)

LONG_ANNOTATIONS_CSV <- file.path(
    OUTPUT_DIR,
    "conos_all_cluster_annotations_long.csv"
)


section <- function(title) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


safe_filename <- function(value) {
    value <- gsub(
        "[^A-Za-z0-9._-]+",
        "_",
        value
    )

    value <- gsub(
        "_+",
        "_",
        value
    )

    value
}


if (!file.exists(INPUT_RDS)) {
    stop(
        "Conos RDS was not found:\n",
        INPUT_RDS
    )
}

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

section("Reading Baituk Conos object")

cat("Input:\n", INPUT_RDS, "\n")

con <- readRDS(INPUT_RDS)

cat("\nObject class:\n")
print(class(con))

cat("\nAvailable R6 fields/methods:\n")
print(ls(con))

sink(INVENTORY_TXT)

section("Baituk Conos object inventory")

cat("Input file:\n")
cat(INPUT_RDS, "\n")

cat("\nObject class:\n")
print(class(con))

cat("\nR6 fields and methods:\n")
print(ls(con))

if (!is.null(con$clusters)) {
    cat("\nAvailable con$clusters entries:\n")
    print(names(con$clusters))
} else {
    cat("\ncon$clusters is NULL.\n")
}

if (!is.null(con$p2.list)) {
    cat("\nNumber of p2 objects:\n")
    print(length(con$p2.list))

    cat("\np2 sample names:\n")
    print(names(con$p2.list))
} else {
    cat("\ncon$p2.list is NULL.\n")
}

sink()

if (is.null(con$clusters)) {
    stop(
        "No con$clusters object was found. Review:\n",
        INVENTORY_TXT
    )
}

cluster_names <- names(con$clusters)

if (is.null(cluster_names) || length(cluster_names) == 0L) {
    stop(
        "con$clusters exists but contains no named entries."
    )
}

section("Extracting cluster and annotation vectors")

candidate_records <- list()
long_records <- list()

record_index <- 0L

for (cluster_name in cluster_names) {

    cat(
        "\nInspecting cluster entry: ",
        cluster_name,
        "\n",
        sep = ""
    )

    cluster_object <- con$clusters[[cluster_name]]

    if (is.null(cluster_object$groups)) {
        cat("  No $groups vector; skipping.\n")
        next
    }

    original_groups <- cluster_object$groups
    cell_ids <- names(original_groups)

    if (is.null(cell_ids)) {
        cat("  $groups has no cell names; skipping.\n")
        next
    }

    group_values <- as.character(original_groups)

    if (length(cell_ids) != length(group_values)) {
        stop(
            "Cell ID and group-vector lengths differ for ",
            cluster_name
        )
    }

    annotation_df <- data.frame(
        conos_cell_id = as.character(cell_ids),
        annotation = group_values,
        stringsAsFactors = FALSE
    )

    annotation_df <- annotation_df[
        !is.na(annotation_df$conos_cell_id) &
        annotation_df$conos_cell_id != "",
        ,
        drop = FALSE
    ]

    annotation_file <- file.path(
        OUTPUT_DIR,
        paste0(
            "conos_annotation_",
            safe_filename(cluster_name),
            ".csv"
        )
    )

    write.csv(
        annotation_df,
        annotation_file,
        row.names = FALSE
    )

    value_counts <- sort(
        table(
            annotation_df$annotation,
            useNA = "ifany"
        ),
        decreasing = TRUE
    )

    count_df <- data.frame(
        annotation = names(value_counts),
        n_cells = as.integer(value_counts),
        stringsAsFactors = FALSE
    )

    count_df$proportion <- (
        count_df$n_cells /
        sum(count_df$n_cells)
    )

    count_file <- file.path(
        OUTPUT_DIR,
        paste0(
            "conos_annotation_",
            safe_filename(cluster_name),
            "_counts.csv"
        )
    )

    write.csv(
        count_df,
        count_file,
        row.names = FALSE
    )

    likely_celltype_name <- grepl(
        "cell.?type|subtype|annotation|annot|class",
        cluster_name,
        ignore.case = TRUE
    )

    record_index <- record_index + 1L

    candidate_records[[record_index]] <- data.frame(
        annotation_name = cluster_name,
        n_cells = nrow(annotation_df),
        n_labels = length(unique(annotation_df$annotation)),
        likely_celltype_name = likely_celltype_name,
        annotation_file = annotation_file,
        count_file = count_file,
        stringsAsFactors = FALSE
    )

    long_records[[record_index]] <- data.frame(
        conos_cell_id = annotation_df$conos_cell_id,
        annotation_name = cluster_name,
        annotation = annotation_df$annotation,
        stringsAsFactors = FALSE
    )

    cat(
        "  Cells: ",
        nrow(annotation_df),
        "\n",
        sep = ""
    )

    cat(
        "  Unique labels: ",
        length(unique(annotation_df$annotation)),
        "\n",
        sep = ""
    )

    cat("  Most common labels:\n")
    print(head(count_df, 20))
}

if (length(candidate_records) == 0L) {
    stop(
        "No named Conos annotation vectors were extracted."
    )
}

candidate_summary <- do.call(
    rbind,
    candidate_records
)

candidate_summary <- candidate_summary[
    order(
        !candidate_summary$likely_celltype_name,
        candidate_summary$n_labels
    ),
    ,
    drop = FALSE
]

write.csv(
    candidate_summary,
    CANDIDATE_SUMMARY_CSV,
    row.names = FALSE
)

all_annotations <- do.call(
    rbind,
    long_records
)

write.csv(
    all_annotations,
    LONG_ANNOTATIONS_CSV,
    row.names = FALSE
)

section("DONE")

cat("Candidate annotation summary:\n")
cat(CANDIDATE_SUMMARY_CSV, "\n")

cat("\nCombined long annotation table:\n")
cat(LONG_ANNOTATIONS_CSV, "\n")

cat("\nObject inventory:\n")
cat(INVENTORY_TXT, "\n")

cat("\nCandidate annotations:\n")
print(candidate_summary)
