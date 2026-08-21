#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT_RDS <- file.path(
    PROJECT_ROOT,
    "data",
    "raw",
    "Baituk",
    "annotations_final.RDS"
)

OUTPUT_DIR <- file.path(
    PROJECT_ROOT,
    "outputs",
    "snrnaseq_branch",
    "04_baituk_celltype_annotation"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


section <- function(title) {
    cat("\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
    cat(title, "\n")
    cat(paste(rep("=", 100), collapse = ""), "\n")
}


safe_filename <- function(x) {
    x <- gsub(
        "[^A-Za-z0-9._-]+",
        "_",
        x
    )

    x <- gsub(
        "_+",
        "_",
        x
    )

    x <- gsub(
        "^_|_$",
        "",
        x
    )

    if (nchar(x) == 0) {
        x <- "annotation"
    }

    x
}


score_candidate <- function(
    candidate_name,
    n_cells,
    n_labels
) {
    name_lower <- tolower(candidate_name)

    score <- 0

    if (grepl("cell.?type", name_lower)) {
        score <- score + 100
    }

    if (grepl("subtype", name_lower)) {
        score <- score + 90
    }

    if (grepl("annotation|annot", name_lower)) {
        score <- score + 80
    }

    if (grepl("class", name_lower)) {
        score <- score + 50
    }

    if (grepl("final", name_lower)) {
        score <- score + 30
    }

    if (grepl("cluster", name_lower)) {
        score <- score + 10
    }

    if (n_cells >= 1000) {
        score <- score + 20
    }

    if (n_labels >= 2 && n_labels <= 100) {
        score <- score + 20
    }

    if (n_labels > 200) {
        score <- score - 30
    }

    score
}


if (!file.exists(INPUT_RDS)) {
    stop(
        "Input annotation file was not found:\n",
        INPUT_RDS
    )
}

section("Reading Baituk annotations")

cat("Input:\n", INPUT_RDS, "\n")

annotations <- readRDS(INPUT_RDS)

cat("Class:", paste(class(annotations), collapse = ", "), "\n")
cat("Type:", typeof(annotations), "\n")

candidates <- list()


add_candidate <- function(
    candidate_name,
    labels,
    cell_ids = NULL
) {
    if (is.null(labels)) {
        return(invisible(NULL))
    }

    if (
        is.list(labels) &&
        !is.factor(labels)
    ) {
        return(invisible(NULL))
    }

    labels <- as.character(labels)

    if (length(labels) < 2) {
        return(invisible(NULL))
    }

    valid <- (
        !is.na(labels) &
        trimws(labels) != ""
    )

    labels <- trimws(labels[valid])

    if (!is.null(cell_ids)) {
        cell_ids <- as.character(cell_ids)
        cell_ids <- cell_ids[valid]
    }

    n_cells <- length(labels)
    n_labels <- length(unique(labels))

    # A cell-type column should normally contain repeated labels.
    # This avoids mistaking unique cell barcodes for annotations.
    if (n_labels < 2) {
        return(invisible(NULL))
    }

    if (n_labels > 500) {
        return(invisible(NULL))
    }

    if (
        n_cells > 100 &&
        (n_labels / n_cells) > 0.50
    ) {
        return(invisible(NULL))
    }

    if (is.null(cell_ids)) {
        cell_ids <- rep(NA_character_, n_cells)
    }

    key <- candidate_name

    if (key %in% names(candidates)) {
        suffix <- 2

        while (
            paste0(key, "_", suffix) %in%
            names(candidates)
        ) {
            suffix <- suffix + 1
        }

        key <- paste0(key, "_", suffix)
    }

    candidates[[key]] <<- data.frame(
        cell_id = cell_ids,
        cell_type = labels,
        stringsAsFactors = FALSE
    )

    invisible(NULL)
}


walk_object <- function(
    object,
    object_name = "annotations",
    depth = 0,
    max_depth = 5
) {
    if (depth > max_depth) {
        return(invisible(NULL))
    }

    if (is.data.frame(object)) {
        row_ids <- rownames(object)

        for (column_name in colnames(object)) {
            column <- object[[column_name]]

            if (
                is.factor(column) ||
                is.character(column) ||
                is.logical(column) ||
                is.numeric(column) ||
                is.integer(column)
            ) {
                add_candidate(
                    candidate_name = paste(
                        object_name,
                        column_name,
                        sep = "$"
                    ),
                    labels = column,
                    cell_ids = row_ids
                )
            }
        }

        return(invisible(NULL))
    }

    if (is.matrix(object)) {
        row_ids <- rownames(object)

        if (!is.null(colnames(object))) {
            for (
                column_index in seq_len(ncol(object))
            ) {
                add_candidate(
                    candidate_name = paste(
                        object_name,
                        colnames(object)[column_index],
                        sep = "$"
                    ),
                    labels = object[, column_index],
                    cell_ids = row_ids
                )
            }
        }

        return(invisible(NULL))
    }

    if (
        is.factor(object) ||
        (
            is.atomic(object) &&
            !is.null(object)
        )
    ) {
        add_candidate(
            candidate_name = object_name,
            labels = object,
            cell_ids = names(object)
        )

        return(invisible(NULL))
    }

    if (is.list(object)) {
        object_names <- names(object)

        if (is.null(object_names)) {
            object_names <- paste0(
                "entry_",
                seq_along(object)
            )
        }

        for (index in seq_along(object)) {
            walk_object(
                object = object[[index]],
                object_name = paste(
                    object_name,
                    object_names[index],
                    sep = "$"
                ),
                depth = depth + 1,
                max_depth = max_depth
            )
        }

        return(invisible(NULL))
    }

    if (isS4(object)) {
        for (slot_name in slotNames(object)) {
            walk_object(
                object = slot(object, slot_name),
                object_name = paste(
                    object_name,
                    slot_name,
                    sep = "@"
                ),
                depth = depth + 1,
                max_depth = max_depth
            )
        }
    }

    invisible(NULL)
}


section("Searching for annotation candidates")

walk_object(
    object = annotations,
    object_name = "annotations_final"
)

if (length(candidates) == 0) {
    stop(
        paste0(
            "No repeated annotation vectors were detected in ",
            "annotations_final.RDS.\n",
            "Review annotations_final_inventory.txt."
        )
    )
}

section("Calculating candidate proportions")

summary_records <- list()

for (candidate_name in names(candidates)) {
    annotation_df <- candidates[[candidate_name]]

    counts <- sort(
        table(
            annotation_df$cell_type,
            useNA = "no"
        ),
        decreasing = TRUE
    )

    proportion_df <- data.frame(
        cell_type = names(counts),
        n_cells = as.integer(counts),
        stringsAsFactors = FALSE
    )

    proportion_df$proportion <- (
        proportion_df$n_cells /
        sum(proportion_df$n_cells)
    )

    proportion_df$percent <- (
        100 * proportion_df$proportion
    )

    candidate_score <- score_candidate(
        candidate_name = candidate_name,
        n_cells = nrow(annotation_df),
        n_labels = nrow(proportion_df)
    )

    output_name <- safe_filename(candidate_name)

    candidate_annotation_file <- file.path(
        OUTPUT_DIR,
        paste0(
            "candidate_",
            output_name,
            "_cell_annotations.csv"
        )
    )

    candidate_proportion_file <- file.path(
        OUTPUT_DIR,
        paste0(
            "candidate_",
            output_name,
            "_proportions.csv"
        )
    )

    write.csv(
        annotation_df,
        candidate_annotation_file,
        row.names = FALSE
    )

    write.csv(
        proportion_df,
        candidate_proportion_file,
        row.names = FALSE
    )

    summary_records[[candidate_name]] <- data.frame(
        candidate_name = candidate_name,
        n_cells = nrow(annotation_df),
        n_cell_types = nrow(proportion_df),
        score = candidate_score,
        annotation_file = candidate_annotation_file,
        proportion_file = candidate_proportion_file,
        stringsAsFactors = FALSE
    )

    cat("\nCandidate:", candidate_name, "\n")
    cat("Cells:", nrow(annotation_df), "\n")
    cat("Labels:", nrow(proportion_df), "\n")
    cat("Score:", candidate_score, "\n")

    print(
        head(
            proportion_df,
            20
        ),
        row.names = FALSE
    )
}

candidate_summary <- do.call(
    rbind,
    summary_records
)

candidate_summary <- candidate_summary[
    order(
        -candidate_summary$score,
        -candidate_summary$n_cells,
        candidate_summary$n_cell_types
    ),
    ,
    drop = FALSE
]

summary_file <- file.path(
    OUTPUT_DIR,
    "baituk_annotation_candidate_summary.csv"
)

write.csv(
    candidate_summary,
    summary_file,
    row.names = FALSE
)


section("Selecting the most likely final annotation")

selected_name <- candidate_summary$candidate_name[1]
selected_annotation <- candidates[[selected_name]]

selected_counts <- sort(
    table(selected_annotation$cell_type),
    decreasing = TRUE
)

selected_proportion <- data.frame(
    cell_type = names(selected_counts),
    n_cells = as.integer(selected_counts),
    stringsAsFactors = FALSE
)

selected_proportion$proportion <- (
    selected_proportion$n_cells /
    sum(selected_proportion$n_cells)
)

selected_proportion$percent <- (
    100 * selected_proportion$proportion
)

selected_annotation_file <- file.path(
    OUTPUT_DIR,
    "baituk_selected_cell_annotations.csv"
)

selected_proportion_file <- file.path(
    OUTPUT_DIR,
    "baituk_selected_celltype_counts_proportions.csv"
)

write.csv(
    selected_annotation,
    selected_annotation_file,
    row.names = FALSE
)

write.csv(
    selected_proportion,
    selected_proportion_file,
    row.names = FALSE
)


section("Creating proportion plot")

plot_height <- max(
    1200,
    65 * nrow(selected_proportion)
)

plot_file <- file.path(
    OUTPUT_DIR,
    "baituk_selected_celltype_proportions.png"
)

png(
    filename = plot_file,
    width = 1800,
    height = plot_height,
    res = 180
)

par(
    mar = c(5, 16, 4, 2) + 0.1
)

plot_data <- selected_proportion[
    order(selected_proportion$percent),
    ,
    drop = FALSE
]

bar_positions <- barplot(
    plot_data$percent,
    names.arg = plot_data$cell_type,
    horiz = TRUE,
    las = 1,
    xlab = "Percentage of annotated nuclei",
    main = "Baituk snRNA-seq cell-type proportions",
    cex.names = 0.85
)

text(
    x = plot_data$percent,
    y = bar_positions,
    labels = sprintf(
        "%.2f%%  (n=%s)",
        plot_data$percent,
        format(
            plot_data$n_cells,
            big.mark = ","
        )
    ),
    pos = 4,
    cex = 0.75
)

dev.off()


section("FINAL RESULT")

cat("Automatically selected annotation:\n")
cat(selected_name, "\n")

cat("\nTotal annotated cells:\n")
cat(
    format(
        sum(selected_proportion$n_cells),
        big.mark = ","
    ),
    "\n"
)

cat("\nNumber of cell types/subtypes:\n")
cat(nrow(selected_proportion), "\n")

cat("\nCell-type counts and proportions:\n")
print(
    selected_proportion,
    row.names = FALSE
)

cat("\nSaved candidate summary:\n")
cat(summary_file, "\n")

cat("\nSaved selected cell annotations:\n")
cat(selected_annotation_file, "\n")

cat("\nSaved selected cell-type proportions:\n")
cat(selected_proportion_file, "\n")

cat("\nSaved plot:\n")
cat(plot_file, "\n")

cat("\nCompleted:", as.character(Sys.time()), "\n")
