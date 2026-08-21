#!/usr/bin/env Rscript

options(
    width = 200,
    stringsAsFactors = FALSE
)

PROJECT_ROOT <- "/beegfs/labs/hulab/projects/mjabin/GeneBridge"

INPUT <- file.path(
    PROJECT_ROOT,
    "data",
    "raw",
    "Baituk",
    "annotations_final.RDS"
)

OUT_DIR <- file.path(
    PROJECT_ROOT,
    "outputs",
    "celltype_comparison",
    "01_annotation_inventory",
    "baituk_published_rds"
)

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(INPUT)) {
    stop(
        paste(
            "File not found:",
            INPUT
        )
    )
}


cat(
    paste(rep("=", 100), collapse = ""),
    "\nBAITUK PUBLISHED ANNOTATION INSPECTION\n",
    paste(rep("=", 100), collapse = ""),
    "\n\n"
)

cat("Input:\n", INPUT, "\n\n")


obj <- readRDS(INPUT)


# =============================================================================
# Overall object structure
# =============================================================================

cat("Object class:\n")
print(class(obj))

cat("\nObject type:\n")
print(typeof(obj))

cat("\nTop-level elements:\n")
print(names(obj))


sink(
    file.path(
        OUT_DIR,
        "annotations_final_structure.txt"
    )
)

cat("INPUT\n")
cat(INPUT, "\n\n")

cat("CLASS\n")
print(class(obj))

cat("\nTYPEOF\n")
print(typeof(obj))

cat("\nTOP-LEVEL NAMES\n")
print(names(obj))

cat("\nSTRUCTURE\n")
str(
    obj,
    max.level = 3,
    list.len = 100
)

sink()


# =============================================================================
# Inspect each annotation vector
# =============================================================================

summary_list <- list()
all_counts <- list()


for (element_name in names(obj)) {

    x <- obj[[element_name]]

    cat(
        "\n",
        paste(rep("=", 100), collapse = ""),
        "\n"
    )

    cat(
        "ELEMENT:",
        element_name,
        "\n"
    )

    cat(
        paste(rep("-", 100), collapse = ""),
        "\n"
    )

    cat(
        "Class:",
        paste(
            class(x),
            collapse = ", "
        ),
        "\n"
    )

    cat(
        "Type:",
        typeof(x),
        "\n"
    )

    cat(
        "Length:",
        length(x),
        "\n"
    )

    cat(
        "Has names:",
        !is.null(names(x)),
        "\n"
    )

    if (!is.null(names(x))) {

        cat(
            "Number of names:",
            length(names(x)),
            "\n"
        )

        cat(
            "First 10 names:\n"
        )

        print(
            head(
                names(x),
                10
            )
        )
    }


    # -------------------------------------------------------------------------
    # Character annotation vectors
    # -------------------------------------------------------------------------

    if (is.atomic(x)) {

        values <- as.character(x)

        values[
            is.na(values)
        ] <- "<NA>"

        cat(
            "Unique labels:",
            length(
                unique(values)
            ),
            "\n"
        )

        cat(
            "\nFirst 20 values:\n"
        )

        print(
            head(
                values,
                20
            )
        )


        counts <- sort(
            table(values),
            decreasing = TRUE
        )

        counts_df <- data.frame(
            annotation_level = element_name,
            label = names(counts),
            n_cells = as.integer(counts),
            stringsAsFactors = FALSE
        )

        counts_df$proportion <- (
            counts_df$n_cells
            / sum(counts_df$n_cells)
        )

        counts_df$percent <- (
            100
            * counts_df$proportion
        )


        cat(
            "\nLABEL COUNTS:\n"
        )

        print(
            counts_df,
            row.names = FALSE
        )


        write.csv(
            counts_df,
            file.path(
                OUT_DIR,
                paste0(
                    "counts__",
                    element_name,
                    ".csv"
                )
            ),
            row.names = FALSE
        )


        # Save cell -> label mapping when vector is named.
        if (!is.null(names(x))) {

            mapping_df <- data.frame(
                cell_id = names(x),
                annotation_level = element_name,
                label = values,
                stringsAsFactors = FALSE
            )

            write.csv(
                mapping_df,
                file.path(
                    OUT_DIR,
                    paste0(
                        "cell_labels__",
                        element_name,
                        ".csv"
                    )
                ),
                row.names = FALSE
            )
        }


        summary_list[[element_name]] <- data.frame(
            annotation_level = element_name,
            n_entries = length(values),
            n_unique_labels = length(
                unique(values)
            ),
            n_missing = sum(
                values == "<NA>"
            ),
            has_cell_names = !is.null(
                names(x)
            ),
            stringsAsFactors = FALSE
        )


        all_counts[[element_name]] <- counts_df

    } else {

        cat(
            "\nNot an atomic vector; additional inspection required.\n"
        )
    }
}


# =============================================================================
# Combined outputs
# =============================================================================

if (length(summary_list) > 0) {

    summary_df <- do.call(
        rbind,
        summary_list
    )

    write.csv(
        summary_df,
        file.path(
            OUT_DIR,
            "annotation_levels_summary.csv"
        ),
        row.names = FALSE
    )

    cat(
        "\n",
        paste(rep("=", 100), collapse = ""),
        "\nANNOTATION LEVEL SUMMARY\n",
        paste(rep("=", 100), collapse = ""),
        "\n"
    )

    print(
        summary_df,
        row.names = FALSE
    )
}


if (length(all_counts) > 0) {

    combined_counts <- do.call(
        rbind,
        all_counts
    )

    write.csv(
        combined_counts,
        file.path(
            OUT_DIR,
            "all_annotation_level_counts.csv"
        ),
        row.names = FALSE
    )
}


# =============================================================================
# Compare whether vectors refer to identical cell sets
# =============================================================================

named_elements <- names(obj)[
    vapply(
        obj,
        function(x) {
            !is.null(names(x))
        },
        logical(1)
    )
]


if (length(named_elements) >= 2) {

    comparison <- matrix(
        NA,
        nrow = length(named_elements),
        ncol = length(named_elements),
        dimnames = list(
            named_elements,
            named_elements
        )
    )

    for (a in named_elements) {
        for (b in named_elements) {

            comparison[a, b] <- identical(
                names(obj[[a]]),
                names(obj[[b]])
            )
        }
    }

    cat(
        "\n",
        paste(rep("=", 100), collapse = ""),
        "\nIDENTICAL CELL NAME / ORDER CHECK\n",
        paste(rep("=", 100), collapse = ""),
        "\n"
    )

    print(comparison)

    write.csv(
        comparison,
        file.path(
            OUT_DIR,
            "annotation_cell_order_comparison.csv"
        )
    )
}


cat(
    "\n",
    paste(rep("=", 100), collapse = ""),
    "\nDONE\n",
    paste(rep("=", 100), collapse = ""),
    "\n"
)

cat(
    "Outputs:\n",
    OUT_DIR,
    "\n"
)
