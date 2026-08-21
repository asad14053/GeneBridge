#!/usr/bin/env Rscript

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

REPORT_FILE <- file.path(
    OUTPUT_DIR,
    "annotations_final_inventory.txt"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(INPUT_RDS)) {
    stop("File not found:\n", INPUT_RDS)
}

annotations <- readRDS(INPUT_RDS)

sink(REPORT_FILE)

cat("================================================================================\n")
cat("BAITUK annotations_final.RDS INVENTORY\n")
cat("================================================================================\n")

cat("\nFile:\n")
cat(INPUT_RDS, "\n")

cat("\nObject class:\n")
print(class(annotations))

cat("\nObject type:\n")
print(typeof(annotations))

cat("\nDimensions:\n")
print(dim(annotations))

cat("\nLength:\n")
print(length(annotations))

cat("\nNames:\n")
print(names(annotations))

cat("\nColumn names, when available:\n")
print(colnames(annotations))

cat("\nRow names preview, when available:\n")
if (!is.null(rownames(annotations))) {
    print(head(rownames(annotations), 20))
} else {
    print(NULL)
}

cat("\nObject structure:\n")
str(
    annotations,
    max.level = 4,
    list.len = 50
)

cat("\nObject preview:\n")

if (is.data.frame(annotations) || is.matrix(annotations)) {
    print(head(annotations, 20))
} else if (is.atomic(annotations)) {
    print(head(annotations, 20))
    cat("\nVector names preview:\n")
    print(head(names(annotations), 20))
} else if (is.list(annotations)) {
    for (nm in head(names(annotations), 30)) {
        cat("\n----------------------------------------\n")
        cat("LIST ENTRY:", nm, "\n")
        cat("----------------------------------------\n")

        item <- annotations[[nm]]

        cat("Class:\n")
        print(class(item))

        cat("Dimensions:\n")
        print(dim(item))

        cat("Length:\n")
        print(length(item))

        cat("Preview:\n")

        if (is.data.frame(item) || is.matrix(item)) {
            print(head(item, 10))
        } else if (is.atomic(item)) {
            print(head(item, 20))
            cat("Names preview:\n")
            print(head(names(item), 20))
        } else {
            str(item, max.level = 2, list.len = 20)
        }
    }
}

sink()

cat("Inspection completed.\n")
cat("Report:\n", REPORT_FILE, "\n", sep = "")
