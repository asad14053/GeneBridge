#!/usr/bin/env Rscript

####################################################################################################
# 01_visualize_matched_N24_visium_layers.R
#
# Task:
#   Visualize Visium layer-wise annotation only for the 24 donors that match Xenium.
#
# Input:
#   data/metadata/fnl_spe_kept_spots_only.rds
#
# Outputs:
#   data/processed/visium/visium_N24_matched_layer_annotated.rds
#
#   outputs/visium_branch/03_annotation/layerwise_matched_N24/tables/
#       visium_N24_matched_layer_annotations.csv
#       visium_N24_matched_layer_counts_by_BrNum.csv
#       visium_N24_matched_layer_proportions_by_BrNum.csv
#       visium_N24_match_report.csv
#
#   outputs/visium_branch/03_annotation/layerwise_matched_N24/figures/
#       visium_N24_matched_layer_maps_panel.png
#       visium_N24_matched_layer_proportions_by_BrNum.png
#       per_sample_layer_maps/*.png
####################################################################################################

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(tidyverse)
  library(ggplot2)
})

PROJECT_ROOT <- getwd()

INPUT_RDS <- file.path(
  PROJECT_ROOT,
  "data",
  "metadata",
  "layer_annotations",
  "fnl_spe_kept_spots_only.rds"
)

OUTPUT_RDS <- file.path(
  PROJECT_ROOT,
  "data",
  "processed",
  "visium",
  "visium_N24_matched_layer_annotated.rds"
)

OUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "visium_branch",
  "03_annotation",
  "layerwise_matched_N24"
)

TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")
PER_SAMPLE_DIR <- file.path(FIGURE_DIR, "per_sample_layer_maps")

dir.create(dirname(OUTPUT_RDS), recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_SAMPLE_DIR, recursive = TRUE, showWarnings = FALSE)

####################################################################################################
# Matched Xenium donors
####################################################################################################

MATCHED_BRNUMS <- c(
  "Br1113", "Br1139", "Br2039", "Br2421", "Br2719", "Br5314",
  "Br5373", "Br5400", "Br5436", "Br5588", "Br5590", "Br5622",
  "Br5639", "Br5746", "Br5931", "Br5973", "Br6032", "Br6389",
  "Br6432", "Br6437", "Br6496", "Br8433", "Br8667", "Br8772"
)

####################################################################################################
# Layer mapping / colors
####################################################################################################

SPD_TO_LAYER <- c(
  "spd07" = "L1/M",
  "spd06" = "L2/3",
  "spd02" = "L3/4",
  "spd05" = "L5",
  "spd03" = "L6",
  "spd01" = "WMtz",
  "spd04" = "WM"
)

LAYER_ORDER <- c("L1/M", "L2/3", "L3/4", "L5", "L6", "WMtz", "WM", "Unknown")

LAYER_COLORS <- c(
  "L1/M" = "#FEAF16",
  "L2/3" = "#3283FE",
  "L3/4" = "#E4E1E3",
  "L5" = "#16FF32",
  "L6" = "#F6222E",
  "WMtz" = "#5A5156",
  "WM" = "#FE00FA",
  "Unknown" = "#999999"
)

####################################################################################################
# Helper functions
####################################################################################################

section <- function(title) {
  cat("\n", paste(rep("=", 100), collapse = ""), "\n", sep = "")
  cat(title, "\n")
  cat(paste(rep("=", 100), collapse = ""), "\n", sep = "")
}

safe_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\:; ,]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-\\.]", "", x)
  x
}

normalize_spd <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x <- gsub("\"", "", x)
  x <- gsub("'", "", x)
  x[x %in% c("", "na", "nan", "none", "null")] <- NA

  idx <- grepl("^spd[0-9]+$", x)
  x[idx] <- sprintf("spd%02d", as.integer(gsub("^spd", "", x[idx])))

  x
}

detect_brnum_column <- function(spe) {
  cd_cols <- colnames(colData(spe))

  candidates <- c(
    "BrNum",
    "brnum",
    "br_num",
    "donor",
    "donor_id",
    "subject",
    "subject_id",
    "sample_id",
    "sample_name",
    "sample_label",
    "Sample"
  )

  found <- candidates[candidates %in% cd_cols]

  if (length(found) > 0) {
    return(found[[1]])
  }

  stop(
    "Could not detect donor/BrNum column. Inspect colData(spe) and set it manually."
  )
}

extract_brnum <- function(x) {
  x <- as.character(x)

  # Extract strings like Br2039 from any sample label.
  br <- stringr::str_extract(x, "Br[0-9]+")
  br
}

detect_layer_column <- function(spe) {
  cd_cols <- colnames(colData(spe))

  cat("\nAvailable colData columns:\n")
  print(cd_cols)

  candidates <- c(
    "layer_annotation",
    "domain_annotations",
    "PRECAST_07",
    "PRECAST_7",
    "spatial_domain",
    "Spatial_Domain",
    "layer_guess",
    "layer",
    "Layer",
    "spatialLIBD"
  )

  found <- candidates[candidates %in% cd_cols]

  if (length(found) > 0) {
    return(found[[1]])
  }

  precast_cols <- cd_cols[grepl("^PRECAST", cd_cols, ignore.case = TRUE)]
  if (length(precast_cols) > 0) {
    return(precast_cols[[1]])
  }

  stop(
    "Could not detect layer annotation column. ",
    "Check colData columns printed above."
  )
}

convert_to_layer_annotation <- function(x) {
  x_chr <- as.character(x)
  x_spd <- normalize_spd(x_chr)

  # If the column is spd/PRECAST labels, map to biological layers.
  if (any(grepl("^spd[0-9]+$", x_spd), na.rm = TRUE)) {
    layer <- unname(SPD_TO_LAYER[x_spd])
    layer[is.na(layer)] <- "Unknown"
    return(layer)
  }

  # Otherwise assume biological layer labels already exist.
  x_chr[is.na(x_chr) | x_chr == ""] <- "Unknown"
  x_chr
}

get_xy_table <- function(spe) {
  coords <- tryCatch(
    spatialCoords(spe),
    error = function(e) NULL
  )

  if (!is.null(coords) && ncol(coords) >= 2) {
    return(data.frame(
      spot_id = colnames(spe),
      x = as.numeric(coords[, 1]),
      y = as.numeric(coords[, 2]),
      stringsAsFactors = FALSE
    ))
  }

  cd <- as.data.frame(colData(spe))

  if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(cd))) {
    return(data.frame(
      spot_id = colnames(spe),
      x = as.numeric(cd$pxl_col_in_fullres),
      y = as.numeric(cd$pxl_row_in_fullres),
      stringsAsFactors = FALSE
    ))
  }

  if (all(c("pxl_col", "pxl_row") %in% colnames(cd))) {
    return(data.frame(
      spot_id = colnames(spe),
      x = as.numeric(cd$pxl_col),
      y = as.numeric(cd$pxl_row),
      stringsAsFactors = FALSE
    ))
  }

  if (all(c("array_col", "array_row") %in% colnames(cd))) {
    return(data.frame(
      spot_id = colnames(spe),
      x = as.numeric(cd$array_col),
      y = as.numeric(cd$array_row),
      stringsAsFactors = FALSE
    ))
  }

  stop("No spatial coordinates found.")
}

####################################################################################################
# Load Visium RDS
####################################################################################################

section("Loading Visium RDS")

if (!file.exists(INPUT_RDS)) {
  stop("Input RDS not found: ", INPUT_RDS)
}

spe <- readRDS(INPUT_RDS)

cat("Object class:\n")
print(class(spe))

cat("\nObject summary:\n")
print(spe)

####################################################################################################
# Detect BrNum and layer columns
####################################################################################################

section("Detecting donor and layer columns")

brnum_col <- detect_brnum_column(spe)
layer_col <- detect_layer_column(spe)

cat("\nDetected donor column:\n")
cat(brnum_col, "\n")

cat("\nDetected layer column:\n")
cat(layer_col, "\n")

raw_brnum <- as.character(colData(spe)[[brnum_col]])
detected_brnum <- extract_brnum(raw_brnum)

# If donor column itself is already BrXXXX, keep it.
detected_brnum[is.na(detected_brnum)] <- raw_brnum[is.na(detected_brnum)]

colData(spe)$BrNum_matched <- detected_brnum

cat("\nDetected BrNum examples:\n")
print(head(unique(colData(spe)$BrNum_matched), 30))

####################################################################################################
# Subset to matched N24 donors
####################################################################################################

section("Subsetting Visium to matched Xenium N24 donors")

available_brnums <- sort(unique(as.character(colData(spe)$BrNum_matched)))
found_brnums <- intersect(MATCHED_BRNUMS, available_brnums)
missing_brnums <- setdiff(MATCHED_BRNUMS, available_brnums)
extra_visium_brnums <- setdiff(available_brnums, MATCHED_BRNUMS)

match_report <- data.frame(
  metric = c(
    "n_expected_xenium_matched_brnums",
    "n_found_in_visium_rds",
    "n_missing_from_visium_rds",
    "n_extra_visium_brnums_not_used"
  ),
  value = c(
    length(MATCHED_BRNUMS),
    length(found_brnums),
    length(missing_brnums),
    length(extra_visium_brnums)
  )
)

write.csv(
  match_report,
  file.path(TABLE_DIR, "visium_N24_match_report.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(BrNum = found_brnums),
  file.path(TABLE_DIR, "visium_N24_found_matched_BrNums.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(BrNum = missing_brnums),
  file.path(TABLE_DIR, "visium_N24_missing_BrNums.csv"),
  row.names = FALSE
)

cat("\nMatch report:\n")
print(match_report)

cat("\nFound matched BrNums:\n")
print(found_brnums)

cat("\nMissing BrNums:\n")
print(missing_brnums)

if (length(found_brnums) == 0) {
  stop("No matched BrNums found. Check donor column and BrNum extraction.")
}

spe_n24 <- spe[, as.character(colData(spe)$BrNum_matched) %in% found_brnums]

cat("\nSubset object:\n")
print(spe_n24)

####################################################################################################
# Add standardized layer annotation
####################################################################################################

section("Adding standardized layer annotation")

colData(spe_n24)$visium_layer_original <- as.character(colData(spe_n24)[[layer_col]])

colData(spe_n24)$visium_layer_annotation <- convert_to_layer_annotation(
  colData(spe_n24)[[layer_col]]
)

colData(spe_n24)$visium_layer_annotation <- factor(
  colData(spe_n24)$visium_layer_annotation,
  levels = LAYER_ORDER
)

colData(spe_n24)$visium_layer_annotation_source <- paste0(
  "matched_N24_from_",
  layer_col
)

cat("\nLayer counts after matching:\n")
print(table(colData(spe_n24)$visium_layer_annotation, useNA = "ifany"))

####################################################################################################
# Save CSV annotation table
####################################################################################################

section("Saving matched N24 Visium layer annotation tables")

xy <- get_xy_table(spe_n24)

annotation_df <- as.data.frame(colData(spe_n24)) %>%
  mutate(spot_id = colnames(spe_n24)) %>%
  left_join(xy, by = "spot_id") %>%
  select(
    spot_id,
    BrNum_matched,
    visium_layer_original,
    visium_layer_annotation,
    visium_layer_annotation_source,
    x,
    y,
    everything()
  )

write.csv(
  annotation_df,
  file.path(TABLE_DIR, "visium_N24_matched_layer_annotations.csv"),
  row.names = FALSE
)

layer_counts <- annotation_df %>%
  count(BrNum_matched, visium_layer_annotation, name = "n_spots")

write.csv(
  layer_counts,
  file.path(TABLE_DIR, "visium_N24_matched_layer_counts_by_BrNum.csv"),
  row.names = FALSE
)

layer_props <- layer_counts %>%
  group_by(BrNum_matched) %>%
  mutate(proportion = n_spots / sum(n_spots)) %>%
  ungroup()

write.csv(
  layer_props,
  file.path(TABLE_DIR, "visium_N24_matched_layer_proportions_by_BrNum.csv"),
  row.names = FALSE
)

cat("Saved tables to:\n")
cat(TABLE_DIR, "\n")

####################################################################################################
# Plot: layer proportions by matched BrNum
####################################################################################################

section("Plotting layer proportions by matched BrNum")

p_props <- layer_props %>%
  mutate(
    BrNum_matched = factor(BrNum_matched, levels = MATCHED_BRNUMS),
    visium_layer_annotation = factor(visium_layer_annotation, levels = LAYER_ORDER)
  ) %>%
  ggplot(
    aes(
      x = BrNum_matched,
      y = proportion,
      fill = visium_layer_annotation
    )
  ) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
  labs(
    title = "Matched N24 Visium layer composition by BrNum",
    x = "Matched Xenium/Visium BrNum",
    y = "Proportion of Visium spots",
    fill = "Layer"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(FIGURE_DIR, "visium_N24_matched_layer_proportions_by_BrNum.png"),
  p_props,
  width = 14,
  height = 6,
  dpi = 300
)

####################################################################################################
# Plot: combined spatial panel for matched N24
####################################################################################################

section("Plotting matched N24 Visium layer map panel")

p_panel <- annotation_df %>%
  mutate(
    BrNum_matched = factor(BrNum_matched, levels = MATCHED_BRNUMS),
    visium_layer_annotation = factor(visium_layer_annotation, levels = LAYER_ORDER)
  ) %>%
  ggplot(aes(x = x, y = y, color = visium_layer_annotation)) +
  geom_point(size = 0.35, alpha = 0.95) +
  scale_color_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
  scale_y_reverse() +
  coord_fixed() +
  facet_wrap(~ BrNum_matched, ncol = 4) +
  labs(
    title = "Matched N24 Visium layer-wise annotation",
    color = "Layer"
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )

ggsave(
  file.path(FIGURE_DIR, "visium_N24_matched_layer_maps_panel.png"),
  p_panel,
  width = 18,
  height = 24,
  dpi = 300
)

####################################################################################################
# Plot: one PNG per matched BrNum
####################################################################################################

section("Plotting one layer map per matched BrNum")

for (br in found_brnums) {
  sub_df <- annotation_df %>%
    filter(BrNum_matched == br) %>%
    mutate(
      visium_layer_annotation = factor(visium_layer_annotation, levels = LAYER_ORDER)
    )

  p <- ggplot(sub_df, aes(x = x, y = y, color = visium_layer_annotation)) +
    geom_point(size = 0.8, alpha = 0.95) +
    scale_color_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = paste0(br, " matched Visium layer annotation"),
      color = "Layer"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

  out_png <- file.path(
    PER_SAMPLE_DIR,
    paste0(safe_filename(br), "_matched_visium_layer_map.png")
  )

  ggsave(
    out_png,
    p,
    width = 7,
    height = 7,
    dpi = 300
  )

  cat("Saved:", out_png, "\n")
}

####################################################################################################
# Save matched N24 RDS
####################################################################################################

section("Saving matched N24 Visium RDS")

saveRDS(spe_n24, OUTPUT_RDS)

cat("Saved RDS:\n")
cat(OUTPUT_RDS, "\n")

section("Done")

cat("Matched BrNums found:", length(found_brnums), "\n")
cat("Matched BrNums missing:", length(missing_brnums), "\n")

cat("\nMain outputs:\n")
cat("RDS:\n")
cat(OUTPUT_RDS, "\n")
cat("Tables:\n")
cat(TABLE_DIR, "\n")
cat("Figures:\n")
cat(FIGURE_DIR, "\n")
cat("Per-sample figures:\n")
cat(PER_SAMPLE_DIR, "\n")