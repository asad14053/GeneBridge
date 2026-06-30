#!/usr/bin/env Rscript

####################################################################################################
# 01_visualize_visium_layerwise_annotation.R
#
# Task:
#   Visualize Visium layer-wise annotation from:
#
#       data/metadata/fnl_spe_kept_spots_only.rds
#
# What this script does:
#   1. Reads Visium SpatialExperiment RDS
#   2. Detects existing layer/PRECAST annotation column
#   3. Creates one combined panel plot across samples
#   4. Creates one separate PNG per sample
#   5. Saves spot-level annotation CSV and layer counts
#
# Expected useful columns:
#   PRECAST_07
#   layer_annotation
#   domain_annotations
#   spatial_domain
#   layer_guess
#
####################################################################################################

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(tidyverse)
  library(escheR)
  library(ggplot2)
  library(here)
})

####################################################################################################
# Paths
####################################################################################################

PROJECT_ROOT <- getwd()

INPUT_RDS <- file.path(
  PROJECT_ROOT,
  "data",
  "metadata",
  "layer_annotations",
  "fnl_spe_kept_spots_only.rds"
)

OUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "visium_branch",
  "03_annotation",
  "layerwise"
)

TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")
PER_SAMPLE_DIR <- file.path(FIGURE_DIR, "per_sample_layer_maps")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PER_SAMPLE_DIR, recursive = TRUE, showWarnings = FALSE)

####################################################################################################
# Layer color logic
####################################################################################################

# If the detected column is already biological layer names, use this order.
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

# If the detected column is PRECAST/spd labels, convert to biological layers.
SPD_TO_LAYER <- c(
  "spd07" = "L1/M",
  "spd06" = "L2/3",
  "spd02" = "L3/4",
  "spd05" = "L5",
  "spd03" = "L6",
  "spd01" = "WMtz",
  "spd04" = "WM"
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

detect_layer_column <- function(spe) {
  cd_cols <- colnames(colData(spe))

  cat("\nAvailable colData columns:\n")
  print(cd_cols)

  # Priority order.
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

  # Fallback: any PRECAST column.
  precast_cols <- cd_cols[grepl("^PRECAST", cd_cols, ignore.case = TRUE)]
  if (length(precast_cols) > 0) {
    return(precast_cols[[1]])
  }

  stop(
    "Could not detect layer annotation column. ",
    "Please inspect colData columns above and set the correct column manually."
  )
}

get_sample_column <- function(spe) {
  cd_cols <- colnames(colData(spe))

  candidates <- c(
    "sample_label",
    "sample_id",
    "Sample",
    "sample",
    "BrNum",
    "brnum",
    "donor",
    "subject",
    "slide_id"
  )

  found <- candidates[candidates %in% cd_cols]

  if (length(found) > 0) {
    return(found[[1]])
  }

  warning("No sample column detected. Creating one sample called all_samples.")
  colData(spe)$sample_label <- "all_samples"
  return("sample_label")
}

make_sample_label <- function(spe, sample_col) {
  cd <- as.data.frame(colData(spe))

  # If sample_label already exists, keep it.
  if ("sample_label" %in% colnames(cd)) {
    return(as.character(cd$sample_label))
  }

  # If both brnum and dx exist, use PI-style label.
  if (all(c("brnum", "dx") %in% colnames(cd))) {
    return(paste0(cd$brnum, "_", toupper(cd$dx)))
  }

  if (all(c("BrNum", "Dx") %in% colnames(cd))) {
    return(paste0(cd$BrNum, "_", toupper(cd$Dx)))
  }

  return(as.character(cd[[sample_col]]))
}

convert_to_layer_annotation <- function(x) {
  x_chr <- as.character(x)

  # If values look like spdXX, map them.
  x_spd <- normalize_spd(x_chr)

  if (any(grepl("^spd[0-9]+$", x_spd), na.rm = TRUE)) {
    layer <- unname(SPD_TO_LAYER[x_spd])
    layer[is.na(layer)] <- "Unknown"
    return(layer)
  }

  # Otherwise assume it is already biological layer labels.
  x_chr[is.na(x_chr) | x_chr == ""] <- "Unknown"
  return(x_chr)
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

  stop("No spatial coordinates found in spatialCoords(spe) or colData pixel columns.")
}

####################################################################################################
# Load RDS
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
# Detect annotation and sample columns
####################################################################################################

section("Detecting layer annotation column")

layer_col <- detect_layer_column(spe)
sample_col <- get_sample_column(spe)

cat("\nDetected layer column:\n")
cat(layer_col, "\n")

cat("\nDetected sample column:\n")
cat(sample_col, "\n")

colData(spe)$sample_label_for_plot <- make_sample_label(spe, sample_col)

colData(spe)$visium_layer_original <- as.character(colData(spe)[[layer_col]])

colData(spe)$visium_layer_annotation <- convert_to_layer_annotation(
  colData(spe)[[layer_col]]
)

colData(spe)$visium_layer_annotation <- factor(
  colData(spe)$visium_layer_annotation,
  levels = unique(c(LAYER_ORDER, sort(unique(colData(spe)$visium_layer_annotation))))
)

cat("\nLayer counts:\n")
print(table(colData(spe)$visium_layer_annotation, useNA = "ifany"))

cat("\nNumber of samples:\n")
print(length(unique(colData(spe)$sample_label_for_plot)))

####################################################################################################
# Save CSV outputs
####################################################################################################

section("Saving CSV tables")

xy <- get_xy_table(spe)

annotation_df <- as.data.frame(colData(spe)) %>%
  mutate(spot_id = colnames(spe)) %>%
  left_join(xy, by = "spot_id") %>%
  select(
    spot_id,
    sample_label_for_plot,
    visium_layer_original,
    visium_layer_annotation,
    x,
    y,
    everything()
  )

write.csv(
  annotation_df,
  file.path(TABLE_DIR, "visium_layerwise_spot_annotations.csv"),
  row.names = FALSE
)

layer_counts <- annotation_df %>%
  count(sample_label_for_plot, visium_layer_annotation, name = "n_spots")

write.csv(
  layer_counts,
  file.path(TABLE_DIR, "visium_layer_counts_by_sample.csv"),
  row.names = FALSE
)

layer_props <- layer_counts %>%
  group_by(sample_label_for_plot) %>%
  mutate(proportion = n_spots / sum(n_spots)) %>%
  ungroup()

write.csv(
  layer_props,
  file.path(TABLE_DIR, "visium_layer_proportions_by_sample.csv"),
  row.names = FALSE
)

cat("Saved tables to:\n")
cat(TABLE_DIR, "\n")

####################################################################################################
# Plot 1: overall layer distribution
####################################################################################################

section("Plotting overall layer distribution")

overall_counts <- annotation_df %>%
  count(visium_layer_annotation, name = "n_spots") %>%
  mutate(
    visium_layer_annotation = factor(
      visium_layer_annotation,
      levels = levels(colData(spe)$visium_layer_annotation)
    )
  )

p_overall <- ggplot(
  overall_counts,
  aes(x = visium_layer_annotation, y = n_spots, fill = visium_layer_annotation)
) +
  geom_col() +
  scale_fill_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
  labs(
    title = "Overall Visium layer-wise annotation",
    x = "Layer",
    y = "Number of spots"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(FIGURE_DIR, "visium_layer_overall_distribution.png"),
  p_overall,
  width = 8,
  height = 5,
  dpi = 300
)

####################################################################################################
# Plot 2: layer proportions by sample
####################################################################################################

section("Plotting layer proportions by sample")

p_props <- layer_props %>%
  mutate(
    visium_layer_annotation = factor(
      visium_layer_annotation,
      levels = levels(colData(spe)$visium_layer_annotation)
    )
  ) %>%
  ggplot(
    aes(
      x = sample_label_for_plot,
      y = proportion,
      fill = visium_layer_annotation
    )
  ) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
  labs(
    title = "Visium layer composition by sample",
    x = "Sample",
    y = "Proportion of spots",
    fill = "Layer"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file.path(FIGURE_DIR, "visium_layer_proportions_by_sample.png"),
  p_props,
  width = 16,
  height = 6,
  dpi = 300
)

####################################################################################################
# Plot 3: combined spatial panel
####################################################################################################

section("Plotting combined spatial panel")

p_panel <- annotation_df %>%
  mutate(
    visium_layer_annotation = factor(
      visium_layer_annotation,
      levels = levels(colData(spe)$visium_layer_annotation)
    )
  ) %>%
  ggplot(aes(x = x, y = y, color = visium_layer_annotation)) +
  geom_point(size = 0.25, alpha = 0.9) +
  scale_color_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
  scale_y_reverse() +
  coord_fixed() +
  facet_wrap(~ sample_label_for_plot, ncol = 7) +
  labs(
    title = "Visium layer-wise annotation by sample",
    color = "Layer"
  ) +
  theme_void() +
  theme(
    strip.text = element_text(size = 6),
    legend.position = "bottom"
  )

ggsave(
  file.path(FIGURE_DIR, "visium_layer_maps_panel.png"),
  p_panel,
  width = 22,
  height = 18,
  dpi = 300
)

####################################################################################################
# Plot 4: one PNG per sample
####################################################################################################

section("Plotting one layer map per sample")

samples <- sort(unique(annotation_df$sample_label_for_plot))

for (smp in samples) {
  sub_df <- annotation_df %>%
    filter(sample_label_for_plot == smp) %>%
    mutate(
      visium_layer_annotation = factor(
        visium_layer_annotation,
        levels = levels(colData(spe)$visium_layer_annotation)
      )
    )

  p <- ggplot(sub_df, aes(x = x, y = y, color = visium_layer_annotation)) +
    geom_point(size = 0.7, alpha = 0.95) +
    scale_color_manual(values = LAYER_COLORS, drop = FALSE, na.value = "#999999") +
    scale_y_reverse() +
    coord_fixed() +
    labs(
      title = paste0(smp, " Visium layer-wise annotation"),
      color = "Layer"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )

  out_png <- file.path(
    PER_SAMPLE_DIR,
    paste0(safe_filename(smp), "_visium_layer_map.png")
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
# Save updated RDS with standardized columns
####################################################################################################

section("Saving standardized Visium layer-annotated RDS")

OUTPUT_RDS <- file.path(
  PROJECT_ROOT,
  "data",
  "processed",
  "visium",
  "visium_layerwise_visualized.rds"
)

dir.create(dirname(OUTPUT_RDS), recursive = TRUE, showWarnings = FALSE)

saveRDS(spe, OUTPUT_RDS)

cat("Saved RDS:\n")
cat(OUTPUT_RDS, "\n")

section("Done")

cat("Main outputs:\n")
cat("Tables:\n")
cat(TABLE_DIR, "\n")
cat("Figures:\n")
cat(FIGURE_DIR, "\n")
cat("Per-sample figures:\n")
cat(PER_SAMPLE_DIR, "\n")