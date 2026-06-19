# Task 1 Meeting Notes: Br2039 Patient Map

## Advisor-requested output

The corrected Task 1 output contains:

1. Xenium cell-level spatial map with broad cell-type annotation.
2. Visium spot-level spatial map with cortical layer-by-layer annotation.

## Key statement

This is the same donor/patient, Br2039, but Xenium and Visium are not necessarily the exact same physical tissue section.

Xenium is cell-level. Each point in the Xenium map is one segmented cell.

Visium is spot-level. Each point in the Visium map is one spot, annotated as a predicted cortical layer.

## Xenium annotation

Xenium is annotated using broad marker-based cell-type scoring.

## Visium layer annotation

Visium is annotated layer-by-layer using first-pass marker-based cortical layer scoring:

- L1
- L2/3
- L4
- L5
- L6
- WM
- Unknown

## Important limitation

This is a first-pass marker-based layer annotation.

For a stronger version, the next step should be to use reference-based cortical layer annotation from LIBD/spatialLIBD, where DLPFC Visium spots have manual layer labels for the six cortical layers plus white matter.

## Generated figures

- outputs/task1_patient_map/Br2039/figures/Br2039_xenium_celltype_map.png
- outputs/task1_patient_map/Br2039/figures/Br2039_visium_layer_map.png
- outputs/task1_patient_map/Br2039/figures/Br2039_combined_xenium_celltype_visium_layer_map.png

## Generated tables

- outputs/task1_patient_map/Br2039/tables/Br2039_xenium_celltype_annotations.csv
- outputs/task1_patient_map/Br2039/tables/Br2039_visium_layer_annotations.csv
- outputs/task1_patient_map/Br2039/tables/Br2039_xenium_celltype_counts.csv
- outputs/task1_patient_map/Br2039/tables/Br2039_visium_layer_counts.csv
