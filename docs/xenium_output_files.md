# Xenium Files

## cell_feature_matrix.h5

Contains:
- Gene expression matrix
- Cells × Genes

Use:
- Main imputation input
- Cell-level expression analysis

---

## cells.parquet.gz

Contains:
- Cell metadata
- Cell coordinates
- Segmentation information

Use:
- Cell-level annotation
- Spatial mapping

---

## cell_boundaries.csv.gz

Contains:
- Cell segmentation boundaries

Use:
- Cell morphology analysis
- Spatial visualization

---

## nucleus_boundaries.csv.gz

Contains:
- Nucleus segmentation boundaries

Use:
- Nucleus morphology analysis
- Segmentation quality control

---

## transcripts.zarr.zip

Contains:
- Individual transcript locations
- Molecule-level coordinates

Use:
- Subcellular spatial analysis
- Transcript localization studies

---

# Xenium Summary

Technology:
- 10x Genomics Xenium

Resolution:
- Single-cell / Subcellular

Measures:
- Targeted gene panel

Output:
- Cell-resolved spatial transcriptomics

Primary Use in This Project:
- Cell-level expression profiling
- Spatial cell-type analysis
- Target of Visium-guided imputation
- Schizophrenia microenvironment characterization

----------------------------------------------------------------

# Visium Files

## matrix.mtx.gz

Contains:
- Raw gene expression count matrix
- Spots × Genes matrix
- Sparse count representation

Use:
- Primary transcriptomics data
- Differential expression analysis
- Spatial gene expression profiling
- Imputation reference

---

## features.tsv.gz

Contains:
- Gene IDs
- Gene symbols
- Feature annotations

Use:
- Gene mapping
- Shared gene discovery with Xenium
- Gene filtering and annotation

---

## barcodes.tsv.gz

Contains:
- Unique spot identifiers

Use:
- Spot indexing
- Linking spots to expression matrix rows

---

## tissue_positions.csv.gz

Contains:
- Spot barcode
- Tissue coordinates
- In-tissue flag
- Pixel coordinates

Use:
- Spatial mapping
- Spot visualization
- Spatial alignment
- Layer-specific analysis

---

## tissue_hires_image.png.gz

Contains:
- High-resolution histology image

Use:
- Tissue visualization
- Histological interpretation
- Future image-based modeling

---

## tissue_lowres_image.png.gz

Contains:
- Low-resolution histology image

Use:
- Quick visualization
- Exploratory analysis

---

## aligned_fiducials.jpg.gz

Contains:
- Fiducial marker alignment image

Use:
- Quality control
- Spatial registration verification

---

## scalefactors_json.json.gz

Contains:
- Image scaling factors
- Coordinate transformation parameters

Use:
- Mapping spots to image coordinates
- Spatial visualization

---

# Visium Summary

Technology:
- 10x Genomics Visium

Resolution:
- Spot-level

Measures:
- Whole-transcriptome gene expression

Output:
- Spatially resolved gene expression matrix

Primary Use in This Project:
- Spatial reference dataset
- Xenium imputation guidance
- Schizophrenia differential expression analysis
- Cortical layer characterization
- Cell-type deconvolution