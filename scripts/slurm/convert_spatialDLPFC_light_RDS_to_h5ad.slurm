#!/bin/bash
#SBATCH --job-name=visium_light_h5ad
#SBATCH --output=logs/visium_light_h5ad_%j.out
#SBATCH --error=logs/visium_light_h5ad_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"

SCRIPT="${PROJECT_ROOT}/src/visium_branch/01_build_datastructure/02_convert_spatialDLPFC_light_RDS_to_h5ad.R"

INPUT_RDS="${PROJECT_ROOT}/data/processed/visium/spatialDLPFC_Visium_sce_light.RDS"

OUTPUT_H5AD="${PROJECT_ROOT}/data/processed/visium/spatialDLPFC_Visium_sce.h5ad"

export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"
export HDF5_USE_FILE_LOCKING=FALSE

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "===================================================================================================="
echo "Convert Huuki Visium light RDS to H5AD"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "Input: ${INPUT_RDS}"
echo "Output: ${OUTPUT_H5AD}"
echo "===================================================================================================="

if [ ! -x "${RSCRIPT}" ]; then
    echo "ERROR: Rscript executable was not found:"
    echo "${RSCRIPT}"
    exit 1
fi

if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: conversion script was not found:"
    echo "${SCRIPT}"
    exit 1
fi

if [ ! -f "${INPUT_RDS}" ]; then
    echo "ERROR: input RDS was not found:"
    echo "${INPUT_RDS}"
    exit 1
fi

echo ""
echo "R version:"
"${RSCRIPT}" --version

echo ""
echo "Checking required R packages..."

"${RSCRIPT}" - <<'RS'
required_packages <- c(
    "anndataR",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors",
    "scuttle"
)

for (package_name in required_packages) {

    if (!requireNamespace(
        package_name,
        quietly = TRUE
    )) {

        stop(
            "Missing required R package: ",
            package_name
        )
    }

    cat(
        package_name,
        " = ",
        as.character(
            packageVersion(
                package_name
            )
        ),
        "\n",
        sep = ""
    )
}
RS

echo ""
echo "Starting conversion..."
echo ""

"${RSCRIPT}" "${SCRIPT}"

echo ""
echo "Checking output..."

if [ ! -f "${OUTPUT_H5AD}" ]; then
    echo "ERROR: expected H5AD output was not created:"
    echo "${OUTPUT_H5AD}"
    exit 1
fi

ls -lh "${OUTPUT_H5AD}"

echo ""
echo "===================================================================================================="
echo "Conversion completed successfully"
echo "Finished: $(date)"
echo "===================================================================================================="
