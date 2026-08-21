#!/bin/bash
#SBATCH --job-name=gene_folds_h5ad
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

INPUT_H5AD="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/spatial_data_xenium_Br8667_vista.h5ad"

OUTPUT_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"

SCRIPT="${PROJECT_ROOT}/src/imputation_beta/01_create_balanced_gene_folds_h5ad.py"

cd "${PROJECT_ROOT}"

mkdir -p logs "${OUTPUT_DIR}"

echo "=================================================="
echo "Balanced gene-fold generation"
echo "=================================================="
echo "Host: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "Python: ${PYTHON}"
echo "Input: ${INPUT_H5AD}"
echo "Output: ${OUTPUT_DIR}"
echo "=================================================="

"${PYTHON}" -u "${SCRIPT}" \
    --input "${INPUT_H5AD}" \
    --output-dir "${OUTPUT_DIR}" \
    --counts-layer auto \
    --spatial-key spatial \
    --expected-genes 300 \
    --n-folds 3 \
    --heldout-per-fold 100 \
    --target-sum 10000 \
    --k-neighbors 15 \
    --n-strata 30 \
    --seed 14053 \
    --n-jobs "${SLURM_CPUS_PER_TASK}"

echo "=================================================="
echo "Finished: $(date)"
echo "=================================================="
