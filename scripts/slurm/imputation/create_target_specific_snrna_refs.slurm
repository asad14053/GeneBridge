#!/bin/bash
#SBATCH --job-name=snrna_refs
#SBATCH --output=logs/snrna_target_references/%A_%a.out
#SBATCH --error=logs/snrna_target_references/%A_%a.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-2%2

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="${PROJECT_ROOT}/../mjabin/conda_envs/genebridge_scanpy/bin/python"

# Explicit canonical Python path
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

REFERENCE="${PROJECT_ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"

OUTPUT_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/snrnaseq_target_references"

DONORS=(
Br6432
Br8667
)

DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID-1))]}"

cd "${PROJECT_ROOT}"

echo "===================================================================================================="
echo "CREATE TARGET-SPECIFIC snRNA REFERENCE"
echo "===================================================================================================="
echo "Job       : ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Donor     : ${DONOR}"
echo "Node      : $(hostname)"
echo "Started   : $(date)"
echo
echo "Source:"
echo "${REFERENCE}"
echo
echo "Output:"
echo "${OUTPUT_ROOT}/${DONOR}"
echo "===================================================================================================="

"${PYTHON}" -u \
    "${PROJECT_ROOT}/src/imputation/reference_preparation/create_target_specific_snrna_references.py" \
    --donor "${DONOR}" \
    --reference "${REFERENCE}" \
    --output-root "${OUTPUT_ROOT}"

echo
echo "===================================================================================================="
echo "SUCCESS"
echo "Donor    : ${DONOR}"
echo "Completed: $(date)"
echo "===================================================================================================="
