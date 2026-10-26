#!/bin/bash
#SBATCH --job-name=aivc_audit
#SBATCH --output=logs/aivc_presplit_%A_%a.out
#SBATCH --error=logs/aivc_presplit_%A_%a.err
#SBATCH --array=0-23%24
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

mkdir -p logs
mkdir -p outputs/AIVC_beta/presplit_audit/per_donor


# ==============================================================================
# Discover all production ENVI donor H5ADs
# ==============================================================================

mapfile -t FILES < <(
    find \
        data/processed/imputation_full/ex1_ntc/envi \
        data/processed/imputation_full/ex2_scz/envi \
        -type f \
        -name 'spatial_data_xenium_*_ENVI_full_transcriptome.h5ad' \
    | sort
)

N=${#FILES[@]}

if [[ "${N}" -ne 24 ]]; then
    echo "ERROR: expected 24 ENVI files; found ${N}"
    printf '%s\n' "${FILES[@]}"
    exit 1
fi

FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

DONOR=$(
    basename "${FILE}" \
    | sed -E \
      's/spatial_data_xenium_(Br[0-9]+)_ENVI_full_transcriptome\.h5ad/\1/'
)

echo "======================================================================"
echo "AIVC PRE-SPLIT DONOR AUDIT"
echo "Array job : ${SLURM_ARRAY_JOB_ID}"
echo "Task      : ${SLURM_ARRAY_TASK_ID}"
echo "Donor     : ${DONOR}"
echo "File      : ${FILE}"
echo "Host      : $(hostname)"
echo "CPUs      : ${SLURM_CPUS_PER_TASK}"
echo "Started   : $(date)"
echo "======================================================================"

"${PYTHON}" -u \
    src/AIVC/data_prep/00_presplit_audit_donor.py \
    --envi "${FILE}" \
    --donor "${DONOR}"

echo "======================================================================"
echo "DONE"
echo "Donor   : ${DONOR}"
echo "Finished: $(date)"
echo "======================================================================"
