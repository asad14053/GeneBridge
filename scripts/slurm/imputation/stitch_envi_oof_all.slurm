#!/bin/bash
#SBATCH --job-name=envi_oof
#SBATCH --output=logs/envi_oof/%x_%A_%a.out
#SBATCH --error=logs/envi_oof/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-28%28

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/ENVI/01_stitch_oof_predictions.py"

MANIFEST="${PROJECT_ROOT}/data/metadata/envi_oof_stitch_manifest.tsv"

cd "${PROJECT_ROOT}"

LINE=$(
    sed -n \
    "$((SLURM_ARRAY_TASK_ID + 1))p" \
    "${MANIFEST}"
)

EXPERIMENT=$(
    echo "${LINE}" |
    cut -f1
)

DONOR=$(
    echo "${LINE}" |
    cut -f2
)

if [[ -z "${EXPERIMENT}" || -z "${DONOR}" ]]; then
    echo "ERROR: manifest lookup failed"
    exit 1
fi

echo "===================================================================================================="
echo "ENVI OOF STITCH"
echo "===================================================================================================="
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Array task  : ${SLURM_ARRAY_TASK_ID}"
echo "Node        : $(hostname)"
echo "Experiment  : ${EXPERIMENT}"
echo "Target donor: ${DONOR}"
echo "Started     : $(date)"
echo "===================================================================================================="

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --donor "${DONOR}"

echo
echo "===================================================================================================="
echo "SLURM STATUS: PASS"
echo "Experiment: ${EXPERIMENT}"
echo "Target donor: ${DONOR}"
echo "Completed: $(date)"
echo "===================================================================================================="
