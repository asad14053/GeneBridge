#!/bin/bash
#SBATCH --job-name=envi_prod
#SBATCH --output=logs/envi_prod_%A_%a.out
#SBATCH --error=logs/envi_prod_%A_%a.err
#SBATCH --array=0-22
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"
MANIFEST="${ROOT}/data/metadata/envi_production_23donors.tsv"

REFERENCE="${ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"

cd "${ROOT}"
mkdir -p logs

# +2 because line 1 is header and array starts at 0.
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "${MANIFEST}")

DONOR=$(echo "${LINE}" | cut -f1)
DX=$(echo "${LINE}" | cut -f2)
EXPERIMENT=$(echo "${LINE}" | cut -f3)
TARGET=$(echo "${LINE}" | cut -f4)

OUTDIR="${ROOT}/data/processed/imputation_full/${EXPERIMENT}/envi/${DONOR}"

mkdir -p "${OUTDIR}"

echo "===================================================================================================="
echo "ENVI PRODUCTION"
echo "Array job : ${SLURM_ARRAY_JOB_ID}"
echo "Task      : ${SLURM_ARRAY_TASK_ID}"
echo "Donor     : ${DONOR}"
echo "Dx        : ${DX}"
echo "Experiment: ${EXPERIMENT}"
echo "Memory    : 200G"
echo "CPUs      : 16"
echo "Started   : $(date)"
echo "===================================================================================================="

# Do not recompute a successfully completed donor.
if [[ -f "${OUTDIR}/complete.flag" ]]; then
    echo "Already complete: ${DONOR}"
    cat "${OUTDIR}/complete.flag"
    exit 0
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

"${PYTHON}" -u \
src/imputation/ENVI/run_envi_production.py \
    --target "${TARGET}" \
    --reference "${REFERENCE}" \
    --output-dir "${OUTDIR}" \
    --donor "${DONOR}" \
    --experiment "${EXPERIMENT}" \
    --training-steps 16000 \
    --batch-size 128 \
    --num-hvg 2048 \
    --num-cov-genes 64 \
    --k-nearest 8 \
    --spatial-dist pois \
    --sc-dist nb \
    --stable-eps 1e-6 \
    --seed "$((8667 + SLURM_ARRAY_TASK_ID))"

echo "===================================================================================================="
echo "SUCCESS: ${DONOR}"
echo "Finished: $(date)"
echo "===================================================================================================="
