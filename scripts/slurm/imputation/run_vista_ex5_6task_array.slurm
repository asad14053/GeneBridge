#!/bin/bash
#SBATCH --job-name=vista_ex5
#SBATCH --output=logs/vista_N24/%x_%A_%a.out
#SBATCH --error=logs/vista_N24/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6%6

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="${PROJECT_ROOT}/../conda_envs/Imputation/bin/python"
SCRIPT="${PROJECT_ROOT}/src/imputation/VISTA/run_vista_fold.py"

FOLD_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/gene_folds_200_100"

REF_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/snrnaseq_target_references"

OUTPUT_BASE="${PROJECT_ROOT}/outputs/imputation_full/VISTA_N24"

TASK="${SLURM_ARRAY_TASK_ID}"

case "${TASK}" in
    1)
        DONOR="Br6432"
        FOLD=1
        ;;
    2)
        DONOR="Br6432"
        FOLD=2
        ;;
    3)
        DONOR="Br6432"
        FOLD=3
        ;;
    4)
        DONOR="Br8667"
        FOLD=1
        ;;
    5)
        DONOR="Br8667"
        FOLD=2
        ;;
    6)
        DONOR="Br8667"
        FOLD=3
        ;;
    *)
        echo "ERROR: invalid task ${TASK}"
        exit 1
        ;;
esac

if [[ "${DONOR}" == "Br6432" ]]; then

    REFERENCE="${REF_ROOT}/Br6432/seq_data_huuki_snrna_Br6432_pseudocount.h5ad"

else

    REFERENCE="${REF_ROOT}/Br8667/seq_data_huuki_snrna_Br8667.h5ad"

fi

EXPERIMENT="ex5"
EXPERIMENT_LABEL="Experiment 5"

FOLD_DIR="${FOLD_ROOT}/${DONOR}"
OUTPUT_ROOT="${OUTPUT_BASE}/${DONOR}"

OBSERVED="${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad"
HELDOUT="${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad"

for FILE in \
    "${REFERENCE}" \
    "${OBSERVED}" \
    "${HELDOUT}"
do
    if [[ ! -f "${FILE}" ]]; then
        echo "ERROR: missing ${FILE}"
        exit 1
    fi
done

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: Python missing: ${PYTHON}"
    exit 1
fi

cd "${PROJECT_ROOT}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg

echo "===================================================================================================="
echo "VISTA EX-5"
echo "===================================================================================================="
echo "Job           : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Target donor  : ${DONOR}"
echo "Fold          : ${FOLD}"
echo "Experiment    : ${EXPERIMENT}"
echo "Reference     : ${REFERENCE}"
echo "Observed      : ${OBSERVED}"
echo "Heldout       : ${HELDOUT}"
echo "Output root   : ${OUTPUT_ROOT}"
echo "CPUs          : ${SLURM_CPUS_PER_TASK}"
echo "Memory        : 200G"
echo "Started       : $(date)"
echo "===================================================================================================="

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --reference-mode benchmark_300 \
    --run-mode full \
    --full-epochs 200 \
    --batch-size 128 \
    --n-latent 32 \
    --neighbor-size 20 \
    --seed 8667 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10

echo
echo "FOLD STATUS: PASS"
echo "Donor     : ${DONOR}"
echo "Experiment: ${EXPERIMENT}"
echo "Fold      : ${FOLD}"
echo
echo "FINAL STATUS: PASS"
