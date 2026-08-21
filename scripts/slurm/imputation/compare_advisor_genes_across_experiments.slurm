#!/bin/bash
#SBATCH --job-name=advisor_gene_compare
#SBATCH --output=logs/advisor_gene_compare_%j.out
#SBATCH --error=logs/advisor_gene_compare_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/common/compare_advisor_genes_across_experiments.py"

MODEL="${1:?Missing model key, for example vista}"
MODEL_LABEL="${2:?Missing model label, for example VISTA}"
COMBINED_DIR="${3:-auto}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "===================================================================================================="
echo "Advisor 11-gene comparison: Ex-5 vs Ex-5.1 vs Ex-5.3"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Memory: 160G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Model: ${MODEL}"
echo "Model label: ${MODEL_LABEL}"
echo "Combined directory selection: ${COMBINED_DIR}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${SCRIPT}" || {
    echo "ERROR: Comparison script not found:"
    echo "${SCRIPT}"
    exit 1
}

test -f \
  "${PROJECT_ROOT}/data/metadata/imputation/advisor_highlight_genes.csv" \
|| {
    echo "ERROR: Advisor highlight-gene file is missing."
    exit 1
}

for EXPERIMENT in ex5 ex5_1 ex5_3
do
    MODEL_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667/${EXPERIMENT}/${MODEL}"

    test -d "${MODEL_ROOT}" || {
        echo "ERROR: Model output directory not found:"
        echo "${MODEL_ROOT}"
        exit 1
    }
done

"${PYTHON}" -m py_compile \
  "${SCRIPT}"

echo "Python syntax check: PASSED"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${SCRIPT}" \
    --model "${MODEL}" \
    --model-label "${MODEL_LABEL}" \
    --combined-dir "${COMBINED_DIR}"

echo
echo "Completed: $(date)"
