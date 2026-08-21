#!/bin/bash
#SBATCH --job-name=vista_ex5_1
#SBATCH --output=logs/vista_ex5_1_%A_%a.out
#SBATCH --error=logs/vista_ex5_1_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --array=1-3%3

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/VISTA/run_vista_fold.py"

REFERENCE="${PROJECT_ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"

OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

PREFLIGHT_MANIFEST="${PROJECT_ROOT}/outputs/imputation_beta/Br8667/experiment_5_1/input_validation/experiment_5_1_manifest.json"

FOLD="${SLURM_ARRAY_TASK_ID}"

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  "${OUTPUT_ROOT}/ex5_1/vista"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

export PYTHONPATH="${PROJECT_ROOT}/src/imputation/VISTA:${PROJECT_ROOT}/src/imputation/common:${PROJECT_ROOT}/src/imputation/common:${PYTHONPATH:-}"

echo "===================================================================================================="
echo "VISTA Experiment 5.1 — all 10 Huuki brains including Br8667"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task/fold: ${FOLD}"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 160G"
echo "Python: ${PYTHON}"
echo "Runner: ${SCRIPT}"
echo "Reference: ${REFERENCE}"
echo "Fold directory: ${FOLD_DIR}"
echo "Output root: ${OUTPUT_ROOT}"
echo "Preflight manifest: ${PREFLIGHT_MANIFEST}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${SCRIPT}" || {
    echo "ERROR: VISTA runner not found:"
    echo "${SCRIPT}"
    exit 1
}

test -f "${REFERENCE}" || {
    echo "ERROR: Experiment 5.1 reference not found:"
    echo "${REFERENCE}"
    exit 1
}

test -f "${PREFLIGHT_MANIFEST}" || {
    echo "ERROR: Experiment 5.1 preflight manifest not found:"
    echo "${PREFLIGHT_MANIFEST}"
    exit 1
}

for REQUIRED_FILE in \
    "${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad" \
    "${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad" \
    "${FOLD_DIR}/gene_metrics_and_fold_assignment.h5ad"
do
    test -f "${REQUIRED_FILE}" || {
        echo "ERROR: Required fold file not found:"
        echo "${REQUIRED_FILE}"
        exit 1
    }
done

"${PYTHON}" - "${PREFLIGHT_MANIFEST}" <<'PY'
import json
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])

with manifest_path.open("r", encoding="utf-8") as handle:
    manifest = json.load(handle)

if manifest.get("overall_pass") is not True:
    raise RuntimeError(
        f"Experiment 5.1 preflight did not pass: {manifest_path}"
    )

print("Experiment 5.1 preflight status: PASS")
PY

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${SCRIPT}"

echo "VISTA runner syntax check: PASSED"
echo "Starting fold ${FOLD}: $(date)"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment ex5_1 \
    --experiment-label "Experiment 5.1" \
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
echo "===================================================================================================="
echo "VISTA Experiment 5.1 fold ${FOLD} finished successfully"
echo "Completed: $(date)"
echo "===================================================================================================="
