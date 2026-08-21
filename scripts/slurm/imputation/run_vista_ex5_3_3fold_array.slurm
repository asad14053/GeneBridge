#!/bin/bash
#SBATCH --job-name=vista_ex5_3
#SBATCH --output=logs/vista_ex5_3_%A_%a.out
#SBATCH --error=logs/vista_ex5_3_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --array=1-3%3

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

RUNNER="${PROJECT_ROOT}/src/imputation/VISTA/run_vista_fold.py"

REFERENCE="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/seq_data_huuki_snrna_9samples_excluding_xenium_overlap.h5ad"

TARGET="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/spatial_data_xenium_Br8667_vista_qc.h5ad"

FOLD_DIR="${PROJECT_ROOT}/data/processed/imputation_beta/Br8667/gene_folds_200_100"

OUTPUT_ROOT="${PROJECT_ROOT}/outputs/imputation_beta/Br8667"

FOLD="${SLURM_ARRAY_TASK_ID}"

cd "${PROJECT_ROOT}"

mkdir -p \
  logs \
  "${OUTPUT_ROOT}/ex5_3/vista"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export HDF5_USE_FILE_LOCKING=FALSE
export MPLBACKEND=Agg

export PYTHONPATH="${PROJECT_ROOT}/src/imputation/VISTA:${PROJECT_ROOT}/src/imputation/common:${PYTHONPATH:-}"

echo "===================================================================================================="
echo "VISTA Experiment 5.3"
echo "Nine Huuki brains excluding all Br8667 cells"
echo "===================================================================================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task / fold: ${FOLD}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo
echo "Memory: 160G"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Python: ${PYTHON}"
echo "Runner: ${RUNNER}"
echo "Reference: ${REFERENCE}"
echo "Target: ${TARGET}"
echo "Fold directory: ${FOLD_DIR}"
echo "Output root: ${OUTPUT_ROOT}"
echo "===================================================================================================="

test -x "${PYTHON}" || {
    echo "ERROR: Python executable not found:"
    echo "${PYTHON}"
    exit 1
}

test -f "${RUNNER}" || {
    echo "ERROR: VISTA runner not found:"
    echo "${RUNNER}"
    exit 1
}

test -f "${REFERENCE}" || {
    echo "ERROR: Experiment 5.3 reference not found:"
    echo "${REFERENCE}"
    exit 1
}

test -f "${TARGET}" || {
    echo "ERROR: Xenium target not found:"
    echo "${TARGET}"
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

echo
echo "Validating Experiment 5.3 reference cohort..."

"${PYTHON}" - "${REFERENCE}" "${TARGET}" <<'PY'
import sys
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse


reference_path = Path(sys.argv[1])
target_path = Path(sys.argv[2])

reference = ad.read_h5ad(reference_path, backed="r")
target = ad.read_h5ad(target_path, backed="r")

try:
    print("Reference shape:", reference.shape)
    print("Target shape:", target.shape)

    if "BrNum" not in reference.obs.columns:
        raise KeyError(
            "Reference obs does not contain the required BrNum column."
        )

    brains = sorted(
        reference.obs["BrNum"]
        .astype(str)
        .unique()
        .tolist()
    )

    print("Reference brains:", brains)
    print("Number of brains:", len(brains))

    if len(brains) != 9:
        raise RuntimeError(
            f"Expected 9 brains, but found {len(brains)}."
        )

    if "Br8667" in brains:
        raise RuntimeError(
            "Br8667 is still present in the Experiment 5.3 reference."
        )

    if "Sample" in reference.obs.columns:
        samples = sorted(
            reference.obs["Sample"]
            .astype(str)
            .unique()
            .tolist()
        )

        br8667_samples = [
            sample
            for sample in samples
            if sample.startswith("Br8667")
        ]

        print("Number of tissue variants:", len(samples))
        print("Br8667 tissue variants:", br8667_samples)

        if br8667_samples:
            raise RuntimeError(
                "A Br8667 tissue variant remains in the reference: "
                + ", ".join(br8667_samples)
            )

    reference_genes = set(
        reference.var_names.astype(str)
    )

    target_genes = set(
        target.var_names.astype(str)
    )

    missing_target_genes = sorted(
        target_genes - reference_genes
    )

    print("Target genes:", len(target_genes))
    print(
        "Target genes missing from reference:",
        len(missing_target_genes),
    )

    if missing_target_genes:
        raise RuntimeError(
            "Target genes are missing from the reference: "
            + ", ".join(missing_target_genes[:20])
        )

    if target.n_vars != 300:
        raise RuntimeError(
            f"Expected 300 target genes, but found {target.n_vars}."
        )

    print("Experiment 5.3 reference validation: PASS")

finally:
    reference.file.close()
    target.file.close()
PY

"${PYTHON}" --version
"${PYTHON}" -m py_compile "${RUNNER}"

echo "VISTA runner syntax check: PASSED"
echo
echo "Starting Experiment 5.3 fold ${FOLD}: $(date)"

srun --cpu-bind=cores \
  "${PYTHON}" -u "${RUNNER}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment ex5_3 \
    --experiment-label "Experiment 5.3" \
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
echo "VISTA Experiment 5.3 fold ${FOLD} completed successfully"
echo "Completed: $(date)"
echo "===================================================================================================="
