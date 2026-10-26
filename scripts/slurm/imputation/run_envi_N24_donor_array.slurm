#!/bin/bash
#SBATCH --job-name=envi_N24
#SBATCH --output=logs/envi_N24/%x_%A_%a.out
#SBATCH --error=logs/envi_N24/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail


# ==================================================================================================
# Usage
#
# Ex-5.1, all 24 donors:
#   sbatch --array=1-24%24 \
#     --job-name=envi_ex5_1_N24 \
#     run_envi_N24_donor_array.slurm ex5_1 "Experiment 5.1"
#
# Ex-5, Br6432 + Br8667:
#   sbatch --array=1-2%2 \
#     --job-name=envi_ex5_matched \
#     run_envi_N24_donor_array.slurm ex5 "Experiment 5"
#
# Ex-5.3, Br6432 + Br8667:
#   sbatch --array=1-2%2 \
#     --job-name=envi_ex5_3_matched \
#     run_envi_N24_donor_array.slurm ex5_3 "Experiment 5.3"
# ==================================================================================================


if [[ $# -lt 2 ]]; then
    echo "Usage:"
    echo "  sbatch --array=... $0 EXPERIMENT EXPERIMENT_LABEL"
    exit 1
fi


PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

SCRIPT="${PROJECT_ROOT}/src/imputation/ENVI/run_envi_fold.py"

FOLD_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/gene_folds_200_100"

OUTPUT_ROOT_BASE="${PROJECT_ROOT}/outputs/imputation_full/ENVI_N24"

REF_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/snrnaseq_target_references"

FULL_REFERENCE="${PROJECT_ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"


EXPERIMENT="$1"
EXPERIMENT_LABEL="$2"


# ==================================================================================================
# DONOR SET + REFERENCE
# ==================================================================================================

ALL_DONORS=(
    Br1113
    Br1139
    Br2039
    Br2421
    Br2719
    Br5314
    Br5373
    Br5400
    Br5436
    Br5588
    Br5590
    Br5622
    Br5639
    Br5746
    Br5931
    Br5973
    Br6032
    Br6389
    Br6432
    Br6437
    Br6496
    Br8433
    Br8667
    Br8772
)

MATCHED_DONORS=(
    Br6432
    Br8667
)


case "${EXPERIMENT}" in

    ex5_1)

        DONORS=("${ALL_DONORS[@]}")

        if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > 24 )); then
            echo "ERROR: Ex-5.1 requires array task 1-24."
            exit 1
        fi

        DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID - 1))]}"

        REFERENCE="${FULL_REFERENCE}"
        ;;


    ex5)

        DONORS=("${MATCHED_DONORS[@]}")

        if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > 2 )); then
            echo "ERROR: Ex-5 requires array task 1-2."
            exit 1
        fi

        DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID - 1))]}"

        if [[ "${DONOR}" == "Br6432" ]]; then
            REFERENCE="${REF_ROOT}/Br6432/seq_data_huuki_snrna_Br6432_pseudocount.h5ad"
        else
            REFERENCE="${REF_ROOT}/${DONOR}/seq_data_huuki_snrna_${DONOR}.h5ad"
        fi
        ;;


    ex5_3)

        DONORS=("${MATCHED_DONORS[@]}")

        if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > 2 )); then
            echo "ERROR: Ex-5.3 requires array task 1-2."
            exit 1
        fi

        DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID - 1))]}"

        REFERENCE="${REF_ROOT}/${DONOR}/seq_data_huuki_snrna_9brains_excluding_${DONOR}.h5ad"
        ;;


    *)

        echo "ERROR: unknown experiment: ${EXPERIMENT}"
        echo "Allowed: ex5, ex5_1, ex5_3"
        exit 1
        ;;
esac


FOLD_DIR="${FOLD_ROOT}/${DONOR}"

OUTPUT_ROOT="${OUTPUT_ROOT_BASE}/${DONOR}"


# ==================================================================================================
# INPUT VALIDATION
# ==================================================================================================

if [[ ! -d "${FOLD_DIR}" ]]; then
    echo "ERROR: fold directory missing:"
    echo "${FOLD_DIR}"
    exit 1
fi


if [[ ! -f "${REFERENCE}" ]]; then
    echo "ERROR: reference missing:"
    echo "${REFERENCE}"
    exit 1
fi


for FOLD in 1 2 3; do

    OBSERVED="${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad"
    HELDOUT="${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad"

    if [[ ! -f "${OBSERVED}" ]]; then
        echo "ERROR: missing observed fold:"
        echo "${OBSERVED}"
        exit 1
    fi

    if [[ ! -f "${HELDOUT}" ]]; then
        echo "ERROR: missing heldout fold:"
        echo "${HELDOUT}"
        exit 1
    fi

done


mkdir -p \
    "${PROJECT_ROOT}/logs/envi_N24" \
    "${OUTPUT_ROOT}"


cd "${PROJECT_ROOT}"


# ==================================================================================================
# CPU / JAX
# ==================================================================================================

export JAX_PLATFORMS=cpu
export JAX_PLATFORM_NAME=cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1


# ==================================================================================================
# FIND ENVI ENVIRONMENT
# ==================================================================================================

PYTHON=""

CANDIDATES=(
    "/beegfs/labs/hulab/projects/mjabin/conda_envs/envi/bin/python"
    "/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"
)


for CANDIDATE in "${CANDIDATES[@]}"; do

    if [[ -x "${CANDIDATE}" ]] && \
       "${CANDIDATE}" -c \
       'import scenvi, jax, scanpy, anndata; assert jax.default_backend().lower() == "cpu"' \
       >/dev/null 2>&1
    then
        PYTHON="${CANDIDATE}"
        break
    fi

done


if [[ -z "${PYTHON}" ]]; then

    echo "ERROR: No working ENVI Python environment found."

    printf 'Checked:\n'
    printf '  %s\n' "${CANDIDATES[@]}"

    exit 1
fi


# ==================================================================================================
# JOB CONFIGURATION
# ==================================================================================================

echo "===================================================================================================="
echo "ENVI N24 — DONOR-PARALLEL / FOLD-SEQUENTIAL"
echo "===================================================================================================="

echo "Job ID          : ${SLURM_JOB_ID}"
echo "Array task      : ${SLURM_ARRAY_TASK_ID}"
echo "Node            : $(hostname)"
echo "Started         : $(date)"

echo
echo "Target donor    : ${DONOR}"
echo "Experiment      : ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"

echo
echo "Reference:"
echo "${REFERENCE}"

echo
echo "Fold directory:"
echo "${FOLD_DIR}"

echo
echo "Output root:"
echo "${OUTPUT_ROOT}"

echo
echo "Python          : ${PYTHON}"

echo
echo "Training steps  : 16000"
echo "Batch size      : 128"
echo "num_HVG         : 2048"
echo "num_cov_genes   : 64"
echo "k_nearest       : 8"
echo "covet_batch_size: 256"
echo "spatial_dist    : pois"
echo "sc_dist         : nb"
echo "stable_eps      : 1e-6"

echo
echo "CPUs            : ${SLURM_CPUS_PER_TASK}"
echo "Memory          : 200G"

echo
echo "Execution policy:"
echo "  donors = parallel SLURM array"
echo "  folds  = sequential 1 -> 2 -> 3"

echo "===================================================================================================="


# ==================================================================================================
# RUN ALL THREE FOLDS IN PARALLEL
# ==================================================================================================

PIDS=()
FOLDS=(1 2 3)

for FOLD in "${FOLDS[@]}"; do

    OBSERVED="${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad"
    HELDOUT="${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad"

    echo
    echo "####################################################################################################"
    echo "# ${DONOR} | ${EXPERIMENT} | FOLD ${FOLD}"
    echo "####################################################################################################"
    echo "Starting fold    : ${FOLD}"
    echo "Started          : $(date)"
    echo "Observed target  : ${OBSERVED}"
    echo "Held-out truth   : ${HELDOUT}"
    echo "Reference        : ${REFERENCE}"

    (
        "${PYTHON}" -u "${SCRIPT}" \
            --project-root "${PROJECT_ROOT}" \
            --experiment "${EXPERIMENT}" \
            --experiment-label "${EXPERIMENT_LABEL}" \
            --fold "${FOLD}" \
            --reference "${REFERENCE}" \
            --fold-dir "${FOLD_DIR}" \
            --output-root "${OUTPUT_ROOT}" \
            --run-mode full \
            --full-training-steps 16000 \
            --batch-size 128 \
            --num-hvg 2048 \
            --num-cov-genes 64 \
            --k-nearest 8 \
            --covet-batch-size 256 \
            --spatial-dist pois \
            --sc-dist nb \
            --stable-eps 1e-6 \
            --seed 8667 \
            --cpus "${SLURM_CPUS_PER_TASK}" \
            --k-spatial 15 \
            --ssim-grid-size 128 \
            --n-clusters 10 \
            --n-pcs 30 \
            --plot-genes 10

        echo
        echo "----------------------------------------------------------------------------------------------------"
        echo "FOLD STATUS: PASS"
        echo "Donor     : ${DONOR}"
        echo "Experiment: ${EXPERIMENT}"
        echo "Fold      : ${FOLD}"
        echo "Completed : $(date)"
        echo "----------------------------------------------------------------------------------------------------"

    ) &

    PIDS+=("$!")

done


# Wait for all three folds.
FAILED=0

for i in "${!PIDS[@]}"; do

    PID="${PIDS[$i]}"
    FOLD="${FOLDS[$i]}"

    if wait "${PID}"; then
        echo "Fold ${FOLD} process completed successfully."
    else
        echo "ERROR: Fold ${FOLD} failed."
        FAILED=1
    fi

done


if [[ "${FAILED}" -ne 0 ]]; then
    echo "ERROR: One or more folds failed."
    exit 1
fi


# ==================================================================================================
# FINAL
# ==================================================================================================

echo
echo "===================================================================================================="
echo "FINAL STATUS: PASS"
echo "===================================================================================================="

echo "Donor      : ${DONOR}"
echo "Experiment : ${EXPERIMENT}"
echo "Folds      : 1, 2, 3"
echo "Completed  : $(date)"

echo "===================================================================================================="
