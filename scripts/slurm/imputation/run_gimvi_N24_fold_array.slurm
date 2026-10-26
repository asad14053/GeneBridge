#!/bin/bash
#SBATCH --output=logs/gimvi_N24/%x_%A_%a.out
#SBATCH --error=logs/gimvi_N24/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail


# ==========================================================================================
# Paths
# ==========================================================================================

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/gimVI/run_gimvi_fold.py"

FOLD_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/gene_folds_200_100"

REF_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/snrnaseq_target_references"

FULL_REFERENCE="${PROJECT_ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"

OUTPUT_BASE="${PROJECT_ROOT}/outputs/imputation_full/GIMVI_N24"


# ==========================================================================================
# Experiment
# ==========================================================================================

EXPERIMENT="${1:?Need experiment: ex5, ex5_1, or ex5_3}"
EXPERIMENT_LABEL="${2:?Need experiment label}"


# ==========================================================================================
# Donors
# ==========================================================================================

DONORS_N24=(
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


# ==========================================================================================
# Array task -> donor × fold
#
# Tasks:
#   1 = donor1 fold1
#   2 = donor1 fold2
#   3 = donor1 fold3
#   4 = donor2 fold1
#   ...
# ==========================================================================================

TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

DONOR_INDEX=$((TASK_INDEX / 3))

FOLD=$((TASK_INDEX % 3 + 1))


# ==========================================================================================
# Experiment-specific reference
# ==========================================================================================

case "${EXPERIMENT}" in


    # ======================================================================================
    # Experiment 5
    #
    # Target donor appears in snRNA reference.
    #
    # Br6432 uses the SAME pseudocount-rescued reference used for ENVI.
    # ======================================================================================

    ex5)

        if (( DONOR_INDEX >= 2 )); then
            echo "ERROR: Ex-5 supports tasks 1-6 only."
            exit 1
        fi

        DONOR="${MATCHED_DONORS[$DONOR_INDEX]}"

        if [[ "${DONOR}" == "Br6432" ]]; then

            REFERENCE="${REF_ROOT}/Br6432/seq_data_huuki_snrna_Br6432_pseudocount.h5ad"

        else

            REFERENCE="${REF_ROOT}/Br8667/seq_data_huuki_snrna_Br8667.h5ad"

        fi
        ;;


    # ======================================================================================
    # Experiment 5.1
    #
    # All N24 Xenium donors.
    # Full 10-brain Huuki snRNA reference.
    # ======================================================================================

    ex5_1)

        if (( DONOR_INDEX >= 24 )); then
            echo "ERROR: Ex-5.1 supports tasks 1-72 only."
            exit 1
        fi

        DONOR="${DONORS_N24[$DONOR_INDEX]}"

        REFERENCE="${FULL_REFERENCE}"
        ;;


    # ======================================================================================
    # Experiment 5.3
    #
    # Target donor excluded from the 9-brain snRNA reference.
    # ======================================================================================

    ex5_3)

        if (( DONOR_INDEX >= 2 )); then
            echo "ERROR: Ex-5.3 supports tasks 1-6 only."
            exit 1
        fi

        DONOR="${MATCHED_DONORS[$DONOR_INDEX]}"

        REFERENCE="${REF_ROOT}/${DONOR}/seq_data_huuki_snrna_9brains_excluding_${DONOR}.h5ad"
        ;;


    *)

        echo "ERROR: Unsupported experiment: ${EXPERIMENT}"
        exit 1
        ;;

esac


# ==========================================================================================
# SAME fold files used for ENVI and VISTA
# ==========================================================================================

FOLD_DIR="${FOLD_ROOT}/${DONOR}"

OBSERVED="${FOLD_DIR}/fold_${FOLD}_observed_genes.h5ad"

HELDOUT="${FOLD_DIR}/fold_${FOLD}_heldout_genes.h5ad"

OUTPUT_ROOT="${OUTPUT_BASE}/${DONOR}"


# ==========================================================================================
# Validate inputs
# ==========================================================================================

for FILE in \
    "${REFERENCE}" \
    "${OBSERVED}" \
    "${HELDOUT}"
do

    if [[ ! -f "${FILE}" ]]; then

        echo "ERROR: Missing required input:"
        echo "${FILE}"

        exit 1

    fi

done


if [[ ! -x "${PYTHON}" ]]; then

    echo "ERROR: Python not found:"
    echo "${PYTHON}"

    exit 1

fi


# ==========================================================================================
# CPU configuration
# ==========================================================================================

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg


cd "${PROJECT_ROOT}"


# ==========================================================================================
# Run information
# ==========================================================================================

echo "===================================================================================================="
echo "gimVI N24 FULLY PARALLEL BENCHMARK"
echo "===================================================================================================="

echo "Job ID          : ${SLURM_JOB_ID}"
echo "Array task      : ${SLURM_ARRAY_TASK_ID}"
echo "Node            : $(hostname)"

echo
echo "Experiment      : ${EXPERIMENT}"
echo "Experiment label: ${EXPERIMENT_LABEL}"

echo
echo "Target donor    : ${DONOR}"
echo "Fold            : ${FOLD}"

echo
echo "Reference       : ${REFERENCE}"
echo "Observed        : ${OBSERVED}"
echo "Heldout         : ${HELDOUT}"

echo
echo "Output root     : ${OUTPUT_ROOT}"

echo
echo "Epochs          : 200"
echo "Batch size      : 128"
echo "Latent          : 32"
echo "Train size      : 0.90"
echo "Learning rate   : 1e-4"
echo "Gradient clip   : 5.0"

echo
echo "CPUs            : ${SLURM_CPUS_PER_TASK}"
echo "Memory          : 200G"

echo
echo "Started         : $(date)"

echo "===================================================================================================="


# ==========================================================================================
# ONE donor × ONE fold
# ==========================================================================================

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
    --train-size 0.90 \
    --learning-rate 1e-4 \
    --gradient-clip-val 5.0 \
    --seed 8667 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --k-spatial 15 \
    --ssim-grid-size 128 \
    --n-clusters 10 \
    --n-pcs 30 \
    --plot-genes 10


echo
echo "===================================================================================================="
echo "FOLD STATUS: PASS"
echo "===================================================================================================="

echo "Donor     : ${DONOR}"
echo "Experiment: ${EXPERIMENT}"
echo "Fold      : ${FOLD}"
echo "Completed : $(date)"

echo
echo "FINAL STATUS: PASS"

