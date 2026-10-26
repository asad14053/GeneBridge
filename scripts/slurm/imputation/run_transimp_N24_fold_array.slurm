#!/bin/bash
#SBATCH --output=logs/transimp_N24/%x_%A_%a.out
#SBATCH --error=logs/transimp_N24/%x_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/TransImp/run_transimpspa_fold.py"

FOLD_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/gene_folds_200_100"

REF_ROOT="${PROJECT_ROOT}/data/processed/imputation_full/snrnaseq_target_references"

FULL_REFERENCE="${PROJECT_ROOT}/data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad"

OUTPUT_BASE="${PROJECT_ROOT}/outputs/imputation_full/TRANSIMP_N24"


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
# One SLURM task = one donor × one fold
# ==========================================================================================

TASK_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

DONOR_INDEX=$((TASK_INDEX / 3))

FOLD=$((TASK_INDEX % 3 + 1))


# ==========================================================================================
# Experiment-specific reference
# ==========================================================================================

case "${EXPERIMENT}" in

    # --------------------------------------------------------------------------------------
    # Ex-5
    # same-donor snRNA reference
    # --------------------------------------------------------------------------------------

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


    # --------------------------------------------------------------------------------------
    # Ex-5.1
    # all 24 Xenium donors + full 10-brain reference
    # --------------------------------------------------------------------------------------

    ex5_1)

        if (( DONOR_INDEX >= 24 )); then
            echo "ERROR: Ex-5.1 supports tasks 1-72 only."
            exit 1
        fi

        DONOR="${DONORS_N24[$DONOR_INDEX]}"

        REFERENCE="${FULL_REFERENCE}"
        ;;


    # --------------------------------------------------------------------------------------
    # Ex-5.3
    # target donor excluded from snRNA reference
    # --------------------------------------------------------------------------------------

    ex5_3)

        if (( DONOR_INDEX >= 2 )); then
            echo "ERROR: Ex-5.3 supports tasks 1-6 only."
            exit 1
        fi

        DONOR="${MATCHED_DONORS[$DONOR_INDEX]}"

        REFERENCE="${REF_ROOT}/${DONOR}/seq_data_huuki_snrna_9brains_excluding_${DONOR}.h5ad"
        ;;


    *)

        echo "ERROR: unsupported experiment: ${EXPERIMENT}"
        exit 1
        ;;

esac


# ==========================================================================================
# SAME fold files used by every other model
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
        echo "ERROR: missing required input:"
        echo "${FILE}"
        exit 1
    fi

done


if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: Python missing:"
    echo "${PYTHON}"
    exit 1
fi


if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: TransImpSpa runner missing:"
    echo "${SCRIPT}"
    exit 1
fi


# ==========================================================================================
# CPU benchmark
#
# The Python runner automatically chooses CUDA if visible.
# Hide CUDA so TransImpSpa uses CPU, consistent with the benchmark runs.
# ==========================================================================================

export CUDA_VISIBLE_DEVICES=""

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg


cd "${PROJECT_ROOT}"


echo "===================================================================================================="
echo "TransImpSpa N24 FULLY PARALLEL BENCHMARK"
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
echo "Mapping mode    : lowrank"
echo "Mapping lowdim  : 256"
echo "Spatial NN      : 6"
echo "Epochs          : 2000"
echo "Learning rate   : 0.01"
echo "Weight decay    : 0.01"
echo "Clip max        : 10"
echo "Spatial weight  : 1.0"

echo
echo "CPUs            : ${SLURM_CPUS_PER_TASK}"
echo "Memory          : 200G"

echo
echo "Started         : $(date)"

echo "===================================================================================================="


"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --experiment "${EXPERIMENT}" \
    --experiment-label "${EXPERIMENT_LABEL}" \
    --fold "${FOLD}" \
    --reference "${REFERENCE}" \
    --fold-dir "${FOLD_DIR}" \
    --output-root "${OUTPUT_ROOT}" \
    --mapping-lowdim 256 \
    --epochs 2000 \
    --lr 0.01 \
    --weight-decay 0.01 \
    --clip-max 10.0 \
    --wt-spa 1.0 \
    --spatial-neighbors 6 \
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

