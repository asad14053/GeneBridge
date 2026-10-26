#!/bin/bash

#SBATCH --job-name=xenv24
#SBATCH --output=logs/download_xenium/%A_%a.out
#SBATCH --error=logs/download_xenium/%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --array=1-24%24

set -euo pipefail


# =============================================================================
# PROJECT
# =============================================================================

PROJECT_ROOT=/beegfs/labs/hulab/projects/mjabin/GeneBridge

PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

cd "${PROJECT_ROOT}"

mkdir -p logs/download_xenium


# =============================================================================
# DONORS
# =============================================================================

DONORS=(
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

DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID - 1))]}"


# =============================================================================
# HEADER
# =============================================================================

echo "===================================================================================================="
echo "GENEBRIDGE: 24-DONOR XENIUM + VISIUM DOWNLOAD"
echo "===================================================================================================="
echo "Job ID       : ${SLURM_JOB_ID}"
echo "Array task   : ${SLURM_ARRAY_TASK_ID}"
echo "Donor        : ${DONOR}"
echo "Node         : $(hostname)"
echo "CPUs         : ${SLURM_CPUS_PER_TASK}"
echo "Memory       : 200G"
echo "Started      : $(date)"
echo "===================================================================================================="


# =============================================================================
# ARGUMENTS
# =============================================================================

ARGS=(
    --donor "${DONOR}"
    --modality both
    --include-heavy
    --max-retries 10
)


# =============================================================================
# Br6432 WAS PREVIOUSLY CORRUPTED.
#
# Force this donor only.
# Existing Br6432 files will be quarantined and downloaded again.
# =============================================================================

if [[ "${DONOR}" == "Br6432" ]]; then

    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "Br6432 is marked as previously corrupted."
    echo "FORCING clean redownload."
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo

    ARGS+=(--force)

fi


# =============================================================================
# RUN
# =============================================================================

"${PYTHON}" -u \
    src/data/download_24_common_xenium_visium.py \
    "${ARGS[@]}"


STATUS=$?


echo
echo "===================================================================================================="
echo "FINISHED"
echo "===================================================================================================="
echo "Donor       : ${DONOR}"
echo "Exit status : ${STATUS}"
echo "Finished    : $(date)"
echo "===================================================================================================="

exit "${STATUS}"
