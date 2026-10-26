#!/bin/bash
#SBATCH --job-name=gimvi_oof
#SBATCH --output=logs/gimvi_oof/gimvi_oof_%A_%a.out
#SBATCH --error=logs/gimvi_oof/gimvi_oof_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-28%28

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

SCRIPT="${PROJECT_ROOT}/src/imputation/gimVI/01_stitch_oof_predictions.py"

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

TASK="${SLURM_ARRAY_TASK_ID}"

# ----------------------------------------------------------------------
# 1-24 : Ex-5.1
# 25   : Ex-5 Br6432
# 26   : Ex-5 Br8667
# 27   : Ex-5.3 Br6432
# 28   : Ex-5.3 Br8667
# ----------------------------------------------------------------------

if (( TASK >= 1 && TASK <= 24 )); then

    DONOR="${DONORS_N24[$((TASK - 1))]}"
    EXPERIMENT="ex5_1"

elif (( TASK == 25 )); then

    DONOR="Br6432"
    EXPERIMENT="ex5"

elif (( TASK == 26 )); then

    DONOR="Br8667"
    EXPERIMENT="ex5"

elif (( TASK == 27 )); then

    DONOR="Br6432"
    EXPERIMENT="ex5_3"

elif (( TASK == 28 )); then

    DONOR="Br8667"
    EXPERIMENT="ex5_3"

else

    echo "ERROR: invalid task ${TASK}"
    exit 1

fi

cd "${PROJECT_ROOT}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PYTHONUNBUFFERED=1

echo "===================================================================================================="
echo "gimVI OOF 300-GENE STITCH"
echo "===================================================================================================="
echo "Job        : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Donor      : ${DONOR}"
echo "Experiment : ${EXPERIMENT}"
echo "Node       : $(hostname)"
echo "Started    : $(date)"
echo "===================================================================================================="

"${PYTHON}" -u "${SCRIPT}" \
    --project-root "${PROJECT_ROOT}" \
    --donor "${DONOR}" \
    --experiment "${EXPERIMENT}"

echo
echo "===================================================================================================="
echo "SLURM STATUS: PASS"
echo "Donor      : ${DONOR}"
echo "Experiment : ${EXPERIMENT}"
echo "Completed  : $(date)"
echo "===================================================================================================="
