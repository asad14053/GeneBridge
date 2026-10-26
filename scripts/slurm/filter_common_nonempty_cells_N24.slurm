#!/bin/bash
#SBATCH --job-name=N24_nonempty
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-24%12
#SBATCH --output=logs/filter_blank_cells_N24/%A_%a.out
#SBATCH --error=logs/filter_blank_cells_N24/%A_%a.err

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

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

DONOR="${DONORS[$((SLURM_ARRAY_TASK_ID-1))]}"

echo "============================================================"
echo "Job ID : ${SLURM_JOB_ID}"
echo "Task   : ${SLURM_ARRAY_TASK_ID}"
echo "Donor  : ${DONOR}"
echo "============================================================"

"${PYTHON}" \
    "${ROOT}/src/imputation/filter_common_nonempty_cells_N24.py" \
    --donor "${DONOR}" \
    --fold-root \
    "${ROOT}/data/processed/imputation_full/gene_folds_200_100"

echo
echo "DONE: ${DONOR}"
