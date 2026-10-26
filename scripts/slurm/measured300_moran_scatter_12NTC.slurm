#!/bin/bash
#SBATCH --job-name=m300Plot
#SBATCH --output=logs/measured300_moran_scatter_12NTC/%A_%a.out
#SBATCH --error=logs/measured300_moran_scatter_12NTC/%A_%a.err
#SBATCH --array=0-11%12
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

DONORS=(
    Br1113
    Br2719
    Br5314
    Br5400
    Br5436
    Br5622
    Br5639
    Br5931
    Br6389
    Br6432
    Br8433
    Br8667
)

export BR="${DONORS[$SLURM_ARRAY_TASK_ID]}"

echo "========================================"
echo "Measured-300 Moran scatter"
echo "BR=${BR}"
echo "TASK=${SLURM_ARRAY_TASK_ID}"
echo "========================================"

"${PYTHON}" \
    src/analysis/13_plot_measured300_patterned_scatter.py

echo
echo "${BR} PASS"
