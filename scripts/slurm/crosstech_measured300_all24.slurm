#!/bin/bash
#SBATCH --job-name=crosstech24
#SBATCH --output=logs/crosstech_all24/%A_%a.out
#SBATCH --error=logs/crosstech_all24/%A_%a.err
#SBATCH --array=0-23%24
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="${ROOT}/../../mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

mapfile -t DONORS < <(find data/metadata/stalign_results_update -maxdepth 1 -type d -name 'Br*_STalign_result' -printf '%f\n' | sed 's/_STalign_result$//' | sort -V)

echo "Total donors detected: ${#DONORS[@]}"

if [[ ${#DONORS[@]} -ne 24 ]]; then
    echo "ERROR: Expected 24 donors, found ${#DONORS[@]}"
    exit 1
fi

BR="${DONORS[$((SLURM_ARRAY_TASK_ID))]}"
export BR

echo "SLURM task: ${SLURM_ARRAY_TASK_ID}"
echo "Running donor: ${BR}"

"${PYTHON}" src/analysis/03_crosstech_measured300.py
