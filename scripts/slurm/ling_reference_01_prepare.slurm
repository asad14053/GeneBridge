#!/bin/bash
#SBATCH --job-name=lingprep
#SBATCH --output=logs/ling_reference_build/lingprep_%j.out
#SBATCH --error=logs/ling_reference_build/lingprep_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

"${PYTHON}" \
scripts/ling_harmonization/06_prepare_ling_array_reference.py
