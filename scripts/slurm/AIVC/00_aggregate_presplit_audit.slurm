#!/bin/bash
#SBATCH --job-name=aivc_aggr
#SBATCH --output=logs/aivc_presplit_aggregate_%j.out
#SBATCH --error=logs/aivc_presplit_aggregate_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

"${PYTHON}" -u \
    src/AIVC/data_prep/00_aggregate_presplit_audit.py
