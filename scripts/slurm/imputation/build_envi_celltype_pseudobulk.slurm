#!/bin/bash
#SBATCH --job-name=envi_ct_pb
#SBATCH --output=logs/envi_ct_pb_%j.out
#SBATCH --error=logs/envi_ct_pb_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

"${PYTHON}" -u \
src/imputation/DE/05_build_envi_celltype_pseudobulk.py
