#!/bin/bash
#SBATCH --job-name=envi_ct_DE
#SBATCH --output=logs/envi_ct_DE_%j.out
#SBATCH --error=logs/envi_ct_DE_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"

cd "${ROOT}"

"${RSCRIPT}" \
src/imputation/DE/06_envi_celltype_limma_de.R
