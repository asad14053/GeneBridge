#!/bin/bash
#SBATCH --job-name=envi_DE
#SBATCH --output=logs/envi_DE_%j.out
#SBATCH --error=logs/envi_DE_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"

cd "${ROOT}"

"${RSCRIPT}" \
src/imputation/DE/03_envi_transcriptome_limma_de.R
