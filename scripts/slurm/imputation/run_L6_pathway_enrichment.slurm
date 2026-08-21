#!/bin/bash
#SBATCH --job-name=L6_enrich
#SBATCH --output=logs/L6_enrich_%j.out
#SBATCH --error=logs/L6_enrich_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
R="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"

cd "${ROOT}"

"${R}" \
src/imputation/DE/09_L6_pathway_enrichment.R
