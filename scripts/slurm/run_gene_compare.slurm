#!/bin/bash
#SBATCH --job-name=gene_compare
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail
mkdir -p logs

source ~/anaconda3/etc/profile.d/conda.sh
conda activate /scratch/mjabin/conda_envs/Imputation

cd /users/mjabin/projects/GeneBridge
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

GENE="${1:-ACTA2}"
RUN_MODE="${2:-full}"

python src/analysis/gene_compare_one_validation_gene.py \
  --project-root /users/mjabin/projects/GeneBridge \
  --gene "$GENE" \
  --run-mode "$RUN_MODE"
