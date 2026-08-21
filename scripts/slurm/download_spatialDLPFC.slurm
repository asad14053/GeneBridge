#!/bin/bash
#SBATCH --job-name=spatialDLPFC
#SBATCH --output=/beegfs/labs/hulab/projects/mjabin/GeneBridge/logs/spatialDLPFC_%j.out
#SBATCH --error=/beegfs/labs/hulab/projects/mjabin/GeneBridge/logs/spatialDLPFC_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

set -euo pipefail

PROJECT=/beegfs/labs/hulab/projects/mjabin/GeneBridge

cd "$PROJECT"
mkdir -p logs

echo "========================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Working directory: $(pwd)"
echo "========================================"

# Use the R environment directly
R_SCRIPT=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript

"$R_SCRIPT" src/data/download_spatialDLPFC.R

echo "========================================"
echo "Finished: $(date)"
echo "========================================"