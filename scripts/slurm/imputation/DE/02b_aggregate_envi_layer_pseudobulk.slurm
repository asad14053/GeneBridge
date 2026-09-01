#!/bin/bash
#SBATCH --job-name=envi_layer_agg
#SBATCH --output=logs/imputation/DE/envi_layer_pseudobulk/aggregate_%j.out
#SBATCH --error=logs/imputation/DE/envi_layer_pseudobulk/aggregate_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=06:00:00

set -euo pipefail

cd /beegfs/labs/hulab/projects/mjabin/GeneBridge

PYTHON=/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

${PYTHON} \
src/imputation/DE/paper_concordance/02b_aggregate_envi_layer_pseudobulk.py
