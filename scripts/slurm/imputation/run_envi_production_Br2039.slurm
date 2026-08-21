#!/bin/bash
#SBATCH --job-name=envi_prod_Br2039
#SBATCH --output=logs/envi_prod_Br2039_%j.out
#SBATCH --error=logs/envi_prod_Br2039_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/Imputation/bin/python"

cd "${ROOT}"

mkdir -p logs
mkdir -p data/processed/imputation_full/ex2_scz/envi/Br2039

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

echo "===================================================================================================="
echo "ENVI PRODUCTION FULL-TRANSCRIPTOME TEST"
echo "===================================================================================================="
echo "Experiment      : Ex2 SCZ"
echo "Donor           : Br2039"
echo "Reference       : All-10 Huuki snRNA"
echo "Measured genes  : 300 Xenium genes"
echo "Training steps  : 16000"
echo "Batch size      : 128"
echo "num_HVG         : 2048"
echo "num_cov_genes   : 64"
echo "k_nearest       : 8"
echo "Memory          : 200G"
echo "CPUs            : 16"
echo "Job ID          : ${SLURM_JOB_ID}"
echo "Host            : $(hostname)"
echo "Started         : $(date)"
echo "===================================================================================================="

"${PYTHON}" -u \
src/imputation/ENVI/run_envi_production.py \
    --target \
data/processed/imputation_full/ex2_scz/targets/spatial_data_xenium_Br2039_ex2_scz.h5ad \
    --reference \
data/processed/snrnaseq/sce_DLPFC_annotated/huuki_snrna_reference_full_allgenes.h5ad \
    --output-dir \
data/processed/imputation_full/ex2_scz/envi/Br2039 \
    --donor Br2039 \
    --experiment ex2_scz \
    --training-steps 16000 \
    --batch-size 128 \
    --num-hvg 2048 \
    --num-cov-genes 64 \
    --k-nearest 8 \
    --spatial-dist pois \
    --sc-dist nb \
    --stable-eps 1e-6 \
    --seed 0

echo "===================================================================================================="
echo "Finished: $(date)"
echo "===================================================================================================="
