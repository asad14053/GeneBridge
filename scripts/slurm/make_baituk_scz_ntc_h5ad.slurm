#!/bin/bash
#SBATCH --job-name=baituk_scz_ntc_rds
#SBATCH --output=logs/baituk_scz_ntc_rds_%j.out
#SBATCH --error=logs/baituk_scz_ntc_rds_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_r/bin/Rscript"

SCRIPT="${PROJECT_ROOT}/src/snrnaseq_branch/03_batiuk_huuki_gene_overlap/03_make_baituk_scz_ntc_h5ad.R"

export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"
export HDF5_USE_FILE_LOCKING=FALSE

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

cd "${PROJECT_ROOT}"
mkdir -p logs

echo "============================================================"
echo "Create separate Batiuk SCZ and NTC H5AD files"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 200G"
echo "============================================================"

if [[ ! -x "${RSCRIPT}" ]]; then
    echo "ERROR: Rscript not found: ${RSCRIPT}"
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: R script not found: ${SCRIPT}"
    exit 1
fi

"${RSCRIPT}" --version

"${RSCRIPT}" -e '
required <- c(
    "Matrix",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors",
    "anndataR"
)

for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing package: ", pkg)
    }

    cat(
        pkg,
        " = ",
        as.character(packageVersion(pkg)),
        "\n",
        sep = ""
    )
}
'

echo
echo "Starting Batiuk SCZ/NTC H5AD creation..."

"${RSCRIPT}" "${SCRIPT}"

echo
echo "============================================================"
echo "Completed successfully"
echo "Finished: $(date)"
echo "============================================================"
