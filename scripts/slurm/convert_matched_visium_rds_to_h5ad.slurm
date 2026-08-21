#!/bin/bash
#SBATCH --job-name=visium_rds_h5ad
#SBATCH --output=logs/visium_rds_h5ad_%j.out
#SBATCH --error=logs/visium_rds_h5ad_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G

set -euo pipefail

PROJECT_ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"

RSCRIPT="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/Rscript"

SCRIPT="${PROJECT_ROOT}/src/visium_branch/03_annotation/02_convert_matched_visium_rds_to_h5ad.R"

export HDF5_USE_FILE_LOCKING=FALSE
export R_LIBS_USER="${PROJECT_ROOT}/software/R_libs"

cd "${PROJECT_ROOT}"

mkdir -p logs

echo "============================================================"
echo "Visium RDS to H5AD conversion"
echo "============================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "Rscript: ${RSCRIPT}"
echo "Script: ${SCRIPT}"
echo "============================================================"

# Check that Rscript exists
if [ ! -x "${RSCRIPT}" ]; then
    echo "ERROR: Rscript not found:"
    echo "${RSCRIPT}"
    exit 1
fi

# Check that conversion script exists
if [ ! -f "${SCRIPT}" ]; then
    echo "ERROR: R conversion script not found:"
    echo "${SCRIPT}"
    exit 1
fi

echo ""
echo "R version:"
"${RSCRIPT}" --version

echo ""
echo "Checking required R packages..."

"${RSCRIPT}" -e '
packages <- c(
    "anndataR",
    "SpatialExperiment",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors",
    "scuttle"
)

for (pkg in packages) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

        stop(
            "Missing required R package: ",
            pkg
        )
    }

    cat(
        pkg,
        " = ",
        as.character(
            packageVersion(pkg)
        ),
        "\n",
        sep = ""
    )
}
'

echo ""
echo "Starting Visium RDS -> H5AD conversion..."
echo ""

"${RSCRIPT}" "${SCRIPT}"

echo ""
echo "============================================================"
echo "Conversion completed successfully"
echo "Finished: $(date)"
echo "============================================================"
