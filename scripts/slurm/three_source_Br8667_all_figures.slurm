#!/bin/bash
#SBATCH --job-name=Br8667_3src
#SBATCH --output=logs/three_source_Br8667/%j.out
#SBATCH --error=logs/three_source_Br8667/%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

echo "============================================================"
echo "BR8667 THREE-SOURCE PI FIGURES"
echo "============================================================"
echo "A = Visium"
echo "B = Measured Xenium"
echo "C = ENVI OOF"
echo
echo "Template:"
echo "a = source 1 map"
echo "b = source 2 map"
echo "c = expression correspondence"
echo "d = Moran's I correspondence"
echo "e = patterned raw-cosine distribution"
echo "f = non-patterned raw-cosine distribution"
echo "============================================================"


# ============================================================
# 1. A vs B + A vs C
# ============================================================

echo
echo ">>> GENERATING A vs B AND A vs C..."

"${PYTHON}" \
    src/analysis/17_plot_Br8667_A_vs_B_A_vs_C_template.py

echo "A vs B + A vs C PASS"


# ============================================================
# 2. B vs C
# ============================================================

echo
echo ">>> GENERATING B vs C..."

"${PYTHON}" \
    src/analysis/16_plot_Br8667_B_vs_C_template.py

echo "B vs C PASS"


# ============================================================
# VERIFY OUTPUTS
# ============================================================

OUT="results/analysis/three_source_validation/Br8667/PI_style_pairwise"

echo
echo "============================================================"
echo "OUTPUT CHECK"
echo "============================================================"

find "${OUT}" \
    -maxdepth 1 \
    -type f \
    -name "fig_Br8667_*template*.png" \
    -printf "%f\n" \
    | sort

N=$(find "${OUT}" \
    -maxdepth 1 \
    -type f \
    -name "fig_Br8667_*template*.png" \
    | wc -l)

echo
echo "Template PNG count: ${N}"

if [ "${N}" -lt 3 ]; then
    echo "ERROR: Expected at least 3 template PNG figures."
    exit 1
fi

echo
echo "============================================================"
echo "ALL THREE BR8667 FIGURES PASS"
echo "============================================================"
