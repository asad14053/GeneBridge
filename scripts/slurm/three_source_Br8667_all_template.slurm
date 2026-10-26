#!/bin/bash
#SBATCH --job-name=Br8667_3srcfig
#SBATCH --output=logs/three_source_Br8667/%j.out
#SBATCH --error=logs/three_source_Br8667/%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
PYTHON="/beegfs/labs/hulab/projects/mjabin/conda_envs/genebridge_scanpy/bin/python"

cd "${ROOT}"

export MPLBACKEND=Agg

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4


echo "============================================================"
echo "BR8667 THREE-SOURCE VALIDATION FIGURES"
echo "============================================================"
echo "A = Visium"
echo "B = Measured Xenium"
echo "C = ENVI OOF"
echo
echo "Generating:"
echo "  A vs B"
echo "  A vs C"
echo "  B vs C"
echo
echo "Each:"
echo "  a = source 1 map"
echo "  b = source 2 map"
echo "  c = normalized-expression correspondence"
echo "  d = Moran's-I correspondence"
echo "  e = patterned raw-cosine distribution"
echo "  f = non-patterned raw-cosine distribution"
echo "============================================================"


"${PYTHON}" \
    src/analysis/18_plot_Br8667_three_source_all_template.py


OUT="${ROOT}/results/analysis/three_source_validation/Br8667/PI_style_pairwise_final"


echo
echo "============================================================"
echo "VERIFY OUTPUTS"
echo "============================================================"


EXPECTED=(
    "fig_Br8667_A_vs_B_PI_template_raw_cosine.png"
    "fig_Br8667_A_vs_C_PI_template_raw_cosine.png"
    "fig_Br8667_B_vs_C_PI_template_raw_cosine.png"
)


for F in "${EXPECTED[@]}"; do

    if [[ ! -s "${OUT}/${F}" ]]; then
        echo "ERROR: missing ${OUT}/${F}"
        exit 1
    fi

    echo "PASS: ${F}"

done


echo
echo "PNG files:"
ls -lh "${OUT}"/*.png


echo
echo "============================================================"
echo "ALL THREE FIGURES PASS"
echo "============================================================"
