#!/bin/bash
set -euo pipefail

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
CACHE="${ROOT}/data/reference/ling_gene_annotation_cache"

mkdir -p "${CACHE}"

download_ensembl () {
    v="$1"

    url="https://ftp.ensembl.org/pub/release-${v}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${v}.gtf.gz"
    out="${CACHE}/Homo_sapiens.GRCh38.${v}.gtf.gz"

    echo "Ensembl ${v}"
    wget -c -nv -O "${out}" "${url}"
}

download_gencode () {
    v="$1"

    url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${v}/gencode.v${v}.annotation.gtf.gz"
    out="${CACHE}/gencode.v${v}.annotation.gtf.gz"

    echo "GENCODE v${v}"
    wget -c -nv -O "${out}" "${url}"
}


for v in 89 90 88 91 92
do
    download_ensembl "${v}"
done


for v in 26 27 29 30 32 35 44
do
    download_gencode "${v}"
done


echo
echo "===== CACHE ====="
ls -lh "${CACHE}"
