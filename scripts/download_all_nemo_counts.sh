#!/bin/bash
set -uo pipefail

BASE="https://data.nemoarchive.org/other/grant/broad/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/counts/"

ROOT="/beegfs/labs/hulab/projects/mjabin/GeneBridge"
DEST="${ROOT}/data/raw/nemo_mccarroll/counts"

MANIFEST="${DEST}/nemo_all_files_urls.txt"
FAILED="${DEST}/nemo_failed_urls.txt"

mkdir -p "${DEST}"
mkdir -p "${ROOT}/logs/nemo_download"

echo "============================================================"
echo "NeMO COMPLETE COUNTS DOWNLOAD"
echo "Destination: ${DEST}"
echo "============================================================"


# =============================================================================
# 1. BUILD MANIFEST
# =============================================================================

if [[ ! -s "${MANIFEST}" ]]; then

    echo
    echo "[1/4] Downloading NeMO directory listing..."

    INDEX=$(mktemp)

    success=0

    for attempt in $(seq 1 10); do

        if wget \
            -q \
            --timeout=120 \
            --tries=1 \
            -O "${INDEX}" \
            "${BASE}"
        then
            success=1
            break
        fi

        sleep_seconds=$((attempt * 20))

        echo "Directory listing retry ${attempt}/10; sleeping ${sleep_seconds}s"

        sleep "${sleep_seconds}"
    done


    if [[ "${success}" -ne 1 ]]; then

        echo "ERROR: Could not retrieve NeMO directory listing."

        rm -f "${INDEX}"

        exit 1
    fi


    python - "${INDEX}" "${BASE}" "${MANIFEST}" <<'PY'
import sys
from html.parser import HTMLParser
from urllib.parse import urljoin

index_file, base, output = sys.argv[1:4]


class Parser(HTMLParser):

    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):

        if tag.lower() != "a":
            return

        href = dict(attrs).get("href")

        if href:
            self.links.append(href)


p = Parser()

with open(
    index_file,
    encoding="utf-8",
    errors="replace"
) as fh:

    p.feed(
        fh.read()
    )


urls = []

for href in p.links:

    if href in (
        "../",
        "./"
    ):
        continue

    if href.startswith("?"):
        continue

    if href.endswith("/"):
        continue

    urls.append(
        urljoin(
            base,
            href
        )
    )


urls = sorted(
    set(
        urls
    )
)


with open(
    output,
    "w"
) as fh:

    for url in urls:
        fh.write(
            url
            +
            "\n"
        )


print(
    "Files discovered:",
    len(urls)
)

print(
    "Manifest:",
    output
)
PY

    rm -f "${INDEX}"

else

    echo
    echo "[1/4] Existing manifest found; reusing:"
    echo "${MANIFEST}"

fi


echo
echo "Manifest entries:"
wc -l "${MANIFEST}"


# =============================================================================
# 2. DOWNLOAD FUNCTION
# =============================================================================

download_one() {

    url="$1"

    filename=$(basename "${url}")

    echo "START  ${filename}"


    for attempt in $(seq 1 12); do

        if (
            cd "${DEST}" &&
            wget \
                -c \
                -nv \
                --timeout=120 \
                --tries=1 \
                "${url}"
        )
        then

            echo "DONE   ${filename}"

            return 0
        fi


        sleep_seconds=$((attempt * 30))

        echo \
            "RETRY ${attempt}/12 ${filename}; sleeping ${sleep_seconds}s"

        sleep "${sleep_seconds}"

    done


    echo "${url}" >> "${FAILED}"

    echo "FAILED ${filename}"

    # Do not kill the whole xargs download.
    return 0
}


export -f download_one

export DEST
export FAILED


# =============================================================================
# 3. FIRST PASS — TWO PARALLEL DOWNLOADS
# =============================================================================

echo
echo "[2/4] Starting downloads with 2 workers..."
echo

: > "${FAILED}"


cat "${MANIFEST}" \
| xargs -P 2 -n 1 bash -c \
    'download_one "$1"' _


# =============================================================================
# 4. RETRY FAILED FILES SEQUENTIALLY
# =============================================================================

echo
echo "[3/4] Parallel pass finished."


if [[ -s "${FAILED}" ]]; then

    FAILED_PASS1="${FAILED}.pass1"

    cp "${FAILED}" "${FAILED_PASS1}"

    : > "${FAILED}"


    echo
    echo "Retrying failures sequentially..."


    while IFS= read -r url; do

        [[ -z "${url}" ]] && continue

        download_one "${url}"

    done < "${FAILED_PASS1}"

fi


# =============================================================================
# SUMMARY
# =============================================================================

echo
echo "[4/4] DOWNLOAD SUMMARY"
echo "============================================================"

echo
echo "Manifest entries:"
wc -l "${MANIFEST}"

echo
echo "H5 files downloaded:"
find "${DEST}" \
    -maxdepth 1 \
    -type f \
    -name "*.h5" \
    | wc -l

echo
echo "All downloaded files:"
find "${DEST}" \
    -maxdepth 1 \
    -type f \
    ! -name "nemo_all_files_urls.txt" \
    ! -name "nemo_failed_urls.txt*" \
    | wc -l

echo
echo "Directory size:"
du -sh "${DEST}"

echo
echo "Remaining failed URLs:"

if [[ -s "${FAILED}" ]]; then

    wc -l "${FAILED}"

    cat "${FAILED}"

else

    echo "0"

fi

echo
echo "============================================================"
echo "DOWNLOAD PASS COMPLETE"
echo "============================================================"
