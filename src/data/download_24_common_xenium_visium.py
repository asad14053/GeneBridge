#!/usr/bin/env python3

"""
Robust downloader for the 24 matched Xenium + Visium DLPFC donors.

Data:
    Xenium: GSE307404
    Visium: GSE307403

Key design:
    - Does NOT scrape GSM HTML pages.
    - Reads GEO series-level filelist.txt.
    - Uses exact file sizes from GEO filelist.txt.
    - Supports one donor per SLURM array task.
    - Restart-safe.
    - Resumes .part files using HTTP Range.
    - Detects incomplete/corrupt files by expected byte size.
    - Can force-redownload known bad files, e.g. Br6432.
    - Writes a separate manifest for each donor to avoid array-task collisions.

Examples:

    # One donor, light Xenium + standard Visium
    python src/data/download_24_common_xenium_visium.py \
        --donor Br1113

    # One donor including huge Xenium files
    python src/data/download_24_common_xenium_visium.py \
        --donor Br1113 \
        --include-heavy

    # Known previously corrupted donor
    python src/data/download_24_common_xenium_visium.py \
        --donor Br6432 \
        --include-heavy \
        --force

    # Dry run
    python src/data/download_24_common_xenium_visium.py \
        --donor Br1113 \
        --include-heavy \
        --dry-run
"""

from __future__ import annotations

from pathlib import Path
from urllib.parse import urlencode
import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
import time

import requests


# =============================================================================
# PATHS
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "data" / "raw"

XENIUM_DIR = RAW_DIR / "xenium"
VISIUM_DIR = RAW_DIR / "visium"

META_DIR = PROJECT_ROOT / "data" / "metadata"
MANIFEST_DIR = META_DIR / "download_manifests_24_common"

XENIUM_DIR.mkdir(parents=True, exist_ok=True)
VISIUM_DIR.mkdir(parents=True, exist_ok=True)
META_DIR.mkdir(parents=True, exist_ok=True)
MANIFEST_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# GEO DATASETS
# =============================================================================

XENIUM_GSE = "GSE307404"
VISIUM_GSE = "GSE307403"

XENIUM_FILELIST_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE307nnn/GSE307404/suppl/filelist.txt"
)

VISIUM_FILELIST_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE307nnn/GSE307403/suppl/filelist.txt"
)

GEO_DOWNLOAD_ENDPOINT = "https://www.ncbi.nlm.nih.gov/geo/download/"


# =============================================================================
# 24 MATCHED DONORS
# =============================================================================

PAIRS = [
    {"patient_index": "P1",  "patient_id": "Br2039", "xenium_gsm": "GSM9223468", "visium_gsm": "GSM9223410"},
    {"patient_index": "P2",  "patient_id": "Br2719", "xenium_gsm": "GSM9223469", "visium_gsm": "GSM9223405"},
    {"patient_index": "P3",  "patient_id": "Br6432", "xenium_gsm": "GSM9223470", "visium_gsm": "GSM9223415"},
    {"patient_index": "P4",  "patient_id": "Br1113", "xenium_gsm": "GSM9223471", "visium_gsm": "GSM9223416"},
    {"patient_index": "P5",  "patient_id": "Br5373", "xenium_gsm": "GSM9223472", "visium_gsm": "GSM9223412"},
    {"patient_index": "P6",  "patient_id": "Br5590", "xenium_gsm": "GSM9223473", "visium_gsm": "GSM9223420"},
    {"patient_index": "P7",  "patient_id": "Br5400", "xenium_gsm": "GSM9223474", "visium_gsm": "GSM9223424"},
    {"patient_index": "P8",  "patient_id": "Br5622", "xenium_gsm": "GSM9223475", "visium_gsm": "GSM9223432"},
    {"patient_index": "P9",  "patient_id": "Br6437", "xenium_gsm": "GSM9223476", "visium_gsm": "GSM9223429"},
    {"patient_index": "P10", "patient_id": "Br5314", "xenium_gsm": "GSM9223477", "visium_gsm": "GSM9223438"},
    {"patient_index": "P11", "patient_id": "Br5588", "xenium_gsm": "GSM9223478", "visium_gsm": "GSM9223434"},
    {"patient_index": "P12", "patient_id": "Br5746", "xenium_gsm": "GSM9223479", "visium_gsm": "GSM9223436"},
    {"patient_index": "P13", "patient_id": "Br5639", "xenium_gsm": "GSM9223480", "visium_gsm": "GSM9223439"},
    {"patient_index": "P14", "patient_id": "Br8433", "xenium_gsm": "GSM9223481", "visium_gsm": "GSM9223454"},
    {"patient_index": "P15", "patient_id": "Br8772", "xenium_gsm": "GSM9223482", "visium_gsm": "GSM9223443"},
    {"patient_index": "P16", "patient_id": "Br1139", "xenium_gsm": "GSM9223483", "visium_gsm": "GSM9223444"},
    {"patient_index": "P17", "patient_id": "Br5973", "xenium_gsm": "GSM9223484", "visium_gsm": "GSM9223459"},
    {"patient_index": "P18", "patient_id": "Br8667", "xenium_gsm": "GSM9223485", "visium_gsm": "GSM9223455"},
    {"patient_index": "P19", "patient_id": "Br5436", "xenium_gsm": "GSM9223486", "visium_gsm": "GSM9223431"},
    {"patient_index": "P20", "patient_id": "Br5931", "xenium_gsm": "GSM9223487", "visium_gsm": "GSM9223448"},
    {"patient_index": "P21", "patient_id": "Br6496", "xenium_gsm": "GSM9223488", "visium_gsm": "GSM9223460"},
    {"patient_index": "P22", "patient_id": "Br2421", "xenium_gsm": "GSM9223489", "visium_gsm": "GSM9223467"},
    {"patient_index": "P23", "patient_id": "Br6032", "xenium_gsm": "GSM9223490", "visium_gsm": "GSM9223435"},
    {"patient_index": "P24", "patient_id": "Br6389", "xenium_gsm": "GSM9223491", "visium_gsm": "GSM9223447"},
]

PAIR_BY_DONOR = {x["patient_id"]: x for x in PAIRS}


# =============================================================================
# FILE SELECTION
# =============================================================================

XENIUM_LIGHT_SUFFIXES = (
    "cell_feature_matrix.h5",
    "cells.parquet.gz",
    "cell_boundaries.csv.gz",
    "nucleus_boundaries.csv.gz",
)

XENIUM_HEAVY_SUFFIXES = (
    "transcripts.zarr.zip",
    "transcripts.csv.gz",
    "morphology.ome.tif.gz",
    "morphology_focus.ome.tif.gz",
)

VISIUM_ANALYSIS_SUFFIXES = (
    "matrix.mtx.gz",
    "features.tsv.gz",
    "barcodes.tsv.gz",
    "tissue_positions.csv.gz",
    "scalefactors_json.json.gz",
    "tissue_hires_image.png.gz",
    "tissue_lowres_image.png.gz",
)

VISIUM_EXTRA_IMAGE_SUFFIXES = (
    "aligned_fiducials.jpg.gz",
    "detected_tissue_image.jpg.gz",
)


# =============================================================================
# HELPERS
# =============================================================================

def section(title: str):
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)


def human_bytes(n: int | None) -> str:
    if n is None:
        return "NA"

    x = float(n)

    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.2f} {unit}"
        x /= 1024

    return str(n)


def fetch_text(url: str, max_retries: int = 10) -> str:
    headers = {
        "User-Agent": "GeneBridge-GEO-downloader/1.0",
        "Accept-Encoding": "identity",
    }

    last_error = None

    for attempt in range(1, max_retries + 1):
        try:
            r = requests.get(
                url,
                timeout=(30, 120),
                headers=headers,
            )
            r.raise_for_status()
            return r.text

        except Exception as e:
            last_error = e
            print(
                f"WARNING: filelist request failed "
                f"{attempt}/{max_retries}: {repr(e)}"
            )

            if attempt < max_retries:
                sleep_time = min(60, attempt * 5)
                time.sleep(sleep_time)

    raise RuntimeError(
        f"Could not retrieve GEO filelist: {url}\n"
        f"Last error: {repr(last_error)}"
    )


def parse_geo_filelist(text: str):
    """
    Parse GEO series filelist.txt.

    Example line:

    File GSM9223471_Br1113-cell_feature_matrix.h5
         09/04/2025 11:32:47 7153434 H5
    """

    entries = []

    pattern = re.compile(
        r"^File\s+"
        r"(\S+)\s+"
        r"(\d{2}/\d{2}/\d{4})\s+"
        r"(\d{2}:\d{2}:\d{2})\s+"
        r"(\d+)\s+"
        r"(\S+)\s*$"
    )

    for raw_line in text.splitlines():

        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        m = pattern.match(line)

        if not m:
            continue

        filename, date, clock, size, filetype = m.groups()

        entries.append(
            {
                "filename": filename,
                "date": date,
                "time": clock,
                "size": int(size),
                "type": filetype,
            }
        )

    return entries


def select_gsm_files(
    entries,
    gsm: str,
    modality: str,
    include_heavy: bool,
    include_visium_images: bool,
):

    prefix = f"{gsm}_"

    gsm_entries = [
        x for x in entries
        if x["filename"].startswith(prefix)
    ]

    selected = []

    for entry in gsm_entries:

        fname = entry["filename"].lower()

        if modality == "xenium":

            light = any(
                fname.endswith(x.lower())
                for x in XENIUM_LIGHT_SUFFIXES
            )

            heavy = any(
                fname.endswith(x.lower())
                for x in XENIUM_HEAVY_SUFFIXES
            )

            if light or (include_heavy and heavy):
                selected.append(entry)

        elif modality == "visium":

            analysis = any(
                fname.endswith(x.lower())
                for x in VISIUM_ANALYSIS_SUFFIXES
            )

            extra_images = any(
                fname.endswith(x.lower())
                for x in VISIUM_EXTRA_IMAGE_SUFFIXES
            )

            if analysis or (
                include_visium_images and extra_images
            ):
                selected.append(entry)

        else:
            raise ValueError(modality)

    return selected


def geo_download_url(gsm: str, filename: str) -> str:

    query = urlencode(
        {
            "acc": gsm,
            "format": "file",
            "file": filename,
        }
    )

    return f"{GEO_DOWNLOAD_ENDPOINT}?{query}"


def md5sum(path: Path, chunk_size: int = 8 * 1024 * 1024):

    h = hashlib.md5()

    with path.open("rb") as f:

        while True:
            chunk = f.read(chunk_size)

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


def quarantine(path: Path, reason: str):

    if not path.exists():
        return None

    stamp = time.strftime("%Y%m%d_%H%M%S")

    new_path = Path(
        str(path) + f".{reason}.{stamp}"
    )

    print(f"Quarantining:")
    print(f"  {path}")
    print(f"  -> {new_path}")

    path.rename(new_path)

    return new_path


# =============================================================================
# DOWNLOAD
# =============================================================================

def stream_download(
    url: str,
    expected_size: int,
    final_path: Path,
    max_retries: int,
    dry_run: bool,
    force: bool,
):

    final_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    part_path = Path(str(final_path) + ".part")

    record = {
        "filename": final_path.name,
        "out_path": str(final_path),
        "expected_size": expected_size,
        "local_size": 0,
        "status": "",
        "md5": "",
    }

    # -------------------------------------------------------------------------
    # Existing FINAL file
    # -------------------------------------------------------------------------

    if final_path.exists():

        local_size = final_path.stat().st_size

        print(f"Existing final file:")
        print(f"  {final_path}")
        print(f"  local:    {human_bytes(local_size)}")
        print(f"  expected: {human_bytes(expected_size)}")

        if force:

            quarantine(
                final_path,
                "forced_redownload",
            )

        elif local_size == expected_size:

            print("Size matches GEO filelist -> SKIP")

            record["local_size"] = local_size
            record["status"] = "skipped_verified_size"

            return record

        else:

            print(
                "WARNING: existing final file size "
                "does NOT match GEO."
            )

            # If smaller, preserve it as a .part so we can attempt resume.
            if (
                local_size < expected_size
                and not part_path.exists()
            ):

                print(
                    "Moving incomplete file to .part "
                    "for resume."
                )

                final_path.rename(part_path)

            else:

                quarantine(
                    final_path,
                    "bad_size",
                )

    # -------------------------------------------------------------------------
    # Existing PART file
    # -------------------------------------------------------------------------

    if part_path.exists():

        part_size = part_path.stat().st_size

        print(f"Existing partial file:")
        print(f"  {part_path}")
        print(f"  size:     {human_bytes(part_size)}")
        print(f"  expected: {human_bytes(expected_size)}")

        if force:

            quarantine(
                part_path,
                "forced_redownload",
            )

        elif part_size == expected_size:

            print(
                ".part file already has exact "
                "expected size -> promote to final."
            )

            part_path.rename(final_path)

            record["local_size"] = expected_size
            record["status"] = "promoted_complete_part"

            return record

        elif part_size > expected_size:

            quarantine(
                part_path,
                "oversize_part",
            )

    if dry_run:

        current = (
            part_path.stat().st_size
            if part_path.exists()
            else 0
        )

        print("[DRY RUN]")
        print(f"  URL:      {url}")
        print(f"  OUT:      {final_path}")
        print(f"  EXPECTED: {expected_size}")
        print(f"  CURRENT:  {current}")

        record["local_size"] = current
        record["status"] = "dry_run"

        return record

    # -------------------------------------------------------------------------
    # Real download / resume
    # -------------------------------------------------------------------------

    headers_base = {
        "User-Agent": "GeneBridge-GEO-downloader/1.0",
        "Accept-Encoding": "identity",
    }

    for attempt in range(1, max_retries + 1):

        resume_from = (
            part_path.stat().st_size
            if part_path.exists()
            else 0
        )

        headers = dict(headers_base)

        if resume_from > 0:
            headers["Range"] = f"bytes={resume_from}-"

        print()
        print(
            f"Download attempt {attempt}/{max_retries}"
        )
        print(f"URL      : {url}")
        print(f"Output   : {final_path}")
        print(f"Expected : {human_bytes(expected_size)}")
        print(f"Resume   : {human_bytes(resume_from)}")

        try:

            with requests.get(
                url,
                stream=True,
                allow_redirects=True,
                timeout=(30, 600),
                headers=headers,
            ) as r:

                r.raise_for_status()

                # Server ignored our Range request.
                if resume_from > 0 and r.status_code == 200:

                    print(
                        "Server ignored HTTP Range. "
                        "Restarting this file from byte 0."
                    )

                    quarantine(
                        part_path,
                        "range_not_supported",
                    )

                    resume_from = 0
                    mode = "wb"

                elif resume_from > 0 and r.status_code == 206:

                    mode = "ab"

                else:

                    resume_from = 0
                    mode = "wb"

                downloaded = resume_from
                last_print = downloaded

                with part_path.open(mode) as f:

                    for chunk in r.iter_content(
                        chunk_size=8 * 1024 * 1024
                    ):

                        if not chunk:
                            continue

                        f.write(chunk)
                        downloaded += len(chunk)

                        if (
                            downloaded - last_print
                            >= 250 * 1024 * 1024
                        ):

                            pct = (
                                100.0
                                * downloaded
                                / expected_size
                            )

                            print(
                                f"  {human_bytes(downloaded)} / "
                                f"{human_bytes(expected_size)} "
                                f"({pct:.1f}%)",
                                flush=True,
                            )

                            last_print = downloaded

            # -----------------------------------------------------------------
            # VERIFY SIZE
            # -----------------------------------------------------------------

            actual_size = part_path.stat().st_size

            print()
            print("Download connection closed.")
            print(
                f"Downloaded size: {actual_size} "
                f"({human_bytes(actual_size)})"
            )
            print(
                f"Expected size  : {expected_size} "
                f"({human_bytes(expected_size)})"
            )

            if actual_size == expected_size:

                os.replace(
                    part_path,
                    final_path,
                )

                print(
                    "PASS: exact GEO byte-size match."
                )

                record["local_size"] = actual_size
                record["status"] = "downloaded_verified_size"

                return record

            if actual_size > expected_size:

                quarantine(
                    part_path,
                    "oversize_after_download",
                )

                raise RuntimeError(
                    f"Downloaded file is larger than expected: "
                    f"{final_path.name}"
                )

            # Smaller: keep .part and retry/resume.
            print(
                "WARNING: download is incomplete. "
                "Keeping .part for resume."
            )

        except Exception as e:

            print(
                f"WARNING: attempt {attempt} failed:"
            )
            print(repr(e))

        if attempt < max_retries:

            sleep_seconds = min(
                120,
                15 * attempt,
            )

            print(
                f"Retrying after {sleep_seconds} seconds..."
            )

            time.sleep(sleep_seconds)

    raise RuntimeError(
        f"FAILED after {max_retries} attempts:\n"
        f"{final_path}"
    )


# =============================================================================
# DONOR DOWNLOAD
# =============================================================================

def download_modality(
    pair,
    modality,
    entries,
    include_heavy,
    include_visium_images,
    dry_run,
    force,
    max_retries,
):

    donor = pair["patient_id"]

    if modality == "xenium":
        gsm = pair["xenium_gsm"]
        out_dir = XENIUM_DIR / donor
    else:
        gsm = pair["visium_gsm"]
        out_dir = VISIUM_DIR / donor

    selected = select_gsm_files(
        entries=entries,
        gsm=gsm,
        modality=modality,
        include_heavy=include_heavy,
        include_visium_images=include_visium_images,
    )

    section(
        f"{donor} | {modality.upper()} | {gsm}"
    )

    print(
        f"Files selected: {len(selected)}"
    )

    if len(selected) == 0:

        raise RuntimeError(
            f"No files selected for "
            f"{donor} / {gsm} / {modality}"
        )

    total_bytes = sum(
        x["size"] for x in selected
    )

    print(
        f"Expected total: {human_bytes(total_bytes)}"
    )

    for x in selected:
        print(
            f"  {human_bytes(x['size']):>12}  "
            f"{x['filename']}"
        )

    records = []

    for entry in selected:

        filename = entry["filename"]
        expected_size = entry["size"]

        url = geo_download_url(
            gsm=gsm,
            filename=filename,
        )

        final_path = out_dir / filename

        result = stream_download(
            url=url,
            expected_size=expected_size,
            final_path=final_path,
            max_retries=max_retries,
            dry_run=dry_run,
            force=force,
        )

        records.append(
            {
                "patient_index": pair["patient_index"],
                "patient_id": donor,
                "modality": modality,
                "gsm": gsm,
                "filename": filename,
                "url": url,
                "out_path": str(final_path),
                "expected_size": expected_size,
                "local_size": result["local_size"],
                "status": result["status"],
                "md5": result["md5"],
            }
        )

    return records


# =============================================================================
# MANIFEST
# =============================================================================

MANIFEST_FIELDS = [
    "patient_index",
    "patient_id",
    "modality",
    "gsm",
    "filename",
    "url",
    "out_path",
    "expected_size",
    "local_size",
    "status",
    "md5",
]


def write_manifest(donor: str, records):

    path = (
        MANIFEST_DIR
        / f"download_manifest_{donor}.csv"
    )

    tmp = Path(str(path) + ".tmp")

    with tmp.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=MANIFEST_FIELDS,
        )

        writer.writeheader()
        writer.writerows(records)

    os.replace(tmp, path)

    print(f"Manifest: {path}")

    return path


# =============================================================================
# MAIN
# =============================================================================

def main():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--donor",
        choices=sorted(PAIR_BY_DONOR),
        help=(
            "Download only one donor. "
            "Recommended for SLURM array jobs."
        ),
    )

    parser.add_argument(
        "--modality",
        choices=["both", "xenium", "visium"],
        default="both",
    )

    parser.add_argument(
        "--include-heavy",
        action="store_true",
        help=(
            "Include Xenium transcripts and morphology. "
            "These files are very large."
        ),
    )

    parser.add_argument(
        "--include-visium-images",
        action="store_true",
        help=(
            "Also download Visium aligned_fiducials "
            "and detected_tissue_image files."
        ),
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Quarantine existing files and redownload. "
            "Useful for known corrupted donor/files."
        ),
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
    )

    parser.add_argument(
        "--max-retries",
        type=int,
        default=10,
    )

    args = parser.parse_args()

    section("GeneBridge GEO downloader")

    print(f"Project root          : {PROJECT_ROOT}")
    print(f"Donor                 : {args.donor or 'ALL'}")
    print(f"Modality              : {args.modality}")
    print(f"Include heavy Xenium  : {args.include_heavy}")
    print(f"Include Visium images : {args.include_visium_images}")
    print(f"Force redownload      : {args.force}")
    print(f"Dry run               : {args.dry_run}")
    print(f"Max retries           : {args.max_retries}")
    print(f"Started               : {time.ctime()}")

    # -------------------------------------------------------------------------
    # Retrieve the two tiny GEO inventories.
    # -------------------------------------------------------------------------

    section("Reading GEO series file lists")

    xenium_text = fetch_text(
        XENIUM_FILELIST_URL,
        max_retries=args.max_retries,
    )

    visium_text = fetch_text(
        VISIUM_FILELIST_URL,
        max_retries=args.max_retries,
    )

    xenium_entries = parse_geo_filelist(
        xenium_text
    )

    visium_entries = parse_geo_filelist(
        visium_text
    )

    print(
        f"Parsed Xenium files: {len(xenium_entries)}"
    )

    print(
        f"Parsed Visium files: {len(visium_entries)}"
    )

    if not xenium_entries:
        raise RuntimeError(
            "Failed to parse Xenium GEO filelist."
        )

    if not visium_entries:
        raise RuntimeError(
            "Failed to parse Visium GEO filelist."
        )

    # -------------------------------------------------------------------------
    # Choose donor(s)
    # -------------------------------------------------------------------------

    if args.donor:

        pairs = [
            PAIR_BY_DONOR[args.donor]
        ]

    else:

        pairs = PAIRS

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    for pair in pairs:

        donor = pair["patient_id"]

        section(
            f"{pair['patient_index']} / {donor}"
        )

        records = []

        if args.modality in (
            "both",
            "xenium",
        ):

            records.extend(
                download_modality(
                    pair=pair,
                    modality="xenium",
                    entries=xenium_entries,
                    include_heavy=args.include_heavy,
                    include_visium_images=args.include_visium_images,
                    dry_run=args.dry_run,
                    force=args.force,
                    max_retries=args.max_retries,
                )
            )

        if args.modality in (
            "both",
            "visium",
        ):

            records.extend(
                download_modality(
                    pair=pair,
                    modality="visium",
                    entries=visium_entries,
                    include_heavy=args.include_heavy,
                    include_visium_images=args.include_visium_images,
                    dry_run=args.dry_run,
                    force=args.force,
                    max_retries=args.max_retries,
                )
            )

        write_manifest(
            donor=donor,
            records=records,
        )

    section("DONE")

    print(
        f"Finished: {time.ctime()}"
    )


if __name__ == "__main__":
    main()
