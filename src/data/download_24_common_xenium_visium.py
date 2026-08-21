#"""
#Download 24 matched/common Xenium + Visium samples from GEO.
#
#Default mode downloads practical analysis files only.
#
#Run from project root:
#
#    python src/data/download_24_common_xenium_visium.py
#
#To include huge Xenium transcript/morphology files:
#
#    python src/data/download_24_common_xenium_visium.py --include-heavy
#
#To only print links without downloading:
#
#    python src/data/download_24_common_xenium_visium.py --dry-run
#"""
#
#from pathlib import Path
#from urllib.parse import urljoin, urlparse, parse_qs, unquote
#import argparse
#import html
#import re
#import time
#
#import requests
#
#
#PROJECT_ROOT = Path(__file__).resolve().parents[2]
#
#RAW_DIR = PROJECT_ROOT / "data" / "raw"
#XENIUM_DIR = RAW_DIR / "xenium"
#VISIUM_DIR = RAW_DIR / "visium"
#
#XENIUM_DIR.mkdir(parents=True, exist_ok=True)
#VISIUM_DIR.mkdir(parents=True, exist_ok=True)
#
#
#PAIRS = [
#    {"patient_index": "P1", "patient_id": "Br2039", "xenium_gsm": "GSM9223468", "visium_gsm": "GSM9223410"},
#    {"patient_index": "P2", "patient_id": "Br2719", "xenium_gsm": "GSM9223469", "visium_gsm": "GSM9223405"},
#    {"patient_index": "P3", "patient_id": "Br6432", "xenium_gsm": "GSM9223470", "visium_gsm": "GSM9223415"},
#    {"patient_index": "P4", "patient_id": "Br1113", "xenium_gsm": "GSM9223471", "visium_gsm": "GSM9223416"},
#    {"patient_index": "P5", "patient_id": "Br5373", "xenium_gsm": "GSM9223472", "visium_gsm": "GSM9223412"},
#    {"patient_index": "P6", "patient_id": "Br5590", "xenium_gsm": "GSM9223473", "visium_gsm": "GSM9223420"},
#    {"patient_index": "P7", "patient_id": "Br5400", "xenium_gsm": "GSM9223474", "visium_gsm": "GSM9223424"},
#    {"patient_index": "P8", "patient_id": "Br5622", "xenium_gsm": "GSM9223475", "visium_gsm": "GSM9223432"},
#    {"patient_index": "P9", "patient_id": "Br6437", "xenium_gsm": "GSM9223476", "visium_gsm": "GSM9223429"},
#    {"patient_index": "P10", "patient_id": "Br5314", "xenium_gsm": "GSM9223477", "visium_gsm": "GSM9223438"},
#    {"patient_index": "P11", "patient_id": "Br5588", "xenium_gsm": "GSM9223478", "visium_gsm": "GSM9223434"},
#    {"patient_index": "P12", "patient_id": "Br5746", "xenium_gsm": "GSM9223479", "visium_gsm": "GSM9223436"},
#    {"patient_index": "P13", "patient_id": "Br5639", "xenium_gsm": "GSM9223480", "visium_gsm": "GSM9223439"},
#    {"patient_index": "P14", "patient_id": "Br8433", "xenium_gsm": "GSM9223481", "visium_gsm": "GSM9223454"},
#    {"patient_index": "P15", "patient_id": "Br8772", "xenium_gsm": "GSM9223482", "visium_gsm": "GSM9223443"},
#    {"patient_index": "P16", "patient_id": "Br1139", "xenium_gsm": "GSM9223483", "visium_gsm": "GSM9223444"},
#    {"patient_index": "P17", "patient_id": "Br5973", "xenium_gsm": "GSM9223484", "visium_gsm": "GSM9223459"},
#    {"patient_index": "P18", "patient_id": "Br8667", "xenium_gsm": "GSM9223485", "visium_gsm": "GSM9223455"},
#    {"patient_index": "P19", "patient_id": "Br5436", "xenium_gsm": "GSM9223486", "visium_gsm": "GSM9223431"},
#    {"patient_index": "P20", "patient_id": "Br5931", "xenium_gsm": "GSM9223487", "visium_gsm": "GSM9223448"},
#    {"patient_index": "P21", "patient_id": "Br6496", "xenium_gsm": "GSM9223488", "visium_gsm": "GSM9223460"},
#    {"patient_index": "P22", "patient_id": "Br2421", "xenium_gsm": "GSM9223489", "visium_gsm": "GSM9223467"},
#    {"patient_index": "P23", "patient_id": "Br6032", "xenium_gsm": "GSM9223490", "visium_gsm": "GSM9223435"},
#    {"patient_index": "P24", "patient_id": "Br6389", "xenium_gsm": "GSM9223491", "visium_gsm": "GSM9223447"},
#]
#
#
#VISIUM_KEEP = [
#    "matrix.mtx.gz",
#    "features.tsv.gz",
#    "barcodes.tsv.gz",
#    "tissue_positions.csv.gz",
#    "scalefactors_json.json.gz",
#    "tissue_hires_image.png.gz",
#    "tissue_lowres_image.png.gz",
#]
#
#XENIUM_KEEP_LIGHT = [
#    "cell_feature_matrix.h5",
#    "cells.parquet.gz",
#    "cell_boundaries.csv.gz",
#    "nucleus_boundaries.csv.gz",
#]
#
#XENIUM_KEEP_HEAVY = [
#    "transcripts.zarr.zip",
#    "transcripts.csv.gz",
#    "morphology.ome.tif.gz",
#    "morphology_focus.ome.tif.gz",
#]
#
#
#def print_section(title: str):
#    print("\n" + "=" * 90)
#    print(title)
#    print("=" * 90)
#
#
#def get_geo_supplementary_links(gsm: str):
#    """
#    Scrape a GEO GSM page and return supplementary download links.
#
#    GEO pages usually expose links like:
#        /geo/download/?acc=GSM...&format=file&file=...
#    """
#    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
#    print(f"Reading GEO page: {url}")
#
#    r = requests.get(url, timeout=60)
#    r.raise_for_status()
#
#    text = html.unescape(r.text)
#
#    hrefs = re.findall(r'href="([^"]+)"', text)
#    links = []
#
#    for href in hrefs:
#        if "geo/download" in href and "format=file" in href:
#            full_url = urljoin("https://www.ncbi.nlm.nih.gov", href)
#            filename = filename_from_geo_url(full_url)
#            links.append({"url": full_url, "filename": filename})
#
#    # Remove duplicates while preserving order
#    seen = set()
#    unique = []
#    for item in links:
#        key = item["url"]
#        if key not in seen:
#            unique.append(item)
#            seen.add(key)
#
#    return unique
#
#
#def filename_from_geo_url(url: str):
#    parsed = urlparse(url)
#    qs = parse_qs(parsed.query)
#
#    if "file" in qs:
#        return unquote(qs["file"][0])
#
#    return Path(parsed.path).name
#
#
#def keep_file(filename: str, modality: str, include_heavy: bool):
#    fname = filename.lower()
#
#    if modality == "visium":
#        return any(fname.endswith(x.lower()) for x in VISIUM_KEEP)
#
#    if modality == "xenium":
#        light = any(fname.endswith(x.lower()) for x in XENIUM_KEEP_LIGHT)
#        heavy = any(fname.endswith(x.lower()) for x in XENIUM_KEEP_HEAVY)
#
#        if include_heavy:
#            return light or heavy
#
#        return light
#
#    return False
#
#
#def download_file(url: str, out_path: Path, dry_run: bool = False):
#    out_path.parent.mkdir(parents=True, exist_ok=True)
#
#    if out_path.exists() and out_path.stat().st_size > 0:
#        print(f"Exists, skipping: {out_path}")
#        return
#
#    if dry_run:
#        print(f"[DRY RUN] {url}")
#        print(f"          -> {out_path}")
#        return
#
#    tmp_path = out_path.with_suffix(out_path.suffix + ".part")
#
#    print(f"Downloading:")
#    print(f"  {url}")
#    print(f"  -> {out_path}")
#
#    with requests.get(url, stream=True, timeout=120) as r:
#        r.raise_for_status()
#
#        total = int(r.headers.get("content-length", 0))
#        downloaded = 0
#        last_print = 0
#
#        with open(tmp_path, "wb") as f:
#            for chunk in r.iter_content(chunk_size=1024 * 1024):
#                if chunk:
#                    f.write(chunk)
#                    downloaded += len(chunk)
#
#                    if downloaded - last_print >= 100 * 1024 * 1024:
#                        if total > 0:
#                            pct = 100 * downloaded / total
#                            print(f"  {downloaded / 1024**2:.1f} MB / {total / 1024**2:.1f} MB ({pct:.1f}%)")
#                        else:
#                            print(f"  {downloaded / 1024**2:.1f} MB")
#                        last_print = downloaded
#
#    tmp_path.rename(out_path)
#    print(f"Saved: {out_path}")
#
#
#def download_sample(gsm: str, patient_id: str, modality: str, include_heavy: bool, dry_run: bool):
#    if modality == "xenium":
#        out_dir = XENIUM_DIR / patient_id
#    elif modality == "visium":
#        out_dir = VISIUM_DIR / patient_id
#    else:
#        raise ValueError(modality)
#
#    out_dir.mkdir(parents=True, exist_ok=True)
#
#    links = get_geo_supplementary_links(gsm)
#
#    if len(links) == 0:
#        print(f"WARNING: no supplementary links found for {gsm}")
#        return []
#
#    selected = [
#        item for item in links
#        if keep_file(item["filename"], modality=modality, include_heavy=include_heavy)
#    ]
#
#    print(f"{gsm} / {patient_id} / {modality}: {len(selected)} selected files")
#
#    downloaded = []
#
#    for item in selected:
#        out_path = out_dir / item["filename"]
#        download_file(item["url"], out_path, dry_run=dry_run)
#        downloaded.append({
#            "gsm": gsm,
#            "patient_id": patient_id,
#            "modality": modality,
#            "filename": item["filename"],
#            "url": item["url"],
#            "out_path": str(out_path),
#        })
#        time.sleep(0.5)
#
#    return downloaded
#
#
#def save_pair_metadata():
#    meta_dir = PROJECT_ROOT / "data" / "metadata"
#    meta_dir.mkdir(parents=True, exist_ok=True)
#
#    out = meta_dir / "patient_xenium_visium_24_common.csv"
#
#    import csv
#    with open(out, "w", newline="", encoding="utf-8") as f:
#        writer = csv.DictWriter(
#            f,
#            fieldnames=["patient_index", "patient_id", "xenium_gsm", "visium_gsm"]
#        )
#        writer.writeheader()
#        writer.writerows(PAIRS)
#
#    print(f"Saved pair metadata: {out}")
#
#
#def main():
#    parser = argparse.ArgumentParser()
#
#    parser.add_argument(
#        "--include-heavy",
#        action="store_true",
#        help="Also download huge Xenium transcripts/morphology files. Not recommended first."
#    )
#
#    parser.add_argument(
#        "--dry-run",
#        action="store_true",
#        help="Print selected URLs and output paths without downloading."
#    )
#
#    args = parser.parse_args()
#
#    print_section("Download 24 common Xenium + Visium samples")
#
#    print(f"Project root: {PROJECT_ROOT}")
#    print(f"Xenium output: {XENIUM_DIR}")
#    print(f"Visium output: {VISIUM_DIR}")
#    print(f"Include heavy Xenium files: {args.include_heavy}")
#    print(f"Dry run: {args.dry_run}")
#
#    save_pair_metadata()
#
#    records = []
#
#    for pair in PAIRS:
#        patient_id = pair["patient_id"]
#
#        print_section(f"{pair['patient_index']} / {patient_id}")
#
#        records.extend(
#            download_sample(
#                gsm=pair["xenium_gsm"],
#                patient_id=patient_id,
#                modality="xenium",
#                include_heavy=args.include_heavy,
#                dry_run=args.dry_run,
#            )
#        )
#
#        records.extend(
#            download_sample(
#                gsm=pair["visium_gsm"],
#                patient_id=patient_id,
#                modality="visium",
#                include_heavy=args.include_heavy,
#                dry_run=args.dry_run,
#            )
#        )
#
#    manifest = PROJECT_ROOT / "data" / "metadata" / "download_manifest_24_common.csv"
#
#    import csv
#    with open(manifest, "w", newline="", encoding="utf-8") as f:
#        fieldnames = ["gsm", "patient_id", "modality", "filename", "url", "out_path"]
#        writer = csv.DictWriter(f, fieldnames=fieldnames)
#        writer.writeheader()
#        writer.writerows(records)
#
#    print_section("Done")
#    print(f"Saved download manifest: {manifest}")
#
#
#if __name__ == "__main__":
#    main()

"""
Download 24 matched/common Xenium + Visium samples from GEO.

This script is restart-safe.

Default mode downloads practical analysis files only.

For each file:
    1. If the final file already exists, compare local size with remote size.
    2. If size matches, skip it.
    3. If file is incomplete, redownload or resume from .part file.
    4. If download is interrupted, keep the .part file for resume.
    5. md5sum is computed locally and saved in the manifest for record keeping.

Run from project root:

    python src/data/download_24_common_xenium_visium.py

To include huge Xenium transcript/morphology files:

    python src/data/download_24_common_xenium_visium.py --include-heavy

To only print links without downloading:

    python src/data/download_24_common_xenium_visium.py --dry-run
"""

from pathlib import Path
from urllib.parse import urljoin, urlparse, parse_qs, unquote
import argparse
import csv
import hashlib
import html
import re
import time

import requests


PROJECT_ROOT = Path(__file__).resolve().parents[2]

RAW_DIR = PROJECT_ROOT / "data" / "raw"
XENIUM_DIR = RAW_DIR / "xenium"
VISIUM_DIR = RAW_DIR / "visium"

XENIUM_DIR.mkdir(parents=True, exist_ok=True)
VISIUM_DIR.mkdir(parents=True, exist_ok=True)


PAIRS = [
    {"patient_index": "P1", "patient_id": "Br2039", "xenium_gsm": "GSM9223468", "visium_gsm": "GSM9223410"},
    {"patient_index": "P2", "patient_id": "Br2719", "xenium_gsm": "GSM9223469", "visium_gsm": "GSM9223405"},
    {"patient_index": "P3", "patient_id": "Br6432", "xenium_gsm": "GSM9223470", "visium_gsm": "GSM9223415"},
    {"patient_index": "P4", "patient_id": "Br1113", "xenium_gsm": "GSM9223471", "visium_gsm": "GSM9223416"},
    {"patient_index": "P5", "patient_id": "Br5373", "xenium_gsm": "GSM9223472", "visium_gsm": "GSM9223412"},
    {"patient_index": "P6", "patient_id": "Br5590", "xenium_gsm": "GSM9223473", "visium_gsm": "GSM9223420"},
    {"patient_index": "P7", "patient_id": "Br5400", "xenium_gsm": "GSM9223474", "visium_gsm": "GSM9223424"},
    {"patient_index": "P8", "patient_id": "Br5622", "xenium_gsm": "GSM9223475", "visium_gsm": "GSM9223432"},
    {"patient_index": "P9", "patient_id": "Br6437", "xenium_gsm": "GSM9223476", "visium_gsm": "GSM9223429"},
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


VISIUM_KEEP = [
    "matrix.mtx.gz",
    "features.tsv.gz",
    "barcodes.tsv.gz",
    "tissue_positions.csv.gz",
    "scalefactors_json.json.gz",
    "tissue_hires_image.png.gz",
    "tissue_lowres_image.png.gz",
]

XENIUM_KEEP_LIGHT = [
    "cell_feature_matrix.h5",
    "cells.parquet.gz",
    "cell_boundaries.csv.gz",
    "nucleus_boundaries.csv.gz",
]

XENIUM_KEEP_HEAVY = [
    "transcripts.zarr.zip",
    "transcripts.csv.gz",
    "morphology.ome.tif.gz",
    "morphology_focus.ome.tif.gz",
]


def print_section(title: str):
    print("\n" + "=" * 90)
    print(title)
    print("=" * 90)


def filename_from_geo_url(url: str):
    parsed = urlparse(url)
    qs = parse_qs(parsed.query)

    if "file" in qs:
        return unquote(qs["file"][0])

    return Path(parsed.path).name


def get_geo_supplementary_links(gsm: str):
    """
    Scrape a GEO GSM page and return supplementary download links.

    GEO pages usually expose links like:
        /geo/download/?acc=GSM...&format=file&file=...
    """
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
    print(f"Reading GEO page: {url}")

    r = requests.get(url, timeout=60)
    r.raise_for_status()

    text = html.unescape(r.text)
    hrefs = re.findall(r'href="([^"]+)"', text)

    links = []

    for href in hrefs:
        if "geo/download" in href and "format=file" in href:
            full_url = urljoin("https://www.ncbi.nlm.nih.gov", href)
            filename = filename_from_geo_url(full_url)
            links.append({"url": full_url, "filename": filename})

    seen = set()
    unique = []

    for item in links:
        key = item["url"]
        if key not in seen:
            unique.append(item)
            seen.add(key)

    return unique


def keep_file(filename: str, modality: str, include_heavy: bool):
    fname = filename.lower()

    if modality == "visium":
        return any(fname.endswith(x.lower()) for x in VISIUM_KEEP)

    if modality == "xenium":
        light = any(fname.endswith(x.lower()) for x in XENIUM_KEEP_LIGHT)
        heavy = any(fname.endswith(x.lower()) for x in XENIUM_KEEP_HEAVY)

        if include_heavy:
            return light or heavy

        return light

    return False


def compute_md5(path: Path, chunk_size: int = 1024 * 1024):
    """
    Compute local md5sum.
    Used for record keeping in the manifest.
    """
    h = hashlib.md5()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)

    return h.hexdigest()


def get_remote_size(url: str):
    """
    Get remote file size.

    First try HEAD.
    If HEAD fails or does not provide content-length, try GET headers.
    """
    headers = {"Accept-Encoding": "identity"}

    try:
        r = requests.head(
            url,
            allow_redirects=True,
            timeout=60,
            headers=headers,
        )

        if r.status_code < 400 and r.headers.get("content-length"):
            return int(r.headers["content-length"])

    except Exception as e:
        print(f"WARNING: HEAD remote size failed: {repr(e)}")

    try:
        with requests.get(
            url,
            stream=True,
            timeout=60,
            headers=headers,
        ) as r:
            r.raise_for_status()

            if r.headers.get("content-length"):
                return int(r.headers["content-length"])

    except Exception as e:
        print(f"WARNING: GET-header remote size failed: {repr(e)}")

    return None


def check_final_file(out_path: Path, url: str):
    """
    Check whether final local file is complete by size.
    """
    result = {
        "exists": out_path.exists(),
        "local_size": 0,
        "remote_size": "",
        "local_md5": "",
        "complete": False,
        "check_method": "",
    }

    if not out_path.exists():
        result["check_method"] = "missing"
        return result

    local_size = out_path.stat().st_size
    result["local_size"] = local_size

    if local_size == 0:
        result["check_method"] = "zero_size"
        return result

    remote_size = get_remote_size(url)

    if remote_size is not None:
        result["remote_size"] = remote_size
        result["check_method"] = "size"

        if local_size == remote_size:
            result["complete"] = True
            result["local_md5"] = compute_md5(out_path)

        return result

    # If GEO does not report size, keep non-empty file.
    result["check_method"] = "nonempty_no_remote_size"
    result["complete"] = True
    result["local_md5"] = compute_md5(out_path)

    return result


def stream_download_with_resume(
    url: str,
    tmp_path: Path,
    max_retries: int = 5,
    chunk_size: int = 1024 * 1024,
):
    """
    Download to tmp_path with resume support.

    If tmp_path exists, continue from the current byte count using HTTP Range.
    If server ignores Range, restart from zero safely.
    """
    headers_base = {"Accept-Encoding": "identity"}

    for attempt in range(1, max_retries + 1):
        resume_from = tmp_path.stat().st_size if tmp_path.exists() else 0

        headers = dict(headers_base)

        if resume_from > 0:
            headers["Range"] = f"bytes={resume_from}-"
            print(f"Resume attempt {attempt}/{max_retries} from {resume_from / 1024**2:.1f} MB")
        else:
            print(f"Download attempt {attempt}/{max_retries} from beginning")

        try:
            with requests.get(
                url,
                stream=True,
                timeout=(30, 300),
                headers=headers,
                allow_redirects=True,
            ) as r:
                r.raise_for_status()

                if resume_from > 0 and r.status_code == 200:
                    print("Server ignored resume request. Restarting download from zero.")
                    tmp_path.unlink()
                    resume_from = 0

                total_header = int(r.headers.get("content-length", 0))

                if r.status_code == 206:
                    total = resume_from + total_header
                else:
                    total = total_header

                downloaded = resume_from
                last_print = downloaded

                mode = "ab" if resume_from > 0 else "wb"

                with open(tmp_path, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)

                            if downloaded - last_print >= 100 * 1024 * 1024:
                                if total > 0:
                                    pct = 100 * downloaded / total
                                    print(
                                        f"  {downloaded / 1024**2:.1f} MB / "
                                        f"{total / 1024**2:.1f} MB ({pct:.1f}%)"
                                    )
                                else:
                                    print(f"  {downloaded / 1024**2:.1f} MB")

                                last_print = downloaded

            return

        except requests.exceptions.RequestException as e:
            print(f"WARNING: download interrupted on attempt {attempt}/{max_retries}")
            print(f"  error: {repr(e)}")

            if attempt == max_retries:
                raise

            sleep_time = 30 * attempt
            print(f"Sleeping {sleep_time} seconds before retry...")
            time.sleep(sleep_time)


def download_file(
    url: str,
    out_path: Path,
    dry_run: bool = False,
    max_retries: int = 5,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    record = {
        "out_path": str(out_path),
        "local_size": 0,
        "remote_size": "",
        "local_md5": "",
        "status": "",
        "check_method": "",
    }

    tmp_path = out_path.with_suffix(out_path.suffix + ".part")

    final_check = check_final_file(out_path, url)

    if final_check["exists"]:
        print(f"Existing final file found: {out_path}")
        print(f"  local size:  {final_check['local_size']}")
        print(f"  remote size: {final_check['remote_size']}")
        print(f"  check:       {final_check['check_method']}")

        if final_check["complete"]:
            print(f"Complete file already exists, skipping: {out_path}")

            record.update({
                "local_size": final_check["local_size"],
                "remote_size": final_check["remote_size"],
                "local_md5": final_check["local_md5"],
                "status": "skipped_complete",
                "check_method": final_check["check_method"],
            })

            return record

        remote_size = final_check["remote_size"]
        local_size = final_check["local_size"]

        print(f"Incomplete final file detected: {out_path}")

        if remote_size != "" and local_size < remote_size and not tmp_path.exists():
            print("Moving incomplete final file to .part so download can resume.")
            out_path.rename(tmp_path)
        else:
            bad_path = out_path.with_suffix(out_path.suffix + ".bad")
            print(f"Moving incomplete/corrupt final file to: {bad_path}")
            out_path.rename(bad_path)

    if dry_run:
        print(f"[DRY RUN] {url}")
        print(f"          -> {out_path}")

        if tmp_path.exists():
            print(f"          existing partial: {tmp_path}")
            print(f"          partial size: {tmp_path.stat().st_size}")

        record.update({
            "status": "dry_run",
            "check_method": final_check["check_method"],
            "local_size": final_check["local_size"],
            "remote_size": final_check["remote_size"],
        })

        return record

    remote_size = get_remote_size(url)

    if tmp_path.exists():
        partial_size = tmp_path.stat().st_size

        print(f"Found partial file: {tmp_path}")
        print(f"  partial size: {partial_size / 1024**2:.1f} MB")

        if remote_size is not None:
            print(f"  remote size:  {remote_size / 1024**2:.1f} MB")

            if partial_size == remote_size:
                print("Partial file already has complete size. Renaming to final file.")
                tmp_path.rename(out_path)

                final_check = check_final_file(out_path, url)

                record.update({
                    "local_size": final_check["local_size"],
                    "remote_size": final_check["remote_size"],
                    "local_md5": final_check["local_md5"],
                    "status": "renamed_complete_part",
                    "check_method": final_check["check_method"],
                })

                return record

            if partial_size > remote_size:
                print("Partial file is larger than remote file. Removing bad partial.")
                tmp_path.unlink()

    print("Downloading:")
    print(f"  {url}")
    print(f"  -> {out_path}")

    stream_download_with_resume(
        url=url,
        tmp_path=tmp_path,
        max_retries=max_retries,
    )

    tmp_path.rename(out_path)

    final_check = check_final_file(out_path, url)

    if not final_check["complete"]:
        print("ERROR: downloaded file failed completeness check.")
        print(f"  file:        {out_path}")
        print(f"  local size:  {final_check['local_size']}")
        print(f"  remote size: {final_check['remote_size']}")

        bad_path = out_path.with_suffix(out_path.suffix + ".bad")
        out_path.rename(bad_path)

        raise RuntimeError(f"Download failed completeness check. Moved to: {bad_path}")

    print(f"Saved and verified complete: {out_path}")
    print(f"  local size: {final_check['local_size']}")
    print(f"  remote size: {final_check['remote_size']}")
    print(f"  md5sum: {final_check['local_md5']}")

    record.update({
        "local_size": final_check["local_size"],
        "remote_size": final_check["remote_size"],
        "local_md5": final_check["local_md5"],
        "status": "downloaded_complete",
        "check_method": final_check["check_method"],
    })

    return record


def download_sample(
    gsm: str,
    patient_id: str,
    modality: str,
    include_heavy: bool,
    dry_run: bool,
    max_retries: int,
):
    if modality == "xenium":
        out_dir = XENIUM_DIR / patient_id
    elif modality == "visium":
        out_dir = VISIUM_DIR / patient_id
    else:
        raise ValueError(modality)

    out_dir.mkdir(parents=True, exist_ok=True)

    links = get_geo_supplementary_links(gsm)

    if len(links) == 0:
        print(f"WARNING: no supplementary links found for {gsm}")
        return []

    selected = [
        item for item in links
        if keep_file(item["filename"], modality=modality, include_heavy=include_heavy)
    ]

    print(f"{gsm} / {patient_id} / {modality}: {len(selected)} selected files")

    downloaded = []

    for item in selected:
        out_path = out_dir / item["filename"]

        result = download_file(
            item["url"],
            out_path,
            dry_run=dry_run,
            max_retries=max_retries,
        )

        downloaded.append({
            "gsm": gsm,
            "patient_id": patient_id,
            "modality": modality,
            "filename": item["filename"],
            "url": item["url"],
            "out_path": str(out_path),
            "local_size": result["local_size"],
            "remote_size": result["remote_size"],
            "local_md5": result["local_md5"],
            "status": result["status"],
            "check_method": result["check_method"],
        })

        time.sleep(0.5)

    return downloaded


def save_pair_metadata():
    meta_dir = PROJECT_ROOT / "data" / "metadata"
    meta_dir.mkdir(parents=True, exist_ok=True)

    out = meta_dir / "patient_xenium_visium_24_common.csv"

    with open(out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["patient_index", "patient_id", "xenium_gsm", "visium_gsm"],
        )
        writer.writeheader()
        writer.writerows(PAIRS)

    print(f"Saved pair metadata: {out}")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--include-heavy",
        action="store_true",
        help="Also download huge Xenium transcripts/morphology files. Not recommended first.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print selected URLs and output paths without downloading.",
    )

    parser.add_argument(
        "--max-retries",
        type=int,
        default=10,
        help="Maximum retry attempts per file. Default: 10.",
    )

    args = parser.parse_args()

    print_section("Download 24 common Xenium + Visium samples")

    print(f"Project root: {PROJECT_ROOT}")
    print(f"Xenium output: {XENIUM_DIR}")
    print(f"Visium output: {VISIUM_DIR}")
    print(f"Include heavy Xenium files: {args.include_heavy}")
    print(f"Dry run: {args.dry_run}")
    print(f"Max retries per file: {args.max_retries}")

    save_pair_metadata()

    records = []

    for pair in PAIRS:
        patient_id = pair["patient_id"]

        print_section(f"{pair['patient_index']} / {patient_id}")

        records.extend(
            download_sample(
                gsm=pair["xenium_gsm"],
                patient_id=patient_id,
                modality="xenium",
                include_heavy=args.include_heavy,
                dry_run=args.dry_run,
                max_retries=args.max_retries,
            )
        )

        records.extend(
            download_sample(
                gsm=pair["visium_gsm"],
                patient_id=patient_id,
                modality="visium",
                include_heavy=args.include_heavy,
                dry_run=args.dry_run,
                max_retries=args.max_retries,
            )
        )

    manifest = PROJECT_ROOT / "data" / "metadata" / "download_manifest_24_common.csv"

    with open(manifest, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "gsm",
            "patient_id",
            "modality",
            "filename",
            "url",
            "out_path",
            "local_size",
            "remote_size",
            "local_md5",
            "status",
            "check_method",
        ]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print_section("Done")
    print(f"Saved download manifest: {manifest}")


if __name__ == "__main__":
    main()