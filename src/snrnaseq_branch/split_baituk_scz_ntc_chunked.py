#!/usr/bin/env python3

import argparse
import gc
import json
import os
import shutil
import time
from pathlib import Path

import anndata as ad
import pandas as pd
from anndata.experimental import concat_on_disk


SCZ_SAMPLES = {
    "MB6", "MB8", "MB8-2", "MB10", "MB12",
    "MB14", "MB22", "MB23", "MB54", "MB56",
}

NTC_SAMPLES = {
    "MB7", "MB9", "MB11", "MB13", "MB15",
    "MB16", "MB17", "MB18-2", "MB19", "MB21",
    "MB51", "MB53", "MB55", "MB57",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Chunked Batiuk SCZ/NTC H5AD splitter."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--scz-output", required=True)
    parser.add_argument("--ntc-output", required=True)
    parser.add_argument("--chunk-dir", required=True)
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def atomic_write_h5ad(adata, output_path):
    output_path = Path(output_path)
    temporary = output_path.with_name(
        output_path.stem + ".tmp.h5ad"
    )

    temporary.unlink(missing_ok=True)

    adata.write_h5ad(
        temporary,
        compression="lzf",
    )

    os.replace(temporary, output_path)


def chunk_is_valid(path, expected_cells, expected_genes):
    path = Path(path)

    if not path.exists() or path.stat().st_size == 0:
        return False

    try:
        obj = ad.read_h5ad(path, backed="r")
        valid = (
            obj.n_obs == expected_cells
            and obj.n_vars == expected_genes
        )
        obj.file.close()
        return valid
    except Exception:
        return False


def load_chunk_with_retry(source, start, stop, retries):
    for attempt in range(1, retries + 1):
        try:
            print(
                f"Loading rows {start:,}:{stop:,} "
                f"(attempt {attempt}/{retries})",
                flush=True,
            )
            return source[start:stop, :].to_memory()

        except Exception as exc:
            print(
                f"Read failed for rows {start:,}:{stop:,}: "
                f"{type(exc).__name__}: {exc}",
                flush=True,
            )

            if attempt == retries:
                raise RuntimeError(
                    f"Could not read source rows {start}:{stop} "
                    f"after {retries} attempts."
                ) from exc

            time.sleep(20)


def concatenate_chunks(chunk_paths, final_output, expected_cells, expected_genes):
    final_output = Path(final_output)
    temporary = final_output.with_name(
        final_output.stem + ".tmp.h5ad"
    )

    temporary.unlink(missing_ok=True)

    print(
        f"Combining {len(chunk_paths)} chunks into:\n"
        f"  {final_output}",
        flush=True,
    )

    concat_on_disk(
        [str(path) for path in chunk_paths],
        str(temporary),
        axis=0,
        join="inner",
        merge="same",
        uns_merge="same",
    )

    check = ad.read_h5ad(temporary, backed="r")
    actual_shape = check.shape
    check.file.close()

    expected_shape = (expected_cells, expected_genes)

    if actual_shape != expected_shape:
        raise RuntimeError(
            f"Final shape mismatch for {final_output}: "
            f"expected {expected_shape}, found {actual_shape}"
        )

    os.replace(temporary, final_output)

    print(
        f"Final output completed: {final_output}\n"
        f"Shape: {expected_shape}",
        flush=True,
    )


def main():
    args = parse_args()

    input_path = Path(args.input).resolve()
    scz_output = Path(args.scz_output).resolve()
    ntc_output = Path(args.ntc_output).resolve()
    chunk_dir = Path(args.chunk_dir).resolve()

    scz_chunk_dir = chunk_dir / "SCZ"
    ntc_chunk_dir = chunk_dir / "NTC"
    manifest_path = chunk_dir / "chunk_manifest.json"

    scz_output.parent.mkdir(parents=True, exist_ok=True)
    scz_chunk_dir.mkdir(parents=True, exist_ok=True)
    ntc_chunk_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 100, flush=True)
    print("CHUNKED BATIUK SCZ/NTC SPLIT", flush=True)
    print("=" * 100, flush=True)
    print(f"Input: {input_path}", flush=True)
    print(f"Chunk size: {args.chunk_size:,}", flush=True)
    print(f"Chunk directory: {chunk_dir}", flush=True)
    print(f"SCZ output: {scz_output}", flush=True)
    print(f"NTC output: {ntc_output}", flush=True)

    source_stat = input_path.stat()

    signature = {
        "input": str(input_path),
        "input_size": source_stat.st_size,
        "input_mtime_ns": source_stat.st_mtime_ns,
        "chunk_size": args.chunk_size,
        "scz_samples": sorted(SCZ_SAMPLES),
        "ntc_samples": sorted(NTC_SAMPLES),
    }

    if manifest_path.exists():
        old_signature = json.loads(manifest_path.read_text())

        if old_signature != signature:
            print(
                "Source or settings changed. Removing stale chunks.",
                flush=True,
            )
            shutil.rmtree(chunk_dir)
            scz_chunk_dir.mkdir(parents=True, exist_ok=True)
            ntc_chunk_dir.mkdir(parents=True, exist_ok=True)

    manifest_path.write_text(
        json.dumps(signature, indent=2)
    )

    print("Opening source in backed mode...", flush=True)
    source = ad.read_h5ad(input_path, backed="r")

    n_cells, n_genes = source.shape

    if "sample" not in source.obs.columns:
        raise KeyError(
            "The source H5AD does not contain obs['sample']."
        )

    sample_ids = (
        source.obs["sample"]
        .astype(str)
        .str.strip()
    )

    detected_samples = set(sample_ids.unique())
    assigned_samples = SCZ_SAMPLES | NTC_SAMPLES
    unknown_samples = sorted(detected_samples - assigned_samples)

    if unknown_samples:
        raise RuntimeError(
            f"Unassigned samples detected: {unknown_samples}"
        )

    sample_counts = (
        sample_ids.value_counts()
        .rename_axis("sample_id")
        .reset_index(name="cell_count")
        .sort_values("sample_id")
    )

    sample_counts["diagnosis"] = sample_counts["sample_id"].map(
        lambda sample: "SCZ" if sample in SCZ_SAMPLES else "NTC"
    )

    counts_output = (
        input_path.parent
        / "baituk_snrna_reference_SCZ_NTC_sample_counts.csv"
    )
    sample_counts.to_csv(counts_output, index=False)

    total_scz = int(sample_ids.isin(SCZ_SAMPLES).sum())
    total_ntc = int(sample_ids.isin(NTC_SAMPLES).sum())

    print(f"Input shape: {source.shape}", flush=True)
    print(f"SCZ cells: {total_scz:,}", flush=True)
    print(f"NTC cells: {total_ntc:,}", flush=True)
    print(f"Genes: {n_genes:,}", flush=True)

    number_of_chunks = (
        n_cells + args.chunk_size - 1
    ) // args.chunk_size

    for chunk_number, start in enumerate(
        range(0, n_cells, args.chunk_size),
        start=1,
    ):
        stop = min(start + args.chunk_size, n_cells)

        chunk_samples = sample_ids.iloc[start:stop]

        expected_scz = int(
            chunk_samples.isin(SCZ_SAMPLES).sum()
        )
        expected_ntc = int(
            chunk_samples.isin(NTC_SAMPLES).sum()
        )

        scz_chunk_path = (
            scz_chunk_dir
            / f"scz_chunk_{chunk_number:05d}.h5ad"
        )
        ntc_chunk_path = (
            ntc_chunk_dir
            / f"ntc_chunk_{chunk_number:05d}.h5ad"
        )

        scz_complete = (
            expected_scz == 0
            or chunk_is_valid(
                scz_chunk_path,
                expected_scz,
                n_genes,
            )
        )

        ntc_complete = (
            expected_ntc == 0
            or chunk_is_valid(
                ntc_chunk_path,
                expected_ntc,
                n_genes,
            )
        )

        print(
            f"\nChunk {chunk_number}/{number_of_chunks}: "
            f"rows {start:,}:{stop:,}; "
            f"SCZ={expected_scz:,}, NTC={expected_ntc:,}",
            flush=True,
        )

        if scz_complete and ntc_complete:
            print("Already completed; skipping.", flush=True)
            continue

        chunk = load_chunk_with_retry(
            source,
            start,
            stop,
            args.retries,
        )

        normalized_samples = (
            chunk.obs["sample"]
            .astype(str)
            .str.strip()
        )

        if expected_scz > 0 and not scz_complete:
            scz_mask = normalized_samples.isin(
                SCZ_SAMPLES
            ).to_numpy()

            scz_chunk = chunk[scz_mask, :].copy()
            scz_chunk.obs["baituk_sample_id"] = (
                normalized_samples.loc[scz_mask].to_numpy()
            )
            scz_chunk.obs["diagnosis"] = "SCZ"

            print(
                f"Writing SCZ chunk: {scz_chunk_path}",
                flush=True,
            )
            atomic_write_h5ad(
                scz_chunk,
                scz_chunk_path,
            )
            del scz_chunk

        if expected_ntc > 0 and not ntc_complete:
            ntc_mask = normalized_samples.isin(
                NTC_SAMPLES
            ).to_numpy()

            ntc_chunk = chunk[ntc_mask, :].copy()
            ntc_chunk.obs["baituk_sample_id"] = (
                normalized_samples.loc[ntc_mask].to_numpy()
            )
            ntc_chunk.obs["diagnosis"] = "NTC"

            print(
                f"Writing NTC chunk: {ntc_chunk_path}",
                flush=True,
            )
            atomic_write_h5ad(
                ntc_chunk,
                ntc_chunk_path,
            )
            del ntc_chunk

        del chunk
        gc.collect()

    source.file.close()

    scz_chunks = sorted(
        scz_chunk_dir.glob("scz_chunk_*.h5ad")
    )
    ntc_chunks = sorted(
        ntc_chunk_dir.glob("ntc_chunk_*.h5ad")
    )

    if not scz_chunks:
        raise RuntimeError("No SCZ chunk files were created.")

    if not ntc_chunks:
        raise RuntimeError("No NTC chunk files were created.")

    concatenate_chunks(
        scz_chunks,
        scz_output,
        total_scz,
        n_genes,
    )

    concatenate_chunks(
        ntc_chunks,
        ntc_output,
        total_ntc,
        n_genes,
    )

    split_manifest = {
        **signature,
        "input_shape": [n_cells, n_genes],
        "scz_shape": [total_scz, n_genes],
        "ntc_shape": [total_ntc, n_genes],
        "scz_output": str(scz_output),
        "ntc_output": str(ntc_output),
        "number_of_scz_chunks": len(scz_chunks),
        "number_of_ntc_chunks": len(ntc_chunks),
    }

    final_manifest = (
        input_path.parent
        / "baituk_snrna_reference_SCZ_NTC_split_manifest.json"
    )
    final_manifest.write_text(
        json.dumps(split_manifest, indent=2)
    )

    print("\n" + "=" * 100, flush=True)
    print("SUCCESS", flush=True)
    print(f"SCZ shape: ({total_scz}, {n_genes})", flush=True)
    print(f"NTC shape: ({total_ntc}, {n_genes})", flush=True)
    print("=" * 100, flush=True)


if __name__ == "__main__":
    main()
