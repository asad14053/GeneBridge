#!/usr/bin/env python3

import json
import os
from pathlib import Path

import anndata as ad
import h5py
from anndata.experimental import concat_on_disk


PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BAITUK_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "snrnaseq"
    / "Baituk"
)

SAMPLE_DIR = BAITUK_DIR / "baituk_rds_sample_h5ad"

SCZ_OUTPUT = (
    BAITUK_DIR
    / "baituk_snrna_reference_SCZ_full_allgenes.h5ad"
)

NTC_OUTPUT = (
    BAITUK_DIR
    / "baituk_snrna_reference_NTC_full_allgenes.h5ad"
)

MANIFEST_OUTPUT = (
    BAITUK_DIR
    / "baituk_snrna_reference_SCZ_NTC_split_manifest.json"
)

SCZ_SAMPLES = [
    "MB6", "MB8", "MB8-2", "MB10", "MB12",
    "MB14", "MB22", "MB23", "MB54", "MB56",
]

NTC_SAMPLES = [
    "MB7", "MB9", "MB11", "MB13", "MB15",
    "MB16", "MB17", "MB18-2", "MB19", "MB21",
    "MB51", "MB53", "MB55", "MB57",
]

EXPECTED_SCZ_SHAPE = (81817, 60617)
EXPECTED_NTC_SHAPE = (127236, 60617)


def inspect_sparse_storage(path: Path) -> dict:
    result = {}

    with h5py.File(path, "r") as handle:
        x = handle["X"]

        result["encoding_type"] = x.attrs.get(
            "encoding-type",
            "unknown",
        )

        result["shape"] = [
            int(value)
            for value in x.attrs.get("shape", [])
        ]

        for dataset_name in ("data", "indices", "indptr"):
            dataset = x[dataset_name]

            result[dataset_name] = {
                "shape": list(dataset.shape),
                "chunks": (
                    list(dataset.chunks)
                    if dataset.chunks is not None
                    else None
                ),
                "compression": dataset.compression,
            }

    return result


def concatenate_group(
    diagnosis: str,
    sample_ids: list[str],
    output_path: Path,
    expected_shape: tuple[int, int],
) -> dict:
    input_files = [
        SAMPLE_DIR / f"{sample_id}.h5ad"
        for sample_id in sample_ids
    ]

    missing = [
        str(path)
        for path in input_files
        if not path.exists()
    ]

    if missing:
        raise FileNotFoundError(
            "Missing sample H5AD files:\n"
            + "\n".join(missing)
        )

    temporary_path = output_path.with_name(
        output_path.stem + ".tmp.h5ad"
    )

    temporary_path.unlink(missing_ok=True)

    print("=" * 100, flush=True)
    print(
        f"Concatenating {diagnosis}: "
        f"{len(input_files)} sample H5AD files",
        flush=True,
    )

    for path in input_files:
        print(f"  {path}", flush=True)

    concat_on_disk(
        [str(path) for path in input_files],
        str(temporary_path),
        axis=0,
        join="inner",
        merge="same",
        uns_merge="same",
        max_loaded_elems=10_000_000,
    )

    result = ad.read_h5ad(
        temporary_path,
        backed="r",
    )

    actual_shape = result.shape
    diagnosis_values = set(
        result.obs["diagnosis"].astype(str).unique()
    )
    result.file.close()

    if actual_shape != expected_shape:
        raise RuntimeError(
            f"{diagnosis} shape mismatch: "
            f"expected {expected_shape}, found {actual_shape}"
        )

    if diagnosis_values != {diagnosis}:
        raise RuntimeError(
            f"{diagnosis} output contains unexpected diagnoses: "
            f"{sorted(diagnosis_values)}"
        )

    os.replace(temporary_path, output_path)

    storage = inspect_sparse_storage(output_path)

    print(
        f"{diagnosis} completed: {output_path}",
        flush=True,
    )
    print(f"Shape: {actual_shape}", flush=True)
    print(
        "X storage:",
        json.dumps(storage, indent=2),
        flush=True,
    )

    return {
        "diagnosis": diagnosis,
        "samples": sample_ids,
        "output": str(output_path),
        "shape": list(actual_shape),
        "storage": storage,
    }


def main() -> None:
    print("=" * 100)
    print("ON-DISK BATIUK SCZ/NTC CONCATENATION")
    print("=" * 100)

    scz_result = concatenate_group(
        diagnosis="SCZ",
        sample_ids=SCZ_SAMPLES,
        output_path=SCZ_OUTPUT,
        expected_shape=EXPECTED_SCZ_SHAPE,
    )

    ntc_result = concatenate_group(
        diagnosis="NTC",
        sample_ids=NTC_SAMPLES,
        output_path=NTC_OUTPUT,
        expected_shape=EXPECTED_NTC_SHAPE,
    )

    manifest = {
        "source": str(SAMPLE_DIR),
        "scz": scz_result,
        "ntc": ntc_result,
    }

    MANIFEST_OUTPUT.write_text(
        json.dumps(manifest, indent=2)
    )

    print("=" * 100)
    print("SUCCESS")
    print(f"SCZ: {SCZ_OUTPUT}")
    print(f"NTC: {NTC_OUTPUT}")
    print(f"Manifest: {MANIFEST_OUTPUT}")
    print("=" * 100)


if __name__ == "__main__":
    main()
