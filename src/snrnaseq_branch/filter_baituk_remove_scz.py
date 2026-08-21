#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


CONTROL_SAMPLES = [
    "MB7",
    "MB9",
    "MB11",
    "MB13",
    "MB15",
    "MB16",
    "MB17",
    "MB19",
    "MB21",
    "MB18-2",
    "MB51",
    "MB53",
    "MB55",
    "MB57",
]

SCZ_SAMPLES = [
    "MB6",
    "MB10",
    "MB12",
    "MB14",
    "MB22",
    "MB23",
    "MB8-2",
    "MB54",
    "MB56",
]

# MB8 is not explicitly assigned in the authors' samplegroups list.
AMBIGUOUS_SAMPLE = "MB8"

OFFICIAL_SAMPLES = set(CONTROL_SAMPLES) | set(SCZ_SAMPLES)

SAMPLE_COLUMN_CANDIDATES = [
    "sample_id",
    "sample",
    "Sample",
    "sample_name",
    "donor",
    "donor_id",
    "brain",
    "brain_id",
    "orig.ident",
    "orig_ident",
    "dataset",
    "batch",
    "individual",
    "individual_id",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove Batiuk schizophrenia samples from a full H5AD."
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--mb8-policy",
        choices=["error", "scz", "control"],
        default="error",
        help=(
            "How to handle an exact MB8 sample. Default is error because "
            "the authors' public samplegroups mapping lists MB8-2, not MB8."
        ),
    )

    return parser.parse_args()


ALL_PARSEABLE_IDS = sorted(
    OFFICIAL_SAMPLES | {AMBIGUOUS_SAMPLE},
    key=len,
    reverse=True,
)

SAMPLE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])("
    + "|".join(re.escape(sample) for sample in ALL_PARSEABLE_IDS)
    + r")(?=$|[-_:/|])"
)


def parse_sample_value(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None

    text = str(value).strip()
    match = SAMPLE_PATTERN.search(text)

    if match is None:
        return None

    return match.group(1)


def parse_series(values, index: pd.Index) -> pd.Series:
    return pd.Series(
        [parse_sample_value(value) for value in values],
        index=index,
        dtype="object",
    )


def detect_sample_ids(adata: ad.AnnData) -> tuple[pd.Series, str]:
    """
    Detect sample IDs from obs metadata or observation names.

    The source with the largest number of recognized Batiuk sample IDs
    is selected.
    """
    candidates: list[tuple[int, int, str, pd.Series]] = []

    column_order = []

    for column in SAMPLE_COLUMN_CANDIDATES:
        if column in adata.obs.columns:
            column_order.append(column)

    for column in adata.obs.columns:
        if column not in column_order:
            dtype = adata.obs[column].dtype

            if (
                pd.api.types.is_object_dtype(dtype)
                or isinstance(dtype, pd.CategoricalDtype)
                or pd.api.types.is_string_dtype(dtype)
            ):
                column_order.append(column)

    for priority, column in enumerate(column_order):
        parsed = parse_series(
            adata.obs[column].to_numpy(),
            adata.obs_names,
        )

        matched = int(parsed.notna().sum())

        if matched:
            candidates.append(
                (
                    matched,
                    -priority,
                    f"obs['{column}']",
                    parsed,
                )
            )

    parsed_names = parse_series(
        adata.obs_names.astype(str),
        adata.obs_names,
    )

    matched_names = int(parsed_names.notna().sum())

    if matched_names:
        candidates.append(
            (
                matched_names,
                -len(column_order),
                "obs_names",
                parsed_names,
            )
        )

    if not candidates:
        raise RuntimeError(
            "Could not detect Batiuk sample IDs from adata.obs or obs_names."
        )

    candidates.sort(
        key=lambda item: (item[0], item[1]),
        reverse=True,
    )

    matched, _, source, sample_ids = candidates[0]

    print("Detected sample source:", source)
    print("Matched observations:", matched, "/", adata.n_obs)

    return sample_ids, source


def main() -> None:
    args = parse_args()

    input_path = args.input.resolve()
    output_path = args.output.resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input H5AD does not exist: {input_path}")

    if input_path == output_path:
        raise ValueError(
            "Input and output paths must differ. "
            "The original full H5AD will not be overwritten."
        )

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 100)
    print("BAITIUK SCZ SAMPLE REMOVAL")
    print("=" * 100)
    print("Input:", input_path)
    print("Output:", output_path)
    print("MB8 policy:", args.mb8_policy)

    source = ad.read_h5ad(
        input_path,
        backed="r",
    )

    print("Input shape:", source.shape)

    sample_ids, sample_source = detect_sample_ids(source)

    unmatched_mask = sample_ids.isna()
    unmatched_count = int(unmatched_mask.sum())

    if unmatched_count:
        examples = source.obs_names[
            unmatched_mask.to_numpy()
        ][:20].astype(str).tolist()

        source.file.close()

        raise RuntimeError(
            f"{unmatched_count} observations could not be assigned to an "
            f"official Batiuk sample ID. Example obs_names: {examples}"
        )

    sample_counts = (
        sample_ids.value_counts()
        .rename_axis("sample_id")
        .reset_index(name="cells_before")
    )

    detected_samples = set(sample_counts["sample_id"])

    print()
    print("Detected samples and cell counts:")
    print(sample_counts.to_string(index=False))

    if AMBIGUOUS_SAMPLE in detected_samples:
        if args.mb8_policy == "error":
            source.file.close()

            raise RuntimeError(
                "The H5AD contains sample MB8. The authors' public "
                "samplegroups mapping explicitly lists MB8-2 as SCZ but "
                "does not assign MB8. No output was written. Confirm MB8 "
                "before rerunning with --mb8-policy scz or control."
            )

        if args.mb8_policy == "scz":
            scz_samples = set(SCZ_SAMPLES) | {AMBIGUOUS_SAMPLE}
            control_samples = set(CONTROL_SAMPLES)
        else:
            scz_samples = set(SCZ_SAMPLES)
            control_samples = set(CONTROL_SAMPLES) | {AMBIGUOUS_SAMPLE}
    else:
        scz_samples = set(SCZ_SAMPLES)
        control_samples = set(CONTROL_SAMPLES)

    unexpected_samples = detected_samples - scz_samples - control_samples

    if unexpected_samples:
        source.file.close()

        raise RuntimeError(
            "Detected unclassified sample IDs: "
            f"{sorted(unexpected_samples)}"
        )

    keep_mask = sample_ids.isin(control_samples).to_numpy()
    remove_mask = sample_ids.isin(scz_samples).to_numpy()

    if np.any(~(keep_mask | remove_mask)):
        source.file.close()
        raise RuntimeError("Some observations were neither retained nor removed.")

    retained_cells = int(keep_mask.sum())
    removed_cells = int(remove_mask.sum())

    if retained_cells == 0:
        source.file.close()
        raise RuntimeError("Filtering would retain zero control cells.")

    if removed_cells == 0:
        source.file.close()
        raise RuntimeError("No SCZ cells were detected for removal.")

    print()
    print("Original cells:", source.n_obs)
    print("Control cells retained:", retained_cells)
    print("SCZ cells removed:", removed_cells)
    print("Genes retained:", source.n_vars)

    # Load only retained control observations into memory.
    filtered = source[keep_mask, :].to_memory()
    source.file.close()

    retained_sample_ids = sample_ids.loc[
        filtered.obs_names
    ].astype(str)

    filtered.obs["baituk_sample_id"] = retained_sample_ids.to_numpy()
    filtered.obs["diagnosis"] = "control"
    filtered.obs["scz_removed"] = False

    filter_metadata = {
        "filter": "remove_scz_samples",
        "diagnosis_source": (
            "Batiuk et al. official Notebook1_preprocessing.md samplegroups"
        ),
        "control_samples": sorted(control_samples),
        "scz_samples_removed": sorted(scz_samples),
        "mb8_policy": args.mb8_policy,
        "sample_id_source": sample_source,
        "input_file": str(input_path),
        "original_shape": [int(source.n_obs), int(source.n_vars)],
        "output_shape": [int(filtered.n_obs), int(filtered.n_vars)],
        "retained_control_cells": retained_cells,
        "removed_scz_cells": removed_cells,
        "created_at": datetime.now().isoformat(),
    }

    filtered.uns["baituk_scz_filter"] = filter_metadata

    output_prefix = output_path.with_suffix("")

    counts_path = Path(
        str(output_prefix) + "_sample_counts.csv"
    )

    manifest_path = Path(
        str(output_prefix) + "_filter_manifest.json"
    )

    diagnosis_map = {
        **{sample: "control" for sample in control_samples},
        **{sample: "SCZ" for sample in scz_samples},
    }

    sample_counts["diagnosis"] = sample_counts["sample_id"].map(
        diagnosis_map
    )

    sample_counts["cells_after"] = np.where(
        sample_counts["diagnosis"].eq("control"),
        sample_counts["cells_before"],
        0,
    )

    sample_counts["removed"] = sample_counts["diagnosis"].eq("SCZ")

    sample_counts.to_csv(
        counts_path,
        index=False,
    )

    temporary_path = output_path.with_name(
        output_path.stem + ".tmp.h5ad"
    )

    if temporary_path.exists():
        temporary_path.unlink()

    print()
    print("Writing filtered H5AD:", output_path)

    filtered.write_h5ad(
        temporary_path,
        compression="gzip",
    )

    os.replace(
        temporary_path,
        output_path,
    )

    manifest_path.write_text(
        json.dumps(filter_metadata, indent=2) + "\n"
    )

    # Final on-disk verification.
    check = ad.read_h5ad(
        output_path,
        backed="r",
    )

    if check.n_obs != retained_cells:
        check.file.close()
        raise RuntimeError(
            f"Output cell count mismatch: {check.n_obs} != {retained_cells}"
        )

    if check.n_vars != filtered.n_vars:
        check.file.close()
        raise RuntimeError("Output gene count differs from the input.")

    remaining_diagnoses = set(
        check.obs["diagnosis"].astype(str).unique()
    )

    if remaining_diagnoses != {"control"}:
        check.file.close()
        raise RuntimeError(
            f"Unexpected diagnoses remain: {remaining_diagnoses}"
        )

    check.file.close()

    print()
    print("=" * 100)
    print("SUCCESS")
    print("=" * 100)
    print("Control-only H5AD:", output_path)
    print("Output shape:", filtered.shape)
    print("SCZ cells removed:", removed_cells)
    print("Sample-count audit:", counts_path)
    print("Manifest:", manifest_path)


if __name__ == "__main__":
    main()
