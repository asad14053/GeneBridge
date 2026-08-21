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


SCZ_SAMPLES = {
    "MB6",
    "MB10",
    "MB12",
    "MB14",
    "MB22",
    "MB23",
    "MB8-2",
    "MB54",
    "MB56",
}

NTC_SAMPLES = {
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
}

AMBIGUOUS_SAMPLE = "MB8"

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
    "batch",
    "individual",
    "individual_id",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Split full Batiuk snRNA-seq H5AD into SCZ and NTC files."
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--scz-output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--ntc-output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--mb8-policy",
        choices=["error", "scz", "ntc"],
        default="error",
        help=(
            "How to classify exact sample MB8. "
            "Default is error because the provided mapping lists MB8-2, "
            "not MB8."
        ),
    )

    return parser.parse_args()


def build_sample_pattern(sample_ids: set[str]) -> re.Pattern:
    ordered = sorted(
        sample_ids,
        key=len,
        reverse=True,
    )

    pattern = "|".join(
        re.escape(sample_id)
        for sample_id in ordered
    )

    return re.compile(
        rf"(?<![A-Za-z0-9])({pattern})(?![A-Za-z0-9])"
    )


def parse_sample_value(
    value: object,
    pattern: re.Pattern,
) -> str | None:
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    match = pattern.search(str(value).strip())

    if match is None:
        return None

    return match.group(1)


def parse_values(
    values,
    index: pd.Index,
    pattern: re.Pattern,
) -> pd.Series:
    return pd.Series(
        [
            parse_sample_value(value, pattern)
            for value in values
        ],
        index=index,
        dtype="object",
    )


def detect_sample_ids(
    adata: ad.AnnData,
    pattern: re.Pattern,
) -> tuple[pd.Series, str]:
    candidates = []

    ordered_columns = []

    for column in SAMPLE_COLUMN_CANDIDATES:
        if column in adata.obs.columns:
            ordered_columns.append(column)

    for column in adata.obs.columns:
        if column in ordered_columns:
            continue

        dtype = adata.obs[column].dtype

        if (
            pd.api.types.is_object_dtype(dtype)
            or pd.api.types.is_string_dtype(dtype)
            or isinstance(dtype, pd.CategoricalDtype)
        ):
            ordered_columns.append(column)

    for priority, column in enumerate(ordered_columns):
        parsed = parse_values(
            adata.obs[column].to_numpy(),
            adata.obs_names,
            pattern,
        )

        matched = int(parsed.notna().sum())

        if matched > 0:
            candidates.append(
                (
                    matched,
                    -priority,
                    f"obs['{column}']",
                    parsed,
                )
            )

    parsed_obs_names = parse_values(
        adata.obs_names.astype(str),
        adata.obs_names,
        pattern,
    )

    matched_obs_names = int(
        parsed_obs_names.notna().sum()
    )

    if matched_obs_names > 0:
        candidates.append(
            (
                matched_obs_names,
                -len(ordered_columns),
                "obs_names",
                parsed_obs_names,
            )
        )

    if not candidates:
        raise RuntimeError(
            "Could not detect MB sample IDs from adata.obs or obs_names."
        )

    candidates.sort(
        key=lambda item: (item[0], item[1]),
        reverse=True,
    )

    matched, _, source, sample_ids = candidates[0]

    print("Selected sample-ID source:", source)
    print("Matched observations:", matched, "/", adata.n_obs)

    return sample_ids, source


def atomic_write_h5ad(
    adata: ad.AnnData,
    destination: Path,
) -> None:
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = destination.with_name(
        destination.stem + ".tmp.h5ad"
    )

    if temporary.exists():
        temporary.unlink()

    adata.write_h5ad(
        temporary,
        compression="lzf",
    )

    os.replace(
        temporary,
        destination,
    )


def validate_output(
    path: Path,
    expected_cells: int,
    expected_genes: int,
    expected_diagnosis: str,
) -> None:
    check = ad.read_h5ad(
        path,
        backed="r",
    )

    try:
        if check.shape != (
            expected_cells,
            expected_genes,
        ):
            raise RuntimeError(
                f"{path.name}: expected "
                f"({expected_cells}, {expected_genes}), "
                f"found {check.shape}."
            )

        diagnoses = set(
            check.obs["diagnosis"]
            .astype(str)
            .unique()
        )

        if diagnoses != {expected_diagnosis}:
            raise RuntimeError(
                f"{path.name}: unexpected diagnoses: "
                f"{sorted(diagnoses)}"
            )

    finally:
        check.file.close()


def main() -> None:
    args = parse_args()

    input_path = args.input.resolve()
    scz_output = args.scz_output.resolve()
    ntc_output = args.ntc_output.resolve()

    if not input_path.exists():
        raise FileNotFoundError(
            f"Missing input H5AD: {input_path}"
        )

    if len({
        input_path,
        scz_output,
        ntc_output,
    }) != 3:
        raise ValueError(
            "Input, SCZ output, and NTC output paths must all differ."
        )

    scz_samples = set(SCZ_SAMPLES)
    ntc_samples = set(NTC_SAMPLES)

    if args.mb8_policy == "scz":
        scz_samples.add(AMBIGUOUS_SAMPLE)
    elif args.mb8_policy == "ntc":
        ntc_samples.add(AMBIGUOUS_SAMPLE)

    recognized_samples = (
        scz_samples
        | ntc_samples
        | {AMBIGUOUS_SAMPLE}
    )

    pattern = build_sample_pattern(
        recognized_samples
    )

    print("=" * 100)
    print("SPLIT BATIUK snRNA-seq INTO SCZ AND NTC")
    print("=" * 100)
    print("Input:", input_path)
    print("SCZ output:", scz_output)
    print("NTC output:", ntc_output)
    print("MB8 policy:", args.mb8_policy)

    print("Loading full input H5AD into memory...")
    source = ad.read_h5ad(input_path)
    print("Full input H5AD loaded.")

    input_cells = int(source.n_obs)
    input_genes = int(source.n_vars)

    print("Input shape:", source.shape)

    sample_ids, sample_id_source = detect_sample_ids(
        source,
        pattern,
    )

    unmatched_mask = sample_ids.isna()

    if unmatched_mask.any():
        unmatched_count = int(
            unmatched_mask.sum()
        )

        examples = source.obs_names[
            unmatched_mask.to_numpy()
        ][:20].astype(str).tolist()

        source.file.close()

        raise RuntimeError(
            f"{unmatched_count} observations could not be assigned "
            f"to a Batiuk sample ID. Examples: {examples}"
        )

    sample_counts = (
        sample_ids
        .value_counts()
        .rename_axis("sample_id")
        .reset_index(name="cell_count")
        .sort_values("sample_id")
        .reset_index(drop=True)
    )

    print()
    print("Detected samples:")
    print(sample_counts.to_string(index=False))

    detected_samples = set(
        sample_counts["sample_id"]
    )

    if (
        AMBIGUOUS_SAMPLE in detected_samples
        and args.mb8_policy == "error"
    ):
        source.file.close()

        raise RuntimeError(
            "Exact sample MB8 was detected. "
            "The supplied mapping explicitly assigns MB8-2 to SCZ "
            "but does not assign MB8. No output was written. "
            "After confirming MB8, rerun with "
            "--mb8-policy scz or --mb8-policy ntc."
        )

    unclassified_samples = (
        detected_samples
        - scz_samples
        - ntc_samples
    )

    if unclassified_samples:
        source.file.close()

        raise RuntimeError(
            "Unclassified sample IDs detected: "
            f"{sorted(unclassified_samples)}"
        )

    scz_mask = sample_ids.isin(
        scz_samples
    ).to_numpy()

    ntc_mask = sample_ids.isin(
        ntc_samples
    ).to_numpy()

    overlap_count = int(
        np.sum(scz_mask & ntc_mask)
    )

    unassigned_count = int(
        np.sum(~(scz_mask | ntc_mask))
    )

    if overlap_count:
        source.file.close()
        raise RuntimeError(
            f"{overlap_count} cells were assigned to both groups."
        )

    if unassigned_count:
        source.file.close()
        raise RuntimeError(
            f"{unassigned_count} cells were not assigned to either group."
        )

    scz_cells = int(
        scz_mask.sum()
    )

    ntc_cells = int(
        ntc_mask.sum()
    )

    if scz_cells == 0:
        source.file.close()
        raise RuntimeError(
            "No SCZ cells were detected."
        )

    if ntc_cells == 0:
        source.file.close()
        raise RuntimeError(
            "No NTC cells were detected."
        )

    if scz_cells + ntc_cells != input_cells:
        source.file.close()
        raise RuntimeError(
            "SCZ and NTC cells do not sum to the input cell count."
        )

    print()
    print("Input cells:", input_cells)
    print("SCZ cells:", scz_cells)
    print("NTC cells:", ntc_cells)
    print("Genes retained in each file:", input_genes)

    scz_adata = source[
        scz_mask,
        :,
    ].to_memory()

    ntc_adata = source[
        ntc_mask,
        :,
    ].to_memory()

    source.file.close()

    scz_ids = sample_ids.loc[
        scz_adata.obs_names
    ].astype(str)

    ntc_ids = sample_ids.loc[
        ntc_adata.obs_names
    ].astype(str)

    scz_adata.obs["baituk_sample_id"] = (
        scz_ids.to_numpy()
    )
    scz_adata.obs["diagnosis"] = "SCZ"

    ntc_adata.obs["baituk_sample_id"] = (
        ntc_ids.to_numpy()
    )
    ntc_adata.obs["diagnosis"] = "NTC"

    created_at = datetime.now().isoformat()

    shared_metadata = {
        "operation": "split_scz_ntc",
        "input_file": str(input_path),
        "input_shape": [
            input_cells,
            input_genes,
        ],
        "sample_id_source": sample_id_source,
        "scz_samples": sorted(scz_samples),
        "ntc_samples": sorted(ntc_samples),
        "mb8_policy": args.mb8_policy,
        "created_at": created_at,
    }

    scz_adata.uns["baituk_diagnosis_split"] = {
        **shared_metadata,
        "group": "SCZ",
        "output_shape": [
            scz_cells,
            input_genes,
        ],
    }

    ntc_adata.uns["baituk_diagnosis_split"] = {
        **shared_metadata,
        "group": "NTC",
        "output_shape": [
            ntc_cells,
            input_genes,
        ],
    }

    diagnosis_map = {
        **{
            sample: "SCZ"
            for sample in scz_samples
        },
        **{
            sample: "NTC"
            for sample in ntc_samples
        },
    }

    sample_counts["diagnosis"] = (
        sample_counts["sample_id"]
        .map(diagnosis_map)
    )

    audit_path = input_path.parent / (
        "baituk_snrna_reference_SCZ_NTC_sample_counts.csv"
    )

    manifest_path = input_path.parent / (
        "baituk_snrna_reference_SCZ_NTC_split_manifest.json"
    )

    sample_counts.to_csv(
        audit_path,
        index=False,
    )

    print()
    print("Writing SCZ H5AD:", scz_output)
    atomic_write_h5ad(
        scz_adata,
        scz_output,
    )

    print("Writing NTC H5AD:", ntc_output)
    atomic_write_h5ad(
        ntc_adata,
        ntc_output,
    )

    manifest = {
        **shared_metadata,
        "scz_output": str(scz_output),
        "scz_shape": [
            scz_cells,
            input_genes,
        ],
        "ntc_output": str(ntc_output),
        "ntc_shape": [
            ntc_cells,
            input_genes,
        ],
        "sample_counts": str(audit_path),
    }

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
        + "\n"
    )

    validate_output(
        scz_output,
        scz_cells,
        input_genes,
        "SCZ",
    )

    validate_output(
        ntc_output,
        ntc_cells,
        input_genes,
        "NTC",
    )

    print()
    print("=" * 100)
    print("SUCCESS")
    print("=" * 100)
    print("SCZ H5AD:", scz_output)
    print("SCZ shape:", scz_adata.shape)
    print("NTC H5AD:", ntc_output)
    print("NTC shape:", ntc_adata.shape)
    print("Total cells:", scz_cells + ntc_cells)
    print("Sample-count audit:", audit_path)
    print("Manifest:", manifest_path)


if __name__ == "__main__":
    main()
