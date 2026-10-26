#!/usr/bin/env python3

from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

COUNTS = (
    ROOT
    / "data/raw/nemo_mccarroll/counts"
)

META = (
    COUNTS
    / "BA46.jointMetadata.txt"
)

OUT = (
    ROOT
    / "results/analysis/ling_reference_build"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


def decode(x):
    return np.asarray([
        v.decode()
        if isinstance(v, (bytes, np.bytes_))
        else str(v)
        for v in x
    ])


def file_prefix(path):

    # Example:
    # BA46_05-07-2019_rxn2_MyelinRemoved.intronic_...
    prefix = path.name.split(
        ".intronic_"
    )[0]

    prefix = re.sub(
        r"_MyelinRemoved$",
        "",
        prefix
    )

    return prefix


print("=" * 100)
print("LOAD LING CELL METADATA")
print("=" * 100)

meta = pd.read_csv(
    META,
    sep="\t",
    low_memory=False
)

print("Rows:", f"{len(meta):,}")
print("Prefixes:", meta["PREFIX"].nunique())
print("Donors:", meta["DONOR"].nunique())


meta["PREFIX"] = (
    meta["PREFIX"]
    .astype(str)
    .str.strip()
)

meta["CELL_BARCODE"] = (
    meta["CELL_BARCODE"]
    .astype(str)
    .str.strip()
)


# ----------------------------------------------------------------------
# Check uniqueness within library
# ----------------------------------------------------------------------

dup = meta.duplicated(
    [
        "PREFIX",
        "CELL_BARCODE"
    ]
)

print(
    "Duplicated PREFIX+BARCODE:",
    int(dup.sum())
)

if dup.any():
    raise RuntimeError(
        "Metadata has duplicated PREFIX + CELL_BARCODE."
    )


# ----------------------------------------------------------------------
# Index metadata by prefix
# ----------------------------------------------------------------------

meta_by_prefix = {
    prefix: df
    for prefix, df
    in meta.groupby(
        "PREFIX",
        sort=False
    )
}


# ----------------------------------------------------------------------
# Find 87 donor H5 files
# ----------------------------------------------------------------------

files = sorted(
    COUNTS.glob(
        "*.gmg.donors.h5"
    )
)

print()
print("=" * 100)
print("H5 LIBRARIES")
print("=" * 100)

print(
    "Files:",
    len(files)
)

if len(files) != 87:
    raise RuntimeError(
        f"Expected 87 H5 files; found {len(files)}"
    )


rows = []

all_cell_ids = set()


for i, path in enumerate(
    files,
    start=1
):

    prefix = file_prefix(
        path
    )


    with h5py.File(
        path,
        "r"
    ) as h:

        barcodes = decode(
            h[
                "matrix/barcodes"
            ][:]
        )


    h5_barcodes = set(
        barcodes
    )


    if prefix not in meta_by_prefix:

        rows.append({
            "file":
                path.name,

            "prefix":
                prefix,

            "h5_cells":
                len(barcodes),

            "metadata_cells_for_prefix":
                0,

            "shared":
                0,

            "h5_missing_in_metadata":
                len(barcodes),

            "metadata_not_in_h5":
                0,

            "coverage_percent":
                0.0,

            "n_donors":
                0,
        })

        print(
            f"{i:02d}/87 FAIL "
            f"{prefix}: metadata prefix missing"
        )

        continue


    m = meta_by_prefix[
        prefix
    ]


    meta_barcodes = set(
        m[
            "CELL_BARCODE"
        ]
    )


    shared = (
        h5_barcodes
        &
        meta_barcodes
    )


    missing = (
        h5_barcodes
        -
        meta_barcodes
    )


    metadata_extra = (
        meta_barcodes
        -
        h5_barcodes
    )


    n_donors = (
        m.loc[
            m[
                "CELL_BARCODE"
            ].isin(
                h5_barcodes
            ),
            "DONOR"
        ]
        .nunique()
    )


    coverage = (
        100
        *
        len(shared)
        /
        len(h5_barcodes)
        if len(h5_barcodes)
        else 0
    )


    rows.append({
        "file":
            path.name,

        "prefix":
            prefix,

        "h5_cells":
            len(barcodes),

        "metadata_cells_for_prefix":
            len(m),

        "shared":
            len(shared),

        "h5_missing_in_metadata":
            len(missing),

        "metadata_not_in_h5":
            len(metadata_extra),

        "coverage_percent":
            coverage,

        "n_donors":
            n_donors,
    })


    # Global nucleus-ID uniqueness
    for barcode in barcodes:

        cell_id = (
            prefix
            + "_"
            + barcode
        )

        if cell_id in all_cell_ids:

            raise RuntimeError(
                f"Duplicate global nucleus ID: {cell_id}"
            )

        all_cell_ids.add(
            cell_id
        )


    print(
        f"{i:02d}/87 "
        f"{prefix:<30} "
        f"H5={len(barcodes):>7,} "
        f"shared={len(shared):>7,} "
        f"missing={len(missing):>5,} "
        f"donors={n_donors:>3}"
    )


qc = pd.DataFrame(
    rows
)


qc.to_csv(
    OUT
    / "ling_h5_metadata_alignment.csv",
    index=False
)


print()
print("=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    "H5 nuclei total:",
    f"{qc['h5_cells'].sum():,}"
)

print(
    "Shared metadata nuclei:",
    f"{qc['shared'].sum():,}"
)

print(
    "H5 nuclei without metadata:",
    f"{qc['h5_missing_in_metadata'].sum():,}"
)

print(
    "Libraries with 100% H5 metadata coverage:",
    int(
        (
            qc[
                "coverage_percent"
            ]
            ==
            100
        ).sum()
    ),
    "/",
    len(qc)
)

print(
    "Minimum H5→metadata coverage:",
    qc[
        "coverage_percent"
    ].min()
)

print(
    "Unique global H5 nucleus IDs:",
    f"{len(all_cell_ids):,}"
)


summary = pd.DataFrame([
    {
        "n_libraries":
            len(qc),

        "total_h5_nuclei":
            int(
                qc[
                    "h5_cells"
                ].sum()
            ),

        "shared_metadata_nuclei":
            int(
                qc[
                    "shared"
                ].sum()
            ),

        "h5_without_metadata":
            int(
                qc[
                    "h5_missing_in_metadata"
                ].sum()
            ),

        "libraries_100pct_metadata":
            int(
                (
                    qc[
                        "coverage_percent"
                    ]
                    ==
                    100
                ).sum()
            ),

        "minimum_coverage_percent":
            float(
                qc[
                    "coverage_percent"
                ].min()
            ),

        "global_unique_cell_ids":
            len(
                all_cell_ids
            ),
    }
])


summary.to_csv(
    OUT
    / "ling_h5_metadata_alignment_summary.csv",
    index=False
)


print()
print(
    summary.to_string(
        index=False
    )
)
