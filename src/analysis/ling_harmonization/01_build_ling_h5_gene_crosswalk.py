#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections
import gzip
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# =============================================================================
# Ling repo annotation strategy
# =============================================================================

PRIMARY = "E89"

FALLBACKS = [
    "E90",
    "G27",
    "E88",
    "E91",
    "E92",
    "G26",
    "G29",
    "G30",
    "G32",
    "G35",
    "G44",
]


GENCODE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/"
    "Gencode_human/release_{v}/"
    "gencode.v{v}.annotation.gtf.gz"
)

ENSEMBL_URL = (
    "https://ftp.ensembl.org/pub/release-{v}/gtf/"
    "homo_sapiens/Homo_sapiens.GRCh38.{v}.gtf.gz"
)


# =============================================================================
# Helpers
# =============================================================================

def split_spec(spec):

    spec = spec.upper()

    if spec.startswith("E"):
        return "ensembl", int(spec[1:])

    return "gencode", int(
        spec.lstrip("GV")
    )


def spec_label(spec):

    source, n = split_spec(spec)

    if source == "ensembl":
        return f"Ensembl {n}"

    return f"GENCODE v{n}"


def decode_array(x):

    return np.asarray([
        v.decode()
        if isinstance(v, (bytes, np.bytes_))
        else str(v)

        for v in x
    ])


# =============================================================================
# Read Ling H5 genes
# =============================================================================

def read_ling_genes(path):

    with h5py.File(
        path,
        "r"
    ) as h:

        key = "matrix/features/name"

        if key not in h:

            raise RuntimeError(
                f"{key} not present in {path}"
            )

        return decode_array(
            h[key][:]
        )


# =============================================================================
# GTF parsing
# =============================================================================

def parse_attributes(s):

    out = {}

    for item in s.strip().strip(";").split(";"):

        item = item.strip()

        if not item:
            continue

        parts = item.split(
            None,
            1
        )

        if len(parts) != 2:
            continue

        key, value = parts

        out[key] = (
            value
            .strip()
            .strip('"')
        )

    return out


def fetch_gtf(
    spec,
    cache_dir
):

    cache_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    source, n = split_spec(
        spec
    )

    if source == "ensembl":

        path = (
            cache_dir
            /
            f"Homo_sapiens.GRCh38.{n}.gtf.gz"
        )

        url = ENSEMBL_URL.format(
            v=n
        )

    else:

        path = (
            cache_dir
            /
            f"gencode.v{n}.annotation.gtf.gz"
        )

        url = GENCODE_URL.format(
            v=n
        )


    if (
        path.exists()
        and
        path.stat().st_size > 1_000_000
    ):

        return path


    print(
        f"Downloading {spec_label(spec)}",
        file=sys.stderr
    )

    subprocess.run(
        [
            "curl",
            "-L",
            "--fail",
            "--retry",
            "3",
            "--max-time",
            "900",
            "-o",
            str(path),
            url,
        ],
        check=True
    )

    return path


def parse_gene_table(
    spec,
    cache_dir
):

    source, n = split_spec(
        spec
    )

    cache_tsv = (
        cache_dir
        /
        f"genes_{source}_{n}.tsv"
    )


    if cache_tsv.exists():

        return pd.read_csv(
            cache_tsv,
            sep="\t",
            dtype=str
        )


    gtf = fetch_gtf(
        spec,
        cache_dir
    )


    rows = []


    with gzip.open(
        gtf,
        "rt"
    ) as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            f = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if (
                len(f) < 9
                or
                f[2] != "gene"
            ):
                continue


            a = parse_attributes(
                f[8]
            )


            gene_id = a.get(
                "gene_id"
            )

            gene_name = a.get(
                "gene_name"
            )


            if (
                gene_id is None
                or
                gene_name is None
            ):
                continue


            gene_type = (
                a.get("gene_type")
                or
                a.get("gene_biotype")
            )


            if source == "ensembl":

                gene_version = a.get(
                    "gene_version"
                )

                if gene_version:

                    versioned = (
                        f"{gene_id}."
                        f"{gene_version}"
                    )

                else:

                    versioned = gene_id

            else:

                versioned = gene_id


            # stable ID
            stable = gene_id.split(
                "."
            )[0]


            # Drop PAR_Y duplicate
            if "_PAR_Y" in stable:
                continue


            rows.append({
                "ensembl_id":
                    stable,

                "ensembl_id_versioned":
                    versioned,

                "gene_name":
                    gene_name,

                "gene_type":
                    gene_type,

                "chrom":
                    f[0],
            })


    df = pd.DataFrame(
        rows
    )


    df = df.drop_duplicates(
        [
            "ensembl_id",
            "gene_name",
        ]
    )


    df.to_csv(
        cache_tsv,
        sep="\t",
        index=False
    )


    return df


# =============================================================================
# Check all 87 Ling H5 files
# =============================================================================

def check_h5_gene_universe(
    files,
    out_dir
):

    reference_file = files[0]

    reference = read_ling_genes(
        reference_file
    )

    reference_set = set(
        reference
    )


    rows = []


    for path in files:

        genes = read_ling_genes(
            path
        )


        same_order = np.array_equal(
            genes,
            reference
        )


        same_set = (
            len(genes)
            ==
            len(reference)
            and
            set(genes)
            ==
            reference_set
        )


        rows.append({
            "file":
                str(path),

            "n_genes":
                len(genes),

            "same_gene_set":
                same_set,

            "same_gene_order":
                same_order,

            "missing_vs_reference":
                len(
                    reference_set
                    -
                    set(genes)
                ),

            "extra_vs_reference":
                len(
                    set(genes)
                    -
                    reference_set
                ),
        })


    qc = pd.DataFrame(
        rows
    )


    qc.to_csv(
        out_dir
        /
        "Ling_H5_gene_universe_QC.csv",
        index=False
    )


    print(
        "\n===== H5 GENE UNIVERSE QC ====="
    )

    print(
        "Files:",
        len(files)
    )

    print(
        "Reference:",
        reference_file.name
    )

    print(
        "Reference genes:",
        len(reference)
    )

    print(
        "Same gene set:",
        int(
            qc[
                "same_gene_set"
            ].sum()
        ),
        "/",
        len(qc)
    )

    print(
        "Same gene order:",
        int(
            qc[
                "same_gene_order"
            ].sum()
        ),
        "/",
        len(qc)
    )


    if not qc[
        "same_gene_set"
    ].all():

        raise RuntimeError(
            "Ling H5 files do not all "
            "have the same gene universe."
        )


    return reference


# =============================================================================
# Build release-specific symbol dictionary
# =============================================================================

def release_symbol_map(
    spec,
    cache_dir
):

    genes = parse_gene_table(
        spec,
        cache_dir
    )


    d = collections.defaultdict(
        list
    )


    for row in genes.itertuples():

        d[
            row.gene_name
        ].append({
            "ensembl_id":
                row.ensembl_id,

            "ensembl_id_versioned":
                row.ensembl_id_versioned,

            "gene_type":
                row.gene_type,

            "chrom":
                row.chrom,
        })


    return d


# =============================================================================
# Map Ling symbols
# =============================================================================

def build_crosswalk(
    symbols,
    cache_dir
):

    releases = [
        PRIMARY
    ] + FALLBACKS


    unresolved = set(
        symbols
    )


    resolved = {}


    # -------------------------------------------------------------------------
    # Metagenes are explicitly NOT assigned one ENSG
    # -------------------------------------------------------------------------

    metagenes = {
        s
        for s in symbols
        if ":" in s
    }


    unresolved -= metagenes


    # -------------------------------------------------------------------------
    # First match E89, then fallback releases
    # -------------------------------------------------------------------------

    for spec in releases:

        if not unresolved:
            break


        print(
            f"\nResolving with {spec_label(spec)} "
            f"| currently unresolved: "
            f"{len(unresolved)}"
        )


        lookup = release_symbol_map(
            spec,
            cache_dir
        )


        hits = (
            unresolved
            &
            set(lookup.keys())
        )


        print(
            "resolved in this release:",
            len(hits)
        )


        for symbol in hits:

            resolved[
                symbol
            ] = {
                "release":
                    spec,

                "records":
                    lookup[
                        symbol
                    ],
            }


        unresolved -= hits


    # -------------------------------------------------------------------------
    # Build output rows
    # -------------------------------------------------------------------------

    rows = []


    for symbol in symbols:

        if symbol in metagenes:

            rows.append({
                "gene_symbol":
                    symbol,

                "ensembl_id":
                    None,

                "ensembl_id_versioned":
                    None,

                "resolved_via":
                    None,

                "gene_type":
                    None,

                "chrom":
                    None,

                "is_metagene":
                    True,

                "n_ensembl_for_symbol":
                    0,
            })

            continue


        if symbol not in resolved:

            rows.append({
                "gene_symbol":
                    symbol,

                "ensembl_id":
                    None,

                "ensembl_id_versioned":
                    None,

                "resolved_via":
                    None,

                "gene_type":
                    None,

                "chrom":
                    None,

                "is_metagene":
                    False,

                "n_ensembl_for_symbol":
                    0,
            })

            continue


        hit = resolved[
            symbol
        ]


        records = hit[
            "records"
        ]


        # remove duplicate ENSG hits
        unique = {}

        for r in records:

            unique[
                r[
                    "ensembl_id"
                ]
            ] = r


        records = list(
            unique.values()
        )


        for r in records:

            rows.append({
                "gene_symbol":
                    symbol,

                "ensembl_id":
                    r[
                        "ensembl_id"
                    ],

                "ensembl_id_versioned":
                    r[
                        "ensembl_id_versioned"
                    ],

                "resolved_via":
                    spec_label(
                        hit[
                            "release"
                        ]
                    ),

                "gene_type":
                    r[
                        "gene_type"
                    ],

                "chrom":
                    r[
                        "chrom"
                    ],

                "is_metagene":
                    False,

                "n_ensembl_for_symbol":
                    len(
                        records
                    ),
            })


    x = pd.DataFrame(
        rows
    )


    # -------------------------------------------------------------------------
    # Preferred row for each ENSG
    #
    # E89 wins over fallback mapping.
    # -------------------------------------------------------------------------

    primary_label = spec_label(
        PRIMARY
    )


    x["_primary_rank"] = np.where(
        x[
            "resolved_via"
        ]
        ==
        primary_label,
        0,
        1
    )


    x = x.sort_values(
        [
            "ensembl_id",
            "_primary_rank",
            "gene_symbol",
        ],
        na_position="last"
    )


    x[
        "preferred_for_ensembl"
    ] = (
        x[
            "ensembl_id"
        ].notna()
        &
        ~x[
            "ensembl_id"
        ].duplicated(
            keep="first"
        )
    )


    x = x.drop(
        columns=[
            "_primary_rank"
        ]
    )


    return (
        x,
        sorted(
            metagenes
        ),
        sorted(
            unresolved
        ),
    )


# =============================================================================
# Main
# =============================================================================

def main():

    ap = argparse.ArgumentParser()


    ap.add_argument(
        "--h5-root",
        required=True,
        type=Path
    )


    ap.add_argument(
        "--cache-dir",
        required=True,
        type=Path
    )


    ap.add_argument(
        "--out-dir",
        required=True,
        type=Path
    )


    args = ap.parse_args()


    args.out_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # -------------------------------------------------------------------------
    # Find the 87 H5 files
    # -------------------------------------------------------------------------

    files = sorted(
        args.h5_root.rglob(
            "*.gmg.donors.h5"
        )
    )


    if not files:

        raise FileNotFoundError(
            f"No Ling H5 files found under "
            f"{args.h5_root}"
        )


    print(
        "Ling H5 files found:",
        len(files)
    )


    if len(files) != 87:

        print(
            "WARNING: expected 87 files, found",
            len(files),
            file=sys.stderr
        )


    # -------------------------------------------------------------------------
    # Verify gene universe
    # -------------------------------------------------------------------------

    genes = check_h5_gene_universe(
        files,
        args.out_dir
    )


    symbols = list(
        dict.fromkeys(
            genes.tolist()
        )
    )


    pd.DataFrame({
        "ling_gene_index":
            np.arange(
                len(genes)
            ),

        "gene_symbol":
            genes,
    }).to_csv(
        args.out_dir
        /
        "Ling_H5_master_gene_symbols.csv",
        index=False
    )


    print(
        "\nUnique Ling symbols:",
        len(symbols)
    )


    # -------------------------------------------------------------------------
    # Crosswalk
    # -------------------------------------------------------------------------

    crosswalk, metagenes, unresolved = (
        build_crosswalk(
            symbols,
            args.cache_dir
        )
    )


    crosswalk.to_csv(
        args.out_dir
        /
        "Ling_H5_symbol_to_ensembl_crosswalk.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Clean unique ENSG crosswalk for downstream merging
    # -------------------------------------------------------------------------

    preferred = crosswalk[
        crosswalk[
            "preferred_for_ensembl"
        ]
    ].copy()


    preferred.to_csv(
        args.out_dir
        /
        "Ling_H5_preferred_unique_ensembl.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Unresolved
    # -------------------------------------------------------------------------

    crosswalk[
        crosswalk[
            "ensembl_id"
        ].isna()
        &
        ~crosswalk[
            "is_metagene"
        ]
    ].to_csv(
        args.out_dir
        /
        "Ling_H5_unresolved_symbols.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Ambiguous
    # -------------------------------------------------------------------------

    crosswalk[
        crosswalk[
            "n_ensembl_for_symbol"
        ]
        >
        1
    ].to_csv(
        args.out_dir
        /
        "Ling_H5_ambiguous_symbol_mappings.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Metagenes
    # -------------------------------------------------------------------------

    pd.DataFrame({
        "metagene":
            metagenes
    }).to_csv(
        args.out_dir
        /
        "Ling_H5_metagenes.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    symbol_level = (
        crosswalk
        .drop_duplicates(
            "gene_symbol"
        )
    )


    nonmeta = symbol_level[
        ~symbol_level[
            "is_metagene"
        ]
    ]


    n_mapped = int(
        nonmeta[
            "ensembl_id"
        ].notna().sum()
    )


    n_nonmeta = len(
        nonmeta
    )


    e89_count = int(
        (
            symbol_level[
                "resolved_via"
            ]
            ==
            spec_label(
                PRIMARY
            )
        ).sum()
    )


    fallback_count = int(
        (
            symbol_level[
                "resolved_via"
            ].notna()
            &
            (
                symbol_level[
                    "resolved_via"
                ]
                !=
                spec_label(
                    PRIMARY
                )
            )
        ).sum()
    )


    summary = pd.DataFrame([
        {
            "n_h5_files":
                len(files),

            "n_genes_per_h5":
                len(genes),

            "n_unique_gene_symbols":
                len(symbols),

            "n_metagenes":
                len(metagenes),

            "n_nonmetagene_symbols":
                n_nonmeta,

            "mapped_nonmetagene_symbols":
                n_mapped,

            "mapping_percent":
                (
                    100
                    *
                    n_mapped
                    /
                    n_nonmeta
                ),

            "resolved_by_Ensembl89":
                e89_count,

            "resolved_by_fallback":
                fallback_count,

            "unresolved_nonmetagene":
                len(unresolved),

            "preferred_unique_ensembl":
                preferred[
                    "ensembl_id"
                ].nunique(),
        }
    ])


    summary.to_csv(
        args.out_dir
        /
        "Ling_H5_gene_crosswalk_summary.csv",
        index=False
    )


    print(
        "\n"
        +
        "=" * 100
    )

    print(
        "LING H5 GENE CROSSWALK SUMMARY"
    )

    print(
        "=" * 100
    )

    print(
        summary.to_string(
            index=False
        )
    )


    print(
        "\nResolution by release:"
    )

    print(
        crosswalk[
            "resolved_via"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )


    print(
        "\nMain downstream merge file:"
    )

    print(
        args.out_dir
        /
        "Ling_H5_preferred_unique_ensembl.csv"
    )


    print(
        "\nIMPORTANT:"
    )

    print(
        "Merge Huuki-Myers / Kwon / Guo "
        "using version-stripped ensembl_id, "
        "not gene_symbol."
    )


if __name__ == "__main__":
    main()
