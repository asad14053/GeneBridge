#!/usr/bin/env python3

from __future__ import annotations

import argparse
import collections
import gzip
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


# =============================================================================
# RELEASE CHAINS
# =============================================================================

# Ling source-repo strategy:
# E89 is primary, other historical releases mop up residual symbols.
LING_RELEASES = [
    "E89",
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


# Huuki/Kwon are more recent.
# If they already carry ENSG, that wins.
# Only use this chain when they contain symbols only.
MODERN_RELEASES = [
    "G44",
    "G35",
    "G32",
    "G30",
    "G29",
    "E92",
    "E91",
    "E90",
    "E89",
    "G27",
    "G26",
    "E88",
]


# =============================================================================
# BASICS
# =============================================================================

ENSG_RE = re.compile(
    r"^ENSG\d+(?:\.\d+)?$"
)


def strip_ensembl_version(x):

    if pd.isna(x):
        return None

    x = str(x).strip()

    if not x:
        return None

    return x.split(".")[0]


def split_spec(spec):

    s = spec.strip().upper()

    if s.startswith("E"):
        return "ensembl", int(s[1:])

    return "gencode", int(
        s.lstrip("GV")
    )


def spec_label(spec):

    source, n = split_spec(
        spec
    )

    if source == "ensembl":
        return f"Ensembl {n}"

    return f"GENCODE v{n}"


def decode(x):

    return np.asarray([
        v.decode()
        if isinstance(v, (bytes, np.bytes_))
        else str(v)

        for v in x
    ])


# =============================================================================
# GTF PARSING
# =============================================================================

def parse_attributes(s):

    out = {}

    for part in s.rstrip(";").split(";"):

        part = part.strip()

        if not part:
            continue

        pieces = part.split(
            None,
            1
        )

        if len(pieces) != 2:
            continue

        key, value = pieces

        out[key] = (
            value
            .strip()
            .strip('"')
        )

    return out


def annotation_path(
    spec,
    cache_dir
):

    source, n = split_spec(
        spec
    )

    if source == "ensembl":

        return (
            cache_dir
            /
            f"Homo_sapiens.GRCh38.{n}.gtf.gz"
        )

    return (
        cache_dir
        /
        f"gencode.v{n}.annotation.gtf.gz"
    )


def parse_annotation(
    spec,
    cache_dir
):

    source, n = split_spec(
        spec
    )

    tsv = (
        cache_dir
        /
        f"parsed_{source}_{n}.tsv"
    )


    if tsv.exists():

        return pd.read_csv(
            tsv,
            sep="\t",
            dtype=str
        )


    gtf = annotation_path(
        spec,
        cache_dir
    )


    if not gtf.exists():

        raise FileNotFoundError(
            f"Missing annotation:\n{gtf}\n"
            "Run scripts/ling_harmonization/"
            "00_download_gene_annotations.sh first."
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


            stable = strip_ensembl_version(
                gene_id
            )


            # Remove pseudoautosomal chrY duplicates.
            if (
                stable is None
                or
                "_PAR_Y" in stable
            ):
                continue


            gene_version = a.get(
                "gene_version"
            )


            if (
                source == "ensembl"
                and
                gene_version
            ):

                versioned = (
                    f"{stable}.{gene_version}"
                )

            else:

                versioned = gene_id


            rows.append({
                "ensembl_id":
                    stable,

                "ensembl_id_versioned":
                    versioned,

                "gene_symbol":
                    gene_name,

                "gene_type":
                    gene_type,

                "chrom":
                    f[0],
            })


    t = pd.DataFrame(
        rows
    )


    t = t.drop_duplicates(
        [
            "ensembl_id",
            "gene_symbol"
        ]
    )


    t.to_csv(
        tsv,
        sep="\t",
        index=False
    )


    return t


# =============================================================================
# SYMBOL LOOKUP
# =============================================================================

def build_lookup(
    releases,
    cache_dir
):

    out = {}


    for spec in releases:

        print(
            "Loading annotation:",
            spec_label(spec),
            flush=True
        )


        t = parse_annotation(
            spec,
            cache_dir
        )


        inv = collections.defaultdict(
            list
        )


        for r in t.itertuples(
            index=False
        ):

            inv[
                r.gene_symbol
            ].append({
                "ensembl_id":
                    r.ensembl_id,

                "gene_type":
                    r.gene_type,

                "chrom":
                    r.chrom,
            })


        # First release in the chain wins.
        for symbol, records in inv.items():

            if symbol not in out:

                # ENSG-deduplicate within release
                uniq = {}

                for r in records:

                    uniq[
                        r["ensembl_id"]
                    ] = r


                out[
                    symbol
                ] = {
                    "resolved_via":
                        spec_label(spec),

                    "records":
                        list(
                            uniq.values()
                        ),
                }


    return out


# =============================================================================
# LING
# =============================================================================

def read_ling_symbols(
    path
):

    with h5py.File(
        path,
        "r"
    ) as h:

        key = "matrix/features/name"

        if key not in h:

            raise RuntimeError(
                f"{path}: missing {key}"
            )

        return decode(
            h[key][:]
        )


def collect_ling_gene_universe(
    root,
    expected_files
):

    files = sorted(
        root.rglob(
            "*.gmg.donors.h5"
        )
    )


    print()
    print("=" * 90)
    print("LING H5 COLLECTION")
    print("=" * 90)

    print(
        "Files found:",
        len(files)
    )


    if len(files) != expected_files:

        raise RuntimeError(
            f"Expected {expected_files} Ling donor H5 files, "
            f"but found {len(files)}.\n"
            "Finish the NeMO download before running the final harmonization."
        )


    presence = collections.Counter()

    library_rows = []


    for i, path in enumerate(
        files,
        1
    ):

        symbols = read_ling_symbols(
            path
        )

        unique = set(
            symbols
        )


        for symbol in unique:
            presence[symbol] += 1


        library_rows.append({
            "file":
                str(path),

            "n_genes":
                len(symbols),

            "n_unique_symbols":
                len(unique),
        })


        if (
            i == 1
            or
            i % 10 == 0
            or
            i == len(files)
        ):

            print(
                f"{i}/{len(files)} "
                f"{path.name}: "
                f"{len(symbols)} genes"
            )


    union = sorted(
        presence
    )


    intersection = sorted([
        gene
        for gene, n
        in presence.items()
        if n == len(files)
    ])


    prevalence = pd.DataFrame({
        "gene_symbol":
            union,

        "n_libraries":
            [
                presence[g]
                for g in union
            ]
    })


    prevalence[
        "fraction_libraries"
    ] = (
        prevalence[
            "n_libraries"
        ]
        /
        len(files)
    )


    print()
    print(
        "Ling union:",
        len(union)
    )

    print(
        "Ling intersection:",
        len(intersection)
    )


    return (
        files,
        union,
        intersection,
        prevalence,
        pd.DataFrame(
            library_rows
        ),
    )


# =============================================================================
# MAP SYMBOLS
# =============================================================================

def map_symbols(
    symbols,
    lookup,
    dataset
):

    rows = []


    for symbol in symbols:

        symbol = str(
            symbol
        ).strip()


        if ":" in symbol:

            rows.append({
                "dataset":
                    dataset,

                "gene_symbol":
                    symbol,

                "ensembl_id":
                    None,

                "resolved_via":
                    None,

                "is_metagene":
                    True,

                "n_ensembl_for_symbol":
                    0,

                "gene_type":
                    None,

                "chrom":
                    None,
            })

            continue


        hit = lookup.get(
            symbol
        )


        if hit is None:

            rows.append({
                "dataset":
                    dataset,

                "gene_symbol":
                    symbol,

                "ensembl_id":
                    None,

                "resolved_via":
                    None,

                "is_metagene":
                    False,

                "n_ensembl_for_symbol":
                    0,

                "gene_type":
                    None,

                "chrom":
                    None,
            })

            continue


        records = hit[
            "records"
        ]


        for r in records:

            rows.append({
                "dataset":
                    dataset,

                "gene_symbol":
                    symbol,

                "ensembl_id":
                    r[
                        "ensembl_id"
                    ],

                "resolved_via":
                    hit[
                        "resolved_via"
                    ],

                "is_metagene":
                    False,

                "n_ensembl_for_symbol":
                    len(records),

                "gene_type":
                    r[
                        "gene_type"
                    ],

                "chrom":
                    r[
                        "chrom"
                    ],
            })


    x = pd.DataFrame(
        rows
    )


    # Deterministically choose one row per stable ENSG.
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


    return x


# =============================================================================
# H5AD GENE EXTRACTION
# =============================================================================

def detect_ensembl_column(
    var,
    var_names
):

    candidates = []


    # var_names can themselves be ENSG
    s = pd.Series(
        var_names,
        dtype=str
    )


    candidates.append(
        (
            "__var_names__",
            float(
                s.str.match(
                    ENSG_RE
                ).mean()
            ),
            s
        )
    )


    for col in var.columns:

        values = (
            var[
                col
            ]
            .astype(str)
        )


        rate = float(
            values.str.match(
                ENSG_RE
            ).mean()
        )


        candidates.append(
            (
                str(col),
                rate,
                values
            )
        )


    candidates.sort(
        key=lambda x: x[1],
        reverse=True
    )


    best_name, best_rate, best_values = (
        candidates[0]
    )


    if best_rate >= 0.80:

        return (
            best_name,
            best_rate,
            best_values
        )


    return (
        None,
        best_rate,
        None
    )


def detect_symbol_column(
    var,
    var_names
):

    preferred = [
        "gene_name",
        "gene_symbol",
        "symbol",
        "gene",
        "external_gene_name",
    ]


    for col in preferred:

        if col in var.columns:

            return (
                col,
                var[
                    col
                ]
                .astype(str)
            )


    return (
        "__var_names__",
        pd.Series(
            var_names,
            dtype=str
        )
    )


def extract_h5ad_gene_table(
    path,
    dataset,
    lookup_if_needed
):

    print()
    print("=" * 90)
    print(
        f"{dataset.upper()} GENE EXTRACTION"
    )
    print("=" * 90)


    a = ad.read_h5ad(
        path,
        backed="r"
    )


    var = a.var.copy()

    var_names = (
        a.var_names
        .astype(str)
        .to_numpy()
    )


    print(
        "Path:",
        path
    )

    print(
        "Shape:",
        a.shape
    )

    print(
        "var columns:",
        list(
            var.columns
        )
    )


    ens_col, ens_rate, ens_values = (
        detect_ensembl_column(
            var,
            var_names
        )
    )


    symbol_col, symbols = (
        detect_symbol_column(
            var,
            var_names
        )
    )


    a.file.close()


    print(
        "Symbol source:",
        symbol_col
    )

    print(
        "Best ENSG source:",
        ens_col
    )

    print(
        "Best ENSG fraction:",
        f"{ens_rate:.4f}"
    )


    symbols = (
        symbols.astype(str)
        .to_numpy()
    )


    # -------------------------------------------------------------------------
    # Direct stable-ID route
    # -------------------------------------------------------------------------

    if ens_col is not None:

        ens = (
            ens_values.astype(str)
            .map(
                strip_ensembl_version
            )
            .to_numpy()
        )


        x = pd.DataFrame({
            "dataset":
                dataset,

            "gene_symbol":
                symbols,

            "ensembl_id":
                ens,

            "resolved_via":
                f"direct:{ens_col}",

            "is_metagene":
                False,

            "n_ensembl_for_symbol":
                1,
        })


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


        return x


    # -------------------------------------------------------------------------
    # Symbol route
    # -------------------------------------------------------------------------

    print(
        f"{dataset}: no sufficiently complete ENSG column; "
        "mapping symbols through annotation releases."
    )


    return map_symbols(
        symbols,
        lookup_if_needed,
        dataset
    )


# =============================================================================
# OVERLAP
# =============================================================================

def preferred_set(
    x
):

    return set(
        x.loc[
            x[
                "preferred_for_ensembl"
            ]
            &
            x[
                "ensembl_id"
            ].notna(),
            "ensembl_id"
        ]
    )


def pair_summary(
    name_a,
    set_a,
    name_b,
    set_b
):

    overlap = (
        set_a
        &
        set_b
    )


    return {
        "comparison":
            f"{name_a}_vs_{name_b}",

        "n_A":
            len(set_a),

        "n_B":
            len(set_b),

        "n_overlap":
            len(overlap),

        "pct_A_recovered_in_B":
            (
                100
                *
                len(overlap)
                /
                len(set_a)
                if set_a
                else np.nan
            ),

        "pct_B_recovered_in_A":
            (
                100
                *
                len(overlap)
                /
                len(set_b)
                if set_b
                else np.nan
            ),
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    p = argparse.ArgumentParser()


    p.add_argument(
        "--ling-h5-root",
        type=Path,
        required=True
    )


    p.add_argument(
        "--huuki-h5ad",
        type=Path,
        required=True
    )


    p.add_argument(
        "--kwon-h5ad",
        type=Path,
        required=True
    )


    p.add_argument(
        "--cache-dir",
        type=Path,
        required=True
    )


    p.add_argument(
        "--out-dir",
        type=Path,
        required=True
    )


    p.add_argument(
        "--expected-ling-files",
        type=int,
        default=87
    )


    args = p.parse_args()


    args.out_dir.mkdir(
        parents=True,
        exist_ok=True
    )


    # -------------------------------------------------------------------------
    # Ling union
    # -------------------------------------------------------------------------

    (
        ling_files,
        ling_union,
        ling_intersection,
        ling_prevalence,
        ling_libraries,
    ) = collect_ling_gene_universe(
        args.ling_h5_root,
        args.expected_ling_files
    )


    ling_libraries.to_csv(
        args.out_dir
        /
        "Ling_library_gene_counts.csv",
        index=False
    )


    ling_prevalence.to_csv(
        args.out_dir
        /
        "Ling_gene_library_prevalence.csv",
        index=False
    )


    pd.DataFrame({
        "gene_symbol":
            ling_intersection
    }).to_csv(
        args.out_dir
        /
        "Ling_all87_intersection_symbols.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Annotation lookup tables
    # -------------------------------------------------------------------------

    print()
    print("=" * 90)
    print("BUILD LING ANNOTATION LOOKUP")
    print("=" * 90)


    ling_lookup = build_lookup(
        LING_RELEASES,
        args.cache_dir
    )


    ling = map_symbols(
        ling_union,
        ling_lookup,
        "Ling"
    )


    ling = ling.merge(
        ling_prevalence,
        on="gene_symbol",
        how="left"
    )


    ling.to_csv(
        args.out_dir
        /
        "Ling_union_symbol_to_ensembl.csv",
        index=False
    )


    print()
    print("=" * 90)
    print("BUILD MODERN ANNOTATION LOOKUP")
    print("=" * 90)


    modern_lookup = build_lookup(
        MODERN_RELEASES,
        args.cache_dir
    )


    # -------------------------------------------------------------------------
    # Huuki
    # -------------------------------------------------------------------------

    huuki = extract_h5ad_gene_table(
        args.huuki_h5ad,
        "Huuki",
        modern_lookup
    )


    huuki.to_csv(
        args.out_dir
        /
        "Huuki_gene_to_ensembl.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Kwon
    # -------------------------------------------------------------------------

    kwon = extract_h5ad_gene_table(
        args.kwon_h5ad,
        "Kwon",
        modern_lookup
    )


    kwon.to_csv(
        args.out_dir
        /
        "Kwon_gene_to_ensembl.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Stable sets
    # -------------------------------------------------------------------------

    L = preferred_set(
        ling
    )

    H = preferred_set(
        huuki
    )

    K = preferred_set(
        kwon
    )


    LH = L & H
    LK = L & K
    HK = H & K

    LHK = (
        L
        &
        H
        &
        K
    )


    # -------------------------------------------------------------------------
    # Pairwise summary
    # -------------------------------------------------------------------------

    summary = pd.DataFrame([
        pair_summary(
            "Ling",
            L,
            "Huuki",
            H
        ),

        pair_summary(
            "Ling",
            L,
            "Kwon",
            K
        ),

        pair_summary(
            "Huuki",
            H,
            "Kwon",
            K
        ),
    ])


    summary.to_csv(
        args.out_dir
        /
        "pairwise_gene_overlap_summary.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Full presence table
    # -------------------------------------------------------------------------

    all_ids = sorted(
        L
        |
        H
        |
        K
    )


    presence = pd.DataFrame({
        "ensembl_id":
            all_ids
    })


    presence[
        "Ling"
    ] = presence[
        "ensembl_id"
    ].isin(
        L
    )


    presence[
        "Huuki"
    ] = presence[
        "ensembl_id"
    ].isin(
        H
    )


    presence[
        "Kwon"
    ] = presence[
        "ensembl_id"
    ].isin(
        K
    )


    presence[
        "n_datasets"
    ] = (
        presence[
            [
                "Ling",
                "Huuki",
                "Kwon"
            ]
        ]
        .sum(
            axis=1
        )
    )


    presence.to_csv(
        args.out_dir
        /
        "Ling_Huuki_Kwon_ensembl_presence_matrix.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Three-way table with symbols
    # -------------------------------------------------------------------------

    def symbol_map(x):

        return (
            x[
                x[
                    "preferred_for_ensembl"
                ]
                &
                x[
                    "ensembl_id"
                ].notna()
            ]
            .drop_duplicates(
                "ensembl_id"
            )
            .set_index(
                "ensembl_id"
            )[
                "gene_symbol"
            ]
            .to_dict()
        )


    lsym = symbol_map(
        ling
    )

    hsym = symbol_map(
        huuki
    )

    ksym = symbol_map(
        kwon
    )


    threeway = pd.DataFrame({
        "ensembl_id":
            sorted(
                LHK
            )
    })


    threeway[
        "Ling_symbol"
    ] = threeway[
        "ensembl_id"
    ].map(
        lsym
    )


    threeway[
        "Huuki_symbol"
    ] = threeway[
        "ensembl_id"
    ].map(
        hsym
    )


    threeway[
        "Kwon_symbol"
    ] = threeway[
        "ensembl_id"
    ].map(
        ksym
    )


    threeway.to_csv(
        args.out_dir
        /
        "Ling_Huuki_Kwon_threeway_common_genes.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Missing lists
    # -------------------------------------------------------------------------

    pd.DataFrame({
        "ensembl_id":
            sorted(
                H - L
            )
    }).to_csv(
        args.out_dir
        /
        "Huuki_not_in_Ling.csv",
        index=False
    )


    pd.DataFrame({
        "ensembl_id":
            sorted(
                K - L
            )
    }).to_csv(
        args.out_dir
        /
        "Kwon_not_in_Ling.csv",
        index=False
    )


    pd.DataFrame({
        "ensembl_id":
            sorted(
                K - H
            )
    }).to_csv(
        args.out_dir
        /
        "Kwon_not_in_Huuki.csv",
        index=False
    )


    # -------------------------------------------------------------------------
    # Overall QC
    # -------------------------------------------------------------------------

    ling_symbol_level = (
        ling
        .drop_duplicates(
            "gene_symbol"
        )
    )


    overall = pd.DataFrame([
        {
            "Ling_H5_files":
                len(
                    ling_files
                ),

            "Ling_union_symbols":
                len(
                    ling_union
                ),

            "Ling_all87_intersection_symbols":
                len(
                    ling_intersection
                ),

            "Ling_mapped_unique_ENSG":
                len(
                    L
                ),

            "Huuki_unique_ENSG":
                len(
                    H
                ),

            "Kwon_unique_ENSG":
                len(
                    K
                ),

            "Ling_Huuki_overlap":
                len(
                    LH
                ),

            "Ling_Kwon_overlap":
                len(
                    LK
                ),

            "Huuki_Kwon_overlap":
                len(
                    HK
                ),

            "threeway_common":
                len(
                    LHK
                ),

            "Kwon_percent_in_Ling":
                (
                    100
                    *
                    len(
                        LK
                    )
                    /
                    len(
                        K
                    )
                ),

            "Kwon_percent_in_Huuki":
                (
                    100
                    *
                    len(
                        HK
                    )
                    /
                    len(
                        K
                    )
                ),

            "Kwon_percent_in_both_Ling_Huuki":
                (
                    100
                    *
                    len(
                        LHK
                    )
                    /
                    len(
                        K
                    )
                ),
        }
    ])


    overall.to_csv(
        args.out_dir
        /
        "gene_harmonization_overall_summary.csv",
        index=False
    )


    print()
    print("=" * 90)
    print("FINAL GENE HARMONIZATION SUMMARY")
    print("=" * 90)

    print(
        overall.to_string(
            index=False
        )
    )


    print()
    print("PAIRWISE:")
    print(
        summary.to_string(
            index=False
        )
    )


    print()
    print(
        "Ling resolution:"
    )

    print(
        ling_symbol_level[
            "resolved_via"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )


    print()
    print(
        "Output:",
        args.out_dir
    )


if __name__ == "__main__":
    main()
