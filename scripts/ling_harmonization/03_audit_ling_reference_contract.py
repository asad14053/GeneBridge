#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import json

import anndata as ad
import h5py
import numpy as np
import pandas as pd


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

LING_DIR = (
    ROOT
    / "data/raw/nemo_mccarroll/counts"
)

HUUKI = (
    ROOT
    / "data/processed/snrnaseq/sce_DLPFC_annotated/"
      "huuki_snrna_reference_full_allgenes.h5ad"
)

HARM = (
    ROOT
    / "results/analysis/"
      "ling_huuki_kwon_gene_harmonization"
)

LING_XWALK = (
    HARM
    / "Ling_union_symbol_to_ensembl.csv"
)

HUUKI_XWALK = (
    HARM
    / "Huuki_gene_to_ensembl.csv"
)

KWON_XWALK = (
    HARM
    / "Kwon_gene_to_ensembl.csv"
)

OUT = (
    ROOT
    / "results/analysis/"
      "ling_reference_build"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


def section(x):
    print(
        "\n" + "=" * 100
    )
    print(x)
    print("=" * 100)


def decode(values):
    return np.asarray([
        x.decode()
        if isinstance(
            x,
            (bytes, np.bytes_)
        )
        else str(x)

        for x in values
    ])


def true_series(x):
    return (
        x.astype(str)
        .str.lower()
        .isin([
            "true",
            "1",
            "yes"
        ])
    )


# =============================================================================
# HUUKI EXACT STRUCTURE
# =============================================================================

section(
    "HUUKI TEMPLATE STRUCTURE"
)

if not HUUKI.exists():
    raise FileNotFoundError(
        HUUKI
    )

huuki = ad.read_h5ad(
    HUUKI,
    backed="r"
)

contract = {
    "file":
        str(HUUKI),

    "n_obs":
        int(
            huuki.n_obs
        ),

    "n_vars":
        int(
            huuki.n_vars
        ),

    "X_dtype":
        str(
            huuki.X.dtype
        ),

    "layers":
        list(
            huuki.layers.keys()
        ),

    "obs_columns":
        list(
            huuki.obs.columns
        ),

    "var_columns":
        list(
            huuki.var.columns
        ),

    "obsm_keys":
        list(
            huuki.obsm.keys()
        ),

    "var_names_unique":
        bool(
            huuki.var_names.is_unique
        ),

    "obs_names_unique":
        bool(
            huuki.obs_names.is_unique
        ),
}


print(
    "Shape:",
    huuki.shape
)

print(
    "X dtype:",
    huuki.X.dtype
)

print(
    "\nLayers:"
)

print(
    list(
        huuki.layers.keys()
    )
)

print(
    "\nOBS columns:"
)

for x in huuki.obs.columns:
    print(
        " ",
        x
    )

print(
    "\nVAR columns:"
)

for x in huuki.var.columns:
    print(
        " ",
        x
    )

print(
    "\nOBSM:"
)

for x in huuki.obsm.keys():

    try:
        print(
            " ",
            x,
            huuki.obsm[
                x
            ].shape
        )

    except Exception:
        print(
            " ",
            x
        )


with open(
    OUT
    / "huuki_reference_contract.json",
    "w"
) as fh:

    json.dump(
        contract,
        fh,
        indent=2
    )


huuki.file.close()


# =============================================================================
# LOAD EXISTING HARMONIZATION
# =============================================================================

section(
    "LOAD GENE HARMONIZATION"
)

for path in [
    LING_XWALK,
    HUUKI_XWALK,
    KWON_XWALK,
]:

    if not path.exists():
        raise FileNotFoundError(
            path
        )


ling = pd.read_csv(
    LING_XWALK
)

huuki_map = pd.read_csv(
    HUUKI_XWALK
)

kwon = pd.read_csv(
    KWON_XWALK
)


# -------------------------------------------------------------------------
# Ling: only unambiguous symbol -> stable ENSG
# -------------------------------------------------------------------------

ling[
    "is_metagene_bool"
] = true_series(
    ling[
        "is_metagene"
    ]
)


ling_clean = ling.loc[
    (~ling[
        "is_metagene_bool"
    ])
    &
    ling[
        "ensembl_id"
    ].notna()
    &
    (
        pd.to_numeric(
            ling[
                "n_ensembl_for_symbol"
            ],
            errors="coerce"
        )
        ==
        1
    )
].copy()


# A symbol should now map to exactly one ENSG.
symbol_n = (
    ling_clean
    .groupby(
        "gene_symbol"
    )[
        "ensembl_id"
    ]
    .nunique()
)


bad_symbols = set(
    symbol_n[
        symbol_n > 1
    ].index
)


if bad_symbols:

    ling_clean = ling_clean.loc[
        ~ling_clean[
            "gene_symbol"
        ].isin(
            bad_symbols
        )
    ].copy()


ling_clean = (
    ling_clean
    .drop_duplicates(
        [
            "gene_symbol",
            "ensembl_id"
        ]
    )
)


symbol_to_ensg = dict(
    zip(
        ling_clean[
            "gene_symbol"
        ].astype(str),

        ling_clean[
            "ensembl_id"
        ].astype(str)
    )
)


# -------------------------------------------------------------------------
# Kwon stable 300-gene universe
# -------------------------------------------------------------------------

kwon_clean = kwon.loc[
    kwon[
        "ensembl_id"
    ].notna()
].copy()


if (
    "preferred_for_ensembl"
    in kwon_clean.columns
):

    pref = true_series(
        kwon_clean[
            "preferred_for_ensembl"
        ]
    )

    if pref.any():

        kwon_clean = (
            kwon_clean.loc[
                pref
            ]
        )


kwon_clean = (
    kwon_clean
    .drop_duplicates(
        "ensembl_id"
    )
)


KWON_ENSG = set(
    kwon_clean[
        "ensembl_id"
    ].astype(str)
)


print(
    "Unambiguous Ling symbols:",
    f"{len(symbol_to_ensg):,}"
)

print(
    "Kwon stable genes:",
    len(
        KWON_ENSG
    )
)


if len(
    KWON_ENSG
) != 300:

    raise RuntimeError(
        "Expected exactly 300 Kwon ENSG IDs."
    )


# =============================================================================
# FIND 87 LING FILES
# =============================================================================

section(
    "LING LIBRARIES"
)

files = sorted(
    LING_DIR.rglob(
        "*.gmg.donors.h5"
    )
)


print(
    "Files:",
    len(files)
)


if len(files) != 87:

    raise RuntimeError(
        f"Expected 87 Ling donor H5 files; "
        f"found {len(files)}"
    )


# =============================================================================
# AUDIT EACH LIBRARY AFTER ENSG HARMONIZATION
# =============================================================================

section(
    "POST-HARMONIZATION LIBRARY AUDIT"
)

presence = Counter()

rows = []

missing_kwon_rows = []


for i, path in enumerate(
    files,
    start=1
):

    with h5py.File(
        path,
        "r"
    ) as h:

        key = (
            "matrix/features/name"
        )

        if key not in h:
            raise KeyError(
                f"{path}: missing {key}"
            )

        symbols = decode(
            h[
                key
            ][:]
        )


    raw_unique = set(
        symbols
    )


    mapped_ids = [
        symbol_to_ensg[
            symbol
        ]

        for symbol in raw_unique

        if symbol
        in symbol_to_ensg
    ]


    mapped_set = set(
        mapped_ids
    )


    for ens in mapped_set:
        presence[
            ens
        ] += 1


    kwon_present = (
        mapped_set
        &
        KWON_ENSG
    )


    kwon_missing = (
        KWON_ENSG
        -
        mapped_set
    )


    duplicate_collapse = (
        len(
            mapped_ids
        )
        -
        len(
            mapped_set
        )
    )


    rows.append({
        "file":
            path.name,

        "raw_genes":
            len(
                symbols
            ),

        "raw_unique_symbols":
            len(
                raw_unique
            ),

        "mapped_symbols":
            len(
                mapped_ids
            ),

        "harmonized_unique_ENSG":
            len(
                mapped_set
            ),

        "symbols_collapsed_by_ENSG":
            duplicate_collapse,

        "unmapped_raw_symbols":
            len(
                raw_unique
            )
            -
            len(
                mapped_ids
            ),

        "kwon_genes_present":
            len(
                kwon_present
            ),

        "kwon_genes_missing":
            len(
                kwon_missing
            ),

        "kwon_coverage_percent":
            (
                100
                *
                len(
                    kwon_present
                )
                /
                300
            ),
    })


    for ens in sorted(
        kwon_missing
    ):

        missing_kwon_rows.append({
            "file":
                path.name,

            "ensembl_id":
                ens
        })


    print(
        f"{i:02d}/87 "
        f"raw={len(raw_unique):,} "
        f"ENSG={len(mapped_set):,} "
        f"Kwon={len(kwon_present)}/300 "
        f"{path.name}"
    )


qc = pd.DataFrame(
    rows
)


qc.to_csv(
    OUT
    / "ling_postharmonization_library_qc.csv",
    index=False
)


pd.DataFrame(
    missing_kwon_rows
).to_csv(
    OUT
    / "ling_kwon_missing_by_library.csv",
    index=False
)


# =============================================================================
# GLOBAL HARMONIZED UNIVERSE
# =============================================================================

section(
    "GLOBAL POST-HARMONIZATION GENE UNIVERSE"
)

n_files = len(
    files
)


prevalence = pd.DataFrame({
    "ensembl_id":
        sorted(
            presence
        )
})


prevalence[
    "n_libraries"
] = prevalence[
    "ensembl_id"
].map(
    presence
)


prevalence[
    "fraction_libraries"
] = (
    prevalence[
        "n_libraries"
    ]
    /
    n_files
)


prevalence[
    "is_kwon_gene"
] = prevalence[
    "ensembl_id"
].isin(
    KWON_ENSG
)


prevalence.to_csv(
    OUT
    / "ling_postharmonization_gene_prevalence.csv",
    index=False
)


def n_at_least(frac):

    return int(
        (
            prevalence[
                "fraction_libraries"
            ]
            >= frac
        ).sum()
    )


union_n = len(
    prevalence
)

intersection_n = int(
    (
        prevalence[
            "n_libraries"
        ]
        ==
        n_files
    ).sum()
)


kwon_prevalence = (
    prevalence.loc[
        prevalence[
            "is_kwon_gene"
        ]
    ]
)


kwon_all87 = int(
    (
        kwon_prevalence[
            "n_libraries"
        ]
        ==
        87
    ).sum()
)


kwon_95pct = int(
    (
        kwon_prevalence[
            "fraction_libraries"
        ]
        >=
        0.95
    ).sum()
)


summary = pd.DataFrame([
    {
        "ling_files":
            n_files,

        "harmonized_union_ENSG":
            union_n,

        "harmonized_intersection_ENSG":
            intersection_n,

        "genes_in_100pct_libraries":
            n_at_least(
                1.00
            ),

        "genes_in_99pct_libraries":
            n_at_least(
                0.99
            ),

        "genes_in_95pct_libraries":
            n_at_least(
                0.95
            ),

        "genes_in_90pct_libraries":
            n_at_least(
                0.90
            ),

        "genes_in_50pct_libraries":
            n_at_least(
                0.50
            ),

        "kwon_total":
            300,

        "kwon_present_in_all_87":
            kwon_all87,

        "kwon_present_in_at_least_95pct":
            kwon_95pct,

        "min_kwon_genes_per_library":
            int(
                qc[
                    "kwon_genes_present"
                ].min()
            ),

        "median_kwon_genes_per_library":
            float(
                qc[
                    "kwon_genes_present"
                ].median()
            ),

        "max_kwon_genes_per_library":
            int(
                qc[
                    "kwon_genes_present"
                ].max()
            ),
    }
])


summary.to_csv(
    OUT
    / "ling_reference_build_audit_summary.csv",
    index=False
)


print()
print(
    summary.to_string(
        index=False
    )
)


# =============================================================================
# KWON GENE PREVALENCE
# =============================================================================

kwon_symbol_map = (
    kwon_clean
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


kwon_report = (
    pd.DataFrame({
        "ensembl_id":
            sorted(
                KWON_ENSG
            )
    })
)


kwon_report[
    "Kwon_symbol"
] = kwon_report[
    "ensembl_id"
].map(
    kwon_symbol_map
)


kwon_report[
    "n_Ling_libraries"
] = kwon_report[
    "ensembl_id"
].map(
    presence
).fillna(
    0
).astype(
    int
)


kwon_report[
    "fraction_Ling_libraries"
] = (
    kwon_report[
        "n_Ling_libraries"
    ]
    /
    87
)


kwon_report = (
    kwon_report.sort_values(
        [
            "n_Ling_libraries",
            "Kwon_symbol"
        ]
    )
)


kwon_report.to_csv(
    OUT
    / "kwon300_presence_across_ling_libraries.csv",
    index=False
)


print()
print(
    "Lowest-prevalence Kwon genes:"
)

print(
    kwon_report.head(
        20
    ).to_string(
        index=False
    )
)


section(
    "FINAL STATUS"
)

print(
    "Audit complete."
)

print(
    "No expression matrices were modified."
)

print(
    "Outputs:",
    OUT
)
