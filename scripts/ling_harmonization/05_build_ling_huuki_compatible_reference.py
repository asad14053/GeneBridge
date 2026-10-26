#!/usr/bin/env python3

import argparse
import gc
import json
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


# =============================================================================
# HELPERS
# =============================================================================

def strip_ensg(x):
    if pd.isna(x):
        return None
    x = str(x).strip()
    if not x:
        return None
    return x.split(".")[0]


def decode_array(x):
    return np.asarray([
        v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v)
        for v in x
    ])


def normalized_name(x):
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def detect_ensg_column(df):
    priority = [
        "ensembl_id",
        "ensembl",
        "ensg",
        "stable_ensembl_id",
        "stable_id",
    ]

    norms = {
        normalized_name(c): c
        for c in df.columns
    }

    for p in priority:
        n = normalized_name(p)
        if n in norms:
            return norms[n]

    for c in df.columns:
        n = normalized_name(c)
        if "ensembl" in n or "ensg" in n:
            return c

    # Last-resort content inspection
    for c in df.columns:
        vals = (
            df[c]
            .dropna()
            .astype(str)
            .head(100)
        )
        if len(vals) and vals.str.startswith("ENSG").mean() > 0.8:
            return c

    raise RuntimeError(
        f"Could not identify Ensembl-ID column. Columns={df.columns.tolist()}"
    )


def detect_symbol_column(df, source):
    source_l = source.lower()

    preferred = [
        f"{source}_symbol",
        f"{source_l}_symbol",
        "gene_symbol",
        "symbol",
        "published_symbol",
        "source_symbol",
        "gene_name",
    ]

    norms = {
        normalized_name(c): c
        for c in df.columns
    }

    for p in preferred:
        n = normalized_name(p)
        if n in norms:
            return norms[n]

    # Prefer a source-specific symbol column.
    for c in df.columns:
        n = normalized_name(c)
        if source_l in n and "symbol" in n:
            return c

    # Then any symbol column.
    for c in df.columns:
        if "symbol" in normalized_name(c):
            return c

    raise RuntimeError(
        f"Could not identify symbol column in {source}. "
        f"Columns={df.columns.tolist()}"
    )


def load_crosswalk(path, source):
    print(f"\nLoading {source} crosswalk:")
    print(path)

    df = pd.read_csv(path)

    ens_col = detect_ensg_column(df)
    sym_col = detect_symbol_column(df, source)

    x = df[[sym_col, ens_col]].copy()
    x.columns = ["symbol", "ensembl_id"]

    x["symbol"] = x["symbol"].astype(str).str.strip()
    x["ensembl_id"] = x["ensembl_id"].map(strip_ensg)

    x = x[
        x["ensembl_id"].notna()
        & x["symbol"].notna()
        & (x["symbol"] != "")
        & (x["symbol"].str.lower() != "nan")
    ].copy()

    # Only retain symbols that map unambiguously to one ENSG.
    n_ens = (
        x.groupby("symbol")["ensembl_id"]
        .nunique()
    )

    good_symbols = set(
        n_ens[n_ens == 1].index
    )

    x = x[
        x["symbol"].isin(good_symbols)
    ].drop_duplicates()

    symbol_to_ens = (
        x.drop_duplicates("symbol")
        .set_index("symbol")["ensembl_id"]
        .to_dict()
    )

    # For ENSG -> symbol choose deterministically.
    ens_to_symbol = {}

    for ens, g in x.groupby("ensembl_id", sort=True):
        symbols = sorted(
            set(g["symbol"].astype(str))
        )
        if symbols:
            ens_to_symbol[ens] = symbols[0]

    print(
        f"  unambiguous symbols: {len(symbol_to_ens):,}"
    )
    print(
        f"  unique ENSG:         {len(ens_to_symbol):,}"
    )

    return x, symbol_to_ens, ens_to_symbol


def find_file(candidates, search_root=None, contains=None):
    for p in candidates:
        p = Path(p)
        if p.exists():
            return p

    if search_root is not None and contains is not None:
        hits = [
            p for p in Path(search_root).rglob("*.csv")
            if contains.lower() in p.name.lower()
        ]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            print("\nCandidate files:")
            for h in hits:
                print(" ", h)
            raise RuntimeError(
                f"Multiple candidate files found for {contains}"
            )

    raise FileNotFoundError(
        "Could not locate expected input file."
    )


def load_gene_prevalence(path, expected_libraries=87):
    print("\nLoading Ling post-harmonization prevalence:")
    print(path)

    df = pd.read_csv(path)

    ens_col = detect_ensg_column(df)

    fraction_col = None
    count_col = None

    for c in df.columns:
        n = normalized_name(c)

        if (
            ("fraction" in n or "prevalence" in n)
            and ("libr" in n or "file" in n)
        ):
            fraction_col = c
            break

    if fraction_col is None:
        for c in df.columns:
            n = normalized_name(c)
            if (
                ("nlibr" in n)
                or ("libraries" in n)
                or ("nfiles" in n)
            ):
                count_col = c
                break

    out = df.copy()
    out["ensembl_id"] = out[ens_col].map(strip_ensg)

    if fraction_col is not None:
        fraction = pd.to_numeric(
            out[fraction_col],
            errors="coerce"
        )

        if fraction.max() > 1.5:
            fraction = fraction / 100.0

        out["library_fraction"] = fraction

    elif count_col is not None:
        out["library_fraction"] = (
            pd.to_numeric(
                out[count_col],
                errors="coerce"
            )
            / float(expected_libraries)
        )

    else:
        raise RuntimeError(
            "Could not identify library prevalence column.\n"
            f"Columns={df.columns.tolist()}"
        )

    out = (
        out.dropna(
            subset=[
                "ensembl_id",
                "library_fraction",
            ]
        )
        .drop_duplicates("ensembl_id")
    )

    return out


def filename_to_prefix(path):
    name = path.name

    prefix = name.split(
        ".intronic_",
        1
    )[0]

    # Metadata omits this suffix.
    prefix = re.sub(
        r"_MyelinRemoved$",
        "",
        prefix
    )

    return prefix


def proportional_stratified_sample(
    meta,
    n_target,
    seed,
):
    """
    Deterministic proportional sampling across DONOR x cellclass.
    """

    eligible = meta[
        meta["DONOR"].notna()
        & meta["cellclass"].notna()
    ].copy()

    # All current jointMetadata doublet values are False.
    # Still enforce if column exists.
    if "doublet" in eligible.columns:
        dbl = (
            eligible["doublet"]
            .astype(str)
            .str.lower()
            .isin(
                ["true", "1", "yes"]
            )
        )
        eligible = eligible[
            ~dbl
        ].copy()

    if len(eligible) < n_target:
        raise RuntimeError(
            f"Only {len(eligible):,} eligible cells, "
            f"but requested {n_target:,}"
        )

    eligible["_stratum"] = (
        eligible["DONOR"].astype(str)
        + "|||"
        + eligible["cellclass"].astype(str)
    )

    counts = (
        eligible["_stratum"]
        .value_counts()
        .sort_index()
    )

    desired = (
        counts.astype(float)
        * float(n_target)
        / float(counts.sum())
    )

    quotas = np.floor(
        desired
    ).astype(int)

    remaining = int(
        n_target - quotas.sum()
    )

    fractional = (
        desired - quotas
    ).sort_values(
        ascending=False,
        kind="mergesort",
    )

    # Largest remainder allocation.
    for stratum in fractional.index:
        if remaining <= 0:
            break

        if quotas.loc[stratum] < counts.loc[stratum]:
            quotas.loc[stratum] += 1
            remaining -= 1

    if remaining != 0:
        raise RuntimeError(
            f"Could not allocate all sample slots; remaining={remaining}"
        )

    rng = np.random.default_rng(
        seed
    )

    chosen_indices = []

    for stratum, group in eligible.groupby(
        "_stratum",
        sort=True,
    ):
        q = int(
            quotas.loc[stratum]
        )

        if q == 0:
            continue

        idx = group.index.to_numpy()

        chosen = rng.choice(
            idx,
            size=q,
            replace=False,
        )

        chosen_indices.extend(
            chosen.tolist()
        )

    selected = meta.loc[
        chosen_indices
    ].copy()

    if len(selected) != n_target:
        raise RuntimeError(
            f"Sampling produced {len(selected):,}, "
            f"expected {n_target:,}"
        )

    selected["_metadata_row"] = (
        selected.index.astype(int)
    )

    # Sorting by library makes H5 reconstruction much more efficient.
    selected = (
        selected.sort_values(
            [
                "PREFIX",
                "_metadata_row",
            ],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )

    selected["cell_id"] = (
        selected["PREFIX"].astype(str)
        + "_"
        + selected["CELL_BARCODE"].astype(str)
    )

    if selected["cell_id"].duplicated().any():
        dup = selected.loc[
            selected["cell_id"].duplicated(),
            "cell_id"
        ].iloc[0]
        raise RuntimeError(
            f"Duplicate sampled cell ID: {dup}"
        )

    selected["reference_row"] = np.arange(
        len(selected),
        dtype=np.int64,
    )

    return selected, eligible, quotas


def load_optional_donor_metadata(path):
    if path is None:
        return None

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(path)

    print("\nLoading donor metadata:")
    print(path)

    # Small metadata table; infer delimiter.
    donor = pd.read_csv(
        path,
        sep=None,
        engine="python",
    )

    donor_col = None

    for c in donor.columns:
        if normalized_name(c) in {
            "donor",
            "donorid",
            "subject",
            "subjectid",
        }:
            donor_col = c
            break

    if donor_col is None:
        raise RuntimeError(
            "Could not identify donor column in donor metadata."
        )

    donor = donor.copy()
    donor["_DONOR_KEY"] = (
        donor[donor_col]
        .astype(str)
        .str.strip()
    )

    donor = donor.drop_duplicates(
        "_DONOR_KEY"
    )

    return donor


def find_column(df, names):
    if df is None:
        return None

    norms = {
        normalized_name(c): c
        for c in df.columns
    }

    for name in names:
        n = normalized_name(name)
        if n in norms:
            return norms[n]

    return None


def normalize_diagnosis(v):
    if pd.isna(v):
        return "NA"

    s = str(v).strip()
    sl = s.lower()

    if sl in {
        "1", "true", "yes",
        "scz", "schizophrenia",
        "case",
    }:
        return "SCZ"

    if sl in {
        "0", "false", "no",
        "control", "ctrl",
        "ntc", "unaffected",
    }:
        return "NTC"

    return s


# =============================================================================
# MAIN
# =============================================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        default="/beegfs/labs/hulab/projects/mjabin/GeneBridge",
    )

    parser.add_argument(
        "--n-cells",
        type=int,
        default=77604,
    )

    parser.add_argument(
        "--gene-prevalence",
        type=float,
        default=0.95,
    )

    parser.add_argument(
        "--expected-genes",
        type=int,
        default=31491,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20261023,
    )

    parser.add_argument(
        "--work-dir",
        required=True,
    )

    parser.add_argument(
        "--donor-meta",
        default=None,
        help=(
            "Optional Ling donor metadata containing Donor, "
            "Age, Sex and Schizophrenia columns."
        ),
    )

    parser.add_argument(
        "--keep-temp",
        action="store_true",
    )

    args = parser.parse_args()

    ROOT = Path(args.root)

    COUNTS_DIR = (
        ROOT
        / "data/raw/nemo_mccarroll/counts"
    )

    CELL_META = (
        COUNTS_DIR
        / "BA46.jointMetadata.txt"
    )

    HARMONIZATION_DIR = (
        ROOT
        / "results/analysis/"
          "ling_huuki_kwon_gene_harmonization"
    )

    BUILD_DIR = (
        ROOT
        / "results/analysis/"
          "ling_reference_build"
    )

    OUTPUT_DIR = (
        ROOT
        / "data/processed/snrnaseq/"
          "ling_compatible"
    )

    WORK_DIR = Path(
        args.work_dir
    )

    BUILD_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    WORK_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    OUTPUT = (
        OUTPUT_DIR
        / "ling_snrna_reference_huuki_matched_77604_allgenes.h5ad"
    )

    print("=" * 100)
    print("LING -> HUUKI-COMPATIBLE REFERENCE BUILD")
    print("=" * 100)

    print("Root:             ", ROOT)
    print("Cell metadata:    ", CELL_META)
    print("Count directory:  ", COUNTS_DIR)
    print("Work directory:   ", WORK_DIR)
    print("Output:           ", OUTPUT)
    print("Target cells:     ", f"{args.n_cells:,}")
    print("Gene prevalence:  ", args.gene_prevalence)
    print("Expected genes:   ", f"{args.expected_genes:,}")
    print("Seed:             ", args.seed)


    # =========================================================================
    # INPUT FILE DISCOVERY
    # =========================================================================

    ling_xwalk_path = find_file(
        [
            HARMONIZATION_DIR
            / "Ling_union_symbol_to_ensembl.csv",
            HARMONIZATION_DIR
            / "ling_union_symbol_to_ensembl.csv",
        ],
        search_root=HARMONIZATION_DIR,
        contains="ling_union_symbol",
    )

    huuki_xwalk_path = find_file(
        [
            HARMONIZATION_DIR
            / "Huuki_gene_to_ensembl.csv",
            HARMONIZATION_DIR
            / "huuki_gene_to_ensembl.csv",
        ],
        search_root=HARMONIZATION_DIR,
        contains="huuki_gene",
    )

    kwon_xwalk_path = find_file(
        [
            HARMONIZATION_DIR
            / "Kwon_gene_to_ensembl.csv",
            HARMONIZATION_DIR
            / "kwon_gene_to_ensembl.csv",
        ],
        search_root=HARMONIZATION_DIR,
        contains="kwon_gene",
    )

    prevalence_path = find_file(
        [
            BUILD_DIR
            / "ling_postharmonization_gene_prevalence.csv",
            HARMONIZATION_DIR
            / "ling_postharmonization_gene_prevalence.csv",
        ],
        search_root=ROOT / "results/analysis",
        contains="postharmonization_gene_prevalence",
    )


    # =========================================================================
    # LOAD CROSSWALKS
    # =========================================================================

    (
        ling_xwalk,
        ling_symbol_to_ens,
        ling_ens_to_symbol,
    ) = load_crosswalk(
        ling_xwalk_path,
        "Ling",
    )

    (
        huuki_xwalk,
        _,
        huuki_ens_to_symbol,
    ) = load_crosswalk(
        huuki_xwalk_path,
        "Huuki",
    )

    (
        kwon_xwalk,
        _,
        kwon_ens_to_symbol,
    ) = load_crosswalk(
        kwon_xwalk_path,
        "Kwon",
    )


    # =========================================================================
    # SELECT >=95% LING GENE UNIVERSE
    # =========================================================================

    prevalence = load_gene_prevalence(
        prevalence_path,
        expected_libraries=87,
    )

    selected_gene_df = prevalence[
        prevalence["library_fraction"]
        >= args.gene_prevalence
    ].copy()

    selected_gene_df = (
        selected_gene_df
        .sort_values("ensembl_id")
        .drop_duplicates("ensembl_id")
        .reset_index(drop=True)
    )

    selected_ensg = (
        selected_gene_df["ensembl_id"]
        .tolist()
    )

    print("\n" + "=" * 100)
    print("GENE UNIVERSE")
    print("=" * 100)

    print(
        "Genes at threshold:",
        f"{len(selected_ensg):,}"
    )

    if (
        args.expected_genes is not None
        and len(selected_ensg)
        != args.expected_genes
    ):
        raise RuntimeError(
            f"Expected {args.expected_genes:,} genes "
            f"at >= {args.gene_prevalence:.2f} prevalence, "
            f"found {len(selected_ensg):,}"
        )

    kwon_ensg = set(
        kwon_ens_to_symbol.keys()
    )

    missing_kwon_ensg = sorted(
        kwon_ensg
        - set(selected_ensg)
    )

    print(
        "Kwon ENSG total:",
        len(kwon_ensg)
    )

    print(
        "Kwon ENSG retained:",
        len(
            kwon_ensg
            & set(selected_ensg)
        )
    )

    if missing_kwon_ensg:
        raise RuntimeError(
            "Some Kwon genes are absent from the "
            ">=95% Ling gene universe:\n"
            + "\n".join(
                missing_kwon_ensg
            )
        )


    # =========================================================================
    # CANONICAL SYMBOLS
    #
    # Priority:
    #   1. Kwon symbol -- required for direct Xenium target matching
    #   2. Huuki symbol
    #   3. Ling symbol
    #   4. ENSG fallback only if needed to maintain uniqueness
    # =========================================================================

    canonical = {}
    canonical_source = {}
    used_symbols = set()

    # Assign Kwon symbols first so none can be displaced.
    for ens in selected_ensg:
        if ens not in kwon_ens_to_symbol:
            continue

        sym = kwon_ens_to_symbol[ens]

        if sym in used_symbols:
            raise RuntimeError(
                f"Duplicate Kwon canonical symbol: {sym}"
            )

        canonical[ens] = sym
        canonical_source[ens] = "Kwon"
        used_symbols.add(sym)


    for ens in selected_ensg:

        if ens in canonical:
            continue

        candidates = []

        if ens in huuki_ens_to_symbol:
            candidates.append(
                (
                    huuki_ens_to_symbol[ens],
                    "Huuki",
                )
            )

        if ens in ling_ens_to_symbol:
            candidates.append(
                (
                    ling_ens_to_symbol[ens],
                    "Ling",
                )
            )

        candidates.append(
            (
                ens,
                "ENSG_fallback",
            )
        )

        chosen = None

        for sym, source in candidates:

            if (
                sym
                and sym not in used_symbols
            ):
                chosen = (
                    sym,
                    source,
                )
                break

        if chosen is None:
            # Should never happen because ENSG itself is unique.
            raise RuntimeError(
                f"Unable to assign unique symbol for {ens}"
            )

        canonical[ens] = chosen[0]
        canonical_source[ens] = chosen[1]
        used_symbols.add(chosen[0])


    var_names = [
        canonical[e]
        for e in selected_ensg
    ]

    if len(set(var_names)) != len(var_names):
        raise RuntimeError(
            "Final canonical var_names are not unique."
        )

    # Explicit Kwon symbol test.
    kwon_symbols = {
        kwon_ens_to_symbol[e]
        for e in kwon_ensg
    }

    missing_kwon_symbols = sorted(
        kwon_symbols
        - set(var_names)
    )

    if missing_kwon_symbols:
        raise RuntimeError(
            "Canonical symbol assignment lost Kwon genes:\n"
            + "\n".join(
                missing_kwon_symbols
            )
        )

    print(
        "Kwon canonical symbols retained:",
        len(kwon_symbols)
    )

    gene_manifest = pd.DataFrame({
        "reference_column":
            np.arange(
                len(selected_ensg)
            ),

        "ensembl_id":
            selected_ensg,

        "canonical_symbol":
            var_names,

        "canonical_source":
            [
                canonical_source[e]
                for e in selected_ensg
            ],
    })

    prevalence_map = (
        selected_gene_df
        .set_index("ensembl_id")[
            "library_fraction"
        ]
        .to_dict()
    )

    gene_manifest[
        "library_fraction"
    ] = [
        prevalence_map[e]
        for e in selected_ensg
    ]

    gene_manifest.to_csv(
        BUILD_DIR
        / "ling_final_reference_gene_manifest.csv",
        index=False,
    )


    # =========================================================================
    # LOAD CELL METADATA
    # =========================================================================

    print("\n" + "=" * 100)
    print("CELL UNIVERSE")
    print("=" * 100)

    meta = pd.read_csv(
        CELL_META,
        sep="\t",
        low_memory=False,
    )

    print(
        "jointMetadata cells:",
        f"{len(meta):,}"
    )

    print(
        "Donors:",
        meta["DONOR"].nunique()
    )

    print(
        "Prefixes:",
        meta["PREFIX"].nunique()
    )

    selected, eligible, quotas = (
        proportional_stratified_sample(
            meta,
            args.n_cells,
            args.seed,
        )
    )

    print(
        "Annotated eligible cells:",
        f"{len(eligible):,}"
    )

    print(
        "Selected cells:",
        f"{len(selected):,}"
    )

    print(
        "Selected donors:",
        selected["DONOR"].nunique()
    )

    print(
        "Selected prefixes:",
        selected["PREFIX"].nunique()
    )

    print("\nSelected broad classes:")
    print(
        selected[
            "cellclass"
        ].value_counts()
    )


    # Save COMPLETE Ling metadata for sampled cells as sidecar.
    selected.to_csv(
        BUILD_DIR
        / "ling_selected_77604_cells_full_metadata.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )


    sampling_summary = (
        selected.groupby(
            [
                "DONOR",
                "cellclass",
            ],
            dropna=False,
        )
        .size()
        .reset_index(
            name="n_selected"
        )
    )

    sampling_summary.to_csv(
        BUILD_DIR
        / "ling_selected_77604_donor_cellclass_counts.csv",
        index=False,
    )


    # =========================================================================
    # FIND 87 H5 COUNT FILES
    # =========================================================================

    h5_files = sorted(
        COUNTS_DIR.glob(
            "*.gmg.donors.h5"
        )
    )

    print("\n" + "=" * 100)
    print("COUNT FILES")
    print("=" * 100)

    print(
        "gmg.donors.h5 files:",
        len(h5_files)
    )

    if len(h5_files) != 87:
        raise RuntimeError(
            f"Expected exactly 87 count H5 files; "
            f"found {len(h5_files)}"
        )

    prefix_to_file = {}

    for path in h5_files:
        prefix = filename_to_prefix(
            path
        )

        if prefix in prefix_to_file:
            raise RuntimeError(
                f"Duplicate normalized library prefix: {prefix}"
            )

        prefix_to_file[prefix] = path


    missing_prefixes = sorted(
        set(selected["PREFIX"])
        - set(prefix_to_file)
    )

    if missing_prefixes:
        raise RuntimeError(
            "Selected metadata prefixes missing H5 files:\n"
            + "\n".join(
                missing_prefixes
            )
        )


    # =========================================================================
    # ALLOCATE DENSE FLOAT64 MEMMAPS
    # =========================================================================

    n_cells = len(selected)
    n_genes = len(selected_ensg)

    x_path = (
        WORK_DIR
        / "ling_X_float64.dat"
    )

    log_path = (
        WORK_DIR
        / "ling_logcounts_float64.dat"
    )

    print("\n" + "=" * 100)
    print("ALLOCATING DENSE MATRICES")
    print("=" * 100)

    gb_per_matrix = (
        n_cells
        * n_genes
        * 8
        / 1024**3
    )

    print(
        f"Matrix shape: {n_cells:,} x {n_genes:,}"
    )

    print(
        f"Approx size per float64 matrix: {gb_per_matrix:.2f} GiB"
    )

    Xmm = np.memmap(
        x_path,
        mode="w+",
        dtype=np.float64,
        shape=(
            n_cells,
            n_genes,
        ),
    )

    # Initialize explicitly.
    Xmm[:] = 0.0
    Xmm.flush()


    # Metrics calculated directly while each library block is built.
    library_size = np.zeros(
        n_cells,
        dtype=np.float64,
    )

    detected = np.zeros(
        n_cells,
        dtype=np.int32,
    )

    mito_sum = np.zeros(
        n_cells,
        dtype=np.float64,
    )

    mito_detected = np.zeros(
        n_cells,
        dtype=np.int32,
    )

    target_ens_to_col = {
        ens: j
        for j, ens
        in enumerate(
            selected_ensg
        )
    }

    mito_cols = np.array(
        [
            j
            for j, sym
            in enumerate(
                var_names
            )
            if str(sym).upper().startswith(
                "MT-"
            )
        ],
        dtype=int,
    )

    print(
        "Mitochondrial target genes:",
        len(mito_cols)
    )


    # =========================================================================
    # RECONSTRUCT COUNTS
    # =========================================================================

    print("\n" + "=" * 100)
    print("RECONSTRUCTING SELECTED LING COUNTS")
    print("=" * 100)

    processed_cells = 0

    for lib_i, (
        prefix,
        group,
    ) in enumerate(
        selected.groupby(
            "PREFIX",
            sort=True,
        ),
        start=1,
    ):

        path = prefix_to_file[
            prefix
        ]

        group = group.sort_values(
            "reference_row"
        )

        row_ids = (
            group["reference_row"]
            .to_numpy(
                dtype=np.int64
            )
        )

        # Because final selected metadata was sorted by PREFIX,
        # each library should occupy one contiguous output block.
        expected_rows = np.arange(
            row_ids.min(),
            row_ids.max() + 1,
        )

        if not np.array_equal(
            row_ids,
            expected_rows,
        ):
            raise RuntimeError(
                f"Output rows are not contiguous for {prefix}"
            )

        with h5py.File(
            path,
            "r",
        ) as h:

            barcodes = decode_array(
                h[
                    "matrix/barcodes"
                ][:]
            )

            genes = decode_array(
                h[
                    "matrix/features/name"
                ][:]
            )

            shape = tuple(
                int(v)
                for v in h[
                    "matrix/shape"
                ][:]
            )

            data = h[
                "matrix/data"
            ][:]

            indices = h[
                "matrix/indices"
            ][:]

            indptr = h[
                "matrix/indptr"
            ][:]

        if (
            shape[0] != len(genes)
            or shape[1] != len(barcodes)
        ):
            raise RuntimeError(
                f"Unexpected matrix orientation for {path.name}: "
                f"shape={shape}, genes={len(genes)}, "
                f"barcodes={len(barcodes)}"
            )

        if len(indptr) != len(barcodes) + 1:
            raise RuntimeError(
                f"Expected CSC cell-column indptr for {path.name}"
            )

        barcode_to_col = {
            b: i
            for i, b
            in enumerate(
                barcodes
            )
        }

        wanted_barcodes = (
            group[
                "CELL_BARCODE"
            ]
            .astype(str)
            .tolist()
        )

        missing = [
            b
            for b in wanted_barcodes
            if b not in barcode_to_col
        ]

        if missing:
            raise RuntimeError(
                f"{prefix}: {len(missing)} selected barcodes "
                f"missing from H5. First={missing[:5]}"
            )

        h5_cols = np.asarray(
            [
                barcode_to_col[b]
                for b in wanted_barcodes
            ],
            dtype=np.int64,
        )

        X_gene_cell = sparse.csc_matrix(
            (
                data,
                indices,
                indptr,
            ),
            shape=shape,
        )

        # Subset CELLS before gene remapping.
        X_sub_gene_cell = (
            X_gene_cell[
                :,
                h5_cols,
            ]
        )


        # -------------------------------------------------------------
        # Source feature -> target ENSG mapping matrix
        # -------------------------------------------------------------

        map_source_rows = []
        map_target_cols = []

        for src_j, symbol in enumerate(
            genes
        ):

            ens = ling_symbol_to_ens.get(
                symbol
            )

            if ens is None:
                continue

            target_j = (
                target_ens_to_col.get(
                    ens
                )
            )

            if target_j is None:
                continue

            map_source_rows.append(
                src_j
            )

            map_target_cols.append(
                target_j
            )


        mapping = sparse.csr_matrix(
            (
                np.ones(
                    len(map_source_rows),
                    dtype=np.float64,
                ),
                (
                    np.asarray(
                        map_source_rows,
                        dtype=np.int64,
                    ),
                    np.asarray(
                        map_target_cols,
                        dtype=np.int64,
                    ),
                ),
            ),
            shape=(
                len(genes),
                n_genes,
            ),
        )


        # Cells x source-genes  @ source-genes x target-genes
        # Multiple source symbols mapping to the same stable ENSG
        # are summed automatically.
        X_block_sparse = (
            X_sub_gene_cell.T
            .tocsr()
            @ mapping
        ).tocsr()


        X_block = (
            X_block_sparse
            .toarray()
            .astype(
                np.float64,
                copy=False,
            )
        )

        start = int(
            row_ids[0]
        )

        end = int(
            row_ids[-1] + 1
        )

        if X_block.shape != (
            end - start,
            n_genes,
        ):
            raise RuntimeError(
                f"Unexpected block shape for {prefix}: "
                f"{X_block.shape}"
            )

        Xmm[
            start:end,
            :
        ] = X_block


        # -------------------------------------------------------------
        # Per-cell QC based on FINAL reference gene universe
        # -------------------------------------------------------------

        library_size[
            start:end
        ] = X_block.sum(
            axis=1
        )

        detected[
            start:end
        ] = np.count_nonzero(
            X_block,
            axis=1,
        )

        if len(mito_cols):
            mito_block = X_block[
                :,
                mito_cols,
            ]

            mito_sum[
                start:end
            ] = mito_block.sum(
                axis=1
            )

            mito_detected[
                start:end
            ] = np.count_nonzero(
                mito_block,
                axis=1,
            )


        processed_cells += len(group)

        print(
            f"[{lib_i:02d}/87] "
            f"{prefix:<30} "
            f"selected={len(group):>5,} "
            f"source_genes={len(genes):>6,} "
            f"mapped_features={len(map_source_rows):>6,} "
            f"processed={processed_cells:>7,}/{n_cells:,}"
        )


        del (
            X_gene_cell,
            X_sub_gene_cell,
            mapping,
            X_block_sparse,
            X_block,
            data,
            indices,
            indptr,
        )

        gc.collect()


    Xmm.flush()


    if processed_cells != n_cells:
        raise RuntimeError(
            f"Processed {processed_cells:,} cells, "
            f"expected {n_cells:,}"
        )


    # =========================================================================
    # HUUKI-IDENTICAL SIZE FACTOR
    # =========================================================================

    print("\n" + "=" * 100)
    print("HUUKI NORMALIZATION")
    print("=" * 100)

    if np.any(
        library_size <= 0
    ):
        zero_cells = np.where(
            library_size <= 0
        )[0]

        raise RuntimeError(
            f"{len(zero_cells)} selected cells have zero "
            f"counts in final gene universe."
        )

    mean_library_size = float(
        library_size.mean()
    )

    size_factor = (
        library_size
        / mean_library_size
    )

    print(
        "Mean library size:",
        mean_library_size
    )

    print(
        "Mean sizeFactor:",
        float(
            size_factor.mean()
        )
    )

    print(
        "sizeFactor min/max:",
        float(
            size_factor.min()
        ),
        float(
            size_factor.max()
        ),
    )


    # =========================================================================
    # LOGCOUNTS
    #
    # Verified against Huuki:
    #
    # logcounts = log2(count / sizeFactor + 1)
    # =========================================================================

    Lmm = np.memmap(
        log_path,
        mode="w+",
        dtype=np.float64,
        shape=(
            n_cells,
            n_genes,
        ),
    )

    normalization_block = 512

    for start in range(
        0,
        n_cells,
        normalization_block,
    ):

        end = min(
            n_cells,
            start
            + normalization_block,
        )

        block = np.asarray(
            Xmm[
                start:end,
                :
            ],
            dtype=np.float64,
        )

        Lmm[
            start:end,
            :
        ] = np.log2(
            block
            / size_factor[
                start:end,
                None
            ]
            + 1.0
        )

        if (
            start == 0
            or end == n_cells
            or start % 5120 == 0
        ):
            print(
                f"logcounts: {end:,}/{n_cells:,}"
            )

    Lmm.flush()


    # =========================================================================
    # OPTIONAL DONOR METADATA
    # =========================================================================

    donor_meta = load_optional_donor_metadata(
        args.donor_meta
    )

    age_map = {}
    sex_map = {}
    diagnosis_map = {}

    if donor_meta is not None:

        age_col = find_column(
            donor_meta,
            [
                "Age",
                "age",
            ],
        )

        sex_col = find_column(
            donor_meta,
            [
                "Sex",
                "sex",
            ],
        )

        dx_col = find_column(
            donor_meta,
            [
                "Schizophrenia",
                "diagnosis",
                "Dx",
            ],
        )

        donor_index = donor_meta.set_index(
            "_DONOR_KEY"
        )

        if age_col:
            age_map = (
                donor_index[
                    age_col
                ].to_dict()
            )

        if sex_col:
            sex_map = (
                donor_index[
                    sex_col
                ].to_dict()
            )

        if dx_col:
            diagnosis_map = {
                str(k):
                    normalize_diagnosis(v)
                for k, v
                in donor_index[
                    dx_col
                ].to_dict().items()
            }


    donors = (
        selected[
            "DONOR"
        ]
        .astype(str)
        .tolist()
    )

    ages = [
        age_map.get(
            d,
            np.nan,
        )
        for d in donors
    ]

    sexes = [
        sex_map.get(
            d,
            "NA",
        )
        for d in donors
    ]

    diagnoses = [
        diagnosis_map.get(
            d,
            "NA",
        )
        for d in donors
    ]


    # =========================================================================
    # EXACT HUUKI OBS COLUMN CONTRACT
    # =========================================================================

    print("\n" + "=" * 100)
    print("BUILDING HUUKI-COMPATIBLE OBS")
    print("=" * 100)

    mito_percent = (
        100.0
        * mito_sum
        / library_size
    )

    broad = (
        selected[
            "cellclass"
        ]
        .fillna("NA")
        .astype(str)
        .to_numpy()
    )

    subclass = (
        selected[
            "subclass"
        ]
        .fillna("NA")
        .astype(str)
        .to_numpy()
    )

    prefixes = (
        selected[
            "PREFIX"
        ]
        .astype(str)
        .to_numpy()
    )

    barcodes = (
        selected[
            "CELL_BARCODE"
        ]
        .astype(str)
        .to_numpy()
    )

    donor_arr = np.asarray(
        donors,
        dtype=object,
    )

    cell_ids = (
        selected[
            "cell_id"
        ]
        .astype(str)
        .to_numpy()
    )

    obs = pd.DataFrame(
        index=pd.Index(
            cell_ids,
            name=None,
        )
    )

    # Same 34 column names and order as Huuki.
    obs["Sample"] = donor_arr
    obs["Barcode"] = barcodes
    obs["key"] = cell_ids
    obs["SAMPLE_ID"] = donor_arr

    obs["pos"] = np.nan
    obs["BrNum"] = donor_arr
    obs["round"] = np.nan
    obs["Position"] = np.nan

    obs["age"] = ages
    obs["sex"] = sexes
    obs["diagnosis"] = diagnoses

    obs["sum"] = library_size
    obs["detected"] = detected.astype(
        np.float64
    )

    obs["subsets_Mito_sum"] = mito_sum
    obs["subsets_Mito_detected"] = (
        mito_detected.astype(
            np.float64
        )
    )

    obs["subsets_Mito_percent"] = (
        mito_percent
    )

    obs["total"] = library_size

    # These Huuki-specific QC flags do not have a direct
    # Ling equivalent. Selected cells come from final jointMetadata,
    # so they are retained cells, but we do NOT fabricate thresholds.
    obs["high_mito"] = np.nan
    obs["low_sum"] = np.nan
    obs["low_detected"] = np.nan
    obs["discard_auto"] = np.nan

    obs["doubletScore"] = np.nan

    obs["prelimCluster"] = "NA"
    obs["collapsedCluster"] = "NA"
    obs["kmeans"] = "NA"

    obs["sizeFactor"] = size_factor

    # Populate both naming conventions with true Ling annotations
    # for downstream compatibility.
    obs["cellType_broad_k"] = broad
    obs["cellType_k"] = subclass

    obs["cellType_broad_hc"] = broad
    obs["cellType_hc"] = subclass

    # Ling subtype is NOT falsely relabeled as cortical layer.
    obs["cellType_layer"] = "NA"
    obs["layer_annotation"] = "NA"

    obs["source_chunk"] = prefixes
    obs["merge_batch"] = prefixes


    HUUKI_OBS_COLUMNS = [
        "Sample",
        "Barcode",
        "key",
        "SAMPLE_ID",
        "pos",
        "BrNum",
        "round",
        "Position",
        "age",
        "sex",
        "diagnosis",
        "sum",
        "detected",
        "subsets_Mito_sum",
        "subsets_Mito_detected",
        "subsets_Mito_percent",
        "total",
        "high_mito",
        "low_sum",
        "low_detected",
        "discard_auto",
        "doubletScore",
        "prelimCluster",
        "collapsedCluster",
        "kmeans",
        "sizeFactor",
        "cellType_broad_k",
        "cellType_k",
        "cellType_broad_hc",
        "cellType_hc",
        "cellType_layer",
        "layer_annotation",
        "source_chunk",
        "merge_batch",
    ]

    obs = obs[
        HUUKI_OBS_COLUMNS
    ]

    if obs.shape != (
        n_cells,
        34,
    ):
        raise RuntimeError(
            f"Unexpected obs shape: {obs.shape}"
        )


    # =========================================================================
    # FINAL ANNDATA
    # =========================================================================

    print("\n" + "=" * 100)
    print("WRITING FINAL H5AD")
    print("=" * 100)

    var = pd.DataFrame(
        index=pd.Index(
            var_names,
            name=None,
        )
    )

    adata = ad.AnnData(
        X=Xmm,
        obs=obs,
        var=var,
    )

    adata.layers[
        "logcounts"
    ] = Lmm

    print("Shape:", adata.shape)
    print("obs:", adata.obs.shape)
    print("var:", adata.var.shape)
    print("layers:", list(adata.layers.keys()))
    print("obsm:", list(adata.obsm.keys()))
    print("X dtype:", adata.X.dtype)
    print(
        "logcounts dtype:",
        adata.layers[
            "logcounts"
        ].dtype,
    )

    adata.write_h5ad(
        OUTPUT,
        compression="gzip",
        compression_opts=4,
    )

    del adata

    Xmm.flush()
    Lmm.flush()

    del Xmm
    del Lmm

    gc.collect()


    # =========================================================================
    # FINAL ON-DISK VALIDATION
    # =========================================================================

    print("\n" + "=" * 100)
    print("FINAL ON-DISK VALIDATION")
    print("=" * 100)

    chk = ad.read_h5ad(
        OUTPUT,
        backed="r",
    )

    print(
        "Shape:",
        chk.shape
    )

    print(
        "OBS:",
        chk.obs.shape
    )

    print(
        "VAR:",
        chk.var.shape
    )

    print(
        "Layers:",
        list(
            chk.layers.keys()
        )
    )

    print(
        "OBSM:",
        list(
            chk.obsm.keys()
        )
    )

    print(
        "var_names unique:",
        chk.var_names.is_unique
    )

    print(
        "Kwon symbols present:",
        len(
            kwon_symbols
            & set(
                chk.var_names
            )
        ),
        "/",
        len(
            kwon_symbols
        ),
    )

    chk.file.close()


    with h5py.File(
        OUTPUT,
        "r",
    ) as h:

        print(
            "HDF5 X:",
            h["X"]
        )

        print(
            "HDF5 logcounts:",
            h[
                "layers/logcounts"
            ]
        )

        test_rows = np.asarray(
            [
                0,
                1,
                10,
                100,
                1000,
            ],
            dtype=int,
        )

        test_rows = test_rows[
            test_rows < n_cells
        ]

        X_test = np.asarray(
            h[
                "X"
            ][
                test_rows,
                :
            ],
            dtype=np.float64,
        )

        L_test = np.asarray(
            h[
                "layers/logcounts"
            ][
                test_rows,
                :
            ],
            dtype=np.float64,
        )

        expected = np.log2(
            X_test
            / size_factor[
                test_rows,
                None
            ]
            + 1.0
        )

        max_error = float(
            np.max(
                np.abs(
                    expected
                    - L_test
                )
            )
        )

        print(
            "Huuki-normalization max abs error:",
            max_error
        )


    summary = {
        "output":
            str(OUTPUT),

        "n_cells":
            int(n_cells),

        "n_genes":
            int(n_genes),

        "source_jointMetadata_cells":
            int(len(meta)),

        "eligible_annotated_cells":
            int(len(eligible)),

        "selected_donors":
            int(
                selected[
                    "DONOR"
                ].nunique()
            ),

        "selected_libraries":
            int(
                selected[
                    "PREFIX"
                ].nunique()
            ),

        "gene_prevalence_threshold":
            float(
                args.gene_prevalence
            ),

        "kwon_genes_present":
            int(
                len(
                    kwon_symbols
                )
            ),

        "mean_library_size":
            float(
                mean_library_size
            ),

        "mean_sizeFactor":
            float(
                size_factor.mean()
            ),

        "normalization":
            "log2(count / sizeFactor + 1)",

        "normalization_max_abs_error":
            max_error,

        "seed":
            int(
                args.seed
            ),
    }

    with open(
        BUILD_DIR
        / "ling_reference_build_summary.json",
        "w",
    ) as f:
        json.dump(
            summary,
            f,
            indent=2,
        )


    # =========================================================================
    # CLEAN TEMP MATRICES ONLY AFTER SUCCESS
    # =========================================================================

    if not args.keep_temp:

        for p in [
            x_path,
            log_path,
        ]:
            if p.exists():
                p.unlink()

        print(
            "\nTemporary dense memmaps removed."
        )

    else:
        print(
            "\nTemporary memmaps retained:",
            WORK_DIR
        )


    print("\n" + "=" * 100)
    print("BUILD COMPLETE")
    print("=" * 100)

    print(OUTPUT)


if __name__ == "__main__":
    main()
