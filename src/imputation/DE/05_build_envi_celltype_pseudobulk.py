#!/usr/bin/env python3

from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = ROOT / "data/metadata/envi_production_23donors.tsv"

ANNOT = (
    ROOT
    / "data/processed/xenium/"
    / "xenium_N24_layer_celltype_annotated.h5ad"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/envi_celltype/pseudobulk"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


PRIMARY_CELL_TYPES = [
    "Oligo",
    "Endo",
    "Ast",
    "L2/3 Ex",
    "Mic",
    "MGE",
    "L4/5 Ex",
    "CGE",
    "L6 Ex",
    "L5 Ex",
]


CHUNK = 512


def section(title):

    print(
        "\n" + "=" * 110
    )

    print(title)

    print(
        "=" * 110
    )


def slugify(x):

    x = x.replace("/", "_")
    x = re.sub(r"[^A-Za-z0-9]+", "_", x)
    x = x.strip("_")

    return x


def chunk_to_array(x):

    if sparse.issparse(x):
        return x.toarray()

    return np.asarray(x)


# =============================================================================
# Manifest
# =============================================================================

section(
    "LOAD 23-DONOR MANIFEST"
)

manifest = pd.read_csv(
    MANIFEST,
    sep="\t"
)

manifest["donor"] = (
    manifest["donor"]
    .astype(str)
    .str.strip()
)

manifest["Dx"] = (
    manifest["Dx"]
    .astype(str)
    .str.strip()
    .str.upper()
)

assert len(manifest) == 23
assert manifest["donor"].nunique() == 23
assert (manifest["Dx"] == "SCZ").sum() == 12
assert (manifest["Dx"] == "NTC").sum() == 11

print(
    manifest[
        ["donor", "Dx", "experiment"]
    ].to_string(index=False)
)


# =============================================================================
# Annotation map
# =============================================================================

section(
    "LOAD CELL-TYPE ANNOTATION"
)

ann = ad.read_h5ad(
    ANNOT,
    backed="r"
)

required = [
    "donor_cell_id",
    "BrNum",
    "cell_type_annotation",
]

missing = [
    c
    for c in required
    if c not in ann.obs.columns
]

if missing:

    raise RuntimeError(
        f"Missing annotation columns: {missing}"
    )


annotation = ann.obs[
    required
].copy()

ann.file.close()


annotation["donor_cell_id"] = (
    annotation["donor_cell_id"]
    .astype(str)
    .str.strip()
)

annotation["BrNum"] = (
    annotation["BrNum"]
    .astype(str)
    .str.strip()
)


# Keep only our 23 usable donors.
annotation = annotation.loc[
    annotation["BrNum"].isin(
        manifest["donor"]
    )
].copy()


if annotation[
    "donor_cell_id"
].duplicated().any():

    raise RuntimeError(
        "donor_cell_id is not unique in annotation source."
    )


annotation = annotation.set_index(
    "donor_cell_id"
)


print(
    "Annotation cells:",
    f"{len(annotation):,}"
)

print("\nPrimary cell types:")

for ct in PRIMARY_CELL_TYPES:

    n = int(
        (
            annotation[
                "cell_type_annotation"
            ].astype(str)
            == ct
        ).sum()
    )

    print(
        f"  {ct:10s}: {n:,}"
    )


# =============================================================================
# Containers
# =============================================================================

pb = {
    ct: []
    for ct in PRIMARY_CELL_TYPES
}

cell_count_rows = []

canonical_genes = None
canonical_source = None

total_unmatched_cells = 0


# =============================================================================
# Donor loop
# =============================================================================

section(
    "BUILD DONOR × CELL-TYPE PSEUDOBULK"
)


for donor_number, r in enumerate(
    manifest.itertuples(index=False),
    start=1,
):

    donor = str(r.donor)
    dx = str(r.Dx)
    experiment = str(r.experiment)


    path = (
        ROOT
        / f"data/processed/imputation_full/{experiment}/envi/{donor}"
        / f"spatial_data_xenium_{donor}_ENVI_full_transcriptome.h5ad"
    )


    if not path.exists():

        raise FileNotFoundError(
            path
        )


    print(
        f"\n[{donor_number:02d}/23] "
        f"{donor} | {dx}",
        flush=True
    )


    a = ad.read_h5ad(
        path,
        backed="r"
    )


    # -------------------------------------------------------------------------
    # Gene universe
    # -------------------------------------------------------------------------

    genes = (
        a.var_names
        .astype(str)
        .tolist()
    )


    if a.n_vars != 34987:

        raise RuntimeError(
            f"{donor}: expected 34,987 genes; "
            f"found {a.n_vars}"
        )


    if "expression_source" not in a.var.columns:

        raise RuntimeError(
            f"{donor}: expression_source missing."
        )


    source = (
        a.var[
            "expression_source"
        ]
        .astype(str)
        .tolist()
    )


    if canonical_genes is None:

        canonical_genes = genes

        canonical_source = source

        reorder_idx = np.arange(
            len(genes)
        )


    else:

        if set(genes) != set(
            canonical_genes
        ):

            missing_genes = sorted(
                set(canonical_genes)
                - set(genes)
            )

            extra_genes = sorted(
                set(genes)
                - set(canonical_genes)
            )

            raise RuntimeError(
                f"{donor}: gene SET differs.\n"
                f"Missing={missing_genes[:20]}\n"
                f"Extra={extra_genes[:20]}"
            )


        donor_index = pd.Index(
            genes
        )

        reorder_idx = (
            donor_index
            .get_indexer(
                canonical_genes
            )
        )


        if np.any(
            reorder_idx < 0
        ):

            raise RuntimeError(
                f"{donor}: gene alignment failed."
            )


        source_aligned = [
            source[j]
            for j in reorder_idx
        ]


        if source_aligned != canonical_source:

            raise RuntimeError(
                f"{donor}: expression_source "
                "differs after gene alignment."
            )


        if genes != canonical_genes:

            print(
                "  gene order differs -> "
                "aligning by gene name",
                flush=True
            )


    # -------------------------------------------------------------------------
    # Cell annotation alignment
    # -------------------------------------------------------------------------

    if "donor_cell_id" not in a.obs.columns:

        raise RuntimeError(
            f"{donor}: donor_cell_id missing."
        )


    cell_ids = (
        a.obs[
            "donor_cell_id"
        ]
        .astype(str)
        .str.strip()
    )


    if cell_ids.duplicated().any():

        raise RuntimeError(
            f"{donor}: duplicate donor_cell_id."
        )


    aligned_annotation = annotation.reindex(
        cell_ids.to_numpy()
    )


    # Confirm matching annotation belongs to same donor.
    matched_record = (
        aligned_annotation[
            "BrNum"
        ].notna()
    )


    wrong_donor = (
        matched_record
        & (
            aligned_annotation[
                "BrNum"
            ].astype(str)
            != donor
        )
    )


    if wrong_donor.any():

        raise RuntimeError(
            f"{donor}: annotation donor mismatch."
        )


    labels = (
        aligned_annotation[
            "cell_type_annotation"
        ]
        .astype("object")
        .to_numpy()
    )


    n_unmatched = int(
        pd.isna(
            labels
        ).sum()
    )


    total_unmatched_cells += (
        n_unmatched
    )


    print(
        f"  cells            : {a.n_obs:,}"
    )

    print(
        f"  unmatched labels : {n_unmatched:,}"
    )


    # -------------------------------------------------------------------------
    # Cell counts
    # -------------------------------------------------------------------------

    for ct in PRIMARY_CELL_TYPES:

        n_ct = int(
            np.sum(
                labels == ct
            )
        )


        if n_ct == 0:

            raise RuntimeError(
                f"{donor}: no cells for {ct}"
            )


        cell_count_rows.append({
            "donor":
                donor,

            "Dx":
                dx,

            "cell_type":
                ct,

            "n_cells":
                n_ct,
        })


        print(
            f"  {ct:10s}: {n_ct:,}"
        )


    # -------------------------------------------------------------------------
    # Chunked pseudobulk
    #
    # Read each 34,987-gene chunk ONCE and aggregate all 10 cell types.
    # -------------------------------------------------------------------------

    if "count_scale" not in a.layers:

        raise RuntimeError(
            f"{donor}: count_scale layer missing."
        )


    matrix = a.layers[
        "count_scale"
    ]


    donor_totals = {
        ct: np.zeros(
            a.n_vars,
            dtype=np.float64
        )
        for ct in PRIMARY_CELL_TYPES
    }


    for start in range(
        0,
        a.n_obs,
        CHUNK
    ):

        stop = min(
            start + CHUNK,
            a.n_obs
        )


        x = chunk_to_array(
            matrix[
                start:stop,
                :
            ]
        )


        if x.ndim != 2:

            raise RuntimeError(
                f"{donor}: unexpected chunk shape {x.shape}"
            )


        chunk_labels = labels[
            start:stop
        ]


        for ct in PRIMARY_CELL_TYPES:

            mask = (
                chunk_labels == ct
            )


            if np.any(mask):

                donor_totals[
                    ct
                ] += np.asarray(
                    x[
                        mask,
                        :
                    ]
                ).sum(
                    axis=0,
                    dtype=np.float64
                )


        if (
            start == 0
            or stop == a.n_obs
            or start % (
                CHUNK * 20
            ) == 0
        ):

            print(
                f"    cells "
                f"{stop:,}/{a.n_obs:,}",
                flush=True
            )


    # -------------------------------------------------------------------------
    # Canonical gene order
    # -------------------------------------------------------------------------

    for ct in PRIMARY_CELL_TYPES:

        total = donor_totals[
            ct
        ][
            reorder_idx
        ]


        if not np.isfinite(
            total
        ).all():

            raise RuntimeError(
                f"{donor} / {ct}: "
                "non-finite pseudobulk values."
            )


        if np.any(
            total < 0
        ):

            raise RuntimeError(
                f"{donor} / {ct}: "
                "negative pseudobulk values."
            )


        pb[
            ct
        ].append(
            pd.Series(
                total,
                index=canonical_genes,
                name=donor,
                dtype=np.float64,
            )
        )


    a.file.close()


# =============================================================================
# Final QC
# =============================================================================

section(
    "WRITE CELL-TYPE PSEUDOBULKS"
)


gene_info = pd.DataFrame({
    "gene":
        canonical_genes,

    "expression_source":
        canonical_source,
})


print(
    "\nGene source:"
)

print(
    gene_info[
        "expression_source"
    ].value_counts()
)


if (
    gene_info[
        "expression_source"
    ]
    .eq(
        "measured_xenium"
    )
    .sum()
    != 300
):

    raise RuntimeError(
        "Expected 300 measured genes."
    )


if (
    gene_info[
        "expression_source"
    ]
    .eq(
        "envi_imputed"
    )
    .sum()
    != 34687
):

    raise RuntimeError(
        "Expected 34,687 imputed genes."
    )


gene_info.to_csv(
    OUT
    / "ENVI_full34987_gene_info.csv",
    index=False
)


mapping_rows = []


for ct in PRIMARY_CELL_TYPES:

    slug = slugify(
        ct
    )

    df = pd.DataFrame(
        pb[
            ct
        ]
    )


    df = df.loc[
        manifest[
            "donor"
        ].tolist()
    ]


    if df.shape != (
        23,
        34987
    ):

        raise RuntimeError(
            f"{ct}: unexpected shape {df.shape}"
        )


    out_file = (
        OUT
        / f"ENVI_{slug}_full34987_"
          f"pseudobulk_countscale_23donors.csv.gz"
    )


    df.to_csv(
        out_file,
        compression="gzip"
    )


    print(
        f"{ct:10s}: "
        f"{df.shape} -> "
        f"{out_file.name}"
    )


    mapping_rows.append({
        "cell_type":
            ct,

        "slug":
            slug,

        "pseudobulk_file":
            str(
                out_file
            ),
    })


cell_counts = pd.DataFrame(
    cell_count_rows
)


cell_counts.to_csv(
    OUT
    / "celltype_donor_cell_counts.csv",
    index=False
)


pd.DataFrame(
    mapping_rows
).to_csv(
    OUT
    / "celltype_file_manifest.csv",
    index=False
)


manifest[
    [
        "donor",
        "Dx",
    ]
].to_csv(
    OUT
    / "donor_diagnosis.csv",
    index=False
)


# =============================================================================
# Summary
# =============================================================================

section(
    "FINAL SUMMARY"
)


print(
    "Donors              :",
    len(
        manifest
    )
)

print(
    "SCZ donors          :",
    int(
        (
            manifest[
                "Dx"
            ] == "SCZ"
        ).sum()
    )
)

print(
    "NTC donors          :",
    int(
        (
            manifest[
                "Dx"
            ] == "NTC"
        ).sum()
    )
)

print(
    "Primary cell types  :",
    len(
        PRIMARY_CELL_TYPES
    )
)

print(
    "Total genes         :",
    len(
        canonical_genes
    )
)

print(
    "Measured genes      :",
    300
)

print(
    "Imputed genes       :",
    34687
)

print(
    "Unmatched cells     :",
    total_unmatched_cells
)


if total_unmatched_cells != 3:

    raise RuntimeError(
        "Expected exactly 3 unmatched cells "
        f"from previous audit; found {total_unmatched_cells}."
    )


print(
    "\nSUCCESS"
)

print(
    "Output:",
    OUT
)
