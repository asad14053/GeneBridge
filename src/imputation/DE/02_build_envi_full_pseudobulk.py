#!/usr/bin/env python3

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = ROOT / "data/metadata/envi_production_23donors.tsv"

OUT = (
    ROOT
    / "outputs/imputation_full/DE/envi_transcriptome"
)
OUT.mkdir(parents=True, exist_ok=True)

CHUNK = 512


def sum_chunk(x):
    if sparse.issparse(x):
        return np.asarray(
            x.sum(axis=0)
        ).ravel()

    return np.asarray(
        x,
        dtype=np.float64
    ).sum(
        axis=0,
        dtype=np.float64
    )


print("=" * 100)
print("BUILD ENVI FULL-TRANSCRIPTOME DONOR PSEUDOBULK")
print("=" * 100)

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
assert (manifest["Dx"] == "NTC").sum() == 11
assert (manifest["Dx"] == "SCZ").sum() == 12


canonical_genes = None
canonical_source = None

rows = []
qc = []


for i, r in enumerate(
    manifest.itertuples(index=False),
    1
):

    donor = r.donor
    experiment = r.experiment

    path = (
        ROOT
        / f"data/processed/imputation_full/{experiment}/envi/{donor}"
        / f"spatial_data_xenium_{donor}_ENVI_full_transcriptome.h5ad"
    )

    if not path.exists():
        raise FileNotFoundError(path)

    print(
        f"\n[{i:02d}/23] {donor} | {r.Dx}",
        flush=True
    )

    a = ad.read_h5ad(
        path,
        backed="r"
    )

    genes = (
        a.var_names
        .astype(str)
        .tolist()
    )

    source = (
        a.var["expression_source"]
        .astype(str)
        .tolist()
    )

    if canonical_genes is None:

        canonical_genes = genes
        canonical_source = source

        # First donor defines canonical gene order.
        reorder_idx = np.arange(len(genes))

    else:

        # Order may differ, but gene identity must be identical.
        if set(genes) != set(canonical_genes):

            missing = sorted(
                set(canonical_genes) - set(genes)
            )

            extra = sorted(
                set(genes) - set(canonical_genes)
            )

            raise RuntimeError(
                f"{donor}: gene SET differs\\n"
                f"Missing ({len(missing)}): {missing[:20]}\\n"
                f"Extra ({len(extra)}): {extra[:20]}"
            )

        donor_index = pd.Index(genes)

        reorder_idx = donor_index.get_indexer(
            canonical_genes
        )

        if np.any(reorder_idx < 0):
            raise RuntimeError(
                f"{donor}: gene alignment failed"
            )

        source_aligned = [
            source[j]
            for j in reorder_idx
        ]

        if source_aligned != canonical_source:
            raise RuntimeError(
                f"{donor}: expression_source mismatch "
                f"after gene alignment"
            )

        if genes != canonical_genes:
            print(
                "    NOTE: gene order differs; "
                "aligning by gene name.",
                flush=True
            )

    if a.n_vars != 34987:
        raise RuntimeError(
            f"{donor}: expected 34,987 genes; "
            f"found {a.n_vars}"
        )

    n_measured = sum(
        x == "measured_xenium"
        for x in source
    )

    n_imputed = sum(
        x == "envi_imputed"
        for x in source
    )

    if n_measured != 300:
        raise RuntimeError(
            f"{donor}: measured={n_measured}"
        )

    if n_imputed != 34687:
        raise RuntimeError(
            f"{donor}: imputed={n_imputed}"
        )

    # Use the explicitly named final count-scale layer.
    matrix = a.layers["count_scale"]

    total = np.zeros(
        a.n_vars,
        dtype=np.float64
    )

    for start in range(
        0,
        a.n_obs,
        CHUNK
    ):

        stop = min(
            start + CHUNK,
            a.n_obs
        )

        x = matrix[
            start:stop,
            :
        ]

        total += sum_chunk(x)

        if (
            start == 0
            or stop == a.n_obs
            or start % (CHUNK * 20) == 0
        ):
            print(
                f"    cells {stop:,}/{a.n_obs:,}",
                flush=True
            )

    # Convert this donor from its native gene order
    # into the canonical Br2039 gene order.
    total = total[reorder_idx]

    if not np.isfinite(total).all():
        raise RuntimeError(
            f"{donor}: non-finite pseudobulk values"
        )

    if np.any(total < 0):
        raise RuntimeError(
            f"{donor}: negative pseudobulk values"
        )

    rows.append(
        pd.Series(
            total,
            index=genes,
            name=donor
        )
    )

    qc.append({
        "donor": donor,
        "Dx": r.Dx,
        "n_cells": a.n_obs,
        "n_genes": a.n_vars,
        "measured_genes": n_measured,
        "imputed_genes": n_imputed,
        "pseudobulk_library": total.sum(),
    })

    a.file.close()


pb = pd.DataFrame(rows)

pb = pb.loc[
    manifest["donor"].tolist()
]

assert pb.shape == (23, 34987)

print("\nFinal pseudobulk shape:", pb.shape)


# Full 23 × 34,987 expected count-scale matrix
pb.to_csv(
    OUT
    / "ENVI_full34987_pseudobulk_countscale_23donors.csv.gz",
    compression="gzip"
)


# Gene annotation
gene_info = pd.DataFrame({
    "gene": canonical_genes,
    "expression_source": canonical_source,
})

gene_info.to_csv(
    OUT
    / "ENVI_full34987_gene_info.csv",
    index=False
)


pd.DataFrame(qc).to_csv(
    OUT
    / "ENVI_full34987_pseudobulk_donor_qc.csv",
    index=False
)


manifest[
    ["donor", "Dx"]
].to_csv(
    OUT
    / "ENVI_full34987_diagnosis.csv",
    index=False
)


print("\nGene source:")
print(
    gene_info[
        "expression_source"
    ].value_counts()
)

print("\nDiagnosis:")
print(
    manifest[
        "Dx"
    ].value_counts()
)

print("\nSUCCESS")
print(
    OUT
    / "ENVI_full34987_pseudobulk_countscale_23donors.csv.gz"
)
