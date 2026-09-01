#!/usr/bin/env python3

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MANIFEST = ROOT / "data/metadata/envi_production_23donors.tsv"
DONOR_META = ROOT / "data/metadata/xenium_DE_metadata_23.csv"

ANNOTATED = (
    ROOT
    / "data/processed/xenium/xenium_N24_layer_celltype_annotated.h5ad"
)

ORIGINAL_META = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "original_layer_pseudobulk/"
      "original_xenium_300gene_donor_layer_pseudobulk_metadata.csv"
)

ORIGINAL_GENES = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "original_layer_pseudobulk/"
      "original_xenium_300gene_order.txt"
)

OUT = (
    ROOT
    / "outputs/imputation_full/DE/paper_concordance/"
      "envi_layer_pseudobulk/per_donor"
)

OUT.mkdir(parents=True, exist_ok=True)

EXPECTED_LAYERS = [
    "L1/M",
    "L2/3",
    "L3/4",
    "L5",
    "L6",
    "WMtz",
    "WM",
]

LAYER_TO_SPD = {
    "L1/M": "spd07",
    "L2/3": "spd06",
    "L3/4": "spd02",
    "L5": "spd05",
    "L6": "spd03",
    "WMtz": "spd01",
    "WM": "spd04",
}

CHUNK = 512
MIN_NCELLS = 10


def sum_block(x):
    if sparse.issparse(x):
        return np.asarray(x.sum(axis=0)).ravel().astype(np.float64)

    return np.asarray(
        x,
        dtype=np.float64,
    ).sum(
        axis=0,
        dtype=np.float64,
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--index",
    required=True,
    type=int,
    help="1-based manifest row / SLURM array index",
)

args = parser.parse_args()

if not 1 <= args.index <= 23:
    raise ValueError("--index must be between 1 and 23")


# =============================================================================
# Resolve donor
# =============================================================================

manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
)

manifest["donor"] = manifest["donor"].astype(str).str.strip()
manifest["Dx"] = manifest["Dx"].astype(str).str.strip().str.upper()

if len(manifest) != 23:
    raise RuntimeError(f"Expected 23 donors; found {len(manifest)}")

r = manifest.iloc[args.index - 1]

donor = r["donor"]
dx = r["Dx"]
experiment = r["experiment"]
target_path = Path(r["target"])


print("=" * 100)
print("ENVI PAPER-CONCORDANCE DONOR × LAYER PSEUDOBULK")
print("=" * 100)
print("Array index:", args.index)
print("Donor:", donor)
print("Dx:", dx)
print("Experiment:", experiment)


# =============================================================================
# Canonical donor metadata
# =============================================================================

dm = pd.read_csv(DONOR_META)
dm["BrNum"] = dm["BrNum"].astype(str).str.strip()

dm = dm.loc[
    dm["BrNum"] == donor
]

if len(dm) != 1:
    raise RuntimeError(
        f"{donor}: canonical metadata row count = {len(dm)}"
    )

dm = dm.iloc[0]

if str(dm["Dx"]).upper() != dx:
    raise RuntimeError(
        f"{donor}: diagnosis mismatch"
    )


# =============================================================================
# Step-2 original pseudobulk structure
# =============================================================================

orig_meta = pd.read_csv(ORIGINAL_META)

orig_donor = (
    orig_meta.loc[
        orig_meta["BrNum"].astype(str) == donor
    ]
    .copy()
)

if len(orig_donor) != 7:
    raise RuntimeError(
        f"{donor}: expected 7 original layer pseudobulks; "
        f"found {len(orig_donor)}"
    )

original_genes = [
    x.strip()
    for x in ORIGINAL_GENES.read_text().splitlines()
    if x.strip()
]

if len(original_genes) != 300:
    raise RuntimeError("Expected 300 original Xenium genes")


# =============================================================================
# Cell annotations
# =============================================================================

print("\nLoading cell annotations...")

ann = ad.read_h5ad(
    ANNOTATED,
    backed="r",
)

mask = (
    ann.obs["BrNum"]
    .astype(str)
    .eq(donor)
)

annotation = (
    ann.obs.loc[
        mask,
        [
            "BrNum",
            "layer_annotation",
            "predictions_smooth",
        ]
    ]
    .copy()
)

ann.file.close()

annotation["BrNum"] = annotation["BrNum"].astype(str)
annotation["layer_annotation"] = (
    annotation["layer_annotation"].astype(str)
)
annotation["predictions_smooth"] = (
    annotation["predictions_smooth"].astype(str)
)

print("Annotated donor cells:", len(annotation))


# =============================================================================
# Open original target + ENVI
# =============================================================================

after_path = (
    ROOT
    / f"data/processed/imputation_full/{experiment}/envi/{donor}"
      f"/spatial_data_xenium_{donor}_ENVI_full_transcriptome.h5ad"
)

if not target_path.exists():
    raise FileNotFoundError(target_path)

if not after_path.exists():
    raise FileNotFoundError(after_path)


target = ad.read_h5ad(
    target_path,
    backed="r",
)

after = ad.read_h5ad(
    after_path,
    backed="r",
)


target_ids = pd.Index(
    target.obs_names.astype(str)
)

after_ids = pd.Index(
    after.obs_names.astype(str)
)


if len(target_ids) != len(after_ids):
    raise RuntimeError(
        f"{donor}: target/ENVI cell counts differ"
    )

if not np.array_equal(
    target_ids.to_numpy(),
    after_ids.to_numpy(),
):
    raise RuntimeError(
        f"{donor}: target/ENVI cell IDs or order differ"
    )


ann_ids = pd.Index(
    annotation.index.astype(str)
)

target_only = target_ids.difference(ann_ids)
annotated_only = ann_ids.difference(target_ids)


print("Target cells:", len(target_ids))
print("Target-only unannotated:", len(target_only))
print("Annotated-only:", len(annotated_only))


if len(annotated_only) > 0:
    raise RuntimeError(
        f"{donor}: {len(annotated_only)} annotated cells "
        "missing from ENVI target"
    )


excluded_rows = [
    {
        "BrNum": donor,
        "cell_id": str(x),
        "reason": "missing_layer_annotation",
    }
    for x in target_only
]


# =============================================================================
# Gene validation
# =============================================================================

genes = after.var_names.astype(str).tolist()

if after.n_vars != 34987:
    raise RuntimeError(
        f"{donor}: expected 34987 genes; found {after.n_vars}"
    )

if "expression_source" not in after.var.columns:
    raise RuntimeError(
        f"{donor}: expression_source missing"
    )

source = (
    after.var["expression_source"]
    .astype(str)
    .tolist()
)

n_measured = sum(x == "measured_xenium" for x in source)
n_imputed = sum(x == "envi_imputed" for x in source)

if n_measured != 300:
    raise RuntimeError(
        f"{donor}: measured genes = {n_measured}"
    )

if n_imputed != 34687:
    raise RuntimeError(
        f"{donor}: imputed genes = {n_imputed}"
    )


gene_index = pd.Index(genes)

measured_idx = gene_index.get_indexer(
    original_genes
)

if np.any(measured_idx < 0):
    missing = np.asarray(original_genes)[
        measured_idx < 0
    ]
    raise RuntimeError(
        f"{donor}: original genes missing: {missing[:20]}"
    )


if not all(
    source[j] == "measured_xenium"
    for j in measured_idx
):
    raise RuntimeError(
        f"{donor}: one or more original genes are not "
        "marked measured_xenium"
    )


# =============================================================================
# Matrix sources
# =============================================================================

# Matches existing AFTER-ENVI measured-300 QC.
measured_matrix = after.X

# Matches existing full-transcriptome ENVI pseudobulk.
if "count_scale" not in after.layers:
    raise RuntimeError(
        f"{donor}: count_scale layer missing"
    )

full_matrix = after.layers["count_scale"]


# =============================================================================
# Allocate layer totals
# =============================================================================

measured_totals = {
    layer: np.zeros(300, dtype=np.float64)
    for layer in EXPECTED_LAYERS
}

full_totals = {
    layer: np.zeros(34987, dtype=np.float64)
    for layer in EXPECTED_LAYERS
}

layer_cells = {
    layer: 0
    for layer in EXPECTED_LAYERS
}

matched_measured_total = np.zeros(
    300,
    dtype=np.float64,
)

matched_full_total = np.zeros(
    34987,
    dtype=np.float64,
)


# =============================================================================
# Chunk through donor cells
# =============================================================================

print("\nAggregating cells...")

for start in range(
    0,
    after.n_obs,
    CHUNK,
):

    stop = min(
        start + CHUNK,
        after.n_obs,
    )

    chunk_ids = target_ids[start:stop]

    keep = chunk_ids.isin(
        ann_ids
    )

    if not keep.any():
        continue

    keep_positions = np.flatnonzero(
        keep
    )

    kept_ids = chunk_ids[
        keep
    ]

    kept_ann = annotation.loc[
        kept_ids,
        [
            "layer_annotation",
            "predictions_smooth",
        ]
    ]


    # ------------------------------------------------------------
    # AFTER-ENVI measured 300
    # ------------------------------------------------------------

    m_all = measured_matrix[
        start:stop,
        measured_idx,
    ]

    if sparse.issparse(m_all):
        m = m_all[
            keep_positions,
            :
        ]
    else:
        m = np.asarray(
            m_all,
            dtype=np.float64,
        )[
            keep_positions,
            :
        ]


    # ------------------------------------------------------------
    # Full count-scale transcriptome
    # ------------------------------------------------------------

    f_all = full_matrix[
        start:stop,
        :
    ]

    if sparse.issparse(f_all):
        f = f_all[
            keep_positions,
            :
        ]
    else:
        f = np.asarray(
            f_all,
            dtype=np.float64,
        )[
            keep_positions,
            :
        ]


    matched_measured_total += sum_block(m)
    matched_full_total += sum_block(f)


    # ------------------------------------------------------------
    # Split into spatial layers
    # ------------------------------------------------------------

    for layer in EXPECTED_LAYERS:

        spd = LAYER_TO_SPD[layer]

        lm = (
            kept_ann["layer_annotation"]
            .eq(layer)
            .to_numpy()
        )

        if not lm.any():
            continue


        if not (
            kept_ann.loc[
                lm,
                "predictions_smooth",
            ]
            == spd
        ).all():
            raise RuntimeError(
                f"{donor} {layer}: SpD/layer mismatch"
            )


        layer_cells[layer] += int(
            lm.sum()
        )

        measured_totals[layer] += sum_block(
            m[lm, :]
        )

        full_totals[layer] += sum_block(
            f[lm, :]
        )


    if (
        start == 0
        or stop == after.n_obs
        or start % (CHUNK * 20) == 0
    ):
        print(
            f"cells {stop:,}/{after.n_obs:,}",
            flush=True,
        )


# =============================================================================
# Validation
# =============================================================================

print("\nLayer counts:")

for layer in EXPECTED_LAYERS:
    print(
        f"{layer:4s}: {layer_cells[layer]:,}"
    )


for layer in EXPECTED_LAYERS:

    if layer_cells[layer] < MIN_NCELLS:
        raise RuntimeError(
            f"{donor} {layer}: "
            f"{layer_cells[layer]} cells < 10"
        )

    expected_n = int(
        orig_donor.loc[
            orig_donor["layer_annotation"] == layer,
            "n_cells",
        ].iloc[0]
    )

    if layer_cells[layer] != expected_n:
        raise RuntimeError(
            f"{donor} {layer}: "
            f"ENVI cells={layer_cells[layer]} "
            f"original cells={expected_n}"
        )


measured_layer_sum = np.stack(
    [
        measured_totals[x]
        for x in EXPECTED_LAYERS
    ]
).sum(axis=0)

full_layer_sum = np.stack(
    [
        full_totals[x]
        for x in EXPECTED_LAYERS
    ]
).sum(axis=0)


if not np.allclose(
    measured_layer_sum,
    matched_measured_total,
    rtol=1e-10,
    atol=1e-6,
):
    raise RuntimeError(
        f"{donor}: measured layer sum mismatch"
    )


if not np.allclose(
    full_layer_sum,
    matched_full_total,
    rtol=1e-10,
    atol=1e-6,
):
    raise RuntimeError(
        f"{donor}: full layer sum mismatch"
    )


# =============================================================================
# Build 7-row donor outputs
# =============================================================================

measured_rows = []
full_rows = []
metadata_rows = []


for layer in EXPECTED_LAYERS:

    orig_row = (
        orig_donor.loc[
            orig_donor["layer_annotation"] == layer
        ]
        .iloc[0]
    )

    pb_id = str(
        orig_row["pseudobulk_id"]
    )

    mv = measured_totals[layer]
    fv = full_totals[layer]


    if not np.isfinite(mv).all():
        raise RuntimeError(
            f"{pb_id}: non-finite measured values"
        )

    if not np.isfinite(fv).all():
        raise RuntimeError(
            f"{pb_id}: non-finite full values"
        )

    if np.any(mv < 0) or np.any(fv < 0):
        raise RuntimeError(
            f"{pb_id}: negative values"
        )


    measured_rows.append(
        pd.Series(
            mv,
            index=original_genes,
            name=pb_id,
        )
    )

    full_rows.append(
        pd.Series(
            fv,
            index=genes,
            name=pb_id,
        )
    )


    metadata_rows.append(
        {
            "pseudobulk_id": pb_id,
            "BrNum": donor,
            "Dx": dx,
            "Age": float(dm["Age"]),
            "Sex": str(dm["Sex"]),
            "slide_id": str(dm["slide_id"]),
            "run_date": str(dm["run_date"]),
            "predictions_smooth":
                LAYER_TO_SPD[layer],
            "layer_annotation":
                layer,
            "n_cells":
                layer_cells[layer],
            "measured300_library":
                float(mv.sum()),
            "full34987_library":
                float(fv.sum()),
        }
    )


measured_df = pd.DataFrame(
    measured_rows
)

full_df = pd.DataFrame(
    full_rows
)

pb_meta = pd.DataFrame(
    metadata_rows
)

gene_info = pd.DataFrame(
    {
        "gene": genes,
        "expression_source": source,
    }
)

excluded = pd.DataFrame(
    excluded_rows,
    columns=[
        "BrNum",
        "cell_id",
        "reason",
    ]
)

qc = pd.DataFrame(
    [
        {
            "BrNum": donor,
            "Dx": dx,
            "target_cells":
                len(target_ids),
            "layer_annotated_cells":
                len(ann_ids),
            "excluded_unannotated_cells":
                len(target_only),
            "measured_genes":
                n_measured,
            "imputed_genes":
                n_imputed,
            "full_genes":
                after.n_vars,
            "measured_layer_sum_matches":
                True,
            "full_layer_sum_matches":
                True,
        }
    ]
)


if measured_df.shape != (7, 300):
    raise RuntimeError(
        f"Measured shape = {measured_df.shape}"
    )

if full_df.shape != (7, 34987):
    raise RuntimeError(
        f"Full shape = {full_df.shape}"
    )


# =============================================================================
# Save donor outputs
# =============================================================================

prefix = OUT / donor

measured_df.to_csv(
    f"{prefix}_ENVI_measured300_layer_pb.csv.gz",
    compression="gzip",
)

full_df.to_csv(
    f"{prefix}_ENVI_full34987_layer_pb.csv.gz",
    compression="gzip",
)

pb_meta.to_csv(
    f"{prefix}_layer_metadata.csv",
    index=False,
)

gene_info.to_csv(
    f"{prefix}_gene_info.csv.gz",
    index=False,
    compression="gzip",
)

qc.to_csv(
    f"{prefix}_qc.csv",
    index=False,
)

excluded.to_csv(
    f"{prefix}_excluded_cells.csv",
    index=False,
)


target.file.close()
after.file.close()


print("\n" + "=" * 100)
print("FINAL SUMMARY")
print("=" * 100)

print("Donor:", donor)
print("Measured matrix:", measured_df.shape)
print("Full matrix:", full_df.shape)
print("Excluded cells:", len(excluded))
print("Measured layer sum validation: PASS")
print("Full layer sum validation: PASS")

print(
    f"\nFINAL STATUS: PASS — {donor}"
)
