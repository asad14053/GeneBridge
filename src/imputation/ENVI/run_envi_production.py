#!/usr/bin/env python3

import gc
import json
import time
import argparse
import inspect
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# ============================================================
# Helpers
# ============================================================

def header(text):
    print("\n" + "=" * 100, flush=True)
    print(text, flush=True)
    print("=" * 100, flush=True)


def dense_float32(x):
    if sparse.issparse(x):
        x = x.toarray()

    return np.asarray(
        x,
        dtype=np.float32,
    )


def row_sums(x):
    if sparse.issparse(x):
        return np.asarray(
            x.sum(axis=1)
        ).ravel()

    return np.asarray(
        x
    ).sum(
        axis=1,
        dtype=np.float64,
    )


def gene_sums(x):
    if sparse.issparse(x):
        return np.asarray(
            x.sum(axis=0)
        ).ravel()

    return np.asarray(
        x
    ).sum(
        axis=0,
        dtype=np.float64,
    )


def assert_finite_nonnegative(
    x,
    name,
    atol=1e-6,
):
    x = np.asarray(x)

    if not np.isfinite(x).all():
        raise FloatingPointError(
            f"{name} contains NaN/Inf."
        )

    minimum = float(
        np.min(x)
    )

    if minimum < -atol:
        raise FloatingPointError(
            f"{name} contains negative values. "
            f"Minimum={minimum}"
        )


def load_reference_counts(path):
    """
    Load ONLY reference .X raw counts + obs/var.

    Do NOT load layers["logcounts"] into memory.
    """

    header("LOAD HUUKI REFERENCE")

    start = time.time()

    backed = ad.read_h5ad(
        path,
        backed="r",
    )

    print(
        "Reference backed shape:",
        backed.shape,
        flush=True,
    )

    print(
        "Reference layers:",
        list(backed.layers.keys()),
        flush=True,
    )

    obs = backed.obs.copy()
    var = backed.var.copy()

    var_names = (
        backed.var_names
        .astype(str)
        .copy()
    )

    # .X has already been audited:
    # raw integer-like Huuki counts.
    X = np.asarray(
        backed.X[:],
        dtype=np.float32,
    )

    backed.file.close()

    ref = ad.AnnData(
        X=X,
        obs=obs,
        var=var,
    )

    ref.var_names = var_names

    print(
        "Loaded reference:",
        ref.shape,
        flush=True,
    )

    print(
        "Reference load seconds:",
        f"{time.time() - start:.2f}",
        flush=True,
    )

    return ref


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "ENVI production full-transcriptome "
        "imputation for one Xenium donor."
    )
)

parser.add_argument(
    "--target",
    type=Path,
    required=True,
)

parser.add_argument(
    "--reference",
    type=Path,
    required=True,
)

parser.add_argument(
    "--output-dir",
    type=Path,
    required=True,
)

parser.add_argument(
    "--donor",
    required=True,
)

parser.add_argument(
    "--experiment",
    required=True,
)

parser.add_argument(
    "--training-steps",
    type=int,
    default=16000,
)

parser.add_argument(
    "--batch-size",
    type=int,
    default=128,
)

parser.add_argument(
    "--num-hvg",
    type=int,
    default=2048,
)

parser.add_argument(
    "--num-cov-genes",
    type=int,
    default=64,
)

parser.add_argument(
    "--k-nearest",
    type=int,
    default=8,
)

parser.add_argument(
    "--spatial-dist",
    default="pois",
)

parser.add_argument(
    "--sc-dist",
    default="nb",
)

parser.add_argument(
    "--stable-eps",
    type=float,
    default=1e-6,
)

parser.add_argument(
    "--seed",
    type=int,
    default=0,
)

args = parser.parse_args()


args.output_dir.mkdir(
    parents=True,
    exist_ok=True,
)

diagnostics_dir = (
    args.output_dir
    / "diagnostics"
)

diagnostics_dir.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# ENVI / JAX
# ============================================================

header("IMPORT ENVI")

import jax
import scenvi

print(
    "ENVI location:",
    inspect.getfile(
        scenvi.ENVI
    ),
)

print(
    "JAX devices:",
    jax.devices(),
)


# ============================================================
# Basic validation
# ============================================================

if not args.target.exists():
    raise FileNotFoundError(
        args.target
    )

if not args.reference.exists():
    raise FileNotFoundError(
        args.reference
    )


header("PRODUCTION ENVI CONFIGURATION")

print("Experiment       :", args.experiment)
print("Donor            :", args.donor)
print("Target           :", args.target)
print("Reference        :", args.reference)
print("Output           :", args.output_dir)
print("Training steps   :", args.training_steps)
print("Batch size       :", args.batch_size)
print("num_HVG          :", args.num_hvg)
print("num_cov_genes    :", args.num_cov_genes)
print("k_nearest        :", args.k_nearest)
print("spatial_dist     :", args.spatial_dist)
print("sc_dist          :", args.sc_dist)
print("stable_eps       :", args.stable_eps)
print("Seed             :", args.seed)


# ============================================================
# Load Xenium target
# ============================================================

header("LOAD XENIUM TARGET")

target = ad.read_h5ad(
    args.target
)

print(
    "Target shape:",
    target.shape,
)

print(
    "Layers:",
    list(target.layers.keys()),
)

print(
    "obsm:",
    list(target.obsm.keys()),
)


if target.n_vars != 300:
    raise RuntimeError(
        f"{args.donor}: expected 300 Xenium genes, "
        f"found {target.n_vars}"
    )


if "spatial" not in target.obsm:
    raise KeyError(
        "Target missing obsm['spatial']."
    )


if "counts" in target.layers:

    spatial_counts = (
        dense_float32(
            target.layers[
                "counts"
            ]
        )
    )

else:

    spatial_counts = (
        dense_float32(
            target.X
        )
    )


assert_finite_nonnegative(
    spatial_counts,
    "Xenium raw counts",
)


coordinates = np.asarray(
    target.obsm[
        "spatial"
    ],
    dtype=np.float32,
)


if coordinates.shape != (
    target.n_obs,
    2,
):
    raise RuntimeError(
        "Unexpected spatial coordinate shape: "
        f"{coordinates.shape}"
    )


if not np.isfinite(
    coordinates
).all():
    raise FloatingPointError(
        "Spatial coordinates contain NaN/Inf."
    )


measured_genes = (
    target.var_names
    .astype(str)
    .tolist()
)

measured_gene_set = set(
    measured_genes
)


observed_library = (
    spatial_counts.sum(
        axis=1,
        dtype=np.float64,
    )
)


bad_spatial_cells = (
    ~np.isfinite(
        observed_library
    )
    | (
        observed_library <= 0
    )
)


if np.any(
    bad_spatial_cells
):
    raise RuntimeError(
        f"{args.donor}: "
        f"{int(bad_spatial_cells.sum())} cells "
        "have zero/invalid library across "
        "the 300 measured Xenium genes."
    )


spatial_gene_library = (
    spatial_counts.sum(
        axis=0,
        dtype=np.float64,
    )
)


zero_spatial_genes = (
    np.asarray(
        measured_genes
    )[
        spatial_gene_library
        <= 0
    ]
    .tolist()
)


if zero_spatial_genes:
    raise RuntimeError(
        "All-zero measured Xenium genes: "
        f"{zero_spatial_genes[:30]}"
    )


print(
    "Cells:",
    f"{target.n_obs:,}",
)

print(
    "Measured genes:",
    len(
        measured_genes
    ),
)

print(
    "Observed 300-gene library range:",
    float(
        observed_library.min()
    ),
    "–",
    float(
        observed_library.max()
    ),
)


# ============================================================
# Load all-10 Huuki reference
# ============================================================

reference = load_reference_counts(
    args.reference
)


if not reference.var_names.is_unique:
    raise RuntimeError(
        "Huuki reference gene names are not unique."
    )


assert_finite_nonnegative(
    reference.X,
    "Huuki reference raw counts",
)


# ============================================================
# Confirm 300/300 overlap
# ============================================================

header("GENE OVERLAP AUDIT")

missing_measured = (
    pd.Index(
        measured_genes
    )
    .difference(
        reference.var_names
    )
)


if len(
    missing_measured
):
    raise RuntimeError(
        "Measured Xenium genes missing "
        "from Huuki reference: "
        f"{missing_measured.tolist()}"
    )


print(
    "Xenium genes:",
    len(
        measured_genes
    ),
)

print(
    "Xenium genes in Huuki:",
    len(
        measured_genes
    ),
)

print(
    "Missing:",
    0,
)


# ============================================================
# Remove all-zero reference genes
# ============================================================

header("REFERENCE GENE QC")

ref_gene_sums = gene_sums(
    reference.X
)


informative_genes = (
    np.isfinite(
        ref_gene_sums
    )
    & (
        ref_gene_sums > 0
    )
)


n_zero_reference_genes = int(
    (~informative_genes).sum()
)


if n_zero_reference_genes:

    zero_names = (
        reference.var_names[
            ~informative_genes
        ]
        .astype(str)
        .tolist()
    )

    pd.DataFrame(
        {
            "gene":
                zero_names
        }
    ).to_csv(
        diagnostics_dir
        / "all_zero_reference_genes.csv",
        index=False,
    )

    print(
        "Removing all-zero reference genes:",
        n_zero_reference_genes,
    )

    reference = (
        reference[
            :,
            informative_genes,
        ]
        .copy()
    )

    gc.collect()


# Make sure no measured Xenium gene disappeared.
missing_after_qc = (
    pd.Index(
        measured_genes
    )
    .difference(
        reference.var_names
    )
)


if len(
    missing_after_qc
):
    raise RuntimeError(
        "Measured Xenium genes are all-zero "
        "in Huuki reference: "
        f"{missing_after_qc.tolist()}"
    )


# ============================================================
# Remove zero-library Huuki cells
# ============================================================

ref_library = row_sums(
    reference.X
)


keep_ref_cells = (
    np.isfinite(
        ref_library
    )
    & (
        ref_library > 0
    )
)


removed_reference_cells = int(
    (~keep_ref_cells).sum()
)


if removed_reference_cells:

    reference = (
        reference[
            keep_ref_cells,
            :
        ]
        .copy()
    )

    gc.collect()


reference_genes = (
    reference.var_names
    .astype(str)
    .tolist()
)


imputation_genes = [
    gene
    for gene
    in reference_genes
    if gene
    not in measured_gene_set
]


print(
    "Reference cells retained:",
    f"{reference.n_obs:,}",
)

print(
    "Reference genes retained:",
    f"{reference.n_vars:,}",
)

print(
    "Reference-only genes:",
    f"{len(imputation_genes):,}",
)

print(
    "Reference cells removed:",
    removed_reference_cells,
)

print(
    "Reference genes removed:",
    n_zero_reference_genes,
)


# ============================================================
# Required metadata
# ============================================================

header("PREPARE ENVI METADATA")

# Spatial batch = target donor.
target_obs = (
    target.obs.copy()
)

target_obs[
    "batch"
] = str(
    args.donor
)


if "cell_type" not in target_obs.columns:
    target_obs[
        "cell_type"
    ] = "unknown"


# Huuki batch by donor.
if "BrNum" in reference.obs.columns:

    reference.obs[
        "batch"
    ] = (
        reference.obs[
            "BrNum"
        ]
        .astype(str)
    )

elif "Sample" in reference.obs.columns:

    reference.obs[
        "batch"
    ] = (
        reference.obs[
            "Sample"
        ]
        .astype(str)
    )

else:

    reference.obs[
        "batch"
    ] = "Huuki"


if "cell_type" not in reference.obs.columns:

    if (
        "cellType_broad_hc"
        in reference.obs.columns
    ):

        reference.obs[
            "cell_type"
        ] = (
            reference.obs[
                "cellType_broad_hc"
            ]
            .astype(str)
        )

    else:

        reference.obs[
            "cell_type"
        ] = "unknown"


# ============================================================
# Spatial AnnData used by ENVI
# ============================================================

spatial_model = ad.AnnData(
    X=(
        spatial_counts.copy()
    ),
    obs=target_obs,
    var=pd.DataFrame(
        index=pd.Index(
            measured_genes,
            name="gene",
        )
    ),
)


spatial_model.obsm[
    "spatial"
] = (
    coordinates.copy()
)


# Explicit float32.
spatial_model.X = (
    dense_float32(
        spatial_model.X
    )
)

reference.X = (
    dense_float32(
        reference.X
    )
)


print(
    "Spatial ENVI input:",
    spatial_model.shape,
)

print(
    "Huuki ENVI input:",
    reference.shape,
)


# ============================================================
# Gene manifest
# ============================================================

pd.DataFrame(
    {
        "gene":
            reference_genes,

        "expression_source":
            [
                (
                    "measured_xenium"
                    if gene
                    in measured_gene_set
                    else
                    "envi_imputed"
                )
                for gene
                in reference_genes
            ],
    }
).to_csv(
    diagnostics_dir
    / "production_gene_manifest.csv",
    index=False,
)


# ============================================================
# ENVI constructor
# ============================================================

header("INITIALIZE ENVI")

requested_kwargs = {
    "spatial_data":
        spatial_model,

    "sc_data":
        reference,

    "spatial_key":
        "spatial",

    "batch_key":
        "batch",

    "k_nearest":
        int(
            args.k_nearest
        ),

    "num_cov_genes":
        int(
            args.num_cov_genes
        ),

    "num_HVG":
        int(
            args.num_hvg
        ),

    # --------------------------------------------------------
    # PRODUCTION:
    # Protect every reference-only gene so ENVI retains
    # the full transcriptome for imputation.
    # --------------------------------------------------------
    "sc_genes":
        imputation_genes,

    "spatial_dist":
        args.spatial_dist,

    "sc_dist":
        args.sc_dist,

    "stable_eps":
        float(
            args.stable_eps
        ),
}


constructor_signature = (
    inspect.signature(
        scenvi.ENVI
    )
)


accepted_args = set(
    constructor_signature.parameters
)


required_args = {
    "spatial_data",
    "sc_data",
    "spatial_key",
    "batch_key",
    "k_nearest",
    "num_cov_genes",
    "num_HVG",
    "sc_genes",
    "spatial_dist",
    "sc_dist",
}


missing_constructor_args = (
    required_args
    - accepted_args
)


if missing_constructor_args:
    raise RuntimeError(
        "Installed ENVI incompatible. "
        f"Missing constructor args: "
        f"{sorted(missing_constructor_args)}"
    )


envi_kwargs = {
    key: value
    for key, value
    in requested_kwargs.items()
    if key
    in accepted_args
}


print(
    "Measured spatial genes:",
    len(
        measured_genes
    ),
)

print(
    "Protected reference-only genes:",
    f"{len(imputation_genes):,}",
)


init_start = time.time()


envi_model = scenvi.ENVI(
    **envi_kwargs
)


initialization_seconds = (
    time.time()
    - init_start
)


print(
    "Initialization seconds:",
    f"{initialization_seconds:.2f}",
)


# ============================================================
# Critical post-init audit
# ============================================================

header("POST-INITIALIZATION AUDIT")

retained_spatial_genes = (
    envi_model
    .spatial_data
    .var_names
    .astype(str)
    .tolist()
)


retained_reference_genes = (
    envi_model
    .sc_data
    .var_names
    .astype(str)
    .tolist()
)


if set(
    retained_spatial_genes
) != set(
    measured_genes
):
    raise RuntimeError(
        "ENVI changed the 300 spatial gene set."
    )


missing_imputation_after_init = (
    pd.Index(
        imputation_genes
    )
    .difference(
        retained_reference_genes
    )
)


if len(
    missing_imputation_after_init
):
    raise RuntimeError(
        "ENVI dropped production genes. "
        f"Missing: "
        f"{len(missing_imputation_after_init):,}. "
        f"First: "
        f"{missing_imputation_after_init[:20].tolist()}"
    )


print(
    "ENVI retained spatial genes:",
    len(
        retained_spatial_genes
    ),
)

print(
    "ENVI retained scRNA genes:",
    f"{len(retained_reference_genes):,}",
)

print(
    "Protected imputation genes retained:",
    f"{len(imputation_genes):,}",
)


output_genes = (
    retained_reference_genes
)


# ============================================================
# Seed
# ============================================================

try:

    jax_key = (
        jax.random.key(
            args.seed
        )
    )

except AttributeError:

    jax_key = (
        jax.random.PRNGKey(
            args.seed
        )
    )


# ============================================================
# Train
# ============================================================

header("TRAIN ENVI")

train_signature = (
    inspect.signature(
        envi_model.train
    )
)


train_kwargs = {
    "training_steps":
        int(
            args.training_steps
        ),

    "batch_size":
        int(
            args.batch_size
        ),
}


if (
    "key"
    in train_signature.parameters
):

    train_kwargs[
        "key"
    ] = jax_key


training_start = time.time()


envi_model.train(
    **train_kwargs
)


training_seconds = (
    time.time()
    - training_start
)


print(
    "Training seconds:",
    f"{training_seconds:.2f}",
)


# ============================================================
# Latent QC
# ============================================================

if (
    "envi_latent"
    not in
    envi_model
    .spatial_data
    .obsm
):
    raise RuntimeError(
        "ENVI did not create "
        "spatial_data.obsm['envi_latent']."
    )


latent = np.asarray(
    envi_model
    .spatial_data
    .obsm[
        "envi_latent"
    ],
    dtype=np.float32,
)


if not np.isfinite(
    latent
).all():
    raise FloatingPointError(
        "ENVI latent contains NaN/Inf."
    )


print(
    "Spatial latent shape:",
    latent.shape,
)


# ============================================================
# Impute full transcriptome
# ============================================================

header("IMPUTE FULL TRANSCRIPTOME")

prediction_start = time.time()


envi_model.impute_genes()


prediction_seconds = (
    time.time()
    - prediction_start
)


print(
    "Prediction seconds:",
    f"{prediction_seconds:.2f}",
)


if (
    "imputation"
    not in
    envi_model
    .spatial_data
    .obsm
):
    raise RuntimeError(
        "ENVI did not create "
        "spatial_data.obsm['imputation']."
    )


imputation = (
    envi_model
    .spatial_data
    .obsm[
        "imputation"
    ]
)


if isinstance(
    imputation,
    pd.DataFrame,
):

    imputation_df = (
        imputation
    )

else:

    imputation_df = pd.DataFrame(
        np.asarray(
            imputation
        ),
        index=(
            envi_model
            .spatial_data
            .obs_names
            .astype(str)
        ),
        columns=(
            envi_model
            .sc_data
            .var_names
            .astype(str)
        ),
    )


imputation_df.index = (
    imputation_df
    .index
    .astype(str)
)

imputation_df.columns = (
    imputation_df
    .columns
    .astype(str)
)


# ============================================================
# Cell alignment
# ============================================================

target_cells = (
    target.obs_names
    .astype(str)
)


if not target_cells.equals(
    pd.Index(
        imputation_df.index
    )
):

    missing_cells = (
        target_cells
        .difference(
            imputation_df.index
        )
    )

    if len(
        missing_cells
    ):
        raise RuntimeError(
            "ENVI output missing "
            f"{len(missing_cells)} Xenium cells."
        )

    imputation_df = (
        imputation_df.loc[
            target_cells
        ]
    )


# ============================================================
# Gene alignment
# ============================================================

missing_output_genes = (
    pd.Index(
        output_genes
    )
    .difference(
        imputation_df.columns
    )
)


if len(
    missing_output_genes
):
    raise RuntimeError(
        "ENVI output missing "
        f"{len(missing_output_genes)} genes."
    )


imputation_df = (
    imputation_df.loc[
        :,
        output_genes,
    ]
)


# ============================================================
# Native ENVI output
# ============================================================

header("NATIVE ENVI OUTPUT")

native = (
    imputation_df
    .to_numpy(
        dtype=np.float32,
        copy=False,
    )
)


assert_finite_nonnegative(
    native,
    "ENVI native full-transcriptome output",
)


np.maximum(
    native,
    0,
    out=native,
)


# ------------------------------------------------------------
# LOCKED OUTPUT:
# preserve original ENVI prediction BEFORE calibration.
# ------------------------------------------------------------

native_unscaled = (
    native.copy()
)


print(
    "Native output shape:",
    native_unscaled.shape,
)


# ============================================================
# Locate the measured 300 genes in ENVI output
# ============================================================

gene_to_index = {
    gene: index
    for index, gene
    in enumerate(
        output_genes
    )
}


missing_measured_output = [
    gene
    for gene
    in measured_genes
    if gene
    not in gene_to_index
]


if missing_measured_output:
    raise RuntimeError(
        "ENVI full output missing measured genes: "
        f"{missing_measured_output}"
    )


measured_indices = np.asarray(
    [
        gene_to_index[
            gene
        ]
        for gene
        in measured_genes
    ],
    dtype=np.int64,
)


# ============================================================
# Count-scale calibration
# ============================================================

header("COUNT-SCALE CALIBRATION")

native_measured = (
    native[
        :,
        measured_indices,
    ]
)


native_measured_library = (
    native_measured.sum(
        axis=1,
        dtype=np.float64,
    )
)


invalid_native_library = (
    ~np.isfinite(
        native_measured_library
    )
    | (
        native_measured_library
        <= 1e-12
    )
)


if np.any(
    invalid_native_library
):

    pd.DataFrame(
        {
            "cell_id":
                target_cells[
                    invalid_native_library
                ],

            "observed_xenium_300_library":
                observed_library[
                    invalid_native_library
                ],

            "native_envi_300_library":
                native_measured_library[
                    invalid_native_library
                ],
        }
    ).to_csv(
        diagnostics_dir
        / "invalid_count_scale_cells.csv",
        index=False,
    )

    raise FloatingPointError(
        "ENVI native prediction has zero/invalid "
        "300-gene library for "
        f"{int(invalid_native_library.sum())} cells."
    )


scale_factor = (
    observed_library
    / native_measured_library
)


if (
    not np.isfinite(
        scale_factor
    ).all()
    or np.any(
        scale_factor <= 0
    )
):
    raise FloatingPointError(
        "Count-scale factors contain invalid values."
    )


quantile_probabilities = (
    np.asarray(
        [
            0.00,
            0.01,
            0.25,
            0.50,
            0.75,
            0.99,
            1.00,
        ]
    )
)


scale_quantiles = (
    np.quantile(
        scale_factor,
        quantile_probabilities,
    )
)


pd.DataFrame(
    {
        "quantile":
            quantile_probabilities,

        "scale_factor":
            scale_quantiles,
    }
).to_csv(
    diagnostics_dir
    / "count_scale_factor_quantiles.csv",
    index=False,
)


print(
    "Scale-factor median:",
    float(
        np.median(
            scale_factor
        )
    ),
)

print(
    "Scale-factor range:",
    float(
        scale_factor.min()
    ),
    "–",
    float(
        scale_factor.max()
    ),
)


# ============================================================
# Calibrate full transcriptome IN PLACE
# ============================================================

native *= (
    scale_factor
    .astype(
        np.float32
    )[
        :,
        None,
    ]
)


assert_finite_nonnegative(
    native,
    "ENVI calibrated full transcriptome",
)


# ============================================================
# Restore measured Xenium counts
# ============================================================

header("RESTORE MEASURED XENIUM COUNTS")

native[
    :,
    measured_indices,
] = spatial_counts


if not np.array_equal(
    native[
        :,
        measured_indices,
    ],
    spatial_counts,
):
    raise RuntimeError(
        "Measured Xenium counts were not "
        "restored exactly."
    )


# This is the final canonical count-scale matrix.
count_scale = native


print(
    "Measured Xenium genes restored:",
    len(
        measured_genes
    ),
)

print(
    "Imputed genes:",
    (
        len(
            output_genes
        )
        - len(
            measured_genes
        )
    ),
)

print(
    "Final count-scale shape:",
    count_scale.shape,
)


# ============================================================
# log1p layer
# ============================================================

header("CREATE log1p LAYER")

log1p_matrix = (
    np.log1p(
        count_scale
    )
    .astype(
        np.float32,
        copy=False,
    )
)


if not np.isfinite(
    log1p_matrix
).all():
    raise FloatingPointError(
        "log1p matrix contains NaN/Inf."
    )


print(
    "log1p shape:",
    log1p_matrix.shape,
)


# ============================================================
# Final metadata
# ============================================================

final_obs = (
    target.obs.copy()
)


final_obs[
    "envi_count_scale_factor"
] = (
    scale_factor.astype(
        np.float32
    )
)


reference_var = (
    reference.var.copy()
)


reference_var.index = (
    reference.var_names
    .astype(str)
)


final_var = (
    reference_var
    .reindex(
        output_genes
    )
    .copy()
)


final_var[
    "expression_source"
] = [
    (
        "measured_xenium"
        if gene
        in measured_gene_set
        else
        "envi_imputed"
    )
    for gene
    in output_genes
]


final_var[
    "measured_in_xenium"
] = [
    gene
    in measured_gene_set
    for gene
    in output_genes
]


# ============================================================
# Manifest
# ============================================================

manifest = {
    "experiment":
        args.experiment,

    "donor":
        args.donor,

    "model":
        "ENVI",

    "reference":
        str(
            args.reference.resolve()
        ),

    "target":
        str(
            args.target.resolve()
        ),

    "reference_cells":
        int(
            reference.n_obs
        ),

    "output_cells":
        int(
            target.n_obs
        ),

    "output_genes":
        int(
            len(
                output_genes
            )
        ),

    "measured_xenium_genes":
        int(
            len(
                measured_genes
            )
        ),

    "envi_imputed_genes":
        int(
            len(
                output_genes
            )
            - len(
                measured_genes
            )
        ),

    "reference_all_zero_genes_removed":
        int(
            n_zero_reference_genes
        ),

    "reference_zero_cells_removed":
        int(
            removed_reference_cells
        ),

    "training_steps":
        int(
            args.training_steps
        ),

    "batch_size":
        int(
            args.batch_size
        ),

    "num_HVG":
        int(
            args.num_hvg
        ),

    "num_cov_genes":
        int(
            args.num_cov_genes
        ),

    "k_nearest":
        int(
            args.k_nearest
        ),

    "spatial_dist":
        args.spatial_dist,

    "sc_dist":
        args.sc_dist,

    "stable_eps":
        float(
            args.stable_eps
        ),

    "seed":
        int(
            args.seed
        ),

    "initialization_seconds":
        float(
            initialization_seconds
        ),

    "training_seconds":
        float(
            training_seconds
        ),

    "prediction_seconds":
        float(
            prediction_seconds
        ),

    "count_scale_calibration":
        (
            "Per-cell actual Xenium library across "
            "all 300 measured genes divided by ENVI "
            "native prediction across the same "
            "300 measured genes."
        ),

    "measured_gene_policy":
        (
            "After count-scale calibration, the "
            "300 measured Xenium genes are replaced "
            "with their original raw Xenium counts."
        ),

    "matrix_definitions": {
        "X":
            (
                "Final count-scale expression. "
                "Measured genes are original Xenium "
                "counts; non-measured genes are "
                "calibrated ENVI predictions."
            ),

        "layers/count_scale":
            (
                "Same final count-scale representation "
                "as X."
            ),

        "layers/log1p":
            (
                "numpy.log1p(count_scale)."
            ),

        "layers/native":
            (
                "Original ENVI decoder output before "
                "per-cell count-scale calibration."
            ),
    },
}


# ============================================================
# Build final AnnData
# ============================================================

header("BUILD FINAL H5AD")

final = ad.AnnData(
    X=count_scale,
    obs=final_obs,
    var=final_var,
)


# ------------------------------------------------------------
# LOCKED full production output format
# ------------------------------------------------------------

final.layers[
    "count_scale"
] = count_scale


final.layers[
    "log1p"
] = log1p_matrix


final.layers[
    "native"
] = native_unscaled


final.obsm[
    "spatial"
] = (
    coordinates.copy()
)


final.uns[
    "production_imputation"
] = manifest


# ============================================================
# Final pre-write audit
# ============================================================

header("PRE-WRITE QC")

print(
    "Shape:",
    final.shape,
)

print(
    "Layers:",
    list(
        final.layers.keys()
    ),
)

print(
    "obsm:",
    list(
        final.obsm.keys()
    ),
)

print(
    "Measured genes:",
    int(
        final.var[
            "measured_in_xenium"
        ].sum()
    ),
)

print(
    "Imputed genes:",
    int(
        (
            ~final.var[
                "measured_in_xenium"
            ]
        ).sum()
    ),
)


if (
    "count_scale"
    not in final.layers
    or
    "log1p"
    not in final.layers
    or
    "native"
    not in final.layers
):
    raise RuntimeError(
        "Required output layers are missing."
    )


if not np.array_equal(
    final.X[
        :,
        measured_indices,
    ],
    spatial_counts,
):
    raise RuntimeError(
        ".X measured Xenium values changed."
    )


if not np.array_equal(
    final.layers[
        "count_scale"
    ][
        :,
        measured_indices,
    ],
    spatial_counts,
):
    raise RuntimeError(
        "count_scale measured Xenium values changed."
    )


# ============================================================
# Write
# ============================================================

output_path = (
    args.output_dir
    / (
        f"spatial_data_xenium_"
        f"{args.donor}_ENVI_"
        f"full_transcriptome.h5ad"
    )
)


header("WRITE FULL PRODUCTION H5AD")

print(
    "Output:",
    output_path,
    flush=True,
)


write_start = time.time()


final.write_h5ad(
    output_path,
    compression="lzf",
)


write_seconds = (
    time.time()
    - write_start
)


# ============================================================
# Save manifest separately
# ============================================================

manifest[
    "output"
] = str(
    output_path.resolve()
)

manifest[
    "write_seconds"
] = float(
    write_seconds
)


with open(
    args.output_dir
    / "run_manifest.json",
    "w",
) as handle:

    json.dump(
        manifest,
        handle,
        indent=2,
    )


# ============================================================
# Completion flag
# ============================================================

with open(
    args.output_dir
    / "complete.flag",
    "w",
) as handle:

    handle.write(
        f"{args.donor}\tSUCCESS\n"
    )


# ============================================================
# Final summary
# ============================================================

header("SUCCESS")

print(
    "Experiment:",
    args.experiment,
)

print(
    "Donor:",
    args.donor,
)

print(
    "Cells:",
    f"{final.n_obs:,}",
)

print(
    "Total genes:",
    f"{final.n_vars:,}",
)

print(
    "Measured Xenium genes:",
    f"{len(measured_genes):,}",
)

print(
    "ENVI-imputed genes:",
    f"{final.n_vars - len(measured_genes):,}",
)

print(
    "X:",
    "final count-scale expression",
)

print(
    'layers["count_scale"]:',
    "final count-scale expression",
)

print(
    'layers["log1p"]:',
    "log1p(count_scale)",
)

print(
    'layers["native"]:',
    "original ENVI decoder output",
)

print(
    'obsm["spatial"]:',
    final.obsm[
        "spatial"
    ].shape,
)

print(
    "Write seconds:",
    f"{write_seconds:.2f}",
)

print(
    "Output:",
    output_path,
)

print(
    "DONE",
)
