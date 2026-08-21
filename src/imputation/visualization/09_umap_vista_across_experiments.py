#!/usr/bin/env python3

from pathlib import Path
import json
import gc

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import umap
import inspect
from sklearn.utils import check_array as _sklearn_check_array

# -------------------------------------------------------------------------
# Compatibility shim:
# Newer umap-learn uses `ensure_all_finite`, while older scikit-learn
# versions use `force_all_finite`.
# This patches only UMAP's local check_array reference and does not modify
# the installed environment.
# -------------------------------------------------------------------------

if (
    "ensure_all_finite"
    not in inspect.signature(
        _sklearn_check_array
    ).parameters
):

    _original_umap_check_array = (
        umap.umap_.check_array
    )

    def _umap_check_array_compat(
        *args,
        **kwargs,
    ):

        if (
            "ensure_all_finite"
            in kwargs
        ):

            kwargs[
                "force_all_finite"
            ] = kwargs.pop(
                "ensure_all_finite"
            )

        return (
            _original_umap_check_array(
                *args,
                **kwargs,
            )
        )

    umap.umap_.check_array = (
        _umap_check_array_compat
    )

    print(
        "UMAP/sklearn compatibility shim: ENABLED"
    )

else:

    print(
        "UMAP/sklearn compatibility shim: not needed"
    )

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =============================================================================
# PATHS
# =============================================================================

ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE = (
    ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

TARGET_ANNOTATED = (
    ROOT
    / "data"
    / "processed"
    / "imputation_beta"
    / "Br8667"
    / "spatial_data_xenium_Br8667_vista_qc.h5ad"
)

OUT_DIR = (
    BASE
    / "final_visualizations"
    / "cell_population"
    / "umap"
    / "vista"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


RUNS = [
    ("ex5", "Experiment 5"),
    ("ex5_1", "Experiment 5.1"),
    ("ex5_3", "Experiment 5.3"),
]

METHOD = "vista"

COUNT_LAYER = "count_scale"

CELLTYPE_COLUMN = "cell_type_annotation"

RANDOM_SEED = 8667

N_PCS = 30

N_NEIGHBORS = 30

MIN_DIST = 0.30


# =============================================================================
# DISCOVER H5AD FILES
# =============================================================================

def flatten_json(obj, prefix=""):

    items = []

    if isinstance(obj, dict):

        for key, value in obj.items():

            current = (
                f"{prefix}.{key}"
                if prefix
                else str(key)
            )

            items.extend(
                flatten_json(
                    value,
                    current,
                )
            )

    elif isinstance(obj, list):

        for i, value in enumerate(obj):

            items.extend(
                flatten_json(
                    value,
                    f"{prefix}[{i}]",
                )
            )

    else:

        items.append(
            (
                prefix,
                obj,
            )
        )

    return items


def resolve_manifest_h5ad(
    run_dir,
    kind,
):

    manifest = (
        run_dir
        / "combined_manifest.json"
    )

    if not manifest.is_file():
        return []


    with manifest.open() as handle:
        obj = json.load(handle)


    if kind == "pred":

        tokens = [
            "pred",
            "prediction",
            "imputed",
            "oof_pred",
        ]

    else:

        tokens = [
            "truth",
            "observed",
            "ground_truth",
            "target",
        ]


    hits = []


    for key, value in flatten_json(obj):

        key_lower = key.lower()

        value_string = str(value)

        if (
            value_string.lower().endswith(".h5ad")
            and
            any(
                token in key_lower
                for token in tokens
            )
        ):

            path = Path(value_string)

            candidates = [
                path,
                run_dir / path,
                ROOT / path,
            ]

            for candidate in candidates:

                if candidate.is_file():
                    hits.append(
                        candidate.resolve()
                    )
                    break


    return sorted(
        set(hits)
    )


def score_file(
    path,
    kind,
):

    text = str(
        path.name
    ).lower()


    if kind == "pred":

        positive = [
            "pred",
            "prediction",
            "imputed",
            "oof_pred",
        ]

        negative = [
            "truth",
            "observed",
            "ground_truth",
        ]

    else:

        positive = [
            "truth",
            "observed",
            "ground_truth",
        ]

        negative = [
            "pred",
            "prediction",
            "imputed",
        ]


    score = 0

    for token in positive:

        if token in text:
            score += 10


    for token in negative:

        if token in text:
            score -= 10


    if "300" in text:
        score += 1

    if "combined" in text:
        score += 1

    if "oof" in text:
        score += 1


    return score


def discover_h5ad(
    run_dir,
    kind,
):

    # -------------------------------------------------------------------------
    # First preference: explicit manifest paths
    # -------------------------------------------------------------------------

    manifest_hits = resolve_manifest_h5ad(
        run_dir,
        kind,
    )

    if len(manifest_hits) == 1:

        return manifest_hits[0]


    # -------------------------------------------------------------------------
    # Fall back to filename-based discovery
    # -------------------------------------------------------------------------

    files = sorted(
        set(
            list(
                run_dir.glob("*.h5ad")
            )
            +
            list(
                run_dir.rglob("*.h5ad")
            )
        )
    )


    if not files:

        raise FileNotFoundError(
            f"No H5AD files found in:\n{run_dir}"
        )


    scored = [
        (
            score_file(
                path,
                kind,
            ),
            path,
        )
        for path in files
    ]


    scored.sort(
        key=lambda x: x[0],
        reverse=True,
    )


    best_score = scored[0][0]

    best = [
        path
        for score, path in scored
        if (
            score == best_score
            and score > 0
        )
    ]


    if len(best) == 1:

        return best[0]


    print()
    print(
        f"Could not uniquely identify {kind} H5AD in:"
    )
    print(run_dir)

    print(
        "\nCandidates:"
    )

    for score, path in scored:

        print(
            f"  score={score:3d}  {path}"
        )


    raise RuntimeError(
        f"Ambiguous {kind} H5AD discovery."
    )


# =============================================================================
# LOAD MATRIX
# =============================================================================

def matrix_to_dense_float32(
    matrix,
):

    if sparse.issparse(matrix):

        matrix = matrix.toarray()

    else:

        matrix = np.asarray(
            matrix
        )


    return np.asarray(
        matrix,
        dtype=np.float32,
    )


def load_count_matrix(
    path,
):

    print()
    print(
        "Loading:",
        path
    )


    a = ad.read_h5ad(
        path
    )


    print(
        "  shape:",
        a.shape
    )

    print(
        "  layers:",
        list(
            a.layers.keys()
        )
    )


    if COUNT_LAYER in a.layers:

        print(
            f"  matrix source: layers[{COUNT_LAYER!r}]"
        )

        matrix = a.layers[
            COUNT_LAYER
        ]

    else:

        print(
            "  matrix source: X "
            "(count_scale layer absent)"
        )

        matrix = a.X


    X = matrix_to_dense_float32(
        matrix
    )


    if not np.isfinite(
        X
    ).all():

        raise ValueError(
            f"Non-finite values in {path}"
        )


    if X.min() < 0:

        raise ValueError(
            f"Negative count-scale values in {path}"
        )


    obs_names = (
        a.obs_names
        .astype(str)
        .to_numpy()
    )

    var_names = (
        a.var_names
        .astype(str)
        .to_numpy()
    )


    del a

    gc.collect()


    return (
        X,
        obs_names,
        var_names,
    )


# =============================================================================
# ANNOTATIONS
# =============================================================================

def load_cell_types(
    observed_cells,
):

    if not TARGET_ANNOTATED.is_file():

        raise FileNotFoundError(
            TARGET_ANNOTATED
        )


    target = ad.read_h5ad(
        TARGET_ANNOTATED,
        backed="r",
    )


    if CELLTYPE_COLUMN not in target.obs.columns:

        raise KeyError(
            f"{CELLTYPE_COLUMN!r} not found.\n"
            f"Available obs columns:\n"
            f"{list(target.obs.columns)}"
        )


    target_cells = (
        target.obs_names
        .astype(str)
    )


    missing = (
        pd.Index(
            observed_cells
        )
        .difference(
            target_cells
        )
    )


    if len(missing) > 0:

        raise ValueError(
            f"{len(missing)} observed cells are missing "
            "from annotated Xenium target."
        )


    annotations = (
        target.obs[
            CELLTYPE_COLUMN
        ]
        .astype(str)
        .reindex(
            observed_cells
        )
    )


    target.file.close()


    if annotations.isna().any():

        raise ValueError(
            "Missing cell-type annotations after alignment."
        )


    return annotations


# =============================================================================
# VALIDATION
# =============================================================================

def validate_alignment(
    reference_cells,
    reference_genes,
    cells,
    genes,
    label,
):

    if not np.array_equal(
        reference_cells,
        cells,
    ):

        raise ValueError(
            f"Cell ordering mismatch: {label}"
        )


    if not np.array_equal(
        reference_genes,
        genes,
    ):

        raise ValueError(
            f"Gene ordering mismatch: {label}"
        )


# =============================================================================
# PREPROCESSING
# =============================================================================

def log1p_counts(
    X,
):

    X = X.copy()

    np.log1p(
        X,
        out=X,
    )

    return X


# =============================================================================
# PLOT
# =============================================================================

def plot_umaps(
    coordinates,
    cell_types,
):

    panel_order = [
        (
            "Observed Xenium",
            coordinates[
                "observed"
            ],
        ),
        (
            "VISTA\nExperiment 5",
            coordinates[
                "ex5"
            ],
        ),
        (
            "VISTA\nExperiment 5.1",
            coordinates[
                "ex5_1"
            ],
        ),
        (
            "VISTA\nExperiment 5.3",
            coordinates[
                "ex5_3"
            ],
        ),
    ]


    categories = sorted(
        pd.unique(
            cell_types
        )
    )


    cmap = plt.get_cmap(
        "tab20",
        len(categories),
    )


    category_colors = {
        category:
            cmap(i)
        for i, category
        in enumerate(categories)
    }


    all_xy = np.vstack(
        [
            xy
            for _, xy
            in panel_order
        ]
    )


    x_min = float(
        np.min(
            all_xy[:, 0]
        )
    )

    x_max = float(
        np.max(
            all_xy[:, 0]
        )
    )

    y_min = float(
        np.min(
            all_xy[:, 1]
        )
    )

    y_max = float(
        np.max(
            all_xy[:, 1]
        )
    )


    x_pad = (
        0.03
        * (
            x_max
            - x_min
        )
    )

    y_pad = (
        0.03
        * (
            y_max
            - y_min
        )
    )


    fig, axes = plt.subplots(
        1,
        4,
        figsize=(
            24,
            6.5,
        ),
        squeeze=False,
    )


    axes = axes[0]


    for ax, (
        title,
        xy,
    ) in zip(
        axes,
        panel_order,
    ):

        for category in categories:

            mask = (
                cell_types.to_numpy()
                == category
            )


            ax.scatter(
                xy[
                    mask,
                    0,
                ],
                xy[
                    mask,
                    1,
                ],
                s=1.5,
                alpha=0.55,
                linewidths=0,
                rasterized=True,
                color=category_colors[
                    category
                ],
            )


        ax.set_title(
            title,
            fontsize=13,
            fontweight="bold",
        )


        ax.set_xlim(
            x_min - x_pad,
            x_max + x_pad,
        )

        ax.set_ylim(
            y_min - y_pad,
            y_max + y_pad,
        )


        ax.set_xlabel(
            "UMAP1"
        )

        ax.set_ylabel(
            "UMAP2"
        )


        ax.set_xticks([])
        ax.set_yticks([])

        ax.spines[
            "top"
        ].set_visible(False)

        ax.spines[
            "right"
        ].set_visible(False)

        ax.spines[
            "bottom"
        ].set_visible(False)

        ax.spines[
            "left"
        ].set_visible(False)


    handles = []

    labels = []


    for category in categories:

        handle = plt.Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=6,
            markerfacecolor=category_colors[
                category
            ],
            markeredgecolor="none",
        )

        handles.append(
            handle
        )

        labels.append(
            category
        )


    fig.legend(
        handles,
        labels,
        title="Xenium cell type",
        bbox_to_anchor=(
            1.005,
            0.5,
        ),
        loc="center left",
        frameon=False,
        fontsize=9,
    )


    fig.suptitle(
        (
            "VISTA cell-population preservation across "
            "reference strategies\n"
            "Observed-fitted preprocessing, PCA and UMAP; "
            "imputed data transformed into the same embedding"
        ),
        fontsize=16,
        fontweight="bold",
        y=1.02,
    )


    fig.tight_layout(
        rect=[
            0,
            0,
            0.88,
            0.94,
        ]
    )


    png = (
        OUT_DIR
        / "vista_shared_umap_celltype.png"
    )

    pdf = (
        OUT_DIR
        / "vista_shared_umap_celltype.pdf"
    )


    fig.savefig(
        png,
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


    return (
        png,
        pdf,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print(
        "=" * 100
    )

    print(
        "VISTA SHARED-BASIS UMAP"
    )

    print(
        "=" * 100
    )


    # -------------------------------------------------------------------------
    # Discover all H5ADs
    # -------------------------------------------------------------------------

    discovered = {}


    for experiment, label in RUNS:

        run_dir = (
            BASE
            / experiment
            / METHOD
            / "combined_v2"
        )


        if not run_dir.is_dir():

            raise FileNotFoundError(
                run_dir
            )


        truth_path = discover_h5ad(
            run_dir,
            "truth",
        )

        pred_path = discover_h5ad(
            run_dir,
            "pred",
        )


        discovered[
            experiment
        ] = {
            "truth":
                truth_path,

            "pred":
                pred_path,
        }


        print()
        print(
            f"{label}"
        )

        print(
            "  truth:",
            truth_path
        )

        print(
            "  pred :",
            pred_path
        )


    # -------------------------------------------------------------------------
    # Load observed truth ONCE
    # -------------------------------------------------------------------------

    observed_path = (
        discovered[
            "ex5"
        ][
            "truth"
        ]
    )


    (
        X_obs,
        obs_cells,
        obs_genes,
    ) = load_count_matrix(
        observed_path
    )


    print()
    print(
        f"Observed matrix: {X_obs.shape}"
    )


    if X_obs.shape != (
        66148,
        300,
    ):

        print(
            "WARNING: expected audit shape "
            "(66148, 300), found",
            X_obs.shape,
        )


    # -------------------------------------------------------------------------
    # Cell-type annotations
    # -------------------------------------------------------------------------

    cell_types = load_cell_types(
        obs_cells
    )


    print()
    print(
        "Cell-type counts:"
    )

    print(
        cell_types
        .value_counts()
        .to_string()
    )


    # -------------------------------------------------------------------------
    # Fit preprocessing on OBSERVED ONLY
    # -------------------------------------------------------------------------

    print()
    print(
        "log1p observed counts..."
    )


    X_obs_log = log1p_counts(
        X_obs
    )


    del X_obs

    gc.collect()


    print(
        "Fitting StandardScaler on observed..."
    )


    scaler = StandardScaler(
        copy=True,
    )


    X_obs_scaled = scaler.fit_transform(
        X_obs_log
    ).astype(
        np.float32,
        copy=False,
    )


    del X_obs_log

    gc.collect()


    print(
        f"Fitting PCA ({N_PCS} PCs) on observed..."
    )


    pca = PCA(
        n_components=N_PCS,
        svd_solver="randomized",
        random_state=RANDOM_SEED,
    )


    Z_obs = pca.fit_transform(
        X_obs_scaled
    ).astype(
        np.float32,
        copy=False,
    )


    del X_obs_scaled

    gc.collect()


    variance = float(
        pca.explained_variance_ratio_
        .sum()
    )


    print(
        f"PCA variance explained by "
        f"{N_PCS} PCs: {variance:.4f}"
    )


    # -------------------------------------------------------------------------
    # Fit UMAP on OBSERVED ONLY
    # -------------------------------------------------------------------------

    print()
    print(
        "Fitting UMAP on observed PCA coordinates..."
    )


    reducer = umap.UMAP(
        n_neighbors=N_NEIGHBORS,
        n_components=2,
        min_dist=MIN_DIST,
        metric="euclidean",
        random_state=RANDOM_SEED,
        transform_seed=RANDOM_SEED,
        n_jobs=1,
        verbose=True,
    )


    U_obs = reducer.fit_transform(
        Z_obs
    ).astype(
        np.float32,
        copy=False,
    )


    coordinates = {
        "observed":
            U_obs
    }


    # -------------------------------------------------------------------------
    # Transform each VISTA prediction into OBSERVED basis
    # -------------------------------------------------------------------------

    for experiment, label in RUNS:

        print()
        print(
            "=" * 100
        )

        print(
            f"Transforming {label}"
        )

        print(
            "=" * 100
        )


        (
            X_pred,
            pred_cells,
            pred_genes,
        ) = load_count_matrix(
            discovered[
                experiment
            ][
                "pred"
            ]
        )


        validate_alignment(
            obs_cells,
            obs_genes,
            pred_cells,
            pred_genes,
            label,
        )


        X_pred_log = log1p_counts(
            X_pred
        )


        del X_pred

        gc.collect()


        X_pred_scaled = scaler.transform(
            X_pred_log
        ).astype(
            np.float32,
            copy=False,
        )


        del X_pred_log

        gc.collect()


        Z_pred = pca.transform(
            X_pred_scaled
        ).astype(
            np.float32,
            copy=False,
        )


        del X_pred_scaled

        gc.collect()


        U_pred = reducer.transform(
            Z_pred
        ).astype(
            np.float32,
            copy=False,
        )


        coordinates[
            experiment
        ] = U_pred


        del Z_pred

        gc.collect()


    # -------------------------------------------------------------------------
    # Save coordinates
    # -------------------------------------------------------------------------

    coord_df = pd.DataFrame(
        {
            "cell":
                obs_cells,

            "cell_type_annotation":
                cell_types.to_numpy(),

            "observed_UMAP1":
                coordinates[
                    "observed"
                ][
                    :,
                    0,
                ],

            "observed_UMAP2":
                coordinates[
                    "observed"
                ][
                    :,
                    1,
                ],

            "ex5_UMAP1":
                coordinates[
                    "ex5"
                ][
                    :,
                    0,
                ],

            "ex5_UMAP2":
                coordinates[
                    "ex5"
                ][
                    :,
                    1,
                ],

            "ex5_1_UMAP1":
                coordinates[
                    "ex5_1"
                ][
                    :,
                    0,
                ],

            "ex5_1_UMAP2":
                coordinates[
                    "ex5_1"
                ][
                    :,
                    1,
                ],

            "ex5_3_UMAP1":
                coordinates[
                    "ex5_3"
                ][
                    :,
                    0,
                ],

            "ex5_3_UMAP2":
                coordinates[
                    "ex5_3"
                ][
                    :,
                    1,
                ],
        }
    )


    coord_file = (
        OUT_DIR
        / "vista_shared_umap_coordinates.csv"
    )


    coord_df.to_csv(
        coord_file,
        index=False,
    )


    # -------------------------------------------------------------------------
    # Descriptive displacement summary
    # -------------------------------------------------------------------------

    displacement_rows = []


    for experiment, label in RUNS:

        distance = np.sqrt(
            np.sum(
                (
                    coordinates[
                        experiment
                    ]
                    -
                    coordinates[
                        "observed"
                    ]
                )
                ** 2,
                axis=1,
            )
        )


        displacement_rows.append(
            {
                "experiment":
                    experiment,

                "experiment_label":
                    label,

                "n_cells":
                    len(distance),

                "median_umap_displacement":
                    float(
                        np.median(
                            distance
                        )
                    ),

                "q75_umap_displacement":
                    float(
                        np.quantile(
                            distance,
                            0.75,
                        )
                    ),

                "q90_umap_displacement":
                    float(
                        np.quantile(
                            distance,
                            0.90,
                        )
                    ),

                "q95_umap_displacement":
                    float(
                        np.quantile(
                            distance,
                            0.95,
                        )
                    ),
            }
        )


    displacement_df = pd.DataFrame(
        displacement_rows
    )


    displacement_df.to_csv(
        OUT_DIR
        / "vista_umap_displacement_summary.csv",
        index=False,
    )


    # -------------------------------------------------------------------------
    # Figure
    # -------------------------------------------------------------------------

    png, pdf = plot_umaps(
        coordinates,
        cell_types,
    )


    manifest = {
        "method":
            "VISTA",

        "observed_basis":
            str(
                observed_path
            ),

        "annotation_source":
            str(
                TARGET_ANNOTATED
            ),

        "annotation_column":
            CELLTYPE_COLUMN,

        "count_layer":
            COUNT_LAYER,

        "preprocessing": {
            "log_transform":
                "log1p(count_scale)",

            "scaler":
                "StandardScaler fit on observed only",

            "pca":
                (
                    f"{N_PCS} components; "
                    "fit on observed only"
                ),

            "umap":
                (
                    "fit on observed PCA only; "
                    "predictions transformed into "
                    "the observed embedding"
                ),

            "n_neighbors":
                N_NEIGHBORS,

            "min_dist":
                MIN_DIST,

            "random_seed":
                RANDOM_SEED,
        },

        "figure_png":
            str(
                png
            ),

        "figure_pdf":
            str(
                pdf
            ),

        "coordinates_csv":
            str(
                coord_file
            ),
    }


    with (
        OUT_DIR
        / "manifest.json"
    ).open(
        "w"
    ) as handle:

        json.dump(
            manifest,
            handle,
            indent=2,
        )


    (
        OUT_DIR
        / "plot_complete.flag"
    ).write_text(
        "PASS\n"
    )


    print()
    print(
        "=" * 100
    )

    print(
        "UMAP COMPLETE"
    )

    print(
        "=" * 100
    )

    print(
        displacement_df.to_string(
            index=False
        )
    )

    print()
    print(
        "Figure:",
        png
    )


if __name__ == "__main__":
    main()
