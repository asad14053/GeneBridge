#!/usr/bin/env python3

from pathlib import Path
import json
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

OUT_DIR = (
    BASE_DIR
    / "final_visualizations"
    / "00_audit"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


EXPERIMENTS = {
    "ex5": "Experiment 5",
    "ex5_1": "Experiment 5.1",
    "ex5_3": "Experiment 5.3",
}


MODELS = {
    "vista": {
        "label": "VISTA",
        "combined_dir": "combined_v2",
    },

    "gimvi": {
        "label": "gimVI",
        "combined_dir": "combined_v2",
    },

    "tangram": {
        "label": "Tangram",
        "combined_dir": "combined_v2",
    },

    "envi": {
        "label": "ENVI",
        "combined_dir": "combined_v2",
    },

    "spage": {
        "label": "SpaGE",
        "combined_dir": "combined_v2",
    },

    "transimpspa": {
        "label": "TransImpSpa",
        "combined_dir": "combined_v2",
    },
}


REQUIRED_FILES = [
    "oof_predictions_300genes.h5ad",
    "oof_ground_truth_300genes.h5ad",
    "gene_level_metrics_300genes.csv",
]


EXPECTED_N_GENES = 300


# =============================================================================
# HELPERS
# =============================================================================

def matrix_min_max_finite(matrix):

    if sparse.issparse(matrix):

        data = matrix.data

        if data.size == 0:
            return {
                "finite": True,
                "min": 0.0,
                "max": 0.0,
            }

        finite = bool(
            np.isfinite(data).all()
        )

        min_val = float(
            data.min()
        )

        max_val = float(
            data.max()
        )

    else:

        arr = np.asarray(
            matrix
        )

        finite = bool(
            np.isfinite(arr).all()
        )

        min_val = float(
            np.nanmin(arr)
        )

        max_val = float(
            np.nanmax(arr)
        )

    return {
        "finite": finite,
        "min": min_val,
        "max": max_val,
    }


def exact_array_equal(a, b):

    if a is None or b is None:
        return False

    if a.shape != b.shape:
        return False

    return bool(
        np.array_equal(
            np.asarray(a),
            np.asarray(b),
        )
    )


def audit_one(
    experiment_key,
    experiment_label,
    model_key,
    model_info,
):

    model_label = (
        model_info["label"]
    )

    combined_dir = (
        model_info["combined_dir"]
    )

    run_dir = (
        BASE_DIR
        / experiment_key
        / model_key
        / combined_dir
    )

    print()
    print("=" * 110)
    print(
        f"{experiment_label} | {model_label}"
    )
    print("=" * 110)

    print(
        "Directory:",
        run_dir,
    )

    result = {
        "experiment":
            experiment_key,

        "experiment_label":
            experiment_label,

        "model":
            model_key,

        "model_label":
            model_label,

        "combined_dir":
            combined_dir,

        "run_dir":
            str(run_dir),

        "directory_exists":
            run_dir.is_dir(),
    }


    # -------------------------------------------------------------------------
    # Required files
    # -------------------------------------------------------------------------

    for filename in REQUIRED_FILES:

        result[
            f"exists__{filename}"
        ] = (
            run_dir
            .joinpath(filename)
            .is_file()
        )


    missing = [
        filename
        for filename in REQUIRED_FILES
        if not (
            run_dir
            / filename
        ).is_file()
    ]


    if missing:

        print(
            "MISSING:",
            missing,
        )

        result["status"] = "FAIL_MISSING_FILES"

        return result, None


    pred_file = (
        run_dir
        / "oof_predictions_300genes.h5ad"
    )

    truth_file = (
        run_dir
        / "oof_ground_truth_300genes.h5ad"
    )

    metrics_file = (
        run_dir
        / "gene_level_metrics_300genes.csv"
    )


    # -------------------------------------------------------------------------
    # Read H5ADs backed
    # -------------------------------------------------------------------------

    pred = ad.read_h5ad(
        pred_file,
        backed="r",
    )

    truth = ad.read_h5ad(
        truth_file,
        backed="r",
    )


    result["pred_n_obs"] = int(
        pred.n_obs
    )

    result["pred_n_vars"] = int(
        pred.n_vars
    )

    result["truth_n_obs"] = int(
        truth.n_obs
    )

    result["truth_n_vars"] = int(
        truth.n_vars
    )


    print(
        "Prediction shape:",
        pred.shape,
    )

    print(
        "Truth shape:",
        truth.shape,
    )


    # -------------------------------------------------------------------------
    # Alignment within this run
    # -------------------------------------------------------------------------

    result[
        "obs_names_match_pred_truth"
    ] = bool(
        pred.obs_names.equals(
            truth.obs_names
        )
    )

    result[
        "var_names_match_pred_truth"
    ] = bool(
        pred.var_names.equals(
            truth.var_names
        )
    )

    result[
        "expected_300_genes"
    ] = bool(
        pred.n_vars
        == EXPECTED_N_GENES
        and
        truth.n_vars
        == EXPECTED_N_GENES
    )


    # -------------------------------------------------------------------------
    # Spatial coordinates
    # -------------------------------------------------------------------------

    result[
        "pred_has_spatial"
    ] = (
        "spatial"
        in pred.obsm
    )

    result[
        "truth_has_spatial"
    ] = (
        "spatial"
        in truth.obsm
    )


    if (
        result["pred_has_spatial"]
        and
        result["truth_has_spatial"]
    ):

        pred_spatial = np.asarray(
            pred.obsm["spatial"]
        )

        truth_spatial = np.asarray(
            truth.obsm["spatial"]
        )

        result[
            "spatial_match_pred_truth"
        ] = exact_array_equal(
            pred_spatial,
            truth_spatial,
        )

    else:

        result[
            "spatial_match_pred_truth"
        ] = False


    # -------------------------------------------------------------------------
    # Prediction layers
    # -------------------------------------------------------------------------

    result[
        "prediction_layers"
    ] = ";".join(
        list(
            pred.layers.keys()
        )
    )

    result[
        "has_count_scale_layer"
    ] = (
        "count_scale"
        in pred.layers
    )

    result[
        "has_log1p_layer"
    ] = (
        "log1p"
        in pred.layers
    )

    result[
        "has_native_layer"
    ] = (
        "native"
        in pred.layers
    )


    # -------------------------------------------------------------------------
    # Basic prediction matrix diagnostics
    # -------------------------------------------------------------------------

    if (
        "count_scale"
        in pred.layers
    ):

        pred_matrix = (
            pred.layers[
                "count_scale"
            ]
        )

        result[
            "prediction_matrix_used"
        ] = "layers[count_scale]"

    else:

        pred_matrix = pred.X

        result[
            "prediction_matrix_used"
        ] = "X"


    pred_stats = (
        matrix_min_max_finite(
            pred_matrix
        )
    )

    result[
        "prediction_finite"
    ] = pred_stats[
        "finite"
    ]

    result[
        "prediction_min"
    ] = pred_stats[
        "min"
    ]

    result[
        "prediction_max"
    ] = pred_stats[
        "max"
    ]

    result[
        "prediction_nonnegative"
    ] = bool(
        pred_stats[
            "min"
        ]
        >= 0
    )


    # -------------------------------------------------------------------------
    # Ground truth diagnostics
    # -------------------------------------------------------------------------

    truth_stats = (
        matrix_min_max_finite(
            truth.X
        )
    )

    result[
        "truth_finite"
    ] = truth_stats[
        "finite"
    ]

    result[
        "truth_min"
    ] = truth_stats[
        "min"
    ]

    result[
        "truth_max"
    ] = truth_stats[
        "max"
    ]

    result[
        "truth_nonnegative"
    ] = bool(
        truth_stats[
            "min"
        ]
        >= 0
    )


    # -------------------------------------------------------------------------
    # Metrics
    # -------------------------------------------------------------------------

    metrics = pd.read_csv(
        metrics_file
    )

    result[
        "metrics_n_rows"
    ] = int(
        len(metrics)
    )

    result[
        "metrics_columns"
    ] = ";".join(
        metrics.columns.astype(str)
    )


    expected_metric_names = [
        "SCC",
        "SSIM",
        "RMSE",
        "MAE",
        "JSD",
        "Moran_error",
    ]


    lower_lookup = {
        str(col).lower():
            str(col)
        for col in metrics.columns
    }


    for metric in expected_metric_names:

        result[
            f"metric_present__{metric}"
        ] = (
            metric.lower()
            in lower_lookup
        )


    # -------------------------------------------------------------------------
    # Completion flag
    # -------------------------------------------------------------------------

    result[
        "aggregate_complete_flag"
    ] = (
        run_dir
        .joinpath(
            "aggregate_complete.flag"
        )
        .is_file()
    )


    # -------------------------------------------------------------------------
    # Overall status
    # -------------------------------------------------------------------------

    critical_checks = [
        result[
            "obs_names_match_pred_truth"
        ],

        result[
            "var_names_match_pred_truth"
        ],

        result[
            "expected_300_genes"
        ],

        result[
            "prediction_finite"
        ],

        result[
            "truth_finite"
        ],

        result[
            "prediction_nonnegative"
        ],

        result[
            "truth_nonnegative"
        ],
    ]


    if all(
        critical_checks
    ):

        result["status"] = "PASS"

    else:

        result["status"] = "FAIL_CHECK"


    print(
        "Status:",
        result["status"],
    )


    reference_info = {
        "obs_names":
            pred.obs_names.copy(),

        "var_names":
            pred.var_names.copy(),

        "truth_obs_names":
            truth.obs_names.copy(),

        "truth_var_names":
            truth.var_names.copy(),

        "spatial":
            (
                np.asarray(
                    truth.obsm[
                        "spatial"
                    ]
                ).copy()
                if "spatial"
                in truth.obsm
                else None
            ),
    }


    try:
        pred.file.close()
    except Exception:
        pass

    try:
        truth.file.close()
    except Exception:
        pass


    return result, reference_info


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 110)
    print(
        "GENEBRIDGE FINAL VISUALIZATION INPUT AUDIT"
    )
    print("=" * 110)

    print()
    print(
        "Expected combinations:",
        len(EXPERIMENTS)
        * len(MODELS),
    )


    rows = []

    references = {}


    for experiment_key, experiment_label in EXPERIMENTS.items():

        for model_key, model_info in MODELS.items():

            result, reference_info = audit_one(
                experiment_key,
                experiment_label,
                model_key,
                model_info,
            )

            rows.append(
                result
            )

            if reference_info is not None:

                references[
                    (
                        experiment_key,
                        model_key,
                    )
                ] = (
                    reference_info
                )


    # =========================================================================
    # Save raw audit table
    # =========================================================================

    audit = pd.DataFrame(
        rows
    )


    audit_file = (
        OUT_DIR
        / "visualization_input_audit.csv"
    )


    audit.to_csv(
        audit_file,
        index=False,
    )


    # =========================================================================
    # Cross-run alignment
    # =========================================================================

    alignment_rows = []


    if references:

        first_key = next(
            iter(
                references
            )
        )

        canonical = (
            references[
                first_key
            ]
        )


        print()
        print("=" * 110)
        print(
            "CROSS-RUN ALIGNMENT"
        )
        print("=" * 110)

        print(
            "Canonical reference:",
            first_key,
        )


        for (
            experiment_key,
            model_key
        ), ref in references.items():

            obs_match = bool(
                ref[
                    "obs_names"
                ].equals(
                    canonical[
                        "obs_names"
                    ]
                )
            )

            var_match = bool(
                ref[
                    "var_names"
                ].equals(
                    canonical[
                        "var_names"
                    ]
                )
            )

            truth_obs_match = bool(
                ref[
                    "truth_obs_names"
                ].equals(
                    canonical[
                        "truth_obs_names"
                    ]
                )
            )

            truth_var_match = bool(
                ref[
                    "truth_var_names"
                ].equals(
                    canonical[
                        "truth_var_names"
                    ]
                )
            )


            if (
                canonical[
                    "spatial"
                ] is None
                or
                ref[
                    "spatial"
                ] is None
            ):

                spatial_match = False

            else:

                spatial_match = (
                    exact_array_equal(
                        ref[
                            "spatial"
                        ],
                        canonical[
                            "spatial"
                        ],
                    )
                )


            row = {
                "experiment":
                    experiment_key,

                "model":
                    model_key,

                "obs_names_match":
                    obs_match,

                "var_names_match":
                    var_match,

                "truth_obs_names_match":
                    truth_obs_match,

                "truth_var_names_match":
                    truth_var_match,

                "spatial_match":
                    spatial_match,
            }


            row[
                "alignment_status"
            ] = (
                "PASS"
                if all(
                    [
                        obs_match,
                        var_match,
                        truth_obs_match,
                        truth_var_match,
                    ]
                )
                else
                "FAIL"
            )


            alignment_rows.append(
                row
            )


            print(
                f"{experiment_key:<8} "
                f"{model_key:<14} "
                f"{row['alignment_status']}"
            )


    alignment = pd.DataFrame(
        alignment_rows
    )


    alignment_file = (
        OUT_DIR
        / "cross_run_alignment.csv"
    )


    alignment.to_csv(
        alignment_file,
        index=False,
    )


    # =========================================================================
    # Compact summary
    # =========================================================================

    summary_columns = [
        "experiment",
        "model",
        "pred_n_obs",
        "pred_n_vars",
        "truth_n_obs",
        "truth_n_vars",
        "has_count_scale_layer",
        "aggregate_complete_flag",
        "status",
    ]


    available_summary_columns = [
        x
        for x in summary_columns
        if x in audit.columns
    ]


    compact = (
        audit[
            available_summary_columns
        ]
        .copy()
    )


    compact_file = (
        OUT_DIR
        / "visualization_input_audit_compact.csv"
    )


    compact.to_csv(
        compact_file,
        index=False,
    )


    print()
    print("=" * 110)
    print(
        "AUDIT SUMMARY"
    )
    print("=" * 110)

    print(
        compact.to_string(
            index=False
        )
    )


    n_pass = int(
        (
            audit[
                "status"
            ]
            == "PASS"
        ).sum()
    )


    print()
    print(
        f"PASS: {n_pass}/{len(audit)}"
    )


    if len(alignment) > 0:

        n_alignment_pass = int(
            (
                alignment[
                    "alignment_status"
                ]
                == "PASS"
            ).sum()
        )

        print(
            "Cross-run alignment PASS:",
            f"{n_alignment_pass}/{len(alignment)}",
        )


    # =========================================================================
    # Manifest
    # =========================================================================

    manifest = {
        "experiments":
            EXPERIMENTS,

        "models":
            MODELS,

        "expected_combinations":
            len(EXPERIMENTS)
            * len(MODELS),

        "audit_pass":
            n_pass,

        "audit_file":
            str(
                audit_file
            ),

        "alignment_file":
            str(
                alignment_file
            ),
    }


    with (
        OUT_DIR
        / "audit_manifest.json"
    ).open(
        "w"
    ) as handle:

        json.dump(
            manifest,
            handle,
            indent=2,
        )


    print()
    print(
        "Outputs:",
        OUT_DIR,
    )


if __name__ == "__main__":
    main()
