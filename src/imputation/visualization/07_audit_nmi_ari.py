#!/usr/bin/env python3

from pathlib import Path
import json
import pandas as pd


ROOT = Path(
    "/beegfs/labs/hulab/projects/mjabin/GeneBridge"
)

BASE = (
    ROOT
    / "outputs"
    / "imputation_beta"
    / "Br8667"
)

OUT_DIR = (
    BASE
    / "final_visualizations"
    / "cell_population"
    / "nmi_ari"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


EXPERIMENTS = [
    "ex5",
    "ex5_1",
    "ex5_3",
]

METHODS = [
    "vista",
    "gimvi",
    "tangram",
    "envi",
    "spage",
    "transimpspa",
]


SEARCH_SUFFIXES = {
    ".csv",
    ".tsv",
    ".json",
}


def inspect_csv(path):

    try:
        if path.suffix.lower() == ".tsv":
            df = pd.read_csv(
                path,
                sep="\t",
            )
        else:
            df = pd.read_csv(
                path
            )

    except Exception as exc:

        return [
            {
                "source_file": str(path),
                "source_type": "table",
                "status": f"READ_ERROR: {exc}",
                "metric": "",
                "column_or_key": "",
                "n_rows": "",
                "values_preview": "",
            }
        ]


    matches = []

    for column in df.columns:

        normalized = str(
            column
        ).lower()

        if (
            "nmi" in normalized
            or
            "ari" in normalized
            or
            "adjusted_rand" in normalized
            or
            "mutual_info" in normalized
        ):

            vals = (
                df[column]
                .dropna()
                .astype(str)
                .head(10)
                .tolist()
            )

            metric = (
                "NMI"
                if (
                    "nmi" in normalized
                    or
                    "mutual_info" in normalized
                )
                else
                "ARI"
            )

            matches.append(
                {
                    "source_file": str(path),
                    "source_type": "table",
                    "status": "FOUND",
                    "metric": metric,
                    "column_or_key": str(column),
                    "n_rows": len(df),
                    "values_preview": " | ".join(vals),
                }
            )

    return matches


def walk_json(
    obj,
    prefix="",
):

    hits = []

    if isinstance(
        obj,
        dict,
    ):

        for key, value in obj.items():

            current = (
                f"{prefix}.{key}"
                if prefix
                else str(key)
            )

            normalized = str(
                key
            ).lower()

            if (
                "nmi" in normalized
                or
                "ari" in normalized
                or
                "adjusted_rand" in normalized
                or
                "mutual_info" in normalized
            ):

                metric = (
                    "NMI"
                    if (
                        "nmi" in normalized
                        or
                        "mutual_info" in normalized
                    )
                    else
                    "ARI"
                )

                hits.append(
                    (
                        metric,
                        current,
                        value,
                    )
                )

            hits.extend(
                walk_json(
                    value,
                    current,
                )
            )


    elif isinstance(
        obj,
        list,
    ):

        for idx, value in enumerate(
            obj
        ):

            hits.extend(
                walk_json(
                    value,
                    f"{prefix}[{idx}]",
                )
            )

    return hits


def inspect_json(path):

    try:

        with path.open() as handle:
            obj = json.load(
                handle
            )

    except Exception as exc:

        return [
            {
                "source_file": str(path),
                "source_type": "json",
                "status": f"READ_ERROR: {exc}",
                "metric": "",
                "column_or_key": "",
                "n_rows": "",
                "values_preview": "",
            }
        ]


    matches = []

    for (
        metric,
        key,
        value,
    ) in walk_json(obj):

        preview = str(
            value
        )

        if len(preview) > 300:
            preview = (
                preview[:300]
                + "..."
            )

        matches.append(
            {
                "source_file": str(path),
                "source_type": "json",
                "status": "FOUND",
                "metric": metric,
                "column_or_key": key,
                "n_rows": "",
                "values_preview": preview,
            }
        )

    return matches


def main():

    results = []


    print(
        "=" * 100
    )
    print(
        "AUDIT NMI / ARI OUTPUTS"
    )
    print(
        "=" * 100
    )


    for experiment in EXPERIMENTS:

        for method in METHODS:

            run_dir = (
                BASE
                / experiment
                / method
                / "combined_v2"
            )


            print()
            print(
                f"{experiment:7s} | "
                f"{method:12s} | "
                f"{run_dir}"
            )


            if not run_dir.is_dir():

                print(
                    "  ERROR: combined_v2 missing"
                )

                results.append(
                    {
                        "experiment": experiment,
                        "method": method,
                        "source_file": "",
                        "source_type": "",
                        "status": "RUN_DIR_MISSING",
                        "metric": "",
                        "column_or_key": "",
                        "n_rows": "",
                        "values_preview": "",
                    }
                )

                continue


            files = sorted(
                path
                for path in run_dir.rglob("*")
                if (
                    path.is_file()
                    and
                    path.suffix.lower()
                    in SEARCH_SUFFIXES
                )
            )


            run_hits = []


            for path in files:

                if path.suffix.lower() in {
                    ".csv",
                    ".tsv",
                }:

                    hits = inspect_csv(
                        path
                    )

                else:

                    hits = inspect_json(
                        path
                    )


                for hit in hits:

                    if (
                        hit[
                            "status"
                        ]
                        == "FOUND"
                    ):

                        run_hits.append(
                            hit
                        )


            if not run_hits:

                print(
                    "  NMI/ARI: NOT FOUND"
                )

                results.append(
                    {
                        "experiment": experiment,
                        "method": method,
                        "source_file": "",
                        "source_type": "",
                        "status": "NOT_FOUND",
                        "metric": "",
                        "column_or_key": "",
                        "n_rows": "",
                        "values_preview": "",
                    }
                )

                continue


            for hit in run_hits:

                print(
                    f"  {hit['metric']:3s} | "
                    f"{Path(hit['source_file']).name} | "
                    f"{hit['column_or_key']} | "
                    f"rows={hit['n_rows']} | "
                    f"{hit['values_preview'][:120]}"
                )


                results.append(
                    {
                        "experiment": experiment,
                        "method": method,
                        **hit,
                    }
                )


    result_df = pd.DataFrame(
        results
    )


    output = (
        OUT_DIR
        / "nmi_ari_source_audit.csv"
    )

    result_df.to_csv(
        output,
        index=False,
    )


    print()
    print(
        "=" * 100
    )
    print(
        "SUMMARY"
    )
    print(
        "=" * 100
    )


    summary = (
        result_df.groupby(
            [
                "experiment",
                "method",
                "metric",
                "status",
            ],
            dropna=False,
        )
        .size()
        .reset_index(
            name="n_matches"
        )
    )


    print(
        summary.to_string(
            index=False
        )
    )


    print()
    print(
        "Saved:"
    )
    print(
        output
    )


if __name__ == "__main__":
    main()
