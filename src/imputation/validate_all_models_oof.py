#!/usr/bin/env python3

from pathlib import Path
import json
import anndata as ad

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

MODELS = {
    "ENVI":     ("ENVI_N24",     "envi"),
    "VISTA":    ("VISTA_N24",    "vista"),
    "gimVI":    ("GIMVI_N24",    "gimvi"),
    "Tangram":  ("TANGRAM_N24",  "tangram"),
    "SpaGE":    ("SPAGE_N24",    "spage"),
    "TransImp": ("TRANSIMP_N24", "transimpspa"),
}

DONORS = [
    "Br1113","Br1139","Br2039","Br2421","Br2719","Br5314",
    "Br5373","Br5400","Br5436","Br5588","Br5590","Br5622",
    "Br5639","Br5746","Br5931","Br5973","Br6032","Br6389",
    "Br6432","Br6437","Br6496","Br8433","Br8667","Br8772",
]

EXPECTED = [(d, "ex5_1") for d in DONORS]
EXPECTED += [
    ("Br6432", "ex5"),
    ("Br8667", "ex5"),
    ("Br6432", "ex5_3"),
    ("Br8667", "ex5_3"),
]

REQUIRED_LAYERS = {
    "count_scale",
    "log1p",
    "native",
    "ground_truth_counts",
    "ground_truth_log1p_norm",
}

# Canonical 300-gene order
canonical_path = (
    ROOT / "data/processed/xenium/"
    "xenium_N24_layer_celltype_annotated.h5ad"
)

canonical = ad.read_h5ad(canonical_path, backed="r")
canonical_genes = list(map(str, canonical.var_names))
canonical.file.close()

print("=" * 120)
print("ALL-MODEL OOF-300 VALIDATION")
print("=" * 120)

grand_pass = 0
grand_expected = len(MODELS) * len(EXPECTED)

model_results = {}

for model, (root_name, model_dir) in MODELS.items():

    print()
    print("=" * 120)
    print(model)
    print("=" * 120)

    model_pass = 0

    for donor, experiment in EXPECTED:

        oof_dir = (
            ROOT
            / "outputs/imputation_full"
            / root_name
            / donor
            / experiment
            / model_dir
            / "oof_300"
        )

        h5ad_path = oof_dir / "oof_predicted_300genes.h5ad"
        summary_path = oof_dir / "oof_stitch_summary.json"

        problems = []

        if not h5ad_path.exists():
            problems.append("MISSING_H5AD")

        if not summary_path.exists():
            problems.append("MISSING_SUMMARY")

        if not problems:
            try:
                a = ad.read_h5ad(h5ad_path, backed="r")

                if a.n_vars != 300:
                    problems.append(f"NGENES={a.n_vars}")

                if list(map(str, a.var_names)) != canonical_genes:
                    problems.append("GENE_ORDER")

                missing_layers = REQUIRED_LAYERS - set(a.layers.keys())
                if missing_layers:
                    problems.append(
                        "LAYERS=" + ",".join(sorted(missing_layers))
                    )

                if "oof_fold" not in a.var.columns:
                    problems.append("NO_OOF_FOLD")
                else:
                    fold_counts = (
                        a.var["oof_fold"]
                        .value_counts()
                        .sort_index()
                        .to_dict()
                    )

                    if fold_counts != {
                        1: 100,
                        2: 100,
                        3: 100,
                    }:
                        problems.append(f"FOLDS={fold_counts}")

                n_cells = a.n_obs
                a.file.close()

                with open(summary_path) as handle:
                    summary = json.load(handle)

                if summary.get("status") != "PASS":
                    problems.append(
                        f"SUMMARY={summary.get('status')}"
                    )

            except Exception as e:
                problems.append(
                    f"ERROR={type(e).__name__}:{e}"
                )

        if problems:
            status = "FAIL"
            detail = " | ".join(problems)
        else:
            status = "PASS"
            detail = f"{n_cells:,} cells × 300 genes"
            model_pass += 1
            grand_pass += 1

        print(
            f"{donor:8s} "
            f"{experiment:6s} "
            f"{status:4s}  "
            f"{detail}"
        )

    model_results[model] = model_pass

    print()
    print(f"{model}: {model_pass}/28 PASS")


print()
print("=" * 120)
print("FINAL SUMMARY")
print("=" * 120)

for model, count in model_results.items():
    marker = "PASS" if count == 28 else "FAIL"
    print(
        f"{model:10s} : "
        f"{count:2d}/28  {marker}"
    )

print("-" * 120)
print(
    f"TOTAL       : "
    f"{grand_pass}/{grand_expected}"
)

if grand_pass == grand_expected:
    print(
        "FINAL STATUS: PASS — "
        "ALL 6 MODELS × 28 OOF MATRICES READY"
    )
else:
    print(
        "FINAL STATUS: FAIL — "
        f"{grand_expected - grand_pass} OOF MATRIX/MATRICES NOT READY"
    )

print("=" * 120)
