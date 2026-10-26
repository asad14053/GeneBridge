#!/usr/bin/env python3

import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")

BUILDER = ROOT / "scripts/ling_harmonization/05_build_ling_huuki_compatible_reference.py"

spec = importlib.util.spec_from_file_location("lingbuilder", BUILDER)
B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)

COUNTS = ROOT / "data/raw/nemo_mccarroll/counts"
HARM = ROOT / "results/analysis/ling_huuki_kwon_gene_harmonization"
OUT = ROOT / "results/analysis/ling_reference_build"

META = COUNTS / "BA46.jointMetadata.txt"

OUT.mkdir(parents=True, exist_ok=True)

ling_path = B.find_file(
    [HARM / "Ling_union_symbol_to_ensembl.csv",
     HARM / "ling_union_symbol_to_ensembl.csv"],
    search_root=HARM,
    contains="ling_union_symbol",
)

huuki_path = B.find_file(
    [HARM / "Huuki_gene_to_ensembl.csv",
     HARM / "huuki_gene_to_ensembl.csv"],
    search_root=HARM,
    contains="huuki_gene",
)

kwon_path = B.find_file(
    [HARM / "Kwon_gene_to_ensembl.csv",
     HARM / "kwon_gene_to_ensembl.csv"],
    search_root=HARM,
    contains="kwon_gene",
)

prev_path = B.find_file(
    [
        OUT / "ling_postharmonization_gene_prevalence.csv",
        HARM / "ling_postharmonization_gene_prevalence.csv",
    ],
    search_root=ROOT / "results/analysis",
    contains="postharmonization_gene_prevalence",
)

_, _, ling_ens_to_symbol = B.load_crosswalk(ling_path, "Ling")
_, _, huuki_ens_to_symbol = B.load_crosswalk(huuki_path, "Huuki")
_, _, kwon_ens_to_symbol = B.load_crosswalk(kwon_path, "Kwon")

prev = B.load_gene_prevalence(prev_path, expected_libraries=87)

genes = (
    prev.loc[
        prev["library_fraction"] >= 0.95
    ]
    .sort_values("ensembl_id")
    .drop_duplicates("ensembl_id")
    .reset_index(drop=True)
)

if len(genes) != 31491:
    raise RuntimeError(
        f"Expected 31,491 genes; got {len(genes):,}"
    )

ensgs = genes["ensembl_id"].tolist()

# Canonical symbol priority:
# Kwon -> Huuki -> Ling -> ENSG
canonical = {}
source = {}
used = set()

for ens in ensgs:
    if ens in kwon_ens_to_symbol:
        sym = kwon_ens_to_symbol[ens]
        if sym in used:
            raise RuntimeError(f"Duplicate Kwon symbol: {sym}")
        canonical[ens] = sym
        source[ens] = "Kwon"
        used.add(sym)

for ens in ensgs:
    if ens in canonical:
        continue

    candidates = []

    if ens in huuki_ens_to_symbol:
        candidates.append((huuki_ens_to_symbol[ens], "Huuki"))

    if ens in ling_ens_to_symbol:
        candidates.append((ling_ens_to_symbol[ens], "Ling"))

    candidates.append((ens, "ENSG_fallback"))

    for sym, src in candidates:
        if sym not in used:
            canonical[ens] = sym
            source[ens] = src
            used.add(sym)
            break

gene_manifest = pd.DataFrame({
    "reference_column": np.arange(len(ensgs)),
    "ensembl_id": ensgs,
    "canonical_symbol": [canonical[e] for e in ensgs],
    "canonical_source": [source[e] for e in ensgs],
    "library_fraction": genes.set_index("ensembl_id")
                             .loc[ensgs, "library_fraction"]
                             .to_numpy(),
})

if gene_manifest["canonical_symbol"].duplicated().any():
    raise RuntimeError("Final canonical symbols are not unique.")

kwon_symbols = set(kwon_ens_to_symbol.values())

missing = kwon_symbols - set(gene_manifest["canonical_symbol"])

if missing:
    raise RuntimeError(
        f"Missing Kwon symbols: {sorted(missing)}"
    )

gene_manifest.to_csv(
    OUT / "array_gene_manifest.csv",
    index=False,
)

print("Genes:", len(gene_manifest))
print("Kwon:", len(kwon_symbols), "/", len(kwon_symbols))


# ----------------------------------------------------------------------
# Cells
# ----------------------------------------------------------------------

meta = pd.read_csv(
    META,
    sep="\t",
    low_memory=False,
)

selected, eligible, quotas = B.proportional_stratified_sample(
    meta,
    n_target=77604,
    seed=20261023,
)

selected.to_csv(
    OUT / "array_selected_cells.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

print("Metadata cells:", len(meta))
print("Eligible:", len(eligible))
print("Selected:", len(selected))
print("Selected donors:", selected["DONOR"].nunique())
print("Selected libraries:", selected["PREFIX"].nunique())


# ----------------------------------------------------------------------
# 87-library manifest
# ----------------------------------------------------------------------

h5s = sorted(COUNTS.glob("*.gmg.donors.h5"))

if len(h5s) != 87:
    raise RuntimeError(
        f"Expected 87 H5s; found {len(h5s)}"
    )

rows = []

for task_id, path in enumerate(h5s):

    prefix = B.filename_to_prefix(path)

    g = selected[
        selected["PREFIX"] == prefix
    ]

    if len(g) == 0:
        raise RuntimeError(
            f"No selected cells for library {prefix}"
        )

    rows.append({
        "task_id": task_id,
        "prefix": prefix,
        "h5_path": str(path),
        "n_selected": len(g),
        "first_reference_row": int(g["reference_row"].min()),
        "last_reference_row": int(g["reference_row"].max()),
    })

libs = pd.DataFrame(rows)

if len(libs) != 87:
    raise RuntimeError("Library manifest != 87")

libs.to_csv(
    OUT / "array_libraries.tsv",
    sep="\t",
    index=False,
)

print()
print(libs.to_string(index=False))

print("\nPREPARATION PASS")
