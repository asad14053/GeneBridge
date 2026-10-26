from pathlib import Path
from collections import defaultdict
import h5py
import pandas as pd

ROOT = Path("/beegfs/labs/hulab/projects/mjabin/GeneBridge")
DIR = ROOT / "data/raw/nemo_mccarroll/counts"
OUT = ROOT / "results/analysis/nemo_mccarroll/inspection"
OUT.mkdir(parents=True, exist_ok=True)

donor_files = sorted(DIR.glob("*.gmg.donors.h5"))[:2]
all_rows = []

def read_barcodes(path):
    with h5py.File(path, "r") as f:
        return list(f["matrix/barcodes"].asstr()[:])

def resolve(labels):
    if not labels:
        return pd.NA
    if len(labels) == 1:
        return next(iter(labels))
    return "AMBIGUOUS:" + ";".join(sorted(labels))

for donor_file in donor_files:

    prefix = donor_file.name.removesuffix(".h5")
    batch = donor_file.name.split(".intronic_")[0]

    print("\nBATCH:", batch)

    barcodes = read_barcodes(donor_file)
    barcode_set = set(barcodes)

    broad = defaultdict(set)
    fine = defaultdict(set)

    subtype_files = sorted(
        DIR.glob(prefix + ".digital_expression.*.counts.h5")
    )

    print("Matching subtype files:", len(subtype_files))

    for path in subtype_files:

        label = path.name.removeprefix(
            prefix + ".digital_expression."
        ).removesuffix(".counts.h5")

        parts = label.split(".", 1)
        if len(parts) != 2:
            print("Skipped unexpected label:", label)
            continue

        cell_class, subtype = parts

        try:
            subset_barcodes = set(read_barcodes(path))
        except (KeyError, OSError) as e:
            print("Could not read:", path.name, e)
            continue

        matched = barcode_set & subset_barcodes

        print(
            f"{label:35s} "
            f"subset={len(subset_barcodes):6d} "
            f"matched={len(matched):6d}"
        )

        for barcode in matched:
            broad[barcode].add(cell_class)

            if subtype != "All":
                fine[barcode].add(label)

    for barcode in barcodes:
        all_rows.append({
            "nucleus_id": batch + "_" + barcode,
            "barcode": barcode,
            "batch": batch,
            "cell_type": resolve(broad[barcode]),
            "cell_subtype": resolve(fine[barcode]),
            "donor_id": pd.NA,
            "diagnosis": pd.NA,
        })

obs = pd.DataFrame(all_rows).set_index("nucleus_id")

output = OUT / "ling_obs_two_batches.csv"
obs.to_csv(output)

print("\n=== COMBINED OBS SUMMARY ===")
print("Shape:", obs.shape)

print("\nCell types:")
print(obs["cell_type"].value_counts(dropna=False))

print("\nCell subtypes:")
print(obs["cell_subtype"].value_counts(dropna=False))

print("\nFirst 10 observations:")
print(obs.head(10).to_string())

print("\nSaved:", output)
