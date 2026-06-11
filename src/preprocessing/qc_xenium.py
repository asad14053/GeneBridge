import gzip
import shutil
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
XENIUM_DIR = PROJECT_ROOT / "data" / "raw" / "xenium"

gz_file = list(XENIUM_DIR.rglob("*cells.parquet.gz"))[0]
out_file = gz_file.with_suffix("")  # removes .gz -> .parquet

print("GZ file:", gz_file)
print("Unzipped file:", out_file)

with gzip.open(gz_file, "rb") as f_in:
    with open(out_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

cells = pd.read_parquet(out_file)

print(cells.shape)
print(cells.columns.tolist())
print(cells.head())