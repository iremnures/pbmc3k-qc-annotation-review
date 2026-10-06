from pathlib import Path
import hashlib
import json

import scanpy as sc
from celltypist import models

for folder in [
    "data/raw", "data/processed",
    "results/tables", "results/figures", "docs"
]:
    Path(folder).mkdir(parents=True, exist_ok=True)

raw_path = Path("data/raw/pbmc3k_raw.h5ad")

if raw_path.is_file():
    print("Existing PBMC3k input retained.")
else:
    print("Downloading PBMC3k...", flush=True)
    data = sc.datasets.pbmc3k()
    assert data.shape == (2700, 32738)
    data.write_h5ad(raw_path, compression="gzip")
    print("PBMC3k input saved.")

models.download_models(
    model="Immune_All_Low.pkl",
    force_update=False
)

model_path = (
    Path.home() / ".celltypist/data/models/Immune_All_Low.pkl"
)
if not model_path.is_file():
    raise FileNotFoundError("Expected CellTypist model file not found.")

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

record = {
    "dataset": "Scanpy PBMC3k; cell-called, unnormalized counts",
    "input_file": str(raw_path),
    "input_sha256": sha256(raw_path),
    "reference_model": model_path.name,
    "reference_model_sha256": sha256(model_path)
}

Path("docs/input_fingerprints.json").write_text(
    json.dumps(record, indent=2) + "\n",
    encoding="utf-8"
)

print("Inputs ready; file fingerprints saved.")
