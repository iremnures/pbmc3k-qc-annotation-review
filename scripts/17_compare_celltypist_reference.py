from pathlib import Path
from importlib.metadata import version

import anndata as ad
import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models

source = Path(
    "data/processed/pbmc3k_mt8_provisionally_annotated.h5ad"
)
adata = sc.read_h5ad(source)

assert adata.obs_names.is_unique
assert adata.var_names.is_unique
assert {"leiden_r10", "manual_annotation"}.issubset(adata.obs.columns)

model = models.Model.load(model="Immune_All_Low.pkl")

shared_genes = adata.var_names.intersection(model.features)
print("\n--- Reference gene coverage ---")
print("Genes in query:", adata.n_vars)
print("Genes in model:", len(model.features))
print("Shared gene symbols:", len(shared_genes))

# Use log-normalized expression; preserve the original object's layers.
query = ad.AnnData(
    X=adata.X.copy(),
    obs=pd.DataFrame(index=adata.obs_names.copy()),
    var=adata.var.copy()
)

print("\nCell-level prediction is starting...", flush=True)
prediction = celltypist.annotate(
    query,
    model=model,
    mode="best match",
    majority_voting=False
)

labels = prediction.predicted_labels["predicted_labels"].reindex(
    adata.obs_names
)
assert labels.notna().all()

adata.obs["celltypist_prediction"] = labels.astype("category")

tables = Path("results/tables")
tables.mkdir(parents=True, exist_ok=True)

adata.obs[
    ["leiden_r10", "manual_annotation", "celltypist_prediction"]
].to_csv(
    tables / "pbmc3k_celltypist_cell_comparison.csv",
    index_label="cell_barcode"
)

comparison = pd.crosstab(
    adata.obs["leiden_r10"],
    adata.obs["celltypist_prediction"]
)
comparison.to_csv(tables / "pbmc3k_celltypist_cluster_comparison.csv")

adata.uns["celltypist_comparison"] = {
    "model": "Immune_All_Low.pkl",
    "software_version": version("celltypist"),
    "mode": "best match",
    "majority_voting": False,
    "shared_gene_symbols": len(shared_genes),
    "interpretation": "Reference predictions; manual labels retained."
}

adata.write_h5ad(
    "data/processed/pbmc3k_mt8_reference_compared.h5ad",
    compression="gzip"
)

print("\n--- Top reference predictions within each cluster ---")
for cluster in sorted(comparison.index, key=lambda value: int(value)):
    row = comparison.loc[cluster]
    top = row[row > 0].sort_values(ascending=False).head(3)
    total = int(row.sum())

    print(f"\nCluster {cluster}: {total} cells")
    for label, count in top.items():
        print(f"  {label}: {count} ({100 * count / total:.1f}%)")

print("\nComparison saved; manual labels preserved.")
