import inspect
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_clustered.h5ad"
)

parameters = {
    "resolution": 1.0,
    "flavor": "igraph",
    "directed": False,
    "n_iterations": 2,
    "key_added": "leiden_r10",
}

signature = inspect.signature(sc.tl.leiden)
if "rng" in signature.parameters:
    parameters["rng"] = 0
elif "random_state" in signature.parameters:
    parameters["random_state"] = 0

# Aynı grafikte yalnızca kümeleme çözünürlüğünü değiştir.
sc.tl.leiden(adata, **parameters)

cross = pd.crosstab(
    adata.obs["leiden_r05"],
    adata.obs["leiden_r10"]
)

print("\n--- Eski küme 5 ve 6'nın yeni kümelere dağılımı ---")
print(cross.loc[["5", "6"]].to_string())

genes = [
    "PF4", "PPBP", "GP9",
    "LST1", "FCGR3A", "MS4A7",
    "FCER1A", "CD1C", "CLEC10A", "CLEC4C"
]
genes = [gene for gene in genes if gene in adata.var_names]

selected = adata.obs["leiden_r05"].astype(str).isin(["5", "6"])
expression = adata[selected, genes].to_df()

expression["old_cluster"] = (
    adata.obs.loc[expression.index, "leiden_r05"].astype(str)
)
expression["new_cluster"] = (
    adata.obs.loc[expression.index, "leiden_r10"].astype(str)
)

grouped = expression.groupby(["old_cluster", "new_cluster"])
means = grouped[genes].mean()
means.insert(0, "n_cells", grouped.size())

print("\n--- Alt gruplarda ortalama log-normalize ifade ---")
print(means.round(2).to_string())

cross.to_csv(
    "results/tables/pbmc3k_resolution_crosstab.csv"
)
means.to_csv(
    "results/tables/pbmc3k_clusters_5_6_resolution_markers.csv"
)

adata.write_h5ad(
    "data/processed/pbmc3k_mt8_resolution_review.h5ad",
    compression="gzip"
)

print("\nKarşılaştırma kaydedildi; ilk küme etiketleri korundu.")
