import inspect
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_normalized.h5ad"
)

hvg_mask = adata.var["highly_variable"].to_numpy()
hvg_names = adata.var_names[hvg_mask]

# Yalnızca PCA için ayrı bir çalışma nesnesi.
work = ad.AnnData(
    X=adata.X[:, hvg_mask].copy(),
    obs=adata.obs.copy(),
    var=adata.var.loc[hvg_names].copy()
)

print("PCA çalışma boyutu:", work.shape, flush=True)

# Genleri standartlaştır; çok yüksek pozitif değerleri 10 ile sınırla.
sc.pp.scale(work, zero_center=True, max_value=10)

parameters = {
    "n_comps": 50,
    "svd_solver": "arpack",
}

signature = inspect.signature(sc.pp.pca)
if "rng" in signature.parameters:
    parameters["rng"] = 0
elif "random_state" in signature.parameters:
    parameters["random_state"] = 0

print("PCA hesaplanıyor...", flush=True)
sc.pp.pca(work, **parameters)

# Hücre sırasının aynı olduğunu doğrula.
assert adata.obs_names.equals(work.obs_names)

# Sonuçları ana nesneye ekle; X ve counts değişmez.
adata.obsm["X_pca"] = work.obsm["X_pca"].copy()
adata.uns["pca"] = work.uns["pca"].copy()
adata.uns["pca_workflow"] = {
    "input": "standardized log-normalized HVGs",
    "gene_count": work.n_vars,
    "scale_max_value": 10,
    "computed_components": 50,
}

ratios = work.uns["pca"]["variance_ratio"]

summary = pd.DataFrame({
    "component": np.arange(1, len(ratios) + 1),
    "explained_variance_percent": 100 * ratios,
    "cumulative_variance_percent": 100 * np.cumsum(ratios),
})

summary.to_csv(
    "results/tables/pbmc3k_mt8_pca_variance.csv",
    index=False
)

pd.DataFrame(
    work.varm["PCs"],
    index=hvg_names,
    columns=[f"PC{i}" for i in range(1, 51)]
).to_csv(
    "results/tables/pbmc3k_mt8_pca_gene_loadings.csv",
    index_label="gene_symbol"
)

adata.write_h5ad(
    "data/processed/pbmc3k_mt8_pca.h5ad",
    compression="gzip"
)

print("\n--- İlk 10 bileşen ---")
print(summary.head(10).round(2).to_string(index=False))

print("\n--- Veri bölümleri ---")
print("Ana X: log-normalize ifade")
print("counts: ham sayımlar")
print("PCA koordinatları:", adata.obsm["X_pca"].shape)
print("Gen ağırlıkları ayrı CSV dosyasına kaydedildi.")
