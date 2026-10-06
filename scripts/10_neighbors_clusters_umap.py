import inspect
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_pca.h5ad"
)
sc.settings.n_jobs = 1

def seeded_call(function, **kwargs):
    parameters = inspect.signature(function).parameters
    if "rng" in parameters:
        kwargs["rng"] = 0
    elif "random_state" in parameters:
        kwargs["random_state"] = 0
    return function(adata, **kwargs)

print("Komşuluk grafiği hesaplanıyor...", flush=True)
seeded_call(
    sc.pp.neighbors,
    n_neighbors=15,
    n_pcs=20,
    use_rep="X_pca",
    transformer="sklearn",
    metric="euclidean"
)

print("Leiden kümeleri hesaplanıyor...", flush=True)
seeded_call(
    sc.tl.leiden,
    resolution=0.5,
    flavor="igraph",
    directed=False,
    n_iterations=2,
    key_added="leiden_r05"
)

print("UMAP hesaplanıyor...", flush=True)
seeded_call(sc.tl.umap, min_dist=0.5)

print("\n--- Küme büyüklükleri ---")
print(
    adata.obs["leiden_r05"]
    .value_counts()
    .sort_index()
)

print("\n--- Veri bölümleri ---")
print("Hücre koordinatları:", list(adata.obsm.keys()))
print("Hücreler arası matrisler:", list(adata.obsp.keys()))
print("UMAP boyutu:", adata.obsm["X_umap"].shape)

adata.obs[["leiden_r05"]].to_csv(
    "results/tables/pbmc3k_mt8_cluster_assignments.csv",
    index_label="cell_barcode"
)

adata.write_h5ad(
    "data/processed/pbmc3k_mt8_clustered.h5ad",
    compression="gzip"
)

sc.pl.umap(
    adata,
    color="leiden_r05",
    title="PBMC3k: Leiden clusters (resolution 0.5)",
    show=False
)
fig = plt.gcf()
fig.savefig(
    "results/figures/pbmc3k_mt8_umap_clusters.png",
    dpi=160,
    bbox_inches="tight"
)
plt.close(fig)

print("\nKümeler, UMAP ve görsel kaydedildi.")
