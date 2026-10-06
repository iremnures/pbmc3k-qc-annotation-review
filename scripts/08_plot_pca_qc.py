import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_pca.h5ad"
)

coordinates = adata.obsm["X_pca"]
ratios = adata.uns["pca"]["variance_ratio"]

fig, axes = plt.subplots(1, 3, figsize=(13, 4))

metrics = [
    ("total_counts", "Total UMI counts"),
    ("n_genes_by_counts", "Detected genes"),
    ("pct_counts_mt", "Mitochondrial counts (%)"),
]

for ax, (column, title) in zip(axes, metrics):
    points = ax.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=adata.obs[column],
        cmap="viridis",
        s=7,
        alpha=0.7,
        linewidths=0
    )
    ax.set_title(title)
    ax.set_xlabel(f"PC1 ({100 * ratios[0]:.2f}%)")
    ax.set_ylabel(f"PC2 ({100 * ratios[1]:.2f}%)")
    fig.colorbar(points, ax=ax)

fig.suptitle("PBMC3k: PCA colored by QC metrics")
fig.tight_layout()
fig.savefig(
    "results/figures/pbmc3k_mt8_pca_qc.png",
    dpi=160
)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4))
components = np.arange(1, len(ratios) + 1)
ax.plot(components, 100 * ratios, marker=".", color="#287C8E")
ax.set_xlabel("Principal component")
ax.set_ylabel("Explained variance (%)")
ax.set_title("PBMC3k: PCA variance by component")
fig.tight_layout()
fig.savefig(
    "results/figures/pbmc3k_mt8_pca_variance.png",
    dpi=160
)
plt.close(fig)

print("PCA QC and variance figures saved.")
