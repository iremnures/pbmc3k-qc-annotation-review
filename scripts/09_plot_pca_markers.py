import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_pca.h5ad"
)

genes = ["LYZ", "NKG7", "MS4A1", "PPBP"]
missing = [gene for gene in genes if gene not in adata.var_names]

if missing:
    raise RuntimeError(f"Genler bulunamadı: {missing}")

# Yalnızca dört geni tabloya çeviriyoruz.
expression = adata[:, genes].to_df()

coordinates = adata.obsm["X_pca"]
ratios = adata.uns["pca"]["variance_ratio"]

fig, axes = plt.subplots(2, 2, figsize=(10, 8))

for ax, gene, y_component in zip(
    axes.ravel(),
    genes,
    [1, 1, 1, 2]
):
    points = ax.scatter(
        coordinates[:, 0],
        coordinates[:, y_component],
        c=expression[gene],
        cmap="viridis",
        s=7,
        alpha=0.75,
        linewidths=0
    )

    ax.set_title(gene)
    ax.set_xlabel(f"PC1 ({100 * ratios[0]:.2f}%)")
    ax.set_ylabel(
        f"PC{y_component + 1} "
        f"({100 * ratios[y_component]:.2f}%)"
    )
    fig.colorbar(
        points,
        ax=ax,
        label="ln(1 + normalized expression)"
    )

fig.suptitle("PBMC3k: marker expression across PCA")
fig.tight_layout()

fig.savefig(
    "results/figures/pbmc3k_mt8_pca_markers.png",
    dpi=160
)
fig.savefig(
    "results/figures/pbmc3k_mt8_pca_markers.pdf"
)
plt.close(fig)

print("PCA marker figures saved.")
