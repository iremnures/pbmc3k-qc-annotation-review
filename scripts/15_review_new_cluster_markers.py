import matplotlib
matplotlib.use("Agg")

import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_resolution_review.h5ad"
)

panels = {
    "T identity": ["CD3D", "CD3E", "TRAC"],
    "T patterns": ["IL7R", "CCR7", "CD8A", "CD8B"],
    "Cytotoxic": ["NKG7", "GNLY", "PRF1"],
    "B identity": ["MS4A1", "CD79A", "CD79B"],
    "Myeloid": ["LYZ", "LST1", "TYROBP"],
    "Monocyte patterns": ["FCN1", "S100A9", "FCGR3A", "MS4A7"],
    "DC patterns": ["FCER1A", "CD1C", "CLEC10A", "CLEC4C"],
    "Platelet patterns": ["PF4", "PPBP", "GP9"],
}

missing = [
    gene for genes in panels.values()
    for gene in genes if gene not in adata.var_names
]
print("\nEksik marker'lar:", missing)

available = {
    name: [gene for gene in genes if gene in adata.var_names]
    for name, genes in panels.items()
}
available = {
    name: genes for name, genes in available.items() if genes
}

genes = [gene for group in available.values() for gene in group]

# Yalnızca seçilen marker genlerini tabloya çevir.
expression = adata[:, genes].to_df()
clusters = adata.obs["leiden_r10"].astype(str)
clusters.name = "cluster"

means = expression.groupby(clusters).mean()
fractions = expression.gt(0).groupby(clusters).mean() * 100

means.to_csv(
    "results/tables/pbmc3k_r10_marker_mean_expression.csv"
)
fractions.to_csv(
    "results/tables/pbmc3k_r10_marker_detected_percent.csv"
)

plot = sc.pl.dotplot(
    adata,
    var_names=available,
    groupby="leiden_r10",
    use_raw=False,
    expression_cutoff=0,
    mean_only_expressed=False,
    dendrogram=False,
    figsize=(17, 6),
    return_fig=True
)

plot.savefig(
    "results/figures/pbmc3k_r10_marker_dotplot.png",
    dpi=180
)
plot.savefig(
    "results/figures/pbmc3k_r10_marker_dotplot.pdf"
)

print("\nYeni kümelerin marker kanıt tabloları ve görseli kaydedildi.")
print("Henüz kesin hücre tipi etiketi atanmadı.")
