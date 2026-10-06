import matplotlib
matplotlib.use("Agg")

import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_clustered.h5ad"
)

panels = {
    "T-associated": ["CD3D", "CD3E", "TRAC"],
    "T subsets": ["IL7R", "CCR7", "CD8A", "CD8B"],
    "Cytotoxic": ["NKG7", "GNLY", "PRF1"],
    "B-associated": ["MS4A1", "CD79A", "CD79B"],
    "Myeloid": ["LYZ", "LST1", "TYROBP"],
    "Monocyte patterns": ["FCN1", "S100A9", "FCGR3A", "MS4A7"],
    "Platelet-associated": ["PF4", "PPBP", "GP9"],
}

missing = [
    gene
    for genes in panels.values()
    for gene in genes
    if gene not in adata.var_names
]

print("\n--- Eksik marker genleri ---")
print(missing if missing else "Yok")

available = {
    group: [gene for gene in genes if gene in adata.var_names]
    for group, genes in panels.items()
}
available = {
    group: genes for group, genes in available.items() if genes
}

plot = sc.pl.dotplot(
    adata,
    var_names=available,
    groupby="leiden_r05",
    use_raw=False,
    expression_cutoff=0,
    mean_only_expressed=False,
    dendrogram=False,
    figsize=(15, 5),
    return_fig=True
)

plot.savefig(
    "results/figures/pbmc3k_mt8_marker_dotplot.png",
    dpi=180
)
plot.savefig(
    "results/figures/pbmc3k_mt8_marker_dotplot.pdf"
)

grouped = adata.obs.groupby("leiden_r05", observed=True)

summary = grouped[
    ["total_counts", "n_genes_by_counts",
     "pct_counts_mt", "doublet_score"]
].median()

summary.insert(0, "n_cells", grouped.size())

summary.to_csv(
    "results/tables/pbmc3k_mt8_cluster_qc_summary.csv",
    index_label="cluster"
)

print("\n--- Küme büyüklükleri ve QC ortancaları ---")
print(summary.round(3).to_string())
print("\nMarker dotplot ve küme QC tablosu kaydedildi.")
