import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_clustered.h5ad"
)

print("Küme marker'ları sıralanıyor...", flush=True)

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_r05",
    groups=["5", "6"],
    reference="rest",
    method="wilcoxon",
    tie_correct=True,
    use_raw=False,
    n_genes=30,
    key_added="markers_5_6"
)

tables = []

for cluster in ["5", "6"]:
    result = sc.get.rank_genes_groups_df(
        adata,
        group=cluster,
        key="markers_5_6"
    )

    result.insert(0, "cluster", cluster)
    tables.append(result)

    print(f"\n--- Küme {cluster}: ilk 15 marker adayı ---")
    print(
        result[["names", "scores", "logfoldchanges", "pvals_adj"]]
        .head(15)
        .round(4)
        .to_string(index=False)
    )

pd.concat(tables, ignore_index=True).to_csv(
    "results/tables/pbmc3k_clusters_5_6_marker_candidates.csv",
    index=False
)

# Küme 5'te hücre düzeyinde birlikte ifade incelemesi.
cells5 = adata.obs["leiden_r05"].astype(str) == "5"
genes = ["PPBP", "PF4", "LST1", "LYZ", "FCGR3A"]

missing = [gene for gene in genes if gene not in adata.var_names]
if missing:
    raise RuntimeError(f"İncelenecek genler bulunamadı: {missing}")

expression = adata[cells5, genes].to_df()

print("\nKüme 5 hücre sayısı:", len(expression))

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

for ax, gene in zip(axes, ["LST1", "LYZ", "FCGR3A"]):
    points = ax.scatter(
        expression[gene],
        expression["PPBP"],
        c=expression["PF4"],
        cmap="viridis",
        s=20,
        alpha=0.75
    )
    ax.set_xlabel(f"{gene}: log-normalized expression")
    ax.set_ylabel("PPBP: log-normalized expression")
    fig.colorbar(points, ax=ax, label="PF4 expression")

fig.suptitle("Cluster 5: platelet and myeloid expression patterns")
fig.tight_layout()

fig.savefig(
    "results/figures/pbmc3k_cluster5_expression_review.png",
    dpi=180
)
plt.close(fig)

expression.to_csv(
    "results/tables/pbmc3k_cluster5_selected_gene_expression.csv",
    index_label="cell_barcode"
)

print("\nMarker tabloları ve hücre düzeyindeki görsel kaydedildi.")
