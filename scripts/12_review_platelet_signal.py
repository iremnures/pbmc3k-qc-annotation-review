import numpy as np
import pandas as pd
import scanpy as sc

before = sc.read_h5ad(
    "data/processed/pbmc3k_qc_doublet_unfiltered.h5ad"
)
current = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_clustered.h5ad"
)

genes = ["PF4", "PPBP", "GP9"]

missing = [gene for gene in genes if gene not in before.var_names]
if missing:
    raise RuntimeError(f"Genler bulunamadı: {missing}")

# Yalnızca üç genin ham sayımlarını tabloya çevir.
counts = before[:, genes].to_df(layer="counts")
counts.columns = [f"{gene}_counts" for gene in genes]

tab = before.obs[
    ["total_counts", "n_genes_by_counts",
     "pct_counts_mt", "doublet_score", "predicted_doublet"]
].join(counts)

# Farklı toplam sayımları dikkate alan tarama değeri.
totals = np.asarray(
    before.layers["counts"].sum(axis=1)
).ravel()

tab["PPBP_log_normalized"] = np.log1p(
    10000 * tab["PPBP_counts"].to_numpy()
    / np.maximum(totals, 1)
)

tab["retained_mt8"] = tab.index.isin(current.obs_names)
tab["current_cluster"] = (
    current.obs["leiden_r05"]
    .astype(str)
    .reindex(tab.index)
    .fillna("excluded")
)

top = tab.sort_values(
    "PPBP_log_normalized", ascending=False
).head(20)

print("\n--- PPBP sinyali en yüksek 20 başlangıç hücresi ---")
print(top.round(3).to_string())

print("\n--- Bu 20 hücrenin mevcut kümeleri ---")
print(top["current_cluster"].value_counts())

top.to_csv(
    "results/tables/pbmc3k_platelet_signal_review.csv",
    index_label="cell_barcode"
)

print("\nİnceleme tablosu kaydedildi; veri değiştirilmedi.")
