import numpy as np
import scanpy as sc

adata = sc.read_h5ad("data/processed/pbmc3k_mt8.h5ad")

# Ham sayımlarla başla.
adata.X = adata.layers["counts"].copy()

genes_before = adata.n_vars

# En az üç hücrede ölçülen genleri tut.
sc.pp.filter_genes(adata, min_cells=3)

print("\n--- Gen filtresi ---")
print("Önce:", genes_before)
print("Sonra:", adata.n_vars)
print("Hücre sayısı:", adata.n_obs)

# Ham sayımları koruyan counts katmanı da genlerle birlikte daralır.
sc.pp.normalize_total(adata, target_sum=10000)

totals = np.asarray(adata.X.sum(axis=1)).ravel()
print("\n--- Log öncesi normalize toplamlar ---")
print("En düşük:", round(float(totals.min()), 2))
print("En yüksek:", round(float(totals.max()), 2))

sc.pp.log1p(adata)

# Genleri silmeden HVG etiketlerini ekle.
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000,
    subset=False
)

hvg = adata.var["highly_variable"]

print("\n--- Değişken gen seçimi ---")
print("İşaretlenen HVG:", int(hvg.sum()))
print("HVG olan mitokondriyal gen:", int((hvg & adata.var["mt"]).sum()))

print("\n--- Veri bölümleri ---")
print("X: ln(1 + normalize ifade)")
print("counts: tutulan genlerin ham sayımları")
print("X türü:", type(adata.X).__name__)
print("Boyut:", adata.shape)

print("\n--- Ham sayım kontrolü ---")
values = adata.layers["counts"].data
print("Tam sayı niteliğinde:",
      bool(np.allclose(values, np.rint(values))))

adata.var.to_csv(
    "results/tables/pbmc3k_mt8_gene_hvg_summary.csv",
    index_label="gene_symbol"
)

adata.uns["preprocessing_choices"] = {
    "minimum_cells_per_gene": 3,
    "normalization_target_sum": 10000,
    "X_representation": "natural_log1p_normalized",
    "counts_representation": "unnormalized_counts",
    "hvg_flavor": "seurat",
    "requested_hvg_count": 2000,
}

adata.write_h5ad(
    "data/processed/pbmc3k_mt8_normalized.h5ad",
    compression="gzip"
)

print("\nNormalize nesne ve gen seçim tablosu kaydedildi.")
