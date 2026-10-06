from pathlib import Path
import scanpy as sc

adata = sc.read_h5ad("data/raw/pbmc3k_raw.h5ad")

# Başlangıç sayımlarını ayrı bir katmanda koru.
adata.layers["counts"] = adata.X.copy()

# Bu insan verisinde MT- ile başlayan genleri işaretle.
adata.var["mt"] = adata.var_names.str.startswith("MT-")

if not adata.var["mt"].any():
    raise RuntimeError("Mitokondriyal gen bulunamadı; gen isimlerini incele.")

# Ham sayımlardan QC ölçütlerini hesapla.
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    layer="counts",
    percent_top=None,
    log1p=False,
    inplace=True
)

columns = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt"
]

print("\n--- Mitokondriyal genler ---")
print(adata.var_names[adata.var["mt"]].tolist())

print("\n--- İlk beş hücrenin QC ölçütleri ---")
print(adata.obs[columns].head().round(2))

print("\n--- QC dağılım özeti ---")
print(
    adata.obs[columns].describe(
        percentiles=[0.01, 0.05, 0.50, 0.95, 0.99]
    ).round(2)
)

Path("results/tables").mkdir(parents=True, exist_ok=True)
Path("data/processed").mkdir(parents=True, exist_ok=True)

adata.obs[columns].to_csv(
    "results/tables/pbmc3k_qc_before_filtering.csv",
    index_label="cell_barcode"
)

adata.write_h5ad(
    "data/processed/pbmc3k_qc_unfiltered.h5ad",
    compression="gzip"
)

print("\nHücre sayısı:", adata.n_obs)
print("Hiçbir hücre veya gen çıkarılmadı.")
print("QC tablosu ve QC bilgileri eklenmiş nesne kaydedildi.")
