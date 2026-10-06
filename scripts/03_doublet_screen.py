import inspect
import numpy as np
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_qc_unfiltered.h5ad"
)

# Scrublet girdisinin ham sayım olduğundan emin ol.
adata.X = adata.layers["counts"].copy()

# Eğitim amaçlı başlangıç varsayımı:
# Bu veride gerçek doublet oranının %5 olduğunu bilmiyoruz.
parameters = {
    "expected_doublet_rate": 0.05,
    "n_prin_comps": 30,
    "use_approx_neighbors": False,
}

# Scanpy sürümüne uygun rastgelelik parametresi.
signature = inspect.signature(sc.pp.scrublet)
if "rng" in signature.parameters:
    parameters["rng"] = 0
elif "random_state" in signature.parameters:
    parameters["random_state"] = 0

print("Scrublet değerlendirmesi başlıyor...", flush=True)
sc.pp.scrublet(adata, **parameters)

print("\n--- Doublet skorlarının özeti ---")
print(adata.obs["doublet_score"].describe().round(4))

if "predicted_doublet" in adata.obs:
    print("\n--- Tahmin etiketleri ---")
    print(adata.obs["predicted_doublet"].value_counts(dropna=False))

print("\n--- Otomatik skor sınırı ---")
print(adata.uns["scrublet"].get("threshold", "Sınır oluşturulmadı"))

columns = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "doublet_score"
]
if "predicted_doublet" in adata.obs:
    columns.append("predicted_doublet")

print("\n--- Skoru en yüksek 10 hücre ---")
print(
    adata.obs[columns]
    .sort_values("doublet_score", ascending=False)
    .head(10).round(4)
)

adata.obs[columns].to_csv(
    "results/tables/pbmc3k_doublet_screen.csv",
    index_label="cell_barcode"
)

raw_unchanged = (adata.X != adata.layers["counts"]).nnz == 0
print("\nX içindeki ham sayımlar korunmuş:", raw_unchanged)

adata.write_h5ad(
    "data/processed/pbmc3k_qc_doublet_unfiltered.h5ad",
    compression="gzip"
)

print("\nDeğerlendirme kaydedildi; hiçbir hücre çıkarılmadı.")
