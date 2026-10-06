import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_resolution_review.h5ad"
)

records = [
    ("0", "IL7R-associated T cells",
     "CD3D, CD3E, IL7R", "T-cell subtype unresolved"),
    ("1", "Cytotoxic T candidate",
     "CD3D, CD3E, CD8A/B, NKG7", "Provisional subtype"),
    ("2", "B cells",
     "MS4A1, CD79A, CD79B", "Broad identity only"),
    ("3", "FCN1/S100A9-associated monocytes",
     "LYZ, LST1, TYROBP, FCN1, S100A9", "Fine subtype unresolved"),
    ("4", "FCGR3A-associated monocytes",
     "LST1, TYROBP, FCGR3A, MS4A7",
     "FCN1/S100A9 signal also present"),
    ("5", "NK candidate",
     "NKG7, GNLY, PRF1; weaker T markers", "Provisional identity"),
    ("6", "IL7R-associated T cells",
     "CD3D, CD3E, IL7R, CCR7", "T-cell subtype unresolved"),
    ("7", "Dendritic cell candidate",
     "FCER1A, CD1C, CLEC10A", "DC subtype unresolved"),
    ("8", "Platelet-associated",
     "PF4, PPBP, GP9; weaker myeloid pattern",
     "Barcode identity and possible mixtures unresolved"),
    ("9", "Uncertain mixed-marker group",
     "T, cytotoxic and other lineage signals",
     "Only seven cells; no definitive identity"),
]

evidence = pd.DataFrame(
    records,
    columns=["cluster", "candidate_label",
             "supporting_evidence", "uncertainty"]
)

labels = evidence.set_index("cluster")["candidate_label"]
clusters = adata.obs["leiden_r10"].astype(str)

assert set(clusters.unique()) == set(labels.index)

adata.obs["manual_annotation"] = (
    clusters.map(labels).astype("category")
)

sizes = clusters.value_counts()
evidence["n_cells"] = evidence["cluster"].map(sizes)

evidence.to_csv(
    "results/tables/pbmc3k_provisional_annotation_evidence.csv",
    index=False
)

adata.obs[["leiden_r05", "leiden_r10", "manual_annotation"]].to_csv(
    "results/tables/pbmc3k_provisional_cell_annotations.csv",
    index_label="cell_barcode"
)

adata.uns["annotation_status"] = (
    "Provisional manual marker-based annotation; "
    "reference comparison and QC sensitivity review pending."
)

adata.write_h5ad(
    "data/processed/pbmc3k_mt8_provisionally_annotated.h5ad",
    compression="gzip"
)

print("\n--- Geçici anotasyon kapsamı ---")
print(
    adata.obs["manual_annotation"]
    .value_counts()
    .to_string()
)

print("\nKanıt tablosu ve hücre etiketleri kaydedildi.")
print("Etiketler geçici; küme numaraları korundu.")
