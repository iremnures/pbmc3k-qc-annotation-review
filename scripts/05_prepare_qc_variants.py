from pathlib import Path
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_qc_doublet_unfiltered.h5ad"
)

assert adata.obs["predicted_doublet"].notna().all()

doublet = adata.obs["predicted_doublet"].astype(bool)
summary = []

variants = [
    ("mt5", 5),
    ("mt8", 8),
    ("mt10", 10),
]

for name, threshold in variants:
    high_mt = adata.obs["pct_counts_mt"] > threshold
    remove = high_mt | doublet
    keep = ~remove

    # .copy(): seçilen hücreler için bağımsız bir nesne oluşturur.
    selected = adata[keep, :].copy()

    selected.uns["qc_variant"] = {
        "name": name,
        "maximum_mitochondrial_percentage": threshold,
        "exclude_predicted_doublets": True,
        "doublet_scores_recomputed": False,
        "additional_gene_count_filter": False,
    }

    path = f"data/processed/pbmc3k_{name}.h5ad"
    selected.write_h5ad(path, compression="gzip")

    summary.append({
        "variant": name,
        "starting_cells": adata.n_obs,
        "mt_above_threshold": int(high_mt.sum()),
        "predicted_doublets": int(doublet.sum()),
        "both_flags": int((high_mt & doublet).sum()),
        "removed_unique_cells": int(remove.sum()),
        "retained_cells": selected.n_obs,
        "retained_genes": selected.n_vars,
    })

    del selected

tab = pd.DataFrame(summary)
tab.to_csv(
    "results/tables/pbmc3k_qc_variant_summary.csv",
    index=False
)

Path("docs/qc_decisions.md").write_text(
    """# PBMC3k QC decisions

This is a learning analysis of a single-donor dataset.

Working analysis: mitochondrial count percentage <= 8%.
Sensitivity comparisons: <= 5% and <= 10%.

Predicted doublets from the initial Scrublet run are excluded
in all three variants. Scores and the automatic threshold are
held fixed to isolate the mitochondrial-threshold comparison.

These are provisional analysis choices, not universal QC rules.
No additional gene-count or upper-count exclusion is applied.
Percentile flags are retained in the review table.

Ambient RNA correction has not been performed.
The starting matrix already contains selected cell barcodes.
Unfiltered source files are preserved.
""",
    encoding="utf-8"
)

print("\n--- QC comparison variants ---")
print(tab.to_string(index=False))
print("\nWorking analysis: mt8")
print("Three count-based objects saved; no normalization applied.")
