from pathlib import Path
import gc
import inspect

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

sc.settings.n_jobs = 1
out = Path("results/tables/qc_sensitivity")
out.mkdir(parents=True, exist_ok=True)

baseline = sc.read_h5ad(
    "data/processed/pbmc3k_mt8_reference_compared.h5ad"
)

markers = [
    "CD3D", "CD3E", "IL7R", "CCR7", "CD8A", "CD8B",
    "NKG7", "GNLY", "PRF1", "MS4A1", "CD79A", "CD79B",
    "LYZ", "LST1", "TYROBP", "FCN1", "S100A9",
    "FCGR3A", "MS4A7", "FCER1A", "CD1C", "CLEC10A",
    "CLEC4C", "PF4", "PPBP", "GP9"
]

def seeded(function, obj, **kwargs):
    parameters = inspect.signature(function).parameters
    if "rng" in parameters:
        kwargs["rng"] = 0
    elif "random_state" in parameters:
        kwargs["random_state"] = 0
    return function(obj, **kwargs)

summary = []

for variant in ["mt5", "mt10"]:
    print(f"\n--- Starting {variant} ---", flush=True)
    data = sc.read_h5ad(f"data/processed/pbmc3k_{variant}.h5ad")
    data.X = data.layers["counts"].copy()

    sc.pp.filter_genes(data, min_cells=3)
    sc.pp.normalize_total(data, target_sum=10000)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(
        data, flavor="seurat", n_top_genes=2000, subset=False
    )

    work = ad.AnnData(
        X=data[:, data.var["highly_variable"]].X.copy(),
        obs=pd.DataFrame(index=data.obs_names.copy())
    )
    sc.pp.scale(work, zero_center=True, max_value=10)
    seeded(sc.pp.pca, work, n_comps=50, svd_solver="arpack")
    data.obsm["X_pca"] = work.obsm["X_pca"].copy()
    del work
    gc.collect()

    seeded(
        sc.pp.neighbors, data,
        n_neighbors=15, n_pcs=20, use_rep="X_pca",
        transformer="sklearn", metric="euclidean"
    )
    seeded(
        sc.tl.leiden, data,
        resolution=1.0, flavor="igraph", directed=False,
        n_iterations=2, key_added="qc_cluster"
    )

    common = baseline.obs_names.intersection(data.obs_names)
    old = baseline.obs.loc[common, "leiden_r10"].astype(str)
    new = data.obs.loc[common, "qc_cluster"].astype(str)
    ari = adjusted_rand_score(old, new)

    cross = pd.crosstab(old, new)
    cross.to_csv(out / f"{variant}_cluster_overlap.csv")

    coverage = baseline.obs[["manual_annotation"]].copy()
    coverage["retained"] = coverage.index.isin(data.obs_names)
    coverage = coverage.groupby(
        "manual_annotation", observed=True
    ).agg(
        baseline_cells=("retained", "size"),
        retained_cells=("retained", "sum")
    )
    coverage["retained_percent"] = (
        100 * coverage["retained_cells"] / coverage["baseline_cells"]
    )
    coverage.to_csv(out / f"{variant}_baseline_label_retention.csv")

    present = [gene for gene in markers if gene in data.var_names]
    means = pd.DataFrame.from_dict(
        {
            cluster: np.asarray(
                data[data.obs["qc_cluster"] == cluster, present]
                .X.mean(axis=0)
            ).ravel()
            for cluster in data.obs["qc_cluster"].cat.categories
        },
        orient="index",
        columns=present
    )
    means.index.name = "qc_cluster"
    means.to_csv(out / f"{variant}_marker_means.csv")

    # Keep baseline labels as comparison information, not new annotations.
    data.obs["baseline_manual_annotation"] = (
        baseline.obs["manual_annotation"]
        .astype(str).reindex(data.obs_names).fillna("Not in mt8 baseline")
        .astype("category")
    )
    data.write_h5ad(
        f"data/processed/pbmc3k_{variant}_stability.h5ad",
        compression="gzip"
    )

    summary.append({
        "variant": variant,
        "cells": data.n_obs,
        "genes": data.n_vars,
        "clusters": data.obs["qc_cluster"].nunique(),
        "shared_cells_with_mt8": len(common),
        "ARI_vs_mt8": ari
    })

    print("\nBaseline cluster x new cluster:")
    print(cross.to_string())
    print("\nBaseline annotation coverage:")
    print(coverage.round(1).to_string())

    del data
    gc.collect()

summary = pd.DataFrame(summary)
summary.to_csv(out / "qc_clustering_stability_summary.csv", index=False)

print("\n--- QC clustering comparison ---")
print(summary.round(4).to_string(index=False))
print("\nSaved. Original mt8 results were not overwritten.")
