import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

adata = sc.read_h5ad(
    "data/processed/pbmc3k_qc_doublet_unfiltered.h5ad"
)

observed = adata.obs["doublet_score"].to_numpy()
simulated = np.asarray(
    adata.uns["scrublet"]["doublet_scores_sim"]
)
threshold = adata.uns["scrublet"].get("threshold")

# Aynı aralıklar, iki dağılımı karşılaştırmayı kolaylaştırır.
bins = np.linspace(
    0, max(observed.max(), simulated.max()), 45
)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

for ax, values, title, color in [
    (axes[0], observed, "Observed cells", "#287C8E"),
    (axes[1], simulated, "Simulated doublets", "#C46A45")
]:
    ax.hist(values, bins=bins, color=color)
    if threshold is not None:
        ax.axvline(
            threshold,
            color="black",
            linestyle="--",
            label=f"Automatic threshold: {threshold:.3f}"
        )
        ax.legend(fontsize=8)

    ax.set_title(title)
    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("Number of profiles")

fig.suptitle("PBMC3k: doublet score distributions")
fig.tight_layout()

fig.savefig(
    "results/figures/pbmc3k_doublet_scores.png",
    dpi=160
)
fig.savefig(
    "results/figures/pbmc3k_doublet_scores.pdf"
)
plt.close(fig)

print("Doublet score figures saved.")
