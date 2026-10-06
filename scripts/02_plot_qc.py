from pathlib import Path
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

tab = pd.read_csv(
    "results/tables/pbmc3k_qc_before_filtering.csv",
    index_col=0
)

out = Path("results/figures")
out.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(13, 8))

metrics = [
    ("total_counts", "Total UMI counts"),
    ("n_genes_by_counts", "Detected genes"),
    ("pct_counts_mt", "Mitochondrial counts (%)")
]

for ax, (column, label) in zip(axes[0], metrics):
    ax.hist(tab[column], bins=50, color="#287C8E")
    ax.set_xlabel(label)
    ax.set_ylabel("Number of cells")

points = axes[1, 0].scatter(
    tab["total_counts"],
    tab["n_genes_by_counts"],
    c=tab["pct_counts_mt"],
    cmap="viridis",
    s=7,
    alpha=0.7
)
axes[1, 0].set_xlabel("Total UMI counts")
axes[1, 0].set_ylabel("Detected genes")
fig.colorbar(
    points, ax=axes[1, 0],
    label="Mitochondrial counts (%)"
)

axes[1, 1].scatter(
    tab["total_counts"], tab["pct_counts_mt"],
    s=7, alpha=0.5, color="#287C8E"
)
axes[1, 1].set_xlabel("Total UMI counts")
axes[1, 1].set_ylabel("Mitochondrial counts (%)")

axes[1, 2].scatter(
    tab["n_genes_by_counts"], tab["pct_counts_mt"],
    s=7, alpha=0.5, color="#287C8E"
)
axes[1, 2].set_xlabel("Detected genes")
axes[1, 2].set_ylabel("Mitochondrial counts (%)")

fig.suptitle("PBMC3k: QC before filtering")
fig.tight_layout()

fig.savefig(out / "pbmc3k_qc_before_filtering.png", dpi=160)
fig.savefig(out / "pbmc3k_qc_before_filtering.pdf")
plt.close(fig)

print("QC figures saved.")
