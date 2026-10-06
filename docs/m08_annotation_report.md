# M08: QC and annotation review

## Dataset and data representation

PBMC3k contains 2,700 cell-associated barcodes from a single donor.
The downloaded matrix is already cell-called and does not include
the full unfiltered droplet population.

In the analysis objects:
- `layers["counts"]` stores raw counts for the retained genes.
- `X` stores log1p expression after normalization to 10,000 counts per cell.
- `obsm["X_pca"]` and `obsm["X_umap"]` store reduced representations.

QC metrics were calculated before gene filtering.

## QC decisions

Scrublet identified 32 candidate doublets; these predictions are
not definitive doublet diagnoses. The same candidates were excluded
from all three mitochondrial-threshold variants.

Retained cells:
- mt5: 2,611
- mt8: 2,660
- mt10: 2,662

The mt8 variant was the working baseline, not a universal optimal
threshold. Low gene counts alone were not used to exclude cells.
Ambient RNA correction was not performed.

## Selected annotation: platelet-associated group

Baseline Leiden resolution 1.0 cluster 8 contained 15 cells.

Supporting evidence:
- Coordinated PF4, PPBP and GP9 expression.
- Weaker LST1, FCGR3A and MS4A7 expression than the monocyte group.
- CellTypist predicted Megakaryocytes/platelets for 11 of 15 cells.

At resolution 0.5, these cells were grouped with predominantly
monocyte-associated cells. Increasing resolution separated the
platelet-associated expression pattern. Cluster membership therefore
does not necessarily imply a single cell identity.

Low RNA complexity alone was insufficient grounds for exclusion.

## QC sensitivity

Normalization, variable-gene selection, PCA and Leiden clustering
were recomputed for mt5 and mt10 with the same parameter settings.

All 15 baseline platelet-associated cells remained together in
both variants. PF4, PPBP and GP9 expression patterns were preserved.
Their numerical cluster label changed from 8 in mt5 to 9 in mt10;
this is not evidence of an identity change.

Adjusted Rand index on shared cells:
- mt5 versus mt8: 0.7568
- mt10 versus mt8: 0.8791

Broad marker patterns were largely preserved, while some T-cell
and monocyte cluster boundaries changed. ARI measures partition
agreement, not annotation accuracy.

## Remaining uncertainty and interpretation limits

Platelet-associated is a supported descriptive annotation.
It does not exclude ambient RNA or cell-mixture contributions in
individual barcodes. Cluster averages cannot establish marker
coexpression in every cell.

Dendritic-cell and cytotoxic T-cell subtypes remain provisional.
The seven-cell mixed-marker group remains uncertain.

CellTypist predictions provide reference-model evidence, not
experimental validation. Best-match prediction assigns a label
even when the biological identity is uncertain.

Cluster-marker comparisons are exploratory identity evidence,
not donor-level differential-expression tests between conditions.
This single-donor dataset cannot support population-level
condition comparisons.

UMAP separation alone does not establish identity, function or
mechanism.
