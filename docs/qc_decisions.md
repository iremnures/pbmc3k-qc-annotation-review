# PBMC3k QC decisions

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
