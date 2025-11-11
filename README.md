# scIso
An R package for quantifying and visualizing **isoform-level heterogeneity** in single-cell and spatial long-read RNA-seq data.  
It provides a reproducible framework to compute isoform Percent Inclusion (PI) and intercellular variability (φ),  
and to visualize isoform usage across cells, conditions, and pseudotime trajectories.

---

## Key features

- **Quantification**
  - `computePIMatrix()` — compute isoform Percent Inclusion (PI) from Seurat assays.  
  - `computePHIMatrix()` — estimate isoform-level variability (φ) using moment-matched Beta parameters.

- **Visualization**
  - `FeaturePlotPI()` — project isoform PI or ratios on embeddings (UMAP, t-SNE, etc.).  
  - `IsoformBoxPlot()` — compare isoform PI/expression across conditions with optional Wilcoxon p-values.  
  - `IsoformPseudotimePlot()` — plot isoform dynamics along pseudotime with annotated metadata bands.

- **Integration**
  - Compatible with `Seurat` v5 objects.

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("alihamraoui/scIso")
```

> Dependencies:  
> `Seurat`, `ggplot2`, `Matrix`, `dplyr`, `tidyr`, `stringr`, `paletteer`, `patchwork`, `RColorBrewer`, `scales`, `rlang`.  
> Optional: `ggpubr (>= 0.6.2)` for statistical comparisons in `IsoformBoxPlot()`.
---

## Development

- Written in **R (≥ 4.3)** using `roxygen2`, `devtools`, and `testthat`.
- All functions include built-in checks for reproducibility and Seurat compatibility.
- Tested with Seurat v5 .
