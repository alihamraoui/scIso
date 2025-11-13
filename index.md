# scIso

[![](logo.svg)](https://alihamraoui.github.io/scIso/images/schema.png)

An R package for quantifying and visualizing **isoform-level
heterogeneity** in single-cell and spatial long-read RNA-seq data.  
It provides a reproducible framework to compute isoform usage
intercellular variability and to visualize it across conditions, and
pseudotime trajectories.

------------------------------------------------------------------------

## Key features

- **Quantification**
  - [`computePIMatrix()`](https://alihamraoui.github.io/scIso/reference/computePIMatrix.md)
    — compute isoform Percent Inclusion (PI) from Seurat assays.  
  - `computePHIMatrix()` — estimate isoform-level variability (φ) using
    moment-matched Beta parameters.
- **Visualization**
  - [`FeaturePlotPI()`](https://alihamraoui.github.io/scIso/reference/FeaturePlotPI.md)
    — project isoform PI or ratios on embeddings (UMAP, t-SNE, etc.).  
  - [`IsoformBoxPlot()`](https://alihamraoui.github.io/scIso/reference/IsoformBoxPlot.md)
    — compare isoform PI/expression across conditions with optional
    Wilcoxon p-values.  
  - [`IsoformPseudotimePlot()`](https://alihamraoui.github.io/scIso/reference/IsoformPseudotimePlot.md)
    — plot isoform dynamics along pseudotime with annotated metadata
    bands.
- **Integration**
  - Use `Seurat` v5 objects and compatible with most scRNA-seq long-read
    pipelines (Sicelore, wf-single-cell…).

------------------------------------------------------------------------

## Installation

You can install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("alihamraoui/scIso")
```

> Dependencies:  
> `Seurat`, `ggplot2`, `Matrix`, `dplyr`, `tidyr`, `stringr`,
> `paletteer`, `patchwork`, `RColorBrewer`, `scales`, `rlang`.  
> Optional: `ggpubr (>= 0.6.2)` for statistical comparisons in
> [`IsoformBoxPlot()`](https://alihamraoui.github.io/scIso/reference/IsoformBoxPlot.md).
> —

## Development

- Written in **R (≥ 4.3)** using `roxygen2`, `devtools`, and `testthat`.
- All functions include built-in checks for reproducibility and Seurat
  compatibility.
- Tested with Seurat v5 .
