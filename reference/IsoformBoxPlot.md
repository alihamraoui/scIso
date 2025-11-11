# Boxplot of isoform expression across conditions

This function visualizes isoform expression (or PI if provided upstream)
across a categorical variable (e.g. cell cycle phase), optionally
faceted by another metadata column, and performs per-panel statistical
comparison between isoforms using Wilcoxon tests.

## Usage

``` r
IsoformBoxPlot(
  seurat_obj,
  isoforms,
  group_col,
  facet_col = NULL,
  gene_expr_col = NULL,
  assay = "ISO",
  slot = "data",
  filter_detected = TRUE,
  add_pvalues = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- isoforms:

  Character vector of isoform feature names to include (e.g.
  c("Cdkn2c..ENSMUST00000097921", "Cdkn2c..ENSMUST00000063531")).

- group_col:

  Metadata column used on the x-axis (e.g. "Phase").

- facet_col:

  Optional metadata column used to facet the plot (e.g.
  "tumor_cell_type"). If `NULL`, no faceting is applied.

- gene_expr_col:

  Optional gene-level expression column name (e.g. "Cdkn2c"). If not
  `NULL`, only cells with non-missing values are kept.

- assay:

  Assay from which to extract isoform data. Default: "ISO".

- slot:

  Slot to use ("data", "counts", etc.). Default: "data".

- filter_detected:

  Logical; if TRUE, keep only cells where all selected isoforms have
  expression \> 0. Default: TRUE.

- add_pvalues:

  Logical; if TRUE, add p-values (Wilcoxon test) using
  ggpubr::stat_compare_means() if ggpubr \>= 0.6.2 is available.
  Default: TRUE.

## Value

A `ggplot` object with grouped boxplots, optional facets, and optional
significance annotations.

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(Matrix)

set.seed(42)

# --- simulate Seurat object ---
iso_names <- c("GeneX..Tx1", "GeneX..Tx2")
mat <- Matrix(
  rpois(200, lambda = 10),
  nrow = 2,
  sparse = TRUE,
  dimnames = list(iso_names, paste0("Cell", 1:100))
)

obj <- CreateSeuratObject(counts = mat, assay = "ISO")

# --- add metadata ---
obj$Phase <- sample(c("G1", "S", "G2M"), 100, TRUE)
obj$tumor_cell_type <- sample(c("Schwann", "Fibroblast", "Immune"), 100, TRUE)

# --- test function ---
IsoformBoxPlot(
  seurat_obj  = obj,
  isoforms    = iso_names,
  group_col   = "Phase",
  facet_col   = "tumor_cell_type",
  assay       = "ISO",
  slot        = "counts",
  filter_detected = TRUE,
  add_pvalues = TRUE
)
} # }
```
