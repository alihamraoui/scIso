# Plot isoform PI or expression along pseudotime with annotation bands

Visualizes isoform-level usage or expression for a given gene along
pseudotime, together with one or multiple categorical annotations
displayed as horizontal bands (e.g. cell type, cell cycle phase). The
two most abundant isoforms of the gene are shown.

## Usage

``` r
IsoformPseudotimePlot(
  seurat_obj,
  gene,
  pseudotime_col,
  band_cols = "cell_type",
  assay = "ISO",
  slot = "counts",
  mode = c("PI", "exp"),
  delimiter = "\\.\\.",
  nbins = 100
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- gene:

  Character scalar. Gene name used as prefix for isoform IDs (e.g.
  "GeneA" for rows like "GeneA..Tx1").

- pseudotime_col:

  Name of the metadata column containing pseudotime values.

- band_cols:

  Character vector of metadata column names to display as annotation
  bands (e.g. `c("cell_type", "Phase")`). Each should be a categorical
  variable.

- assay:

  Assay name containing isoform-level counts. Default: `"ISO"`.

- slot:

  Slot of the assay to use (e.g. `"counts"`, `"data"`). Default:
  `"counts"`.

- mode:

  Either `"PI"` (default) to plot isoform Percent Inclusion (per-cell
  proportions across isoforms of the gene) or `"exp"` to plot raw
  isoform expression.

- delimiter:

  Regex delimiter between gene and isoform ID. Default: `"\\.\\."`.

- nbins:

  Integer. Number of bins along pseudotime for aggregating annotation
  bands. Default: 100.

## Value

A `patchwork` object combining the annotation band plot(s) and the
isoform trend plot.

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(Matrix)

set.seed(123)

iso_names <- c("GeneA..Tx1", "GeneA..Tx2", "GeneA..Tx3")
mat <- Matrix(
  rpois(300, lambda = 10),
  nrow = 3,
  sparse = TRUE,
  dimnames = list(iso_names, paste0("Cell", 1:100))
)

obj <- CreateSeuratObject(counts = mat, assay = "ISO")
obj$pseudotime <- sort(runif(100, 0, 1))
obj$cell_type  <- sample(c("Tumor", "Stromal", "Immune"), 100, TRUE)
obj$Phase      <- sample(c("G1", "S", "G2M"), 100, TRUE)

IsoformPseudotimePlot(
  seurat_obj     = obj,
  gene           = "GeneA",
  pseudotime_col = "pseudotime",
  band_cols      = c("cell_type", "Phase"),
  assay          = "ISO",
  mode           = "PI"
)
} # }
```
