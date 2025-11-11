# Visualize isoform percent inclusion (PI), ratios and distributions on embeddings

This function visualizes the relative usage (PI, Percent Inclusion) of
two isoforms across cells on a low-dimensional embedding stored in a
Seurat object. It computes per-cell PI for each isoform, their ratio and
log2-ratio, adds these values to the object metadata, and displays:

## Usage

``` r
FeaturePlotPI(
  obj,
  iso1,
  iso2,
  reduction,
  assay = "ISO",
  slot = "counts",
  name = NULL,
  add_ratio = TRUE,
  min_total = 0,
  col = "PI",
  bins = 30
)
```

## Arguments

- obj:

  A Seurat object.

- iso1:

  Isoform ID for the first isoform (must be a row name in the assay).

- iso2:

  Isoform ID for the second isoform (must be a row name in the assay).

- reduction:

  Name of the dimensional reduction to use (e.g. `"umap"`).

- assay:

  Assay name containing isoform counts. Default: `"ISO"`.

- slot:

  Slot to use in the assay (e.g. `"counts"`, `"data"`). Default:
  `"counts"`.

- name:

  Optional base name used to construct metadata column names.

- add_ratio:

  Logical; if `TRUE`, add ratio and log2-ratio to metadata. Default:
  `TRUE`.

- min_total:

  Numeric threshold: cells with `iso1 + iso2 < min_total` are set to
  `NA` for PI. Default: `0`.

- col:

  Suffix used for PI-based metadata/label naming. Default: `"PI"`.

- bins:

  Number of bins for the PI histogram. Default: `30`.

## Value

A `patchwork` object combining three feature plots and one histogram.

## Details

- PI of isoform 1 on the embedding,

- PI of isoform 2 on the embedding,

- log2(ratio iso1 / iso2) on the embedding,

- histogram of the PI distribution for isoform 1.

## Examples

``` r
if (FALSE) { # \dontrun{
library(Seurat)
library(Matrix)

set.seed(42)

mat <- Matrix(
  rpois(200, lambda = 10),
  nrow = 2,
  sparse = TRUE,
  dimnames = list(
    c("GeneA..Tx1", "GeneA..Tx2"),
    paste0("Cell", 1:100)
  )
)

obj <- CreateSeuratObject(counts = mat, assay = "ISO")

emb <- matrix(
  rnorm(200),
  ncol = 2,
  dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2"))
)

obj[["umap"]] <- CreateDimReducObject(
  embeddings = emb,
  key = "UMAP_",
  assay = "ISO"
)

FeaturePlotPI(obj, "GeneA..Tx1", "GeneA..Tx2", reduction = "umap")
} # }
```
