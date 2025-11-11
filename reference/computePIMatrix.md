# Compute PI (Percent Isoform) Matrix from Isoform Counts

This function calculates per-cell isoform usage (PI) values from an
isoform-level count matrix stored in a Seurat object. For each gene,
isoform counts are divided by the total counts of all isoforms belonging
to that gene within each cell.

## Usage

``` r
computePIMatrix(seurat_obj, assay = "ISO", slot = "counts")
```

## Arguments

- seurat_obj:

  A Seurat object containing isoform-level counts in one assay.

- assay:

  Character string specifying the assay that contains isoform counts.
  Default is `"ISO"`.

- slot:

  Character string indicating which slot of the assay to use (e.g.,
  `"counts"` or `"data"`). Default is `"counts"`.

## Value

A matrix of PI values (isoforms - cells) with `NA` where total gene
counts are zero.

## Details

The function expects rownames in the form `GeneName..TranscriptID`. The
part before the first double dot (`..`) is interpreted as the gene
identifier. PI values are computed per isoform and per cell as: \$\$PI =
count(isoform) / sum(counts of all isoforms from the same gene)\$\$

Cells with zero total counts for a given gene return `NA` PI values for
its isoforms.

## Examples

``` r
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: ‘SeuratObject’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t
iso_names <- c("G1..Tx1", "G1..Tx2", "G2..Tx1", "G2..Tx2")
counts <- matrix(c(10, 5, 3, 7, 0, 4, 2, 6),
                 nrow = 4,
                 dimnames = list(iso_names, paste0("Cell", 1:2)))
seurat_obj <- CreateSeuratObject(counts = counts, assay = "ISO")
#> Warning: Data is of class matrix. Coercing to dgCMatrix.

pi <- computePIMatrix(seurat_obj)
#> Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
#> ℹ Please use the `layer` argument instead.
#> ℹ The deprecated feature was likely used in the scIso package.
#>   Please report the issue at <https://github.com/alihamraoui/scIso/issues>.
#> [computePIMatrix] Extracted assay 'ISO' (4 isoforms - 2 cells).
#> [computePIMatrix] PI matrix computed successfully (dim: 4 - 2).
head(pi)
#> 4 x 2 Matrix of class "dgeMatrix"
#>             Cell1 Cell2
#> G1..Tx1 0.6666667  0.00
#> G1..Tx2 0.3333333  1.00
#> G2..Tx1 0.3000000  0.25
#> G2..Tx2 0.7000000  0.75
```
