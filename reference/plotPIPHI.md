# Plot mean PI vs Phi for isoforms

This function summarizes isoform-level PI values and visualizes the
relationship between mean PI and Phi across isoforms. Isoforms can be
filtered based on their mean PI and total gene-level support before
plotting.

## Usage

``` r
plotPIPHI(pi_mat, iso_mat, pi_min = 0.3, pi_max = 0.7, min_gene_count = 1000)
```

## Arguments

- pi_mat:

  A numeric or sparse matrix of PI values (isoforms × cells).

- iso_mat:

  A numeric or sparse matrix of isoform-level counts (same dimensions
  and rownames as `pi_mat`), used to derive gene-level coverage.

- pi_min:

  Minimum mean PI threshold to retain isoforms (default: 0.3).

- pi_max:

  Maximum mean PI threshold to retain isoforms (default: 0.7).

- min_gene_count:

  Minimum total gene count required to include an isoform in the plot
  (default: 3500).

## Value

A `ggplot` object visualizing mean PI vs Phi, colored by gene-level
count.

## Details

For each isoform, the function computes:

- mean PI across cells;

- Phi using
  [`computePHI`](https://alihamraoui.github.io/scIso/reference/computePHI.md)
  (inter-cell variability);

- number of cells with non-zero counts;

- total gene-level count (sum across all isoforms of a gene).

Isoforms are filtered according to `pi_min`, `pi_max`, and
`min_gene_count`, and only those passing all filters are plotted.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
pi_mat <- matrix(runif(200, 0.05, 0.95), nrow = 50,
                 dimnames = list(paste0("Gene", rep(1:10, each = 5), "..Tx", 1:5),
                                 paste0("Cell", 1:4)))
iso_mat <- matrix(sample(0:20, 200, replace = TRUE), nrow = 50,
                  dimnames = dimnames(pi_mat))

# Default thresholds
plotPIPHI(pi_mat, iso_mat)

# Custom thresholds
plotPIPHI(pi_mat, iso_mat, pi_min = 0.2, pi_max = 0.8, min_gene_count = 10)
} # }
```
