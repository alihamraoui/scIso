# Compute Phi (ϕ) from a PI matrix

This function estimates the dispersion parameter \\\phi\\ of a Beta
distribution for each isoform (row) in a matrix of PI (Percent
Inclusion) values. The parameter \\\phi\\ reflects the intercellular
variability of isoform usage across cells.

## Usage

``` r
computePHI(pi_mat)
```

## Arguments

- pi_mat:

  A numeric or sparse matrix of PI values (isoforms × cells).

## Value

A named numeric vector of \\\phi\\ values, one per isoform (row of
`pi_mat`). Rows with insufficient data or invalid parameters return
`NA`.

## Details

For each isoform, \\\phi\\ is estimated from the sample mean (\\m\\) and
variance (\\v\\) of PI values across cells, using the method of moments:
\$\$\phi = 1 / (\alpha + \beta + 1)\$\$ where \\\alpha = m((1 - m)/v -
1)\\ and \\\beta = (1 - m)((1 - m)/v - 1)\\.

Isoforms with fewer than three valid observations or with degenerate
mean or variance values return `NA`.

## Examples

``` r
set.seed(123)
pi_mat <- matrix(runif(20, 0.1, 0.9), nrow = 5,
                 dimnames = list(paste0("Iso", 1:5), paste0("Cell", 1:4)))
phi_vec <- computePHI(pi_mat)
phi_vec
#>      Iso1      Iso2      Iso3      Iso4      Iso5 
#> 0.5245752 0.1284286 0.3442660 0.1357083 0.4487447 
```
