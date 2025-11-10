test_that("computePIMatrix computes correct PI values", {
  iso_names <- c("G1..Tx1", "G1..Tx2", "G2..Tx1", "G2..Tx2")
  counts <- matrix(
    c(10, 5, 3, 7,   # Cell1
      0,  4, 2, 6),  # Cell2
    nrow = 4,
    dimnames = list(iso_names, paste0("Cell", 1:2))
  )
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, assay = "ISO")

  pi <- computePIMatrix(seurat_obj)

  expected <- matrix(
    c(
      0.6666667, 0.0,   # G1..Tx1
      0.3333333, 1.0,   # G1..Tx2
      0.3000000, 0.25,  # G2..Tx1
      0.7000000, 0.75   # G2..Tx2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(iso_names, paste0("Cell", 1:2))
  )

  # 1) Values match
  expect_equal(
    round(as.matrix(pi), 6),
    round(expected, 6)
  )

  # 2) For each gene, PI sums to 1 across its isoforms (per cell)
  gene_names <- sub("\\.\\..*$", "", rownames(pi))
  gene_split <- split(seq_along(gene_names), gene_names)

  for (idx in gene_split) {
    sums <- Matrix::colSums(pi[idx, , drop = FALSE])
    expect_true(all(abs(sums - 1) < 1e-6))
  }
})
