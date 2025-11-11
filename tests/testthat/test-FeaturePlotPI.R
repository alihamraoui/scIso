test_that("vis_iso_psi runs correctly and returns patchwork", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("patchwork")

  set.seed(123)

  # simulate small dataset
  mat <- Matrix::Matrix(
    rpois(400, lambda = 5),
    nrow = 4,
    sparse = TRUE,
    dimnames = list(
      c("Gene1..Tx1", "Gene1..Tx2", "Gene2..Tx1", "Gene2..Tx2"),
      paste0("Cell", 1:100)
    )
  )

  obj <- Seurat::CreateSeuratObject(counts = mat, assay = "ISO")

  # add fake UMAP reduction
  emb <- matrix(
    rnorm(200),
    ncol = 2,
    dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2"))
  )
  obj[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = emb,
    key = "UMAP_",
    assay = "ISO"
  )

  # run function silently
  expect_silent({
    p <- FeaturePlotPI(
      obj,
      iso1 = "Gene1..Tx1",
      iso2 = "Gene1..Tx2",
      reduction = "umap"
    )
  })

  # class check
  expect_true(inherits(p, "patchwork") || inherits(p, "gg"))
})
