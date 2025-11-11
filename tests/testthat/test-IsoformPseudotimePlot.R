test_that("IsoformPseudotimePlot runs correctly and returns a patchwork object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("patchwork")

  set.seed(42)

  # --- simulate Seurat object ---
  iso_names <- c("GeneA..Tx1", "GeneA..Tx2", "GeneA..Tx3")
  mat <- Matrix::Matrix(
    rpois(300, lambda = 10),
    nrow = 3,
    sparse = TRUE,
    dimnames = list(iso_names, paste0("Cell", 1:100))
  )

  obj <- Seurat::CreateSeuratObject(counts = mat, assay = "ISO")

  # metadata
  obj$pseudotime <- sort(runif(100))
  obj$cell_type  <- sample(c("Tumor", "Stromal", "Immune"), 100, TRUE)
  obj$Phase      <- sample(c("G1", "S", "G2M"), 100, TRUE)

  # --- test function ---
  result <- IsoformPseudotimePlot(
    seurat_obj     = obj,
    gene           = "GeneA",
    pseudotime_col = "pseudotime",
    band_cols      = c("cell_type", "Phase"),
    assay          = "ISO",
    mode           = "PI"
  )

  # --- check result type ---
  expect_s3_class(result, "patchwork")
})
