test_that("IsoformBoxPlot runs and returns a ggplot object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")

  set.seed(42)

  # --- simulate Seurat object ---
  iso_names <- c("GeneX..Tx1", "GeneX..Tx2")
  mat <- Matrix::Matrix(
    rpois(200, lambda = 10),
    nrow = 2,
    sparse = TRUE,
    dimnames = list(iso_names, paste0("Cell", 1:100))
  )

  obj <- Seurat::CreateSeuratObject(counts = mat, assay = "ISO")

  # --- add metadata ---
  obj$Phase <- sample(c("G1", "S", "G2M"), 100, TRUE)
  obj$tumor_cell_type <- sample(c("Schwann", "Fibroblast", "Immune"), 100, TRUE)

  # --- run function ---
  p <- IsoformBoxPlot(
    seurat_obj  = obj,
    isoforms    = iso_names,
    group_col   = "Phase",
    facet_col   = "tumor_cell_type",
    assay       = "ISO",
    slot        = "count",
    filter_detected = TRUE,
    add_pvalues = TRUE
  )

  # --- check ---
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl("Isoform", names(p$mapping))) || is.list(p$layers))
})
