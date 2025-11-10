test_that("plotPIPHI runs correctly and returns a ggplot object", {
  set.seed(123)

  # Fake data with isoform and cell structure
  pi_mat <- matrix(
    runif(200, 0.05, 0.95),
    nrow = 50,
    dimnames = list(
      paste0("Gene", rep(1:10, each = 5), "..Tx", 1:5),
      paste0("Cell", 1:4)
    )
  )

  iso_mat <- matrix(
    sample(0:20, 200, replace = TRUE),
    nrow = 50,
    dimnames = dimnames(pi_mat)
  )

  # Call the function (should not error)
  expect_silent({
    p <- plotPIPHI(pi_mat, iso_mat, min_gene_count = 10)
  })

  # Check class of output
  expect_s3_class(p, "ggplot")

  # Call with custom thresholds (should still return ggplot)
  expect_silent({
    p2 <- plotPIPHI(pi_mat, iso_mat, pi_min = 0.2, pi_max = 0.8, min_gene_count = 10)
  })
  expect_s3_class(p2, "ggplot")
})
