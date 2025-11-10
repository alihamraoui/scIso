test_that("computePHI computes valid phi values", {
  # Matrix with variation across cells (should give finite phi)
  set.seed(1)
  pi_mat <- matrix(
    runif(20, 0.1, 0.9),
    nrow = 5,
    dimnames = list(paste0("Iso", 1:5), paste0("Cell", 1:4))
  )

  phi_vec <- computePHI(pi_mat)

  # 1) Output type and size
  expect_type(phi_vec, "double")
  expect_length(phi_vec, nrow(pi_mat))
  expect_named(phi_vec)

  # 2) Phi should be between 0 and 1 for valid isoforms
  expect_true(all(phi_vec[!is.na(phi_vec)] > 0))
  expect_true(all(phi_vec[!is.na(phi_vec)] < 1))

  # 3) Constant PI per isoform -> phi = NA
  const_mat <- matrix(0.5, nrow = 3, ncol = ncol(pi_mat),
                      dimnames = list(paste0("Const", 1:3), colnames(pi_mat)))
  phi_const <- computePHI(const_mat)
  expect_true(all(is.na(phi_const)))

  # 4) Mix of variable and constant isoforms
  mix_mat <- rbind(pi_mat, const_mat)
  phi_mix <- computePHI(mix_mat)

  # Expect same length as number of rows
  expect_length(phi_mix, nrow(mix_mat))

  # Expect NAs for constant isoforms only
  expect_true(sum(is.na(phi_mix)) >= nrow(const_mat))
})
