test_that("strsplit1 spluts a string", {
  expect_equal(strsplit1("a,b,c", ","), c("a","b","c"))
})

