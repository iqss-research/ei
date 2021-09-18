test_that("bounds works", {
  check <- bounds1(matproii$x, matproii$t, matproii$n)

  expect_gte(min(check), expected = 0)
  expect_lte(max(check), expected = 1)
  expect_true(all(check[, 1] <= check[, 2]))
  expect_true(all(check[, 3] <= check[, 4]))
})
