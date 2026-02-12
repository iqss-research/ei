test_that("bounds works", {
  check <- bounds1(matproii$x, matproii$t, matproii$n)

  expect_gte(min(check), expected = 0)
  expect_lte(max(check), expected = 1)
  expect_true(all(check[, 1] <= check[, 2]))
  expect_true(all(check[, 3] <= check[, 4]))

  expect_equal(colnames(check), c("LbetaB", "UbetaB", "LbetaW", "UbetaW"))
  expect_equal(nrow(check), nrow(matproii))

  expect_equal(mean(check[, "LbetaB"]), 0.3045, tolerance = 0.01)
  expect_equal(mean(check[, "UbetaB"]), 0.9935, tolerance = 0.01)
  expect_equal(mean(check[, "LbetaW"]), 0.6669, tolerance = 0.01)
  expect_equal(mean(check[, "UbetaW"]), 0.9331, tolerance = 0.01)
})
