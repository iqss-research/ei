test_that("bounds works", {
  check <- bounds1(matproii$x, matproii$t, matproii$n)

  expect_gte(min(check), expected = 0)
  expect_lte(max(check), expected = 1)
  expect_true(all(check[, 1] <= check[, 2]))
  expect_true(all(check[, 3] <= check[, 4]))

  test <- summary(check)
  goal <- structure(c("Min.   :0.0000  ", "1st Qu.:0.0000  ", "Median :0.2361  ",
                      "Mean   :0.3045  ", "3rd Qu.:0.5867  ", "Max.   :0.9574  ", "Min.   :0.4042  ",
                      "1st Qu.:1.0000  ", "Median :1.0000  ", "Mean   :0.9935  ", "3rd Qu.:1.0000  ",
                      "Max.   :1.0000  ", "Min.   :0.0000  ", "1st Qu.:0.5016  ", "Median :0.7196  ",
                      "Mean   :0.6669  ", "3rd Qu.:0.8609  ", "Max.   :0.9916  ", "Min.   :0.4250  ",
                      "1st Qu.:0.9165  ", "Median :1.0000  ", "Mean   :0.9331  ", "3rd Qu.:1.0000  ",
                      "Max.   :1.0000  "), .Dim = c(6L, 4L),
                    .Dimnames = list(c("","", "", "", "", ""),
                                     c("    LbetaB", "    UbetaB", "    LbetaW",
                                       "    UbetaW")), class = "table")
  expect_equal(test, goal)

})
