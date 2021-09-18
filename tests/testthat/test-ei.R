# Test for the `ei()` function

data(matproii)
data(RxCdata)

test_that("`ei()` 2x2", {
  set.seed(225)
  suppressMessages({
    dbuf <- ei(formula = t ~ x, total = "n", data = matproii[1:50, ])
  })
  dbuf_summary <- summary(dbuf)

  expect_equal(dbuf_summary$Erho, 0.5)
  expect_equal(dbuf_summary$Esigma, 0.5)
  expect_equal(dbuf_summary$Ebeta, 0.5)
  expect_equal(dbuf_summary$N, 50)
  expect_equal(dbuf_summary$Resamp, 39)

  expect_equal(dbuf_summary[[6]][1, 1], 1.581521, tolerance = 0.00001)
  expect_equal(dbuf_summary[[6]][2, 5], 1.124315, tolerance = 0.00001)

  expect_equal(dbuf_summary[[7]][2], 1.397145, tolerance = 0.00001)

  expect_equal(dbuf_summary[[8]][3], 0.1889316, tolerance = 0.00001)

  expect_equal(dbuf_summary[[9]][1, 1], 0.2152043, tolerance = 0.00001)
  expect_equal(dbuf_summary[[9]][2, 2], 0.9507709, tolerance = 0.00001)

  expect_equal(dbuf_summary[[10]][1, 1], 0.7029558, tolerance = 0.00001)
  expect_equal(dbuf_summary[[10]][2, 2], 0.01764645, tolerance = 0.00001)
})

test_that("`ei()` RxC", {
  set.seed(225)
  suppressMessages({
    dbuf <- ei(
      formula = cbind(turnout, noturnout) ~ cbind(white, black, hisp),
      data = RxCdata
    )
  })
  out <- dbuf$draws$Beta[, "beta.white.turnout.3"]
  out_summary <- summary(out)

  expect_equal(out_summary$nchain, 1)
  expect_equal(out_summary$thin, 1)
  expect_equal(out_summary$start, 1)
  expect_equal(out_summary$end, 1000)

  expect_equal(as.numeric(out_summary$quantiles[1]), 0.3345407, tolerance = 0.00001)
  expect_equal(as.numeric(out_summary$quantiles[3]), 0.4154687, tolerance = 0.00001)

  expect_equal(as.numeric(out_summary$statistics[1]), 0.4159187, tolerance = 0.00001)
  expect_equal(as.numeric(out_summary$statistics[3]), 0.001258285, tolerance = 0.00001)
})
