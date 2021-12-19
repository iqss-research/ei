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
  expect_equal(dbuf_summary$Resamp, 33, tolerance = 1)

  expect_equal(dbuf_summary[[6]][1, 1], 1.227206, tolerance = 0.001)
  expect_equal(dbuf_summary[[6]][2, 5], 0.6450481, tolerance = 0.01)

  expect_equal(dbuf_summary[[7]][2], 1.762994, tolerance = 0.01)

  expect_equal(dbuf_summary[[8]][3], 0.01229979, tolerance = 0.1)

  expect_equal(dbuf_summary[[9]][1, 1], 0.2152043, tolerance = 0.00001)
  expect_equal(dbuf_summary[[9]][2, 2], 0.9507709, tolerance = 0.00001)

  expect_equal(dbuf_summary[[10]][1, 1], 0.8056159, tolerance = 0.001)
  expect_equal(dbuf_summary[[10]][2, 2], 0.001275387, tolerance = 0.01)
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

test_that("`ei.estimate()` runs", {
  suppressMessages({
    dbuf <- ei(formula = t ~ x, total = "n", data = matproii[1:50, ], simulate = FALSE)
  })

  expect_equal(names(dbuf), c(
    "phi", "hessian", "hessianC", "erho",
    "esigma", "ebeta", "ealphab", "ealphaw", "numb",
    "x", "t", "n", "Zb", "Zw",
    "truth", "precision", "covs", "Rfun", "id"
  ))
})
