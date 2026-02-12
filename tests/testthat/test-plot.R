# Test for the plot functions

data(matproii)
truth <- cbind(matproii$tb, matproii$tw)
dbuf_v <- ei(formula = t ~ x, total = "n", data = matproii, truth = truth)

data(RxCdata)
formula <- cbind(turnout, noturnout) ~ cbind(white, black, hisp)
dbuf <- ei(formula, data = RxCdata)

test_that("`plot_tomog()`", {
  p <- plot_tomog(dbuf_v)
  expect_s3_class(p, "gg")

  p <- plot_tomog(dbuf_v, options = list(contour_ML = TRUE))
  expect_s3_class(p, "gg")

  p <- plot_tomog(dbuf_v,
    options = list(
      category = 5, linecolor = "betab",
      points = TRUE, CI = 0.95, contour_ML = TRUE
    )
  )
  expect_s3_class(p, "gg")

  p <- plot_tomog(dbuf_v, options = list(contour_posterior = TRUE))
  expect_s3_class(p, "gg")
})

test_that("`plot_density()`", {
  p <- plot_density(dbuf_v)
  expect_s3_class(p, "gg")
  p <- plot_density(dbuf_v, options = list(parameter = "betab"))
  expect_s3_class(p, "gg")
  p <- plot_density(dbuf_v, options = list(parameter = "betaw"))
  expect_s3_class(p, "gg")
})

test_that("`plot_xt()`", {
  p <- plot_xt(dbuf_v)
  expect_s3_class(p, "gg")
  p <- plot_xt(dbuf_v, options = list(density = TRUE))
  expect_s3_class(p, "gg")
  p <- plot_xt(dbuf_v, options = list(fit = TRUE))
  expect_s3_class(p, "gg")
  p <- plot_xt(dbuf_v, options = list(density = TRUE, fit = TRUE, CI = 0.95))
  expect_s3_class(p, "gg")
  p <- plot_xt(dbuf_v, options = list(goodman = TRUE))
  expect_s3_class(p, "gg")
  p <- plot_xt(dbuf_v, options = list(
    density = TRUE, fit = TRUE,
    CI = 0.80, goodman = TRUE
  ))
  expect_s3_class(p, "gg")
})

test_that("`plot_sims()`", {
  p <- plot_sims(dbuf_v)
  expect_s3_class(p, "gg")
})

test_that("`plot_bound()`", {
  p <- plot_bound(dbuf_v, options = list(parameter = "betab"))
  expect_s3_class(p, "gg")
  p <- plot_bound(dbuf_v, options = list(parameter = "betaw"))
  expect_s3_class(p, "gg")

  p <- plot_bound(dbuf)
  expect_s3_class(p, "gg")
})

test_that("`plot_tomogRxC()`", {
  data(RxCdata)
  formula <- cbind(turnout, noturnout) ~ cbind(white, black, hisp)
  p <- plot_tomogRxC(formula, RxCdata)
  expect_s3_class(p, "gg")
})
