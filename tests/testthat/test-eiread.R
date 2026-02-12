set.seed(1)
data(sample_ei)
formula <- t ~ x
dbuf <- ei(formula = formula, total = "n", data = sample_ei)

n_obs <- nrow(sample_ei)

# --- Simulation-dependent outputs: structural & approximate checks ---

test_that("eiread betab works", {
  out <- eiread(dbuf, "betab")
  expect_length(out, n_obs)
  expect_true(all(out >= 0 & out <= 1))
  expect_equal(mean(out), 0.20, tolerance = 0.1)
})

test_that("eiread betaw works", {
  out <- eiread(dbuf, "betaw")
  expect_length(out, n_obs)
  expect_true(all(out >= 0 & out <= 1))
  expect_equal(mean(out), 0.72, tolerance = 0.05)
})

test_that("eiread phi works", {
  out <- eiread(dbuf, "phi")
  expect_length(out, 7)
  expect_equal(out[6], 0)
  expect_equal(out[7], 0)
  expect_true(all(is.finite(out)))
})

test_that("eiread sbetab works", {
  out <- eiread(dbuf, "sbetab")
  expect_length(out, n_obs)
  expect_true(all(out >= 0))
  expect_equal(mean(out), 0.08, tolerance = 0.15)
})

test_that("eiread sbetaw works", {
  out <- eiread(dbuf, "sbetaw")
  expect_length(out, n_obs)
  expect_true(all(out >= 0))
  expect_equal(mean(out), 0.07, tolerance = 0.15)
})

test_that("eiread psisims works", {
  out <- eiread(dbuf, "psisims")
  expect_true(is.matrix(out))
  expect_true(ncol(out) >= n_obs)
  expect_true(nrow(out) > 0)
  col_means <- colMeans(out)
  expect_true(all(is.finite(col_means)))
})

test_that("eiread CI80b works", {
  out <- eiread(dbuf, "CI80b")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(n_obs, 2L))
  expect_equal(colnames(out), c("lwr", "upr"))
  expect_true(all(out >= 0 & out <= 1))
  expect_true(all(out[, "lwr"] <= out[, "upr"]))
})

test_that("eiread CI80w works", {
  out <- eiread(dbuf, "CI80w")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(n_obs, 2L))
  expect_equal(colnames(out), c("lwr", "upr"))
  expect_true(all(out >= 0 & out <= 1))
  expect_true(all(out[, "lwr"] <= out[, "upr"]))
})

test_that("eiread aggs works", {
  out <- eiread(dbuf, "aggs")
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 2L)
  expect_equal(colnames(out), c("Bbgg", "Bwgg"))
  expect_true(nrow(out) > 0)
  expect_true(all(out >= 0 & out <= 1))
  expect_equal(mean(out[, "Bbgg"]), 0.20, tolerance = 0.1)
  expect_equal(mean(out[, "Bwgg"]), 0.72, tolerance = 0.1)
})

test_that("eiread maggs works", {
  out <- eiread(dbuf, "maggs")
  expect_length(out, 4)
  expect_equal(out[1], 0.20, tolerance = 0.1)
  expect_equal(out[2], 0.72, tolerance = 0.1)
  expect_true(all(out[3:4] >= 0))
})

test_that("eiread VCaggs works", {
  out <- eiread(dbuf, "VCaggs")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 2L))
  expect_equal(out[1, 2], out[2, 1])
})

# --- Deterministic outputs: bounds, abounds, goodman ---

test_that("eiread bounds works", {
  out <- eiread(dbuf, "bounds")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(n_obs, 4L))
  expect_equal(colnames(out), c("LbetaB", "UbetaB", "LbetaW", "UbetaW"))
  expect_true(all(out >= 0 & out <= 1))
  expect_true(all(out[, "LbetaB"] <= out[, "UbetaB"]))
  expect_true(all(out[, "LbetaW"] <= out[, "UbetaW"]))
})

test_that("eiread abounds works", {
  out <- eiread(dbuf, "abounds")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 2L))
  expect_true(all(out >= 0 & out <= 1))
  expect_true(out[1, 1] <= out[1, 2])
  expect_true(out[2, 1] <= out[2, 2])
  expect_equal(out, structure(c(0.0954, 0.6089, 0.3979, 0.8031), .Dim = c(2L, 2L)))
})

test_that("eiread goodman works", {
  out <- eiread(dbuf, "goodman")
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(2L, 4L))
  expect_equal(rownames(out), c("BetaB", "BetaW"))
  expect_equal(colnames(out), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  expect_equal(
    out,
    structure(c(0.1956, 0.7228, 0.0275, 0.023, 7.0995, 31.3996, 0, 0),
      .Dim = c(2L, 4L),
      .Dimnames = list(
        c("BetaB", "BetaW"),
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      )
    )
  )
})
