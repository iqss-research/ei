test_that("ei.sim() works", {
  set.seed(02138)
  suppressMessages({
    ei_obj <- ei(formula = t ~ x, total = "n", data = matproii[201:268, ], simulate = FALSE)
    sims <- ei.sim(ei_obj)
  })

  expect_equal(sims$phi,
               c(0.832247337021409, 0.71886770438913, -2.21162496298499, -1.85832226505417,
                 0.137461849932152, 0, 0),
               tolerance = 0.001)


  expect_equal(names(sims),
               c("phi", "hessian", "hessianC", "psi", "betab", "betaw", "sbetab",
                 "sbetaw", "betabs", "betaws", "resamp", "erho",
                 "esigma", "ebeta", "ealphab", "ealphaw", "numb",
                 "x", "t", "n", "Zb", "Zw", "truth", "precision", "id"))

  expect_equal(mean(sims$sbetab), 0.1105136, tolerance = 0.001)
  expect_equal(mean(sims$sbetaw), 0.05280262, tolerance = 0.001)

})
