#' Run (tidy) Ecological Inference Estimation and Simulation
#'
#' @param data data where `x`, `t`, `total`, `Zb`, `Zw` are found
#' @param x <[`data-masking`][dplyr_data_masking]> column of subgroup proportions in data
#' @param t <[`data-masking`][dplyr_data_masking]> column of turnout in data
#' @param n <[`data-masking`][dplyr_data_masking]> column of total in data
#' @param Zb <[`data-masking`][dplyr_tidy_select]> columns of covariates in data
#' @param Zw <[`data-masking`][dplyr_tidy_select]> columns of covariates in data
#' @param id <[`data-masking`][dplyr_data_masking]> column of unique ids in data
#' @param erho The standard deviation of the normal prior on \eqn{\phi_5} for
#' the correlation. Numeric vector, used one at a time, in order. Default `c(.5, 3, 5, .1, 10)`.
#' @param esigma The standard deviation of an underlying normal distribution,
#' from which a half normal is constructed as a prior for both
#' \eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}. Default \eqn{= 0.5}
#' @param ebeta Standard deviation of the "flat normal" prior on
#' \eqn{\breve{B}^b} and \eqn{\breve{B}^w}.  The flat normal prior is uniform
#' within the unit square and dropping outside the square according to the
#' normal distribution.  Set to zero for no prior. Setting to positive values
#' probabilistically keeps the estimated mode within the unit square.
#' Default\eqn{=0.5}
#' @param ealphab cols(Zb) x 2 matrix of means (in the first column) and
#' standard deviations (in the second) of an independent normal prior
#' distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you should
#' probably specify a prior, at least with mean zero and some variance (default
#' is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).
#' @param ealphaw cols(Zw) x 2 matrix of means (in the first column) and
#' standard deviations (in the second) of an independent normal prior
#' distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you should
#' probably specify a prior, at least with mean zero and some variance (default
#' is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).
#' @param truth A length(t) x 2 matrix of the true values of the quantities of
#' interest.
#' @param simulate default = TRUE:see documentation in \code{eiPack} for
#' options for RxC ei.
#' @param covariate see documentation in \code{eiPack} for options for RxC ei.
#' @param lambda1 default = 4:see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param lambda2 default = 2:see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param covariate.prior.list see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param tune.list see documentation in \code{eiPack} for options for RxC ei.
#' @param start.list see documentation in \code{eiPack} for options for RxC ei.
#' @param sample default = 1000
#' @param thin default = 1
#' @param burnin default = 1000
#' @param verbose default = 0:see documentation in \code{eiPack} for options
#' for RxC ei.
#' @param ret.beta default = "r": see documentation in \code{eiPack} for
#' options for RxC ei.
#' @param ret.mcmc default = TRUE: see documentation in \code{eiPack} for
#' options for RxC ei.
#' @param usrfun see documentation in \code{eiPack} for options for RxC ei.
#'
#' @concept tidy
#' @return ei_tbl
#' @export
#'
#' @examples
#' data(sample_ei)
#' dbuf <- ei_(sample_ei, x, t, n)
ei_ <- function(data, x, t, n, Zb = NULL, Zw = NULL, id = NA,
                erho = c(.5, 3, 5, .1, 10), esigma = .5, ebeta = .5, ealphab = NA, ealphaw = NA,
                truth = NA, simulate = TRUE, covariate = NULL, lambda1 = 4,
                lambda2 = 2, covariate.prior.list = NULL, tune.list = NULL,
                start.list = NULL, sample = 1000, thin = 1, burnin = 1000,
                verbose = 0, ret.beta = "r", ret.mcmc = TRUE, usrfun = NULL) {
  x_name <- rlang::as_name(rlang::enquo(x))
  t_name <- rlang::as_name(rlang::enquo(t))
  n_name <- rlang::as_name(rlang::enquo(n))
  x <- rlang::eval_tidy(rlang::enquo(x), data)
  t <- rlang::eval_tidy(rlang::enquo(t), data)
  n <- rlang::eval_tidy(rlang::enquo(n), data)

  Zb <- rlang::eval_tidy(rlang::enquo(Zb), data)
  Zw <- rlang::eval_tidy(rlang::enquo(Zw), data)
  if (is.null(Zb)) Zb <- 1
  if (is.null(Zw)) Zw <- 1

  cli::cli_progress_step("Running 2x2 ei")

  dbuf <- NULL
  i <- 1
  while (i <= length(erho) & is.null(dbuf)) {
    try(
      {
        dbuf <- ei.estimate(t, x, n,
                            id = id,
                            data = data, Zb = Zb, Zw = Zw, erho = erho[i],
                            esigma = esigma, ebeta = ebeta,
                            ealphab = ealphab, ealphaw = ealphaw,
                            truth = truth
        )
      },
      silent = TRUE
    )
    i <- i + 1
  }
  if (is.null(dbuf)) {
    cli::cli_abort("{.fn ei.estimate} did not converge. Try a different value of {.arg erho}.")
  }
  cli::cli_progress_done()
  if (simulate) {
    dbuf.sim <- ei.sim(dbuf)
    cli::cli_progress_done()
    return(new_ei_tbl(
      data,
      x = x_name, t = t_name, n = n_name,
      phi = dbuf.sim$phi,
      hessian = dbuf.sim$hessian, hessianC = dbuf.sim$hessianC, psi = dbuf.sim$psi,
      betab = dbuf.sim$betab, betaw = dbuf.sim$betaw, sbetab = dbuf.sim$sbetab,
      sbetaw = dbuf.sim$sbetaw, betabs = dbuf.sim$betabs, betaws = dbuf.sim$betaws,
      resamp = dbuf.sim$resamp,
      erho = dbuf.sim$erho, esigma = dbuf.sim$esigma, ebeta = dbuf.sim$ebeta,
      ealphab = dbuf.sim$ealphab, ealphaw = dbuf.sim$ealphaw, numb = dbuf.sim$numb,
      Zb = dbuf.sim$Zb, Zw = dbuf.sim$Zw, truth = dbuf.sim$truth, precision = dbuf.sim$precision,
      id = dbuf.sim$id
    ))
  } else {
    return(new_ei_tbl(
      data,
      x = x_name, t = t_name, n = n_name,
      phi = dbuf$phi,
      hessian = dbuf$hessian, hessianC = dbuf$hessianC, psi = dbuf$psi,
      betab = dbuf$betab, betaw = dbuf$betaw, sbetab = dbuf$sbetab,
      sbetaw = dbuf$sbetaw, betabs = dbuf$betabs, betaws = dbuf$betaws,
      resamp = dbuf$resamp,
      erho = dbuf$erho, esigma = dbuf$esigma, ebeta = dbuf$ebeta,
      ealphab = dbuf$ealphab, ealphaw = dbuf$ealphaw, numb = dbuf$numb,
      Zb = dbuf$Zb, Zw = dbuf$Zw, truth = dbuf$truth, precision = dbuf$precision,
      id = dbuf$id
    ))
  }
}


#' Run (tidy) Ecological Inference Estimation
#'
#' @param data data where `x`, `t`, `total`, `Zb`, `Zw` are found
#' @param x <[`data-masking`][dplyr_data_masking]> column of subgroup proportions in data
#' @param t <[`data-masking`][dplyr_data_masking]> column of turnout in data
#' @param n <[`data-masking`][dplyr_data_masking]> column of total in data
#' @param Zb <[`data-masking`][dplyr_tidy_select]> columns of covariates in data
#' @param Zw <[`data-masking`][dplyr_tidy_select]> columns of covariates in data
#' @param id <[`data-masking`][dplyr_data_masking]> column of unique ids in data
#' @param erho The standard deviation of the normal prior on \eqn{\phi_5} for
#' the correlation. Numeric vector, used one at a time, in order. Default `c(.5, 3, 5)`.
#' @param esigma The standard deviation of an underlying normal distribution,
#' from which a half normal is constructed as a prior for both
#' \eqn{\breve{\sigma}_b} and \eqn{\breve{\sigma}_w}. Default \eqn{= 0.5}
#' @param ebeta Standard deviation of the "flat normal" prior on
#' \eqn{\breve{B}^b} and \eqn{\breve{B}^w}.  The flat normal prior is uniform
#' within the unit square and dropping outside the square according to the
#' normal distribution.  Set to zero for no prior. Setting to positive values
#' probabilistically keeps the estimated mode within the unit square.
#' Default\eqn{=0.5}
#' @param ealphab cols(Zb) x 2 matrix of means (in the first column) and
#' standard deviations (in the second) of an independent normal prior
#' distribution on elements of \eqn{\alpha^b}.  If you specify Zb, you should
#' probably specify a prior, at least with mean zero and some variance (default
#' is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^b}).
#' @param ealphaw cols(Zw) x 2 matrix of means (in the first column) and
#' standard deviations (in the second) of an independent normal prior
#' distribution on elements of \eqn{\alpha^w}.  If you specify Zw, you should
#' probably specify a prior, at least with mean zero and some variance (default
#' is no prior).  (See Equation 9.2, page 170, to interpret \eqn{\alpha^w}).
#' @param truth A length(t) x 2 matrix of the true values of the quantities of
#' interest.
#'
#' @concept tidy
#' @return ei_tbl
#' @export
#'
#' @examples
#' data(sample_ei)
#' dbuf <- ei_est(sample_ei, x, t, n)
ei_est <- function(data, t, x, n, id = seq_len(nrow(data)), Zb = NULL, Zw = NULL, erho = .5,
                   esigma = .5, ebeta = .5, ealphab = NA, ealphaw = NA,
                   truth = NA) {
  # set unused args for future?
  Rfun <- 2
  precision <- 4

  x_name <- rlang::as_name(rlang::enquo(x))
  t_name <- rlang::as_name(rlang::enquo(t))
  n_name <- rlang::as_name(rlang::enquo(n))

  x <- rlang::eval_tidy(rlang::enquo(x), data)
  t <- rlang::eval_tidy(rlang::enquo(t), data)
  n <- rlang::eval_tidy(rlang::enquo(n), data)

  Zb <- rlang::eval_tidy(rlang::enquo(Zb), data)
  Zw <- rlang::eval_tidy(rlang::enquo(Zw), data)
  if (is.null(Zb)) Zb <- 1
  if (is.null(Zw)) Zw <- 1

  Zb <- as.matrix(Zb)
  Zw <- as.matrix(Zw)

  # If there are no covariates, run simplest R function
  if (dim(Zb)[1] == 1 & Zb[1, 1] == 1 & dim(Zw)[1] == 1 & Zw[1, 1] == 1) {
    Rfun <- 5
  }
  if (dim(Zb)[1] == 1 & Zb[1, 1] == 1) {
    Zb <- as.matrix(rep(1, length(x)))
  }
  if (dim(Zw)[1] == 1 & Zw[1, 1] == 1) {
    Zw <- as.matrix(rep(1, length(x)))
  }

  # Extract the number of covariates
  numb <- dim(Zb)[2]
  numw <- dim(Zw)[2]

  # Starting values
  start <- c(0, 0, -1.2, -1.2, 0, rep(0, numb + numw))

  cli::cli_alert_info("Maximizing likelihood")

  solution <- optim(start, like,
                    y = t, x = x, n = n, Zb = Zb,
                    Zw = Zw, numb = numb, erho = erho, esigma = esigma,
                    ebeta = ebeta, ealphab = ealphab, ealphaw = ealphaw, Rfun = Rfun,
                    hessian = TRUE,
                    method = "BFGS"
  )

  # Find values of the Hessian that are 0 or 1.
  covs <- as.logical(ifelse(diag(solution$hessian) == 0 |
                              diag(solution$hessian) == 1, 0, 1))
  new_ei_tbl(
    data,
    x = x_name, t = t_name, n = n_name,
    phi = solution$par,
    hessian = solution$hessian, hessianC = solution$hessian[covs, covs],
    erho = erho, esigma = esigma, ebeta = ebeta,
    ealphab = ealphab, ealphaw = ealphaw, numb = numb,
    Zb = Zb, Zw = Zw, truth = truth, precision = precision, covs = covs,
    Rfun = Rfun, id = id
  )
}


#' Run Ecological Inference Simulation
#'
#' @param data an `ei_tbl` object from `ei_est()`
#' @param ndraws integer, default 99. The number of draws.
#' @param nsims integer, default 10. The number of simulations with each draw.
#'
#' @concept tidy
#' @return ei_tbl
#' @export
#'
#' @examples
#' data(sample_ei)
#' dbuf <- ei_est(sample_ei, x, t, n) %>% ei_sim()
ei_sim <- function(data, ndraws = 99, nsims = 100) {
  check_ei_types(data)

  hessian <- attr(data, "hessianC")
  erho <- attr(data, "erho")
  esigma <- attr(data, "esigma")
  ebeta <- attr(data, "ebeta")
  ealphab <- attr(data, "ealphab")
  ealphaw <- attr(data, "ealphaw")
  numb <- attr(data, "numb")
  covs <- attr(data, "covs")
  Rfun <- attr(data, "Rfun")
  x <- data[[attr(data, "x")]]
  t <- data[[attr(data, "t")]]
  n <- data[[attr(data, "n")]]
  Zb <- attr(data, "z_b")
  Zw <- attr(data, "z_w")
  truth <- attr(data, "truth")
  id <- attr(data, "id")
  precision <- attr(data, "precision")
  # Begin Importance Sampling
  cli::cli_progress_step("Beginning importance sampling.", spinner = TRUE)

  keep <- matrix(data = NA, nrow = ndraws, ncol = length(attr(data, "phi")))
  resamp <- 0
  cur_row <- 1 # 99 resamples
  while (cur_row <= ndraws) {
    out_samp <- .samp(t, x, n, Zb, Zw, attr(data, "phi"), hessian, nsims, keep,
                      numb = numb, covs, erho, esigma,
                      ebeta, ealphab, ealphaw, Rfun
    )

    if (!is.null(out_samp)) {
      nro <- nrow(out_samp)
      keep[cur_row:min(cur_row + nro - 1, ndraws), ] <- out_samp[1:min(nro, 1 + ndraws - cur_row), ]
      cur_row <- cur_row + nro
    }
    resamp <- resamp + 1
  }

  # Extract values from importance sampling
  mu <- keep[, 1:2]
  sd <- keep[, 3:4]
  rho <- keep[, 5]
  Bb0v <- keep[, 6:(5 + numb)]
  Bw0v <- keep[, (6 + numb):length(attr(data, "phi"))]
  sd[, 1] <- exp(sd[, 1])
  sd[, 2] <- exp(sd[, 2])

  # Reparamterize
  Zb <- as.matrix(Zb)
  Zw <- as.matrix(Zw)

  mu1 <- mu[, 1] * (.25 + sd[, 1]^2) + .5 + t(as.matrix(apply(
    Zb, 2,
    function(x) x - mean(x)
  )) %*% t(Bb0v))
  mu2 <- mu[, 2] * (.25 + sd[, 2]^2) + .5 + t(as.matrix(apply(
    Zw, 2,
    function(x) x - mean(x)
  )) %*% t(Bw0v))

  # phin <- dmvnorm(psi, par, log=T)
  rho <- (exp(2 * rho) - 1) / (exp(2 * rho) + 1)
  psi <- cbind(mu1, mu2, sd, rho)
  bb <- psi[, 1:length(x)]
  bw <- psi[, (length(x) + 1):(length(x) * 2)]
  sb <- psi[, (length(x) * 2 + 1)]
  sw <- psi[, (length(x) * 2 + 2)]
  rho <- psi[, (length(x) * 2 + 3)]
  omx <- 1 - x
  sbw <- rho * sb * sw
  betab <- matrix(nrow = length(x), ncol = nrow(keep))
  betaw <- matrix(nrow = length(x), ncol = nrow(keep))
  homoindx <- ifelse(x == 0, 1, 0)
  homoindx <- ifelse(x == 1, 2, homoindx)
  enumtol <- .0001
  cT0 <- t < enumtol & homoindx == 0
  cT1 <- t > (1 - enumtol) & homoindx == 0
  ok <- ifelse(homoindx == 0 & cT0 == 0 & cT1 == 0, T, F)
  wh <- homoindx == 1
  bl <- homoindx == 2
  for (i in 1:nrow(keep)) {
    sig2 <- sb[i]^2 * x^2 + sw[i]^2 * omx^2 + sbw[i] * 2 * x * omx
    omega <- sb[i]^2 * x + sbw[i] * omx
    eps <- t - (bb[i, ]) * x - (bw[i, ]) * omx
    mbb <- bb[i, ] + omega / sig2 * eps
    vbb <- sb[i]^2 - (omega^2) / sig2
    vbb <- ifelse(vbb < 1 * 10^-32, .0001, vbb)
    s <- ifelse(vbb >= 0 & vbb != Inf & !is.na(vbb), sqrt(vbb), NaN)
    bounds <- bounds1(x, t, n)
    out <- NULL
    for (j in 1:length(x[ok])) {
      out[ok][j] <- rtnorm(1,
                           mean = mbb[ok][j], sd = s[ok][j],
                           lower = bounds[ok, ][j, 1],
                           upper = bounds[ok, ][j, 2]
      )
    }
    out[wh] <- NA
    out[bl] <- t[bl]
    out[cT1] <- bounds[cT1, 1]
    out[cT0] <- bounds[cT0, 1]
    betab[, i] <- out
  }
  omx <- 1 - x
  for (j in 1:length(x[ok])) {
    betabs <- betab[ok, ][j, ]
    betaw[ok, ][j, ] <- t[ok][j] / omx[ok][j] - betabs * x[ok][j] / omx[ok][j]
  }

  if (sum(wh) > 0) {
    betaw[wh, ] <- as.matrix(rep(1, nrow(keep))) %*% t(as.matrix(t[wh]))
  }

  if (sum(bl) > 0) {
    # betaw[bl, ] <- NA
    betaw[bl, ] <- as.matrix(rep(1, nrow(keep))) %*% t(as.matrix(t[bl]))
  }
  if (sum(cT1) > 0) {
    betaw[cT1, ] <-
      as.matrix(rep(1, nrow(keep))) %*% t(as.matrix(bounds[cT1, 3]))
  }
  if (sum(cT0) > 0) {
    betaw[cT0, ] <-
      as.matrix(rep(1, nrow(keep))) %*% t(as.matrix(bounds[cT0, 3]))
  }

  mbetab <- rowMeans(betab)
  mbetaw <- rowMeans(betaw)
  sdbetab <- apply(betab, 1, sd)
  sdbetaw <- apply(betaw, 1, sd)


  data <- data %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      beta_b = mbetab,
      beta_w = mbetaw,
      s_beta_b = sdbetab,
      s_beta_w = sdbetaw,
      beta_bs = betab,
      beta_ws = betaw,
      z_b = Zb,
      z_w = Zw
    )
  cli::cli_progress_done()

  new_ei_tbl(
    data, attr(data, "x"), attr(data, "t"), attr(data, "n"),
    phi = attr(data, "phi"),
    hessian = attr(data, "hessian"),
    hessianC = hessian, psi = psi,
    betab = "beta_b", betaw = "beta_w", sbetab = "s_beta_b",
    sbetaw = "s_beta_w", betabs = "beta_bs", betaws = "beta_ws",
    resamp = resamp, erho = erho, esigma = esigma,
    ebeta = ebeta, ealphab = ealphab,
    ealphaw = ealphaw, numb = numb,
    Zb = "z_b", Zw = "z_w",
    truth = truth,
    precision = precision, id = id
  )
}

check_ei_types <- function(data) {
  if (!is.null(data) && !inherits(data, "ei_tbl")) {
    cli::cli_abort("{.arg data} must be a {.cls ei_tbl}.")
  }
}
