#' Simulate EI Solution via Importance Sampling
#'
#' @param ei.object \code{ei} object
#' @param ndraws integer. The number of draws.
#' @param nsims integer. The number of simulations within each draw.
#'
#' @author Gary King <<email: king@@harvard.edu>> and Molly Roberts <<email:
#' molly.e.roberts@@gmail.com>>
#' @references Gary King (1997). A Solution to the Ecological Inference
#' Problem.  Princeton: Princeton University Press.
#'
#' @export
#' @return TODO
#'
#' @examples
#' data(sample_ei)
#' form <- t ~ x
#' ei_obj <- ei(form, total = "n", data = sample_ei, simulate = FALSE)
#' sims <- ei.sim(ei_obj)
ei.sim <- function(ei.object, ndraws = 99, nsims = 100) {
  hessian <- ei.object$hessianC
  erho <- ei.object$erho
  esigma <- ei.object$esigma
  ebeta <- ei.object$ebeta
  ealphab <- ei.object$ealphab
  ealphaw <- ei.object$ealphaw
  numb <- ei.object$numb
  covs <- ei.object$covs
  Rfun <- ei.object$Rfun
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  Zb <- ei.object$Zb
  Zw <- ei.object$Zw
  truth <- ei.object$truth
  id <- ei.object$id
  precision <- ei.object$precision
  # Begin Importance Sampling
  cli::cli_progress_step("Beginning importance sampling.", spinner = TRUE)

  keep <- matrix(data = NA, nrow = ndraws, ncol = length(ei.object$phi))
  resamp <- 0
  cur_row <- 1 # 99 resamples
  while (cur_row <= ndraws) {
    out_samp <- .samp(t, x, n, Zb, Zw, ei.object$phi, hessian, nsims, keep,
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
  Bw0v <- keep[, (6 + numb):length(ei.object$phi)]
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
    #betaw[bl, ] <- NA
    betaw[bl,] <- as.matrix(rep(1,nrow(keep)))%*%t(as.matrix(t[bl]))
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
  output <- list(
    phi = ei.object$phi,
    hessian = ei.object$hessian,
    hessianC = hessian, psi = psi,
    betab = mbetab, betaw = mbetaw,
    sbetab = sdbetab, sbetaw = sdbetaw,
    betabs = betab, betaws = betaw,
    resamp = resamp, erho = erho, esigma = esigma,
    ebeta = ebeta, ealphab = ealphab,
    ealphaw = ealphaw, numb = numb,
    x = x, t = t, n = n,
    Zb = Zb, Zw = Zw,
    truth = truth,
    precision = precision, id = id
  )

  cli::cli_progress_done()
  class(output) <- "ei"
  return(output)
}
