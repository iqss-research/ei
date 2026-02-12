#' Simulate EI Solution via Importance Sampling
#'
#' @param ei.object \code{ei} object
#' @param ndraws integer. The number of draws. Default is 99.
#' @param nsims integer. The number of simulations within each draw. Default is 100.
#'
#' @author Gary King <<email: king@@harvard.edu>> and Molly Roberts <<email:
#' molly.e.roberts@@gmail.com>>
#' @references Gary King (1997). A Solution to the Ecological Inference
#' Problem.  Princeton: Princeton University Press.
#'
#' @export
#' @return `ei.sim` object
#'
#' @examples
#' data(sample_ei)
#' form <- t ~ x
#' ei_obj <- ei(form, total = "n", data = sample_ei, simulate = FALSE)
#' sims <- ei.sim(ei_obj)
ei.sim <- function(ei.object, ndraws = 99, nsims = 100) {
  # Check the output
  if ("eiRxC" %in% class(ei.object)) {
    cli::cli_abort("`ei.sim()` does not support the RxC case.")
  }

  # Preparation
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
  max_resamp <- ndraws * 50
  cur_row <- 1 # 99 resamples
  while (cur_row <= ndraws) {
    out_samp <- .samp(t, x, n, Zb, Zw, ei.object$phi, hessian, nsims, keep,
      numb = numb, covs, erho, esigma,
      ebeta, ealphab, ealphaw, Rfun
    )

    if (!is.null(out_samp)) {
      nro <- nrow(out_samp)
      if (nro > 0) {
        keep[cur_row:min(cur_row + nro - 1, ndraws), ] <- out_samp[1:min(nro, 1 + ndraws - cur_row), ]
        cur_row <- cur_row + nro
      }
    }
    resamp <- resamp + 1
    if (resamp > max_resamp) {
      cli::cli_abort("Importance sampling did not converge after {max_resamp} attempts.")
    }
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
  np <- length(x)
  nd <- nrow(keep)
  betab <- matrix(nrow = np, ncol = nd)
  betaw <- matrix(nrow = np, ncol = nd)

  # Precompute categories using direct logical indexing (avoid ifelse)
  homoindx <- integer(np)
  homoindx[x == 0] <- 1L
  homoindx[x == 1] <- 2L
  enumtol <- .0001
  hetero <- homoindx == 0L
  cT0 <- hetero & t < enumtol
  cT1 <- hetero & t > (1 - enumtol)
  ok <- hetero & !cT0 & !cT1
  wh <- homoindx == 1L
  bl <- homoindx == 2L

  # Hoist bounds computation out of the loop (constant across iterations)
  bounds <- bounds1(x, t, n)
  ok_idx <- which(ok)
  n_ok <- length(ok_idx)
  bounds_ok_lo <- bounds[ok_idx, 1]
  bounds_ok_hi <- bounds[ok_idx, 2]

  # Precompute fixed values for special categories
  out_template <- numeric(np)
  out_template[wh] <- NA
  out_template[bl] <- t[bl]
  out_template[cT1] <- bounds[cT1, 1]
  out_template[cT0] <- bounds[cT0, 1]

  x2 <- x^2
  omx2 <- omx^2
  x_omx_2 <- 2 * x * omx

  for (i in 1:nd) {
    sig2 <- sb[i]^2 * x2 + sw[i]^2 * omx2 + sbw[i] * x_omx_2
    omega <- sb[i]^2 * x + sbw[i] * omx
    eps <- t - bb[i, ] * x - bw[i, ] * omx
    mbb <- bb[i, ] + omega / sig2 * eps
    vbb_raw <- sb[i]^2 - (omega^2) / sig2
    vbb_raw[vbb_raw < 1e-32] <- .0001
    s_all <- sqrt(vbb_raw)
    s_all[!is.finite(s_all)] <- NaN

    # Vectorized rtnorm call for all ok indices at once
    out <- out_template
    out[ok_idx] <- rtnorm(n_ok,
      mean = mbb[ok_idx], sd = s_all[ok_idx],
      lower = bounds_ok_lo, upper = bounds_ok_hi
    )
    betab[, i] <- out
  }

  # Vectorized betaw computation (no loop needed)
  if (n_ok > 0) {
    t_over_omx <- t[ok_idx] / omx[ok_idx]
    x_over_omx <- x[ok_idx] / omx[ok_idx]
    betaw[ok_idx, ] <- t_over_omx - betab[ok_idx, , drop = FALSE] * x_over_omx
  }

  if (sum(wh) > 0) {
    betaw[wh, ] <- rep(1, nd) %o% t[wh]
  }
  if (sum(bl) > 0) {
    betaw[bl, ] <- rep(1, nd) %o% t[bl]
  }
  if (sum(cT1) > 0) {
    betaw[cT1, ] <- rep(1, nd) %o% bounds[cT1, 3]
  }
  if (sum(cT0) > 0) {
    betaw[cT0, ] <- rep(1, nd) %o% bounds[cT0, 3]
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
  class(output) <- c("ei", "ei2X2", "eisim")
  output
}
