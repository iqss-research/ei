#' Reparameterize
#'
#' @param Bb0 TODO
#' @param Bw0 TODO
#' @param sb0 TODO
#' @param sw0 TODO
#' @param rho0 TODO
#' @param Bb0v TODO
#' @param Bw0v TODO
#' @param Zb TODO
#' @param Zw TODO
#'
#' @return TODO
#' @noRd
repar <- function(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw) {
  sb <- exp(sb0)
  sw <- exp(sw0)
  bb <- Bb0 * (.25 + sb^2) + .5 +
    as.matrix(apply(Zb, 2, function(x) {x - mean(x)})) %*% as.matrix(Bb0v)
  bw <- Bw0 * (.25 + sw^2) + .5 +
    as.matrix(apply(Zw, 2, function(x) {x - mean(x)})) %*% as.matrix(Bw0v)
  rho <- (exp(2 * rho0) - 1) / (exp(2 * rho0) + 1)

  c(t(bb), t(bw), sb, sw, rho)
}
