#' Reparameterize
#' @noRd
repar <- function(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw) {
  sb <- exp(sb0)
  sw <- exp(sw0)
  bb <- Bb0 * (.25 + sb^2) + .5 +
    as.matrix(apply(Zb, 2, function(x) {
      x - mean(x)
    })) %*% as.matrix(Bb0v)
  bw <- Bw0 * (.25 + sw^2) + .5 +
    as.matrix(apply(Zw, 2, function(x) {
      x - mean(x)
    })) %*% as.matrix(Bw0v)
  rho <- (exp(2 * rho0) - 1) / (exp(2 * rho0) + 1)

  c(t(bb), t(bw), sb, sw, rho)
}


check_object <- function(obj, name, msg = "") {
  if (!name %in% names(obj)) {
    if (msg != "") {
      stop(msg)
    } else {
      stop(paste0(name, " is missing."))
    }
  }
}

#' Get data used to create a plot
#'
#' @param x a plot object.
#' @concept visualization
#' @return a data.frame
#' @export
plot_data <- function(x) {
  data <- attr(x, "data")
  return(data)
}
