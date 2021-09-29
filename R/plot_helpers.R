#
# For plot_tomog*()
#

tomog80CI <- function(ei.object) {
  if (!("betabs" %in% names(ei.object))) {
    message("Error: This plot function requires an ei.sim object.")
  }
  if ("betabs" %in% names(ei.object)) {
    # Only consider precincts that are heterogeneous
    ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
    x <- ei.object$x[ok]
    t <- ei.object$t[ok]
    n <- ei.object$n[ok]
    betabs <- ei.object$betabs[ok, ]
    betaws <- ei.object$betaws[ok, ]
    betabcd <- apply(betabs, 1, function(x) quantile(x, probs = c(.1, .9)))
    betawcd <- apply(betaws, 1, function(x) quantile(x, probs = c(.1, .9)))
    n <- dim(betabcd)[2]
    return(list(x = x, t = t, n = n, betabcd = betabcd, betawcd = betawcd))
  }
}


tomog95CI <- function(ei.object) {
  if (!("betabs" %in% names(ei.object))) {
    message("Error: This plot function requires an ei.sim object.")
  }
  if ("betabs" %in% names(ei.object)) {
    ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
    x <- ei.object$x[ok]
    t <- ei.object$t[ok]
    n <- ei.object$n[ok]
    betabs <- ei.object$betabs[ok, ]
    betaws <- ei.object$betaws[ok, ]
    betabcd <- apply(betabs, 1, function(x) {
      quantile(x,
               probs = c(.025, .975)
      )
    })
    betawcd <- apply(betaws, 1, function(x) {
      quantile(x, probs = c(.025, .975))
    })
    n <- dim(betabcd)[2]
    return(list(x = x, t = t, n = n, betabcd = betabcd, betawcd = betawcd))
  }
}


tomogE <- function(ei.object) {
  if (!("betabs" %in% names(ei.object))) {
    message("Error: This plot function requires an ei.sim object.")
  }
  if ("betabs" %in% names(ei.object)) {
    ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
    x <- ei.object$x[ok]
    t <- ei.object$t[ok]
    n <- ei.object$n[ok]
    betabs <- ei.object$betabs[ok, ]
    betaws <- ei.object$betaws[ok, ]
    betabm <- apply(betabs, 1, mean)
    betawm <- apply(betaws, 1, mean)
    return(tibble::tibble(betabm = betabm, betawm = betawm))
  }
}



tomogd <- function(x, t, n, title, lci = T) {
  bounds <- bounds1(x, t, n)
  bbounds <- cbind(bounds[, 1], bounds[, 2])
  wbounds <- cbind(bounds[, 4], bounds[, 3])
  n <- dim(bounds)[1]
  return(bounds)
}

bounds <- function(x, t, n) {
  # Migrate bounds() later here?
  #   Original package has `.bounds()` and `bounds1()`
  return(bounds1(x, t, n))
}
