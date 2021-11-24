#
# For plot_tomog*()
#


# tomogd <- function(x, t, n, title, lci = T) {
#   bounds <- bounds1(x, t, n)
#   bbounds <- cbind(bounds[, 1], bounds[, 2])
#   wbounds <- cbind(bounds[, 4], bounds[, 3])
#   n <- dim(bounds)[1]
#   return(bounds)
# }

bounds <- function(x, t, n) {
  # Migrate bounds() later here?
  #   Original package has `.bounds()` and `bounds1()`
  return(bounds1(x, t, n))
}

calc_ellipse <- function(x, scale = c(1, 1), centre = c(0, 0), level = 0.95,
                         t = sqrt(stats::qchisq(level, 2)), which = c(1, 2), npoints = 350) {
  # From R package `elli[se`
  names <- c("x", "y")
  if (is.matrix(x)) {
    xind <- which[1]
    yind <- which[2]
    r <- x[xind, yind]
    if (missing(scale)) {
      scale <- sqrt(c(x[xind, xind], x[yind, yind]))
      if (scale[1] > 0) r <- r / scale[1]
      if (scale[2] > 0) r <- r / scale[2]
    }
    if (!is.null(dimnames(x)[[1]])) {
      names <- dimnames(x)[[1]][c(xind, yind)]
    }
  } else {
    r <- x
  }
  r <- min(max(r, -1), 1) # clamp to -1..1, in case of rounding errors
  d <- acos(r)
  a <- seq(0, 2 * pi, len = npoints)
  res <- matrix(c(t * scale[1] * cos(a + d / 2) + centre[1], t * scale[2] *
    cos(a - d / 2) + centre[2]), npoints, 2, dimnames = list(
    NULL,
    names
  ))
  return(tibble::as_tibble(res))
}
