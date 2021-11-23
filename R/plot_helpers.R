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
                         t = sqrt(qchisq(level, 2)), which = c(1, 2), npoints = 350) {
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




# plot_tomog <- function(ei.object, title = "Tomography Plot with the Data", lci = TRUE) {
#   plot_tomogd(ei.object$x, ei.object$t, ei.object$n, title, lci)
# }

# plot_tomogd <- function(x, t, n, title, lci = TRUE) {
#   # Take out the bounds
#   bounds <- bounds1(x, t, n)

#   tb <- tibble::tibble(
#     b_bounds = cbind(bounds[, 1], bounds[, 2]),
#     w_bounds = cbind(bounds[, 4], bounds[, 3]),
#     length = sqrt(abs(b_bounds[, 1] - b_bounds[, 2])^2 +
#       abs(w_bounds[, 1] - w_bounds[, 2])^2),
#     scale = ((length - min(length)) / (max(length) - min(length))) * 100
#   )

#   # Plot
#   p <- tb %>%
#     ggplot2::ggplot() +
#     ggplot2::guides(color = "none") +
#     ggplot2::coord_fixed() +
#     ggplot2::labs(x = latex2exp::TeX("$\\beta_B$"), y = latex2exp::TeX("$\\beta_W$")) +
#     theme_ei()

#   if (lci) {
#     p <- p +
#       ggplot2::geom_segment(aes(
#         x = b_bounds[, 1], y = w_bounds[, 1],
#         xend = b_bounds[, 2], yend = w_bounds[, 2],
#         color = hcl(h = 30, c = 100, l = scale)
#       )) +
#       ggplot2::scale_color_manual(values = hcl(h = 30, c = 100, l = scale))
#   } else {
#     p <- p +
#       ggplot2::geom_segment(aes(
#         x = b_bounds[, 1], y = w_bounds[, 1],
#         xend = b_bounds[, 2], yend = w_bounds[, 2]
#       ),
#       color = "yellow"
#       )
#   }

#   return(p)
# }


# plot_tomogl <- function(ei.object, lci = TRUE) {
#   x <- ei.object$x
#   t <- ei.object$t
#   n <- ei.object$n
#   Zb <- ei.object$Zb
#   Zw <- ei.object$Zw
#   phi <- ei.object$phi
#   p <- plot_tomogd(x, t, n, "Tomography Plot with ML Contours", lci = lci)
#   numb <- dim(Zb)[2]
#   numw <- dim(Zw)[2]
#   Bb0 <- phi[1]
#   Bw0 <- phi[2]
#   sb0 <- phi[3]
#   sw0 <- phi[4]
#   rho0 <- phi[5]
#   Bb0v <- phi[6:(5 + numb)]
#   Bw0v <- phi[(6 + numb):length(phi)]
#   vars <- repar(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw)
#   bb <- vars[1:length(x)]
#   bw <- vars[(length(x) + 1):(2 * length(x))]
#   sb <- vars[2 * length(x) + 1]
#   sw <- vars[2 * length(x) + 2]
#   rho <- vars[2 * length(x) + 3]
#   .tomog3 <- function(bb, bw, sb, sw, rho) {
#     lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
#       scale = c(sb, sw),
#       centre = c(mean(bb), mean(bw)), level = .914
#     ), col = "blue", lwd = 4)
#     lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
#       scale = c(sb, sw),
#       centre = c(mean(bb), mean(bw)), level = .35
#     ), col = "red", lwd = 4)
#     points(mean(bb), mean(bw), col = "pink", pch = 15)
#   }

#   .tomog3(bb, bw, sb, sw, rho)
# }

# .tomog3 <- function(bb, bw, sb, sw, rho) {
#   lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
#     scale = c(sb, sw), centre = c(mean(bb), mean(bw)), level = .914
#   ),
#   col = "blue", lwd = 4
#   )
#   lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
#     scale = c(sb, sw),
#     centre = c(mean(bb), mean(bw)), level = .35
#   ), col = "red", lwd = 4)
#   points(mean(bb), mean(bw), col = "pink", pch = 15)
# }
