#' Visualizing EI (xt-plot)
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @return a ggplot object
#' @concept visualization
#' @examples
#' data(matproii)
#' suppressMessages({
#'   ei_res <- ei(formula = t ~ x, total = "n", data = matproii)
#' })
#' plot_xt(ei_res)
#' plot_xt(ei_res, options = list(CI = 0.95, goodman = TRUE))
#' @export
plot_xt <- function(ei.object, options = list(density = TRUE, fit = TRUE, CI = 0.8, goodman = FALSE)) {
  options <- plot_xt_options(options)

  p <- plot_xt_base(ei.object, options)

  if (options$fit) {
    p <- plot_add_fit(p, ei.object, options)
  }

  if (options$goodman) {
    p <- plot_add_goodman(p, ei.object, options)
  }

  return(p)
}

plot_xt_options <- function(options) {
  if (!"density" %in% names(options)) {
    options$density <- TRUE
  }

  if (!options$density %in% c(TRUE, FALSE)) {
    stop("`options$density` takes either TRUE or FALSE.")
  }

  if (!"fit" %in% names(options)) {
    options$fit <- TRUE
  }

  if (!options$fit %in% c(TRUE, FALSE)) {
    stop("`options$fit` takes either TRUE or FALSE.")
  }

  if (!"CI" %in% names(options)) {
    options$CI <- 0.8
  }

  if (!(options$CI > 0 & options$CI < 1)) {
    stop("Use the appropriate value for CI")
  }

  if (!"goodman" %in% names(options)) {
    options$goodman <- FALSE
  }

  if (!options$goodman %in% c(TRUE, FALSE)) {
    stop("`options$goodman` takes either TRUE or FALSE.")
  }

  return(options)
}


plot_xt_base <- function(ei.object, options) {
  points <- tibble::tibble(
    x = ei.object$x,
    t = ei.object$t,
    n = ei.object$n
  )

  minn <- min(points$n)
  maxn <- max(points$n)
  points$scale <- (points$n - minn + 1) / (1 + maxn - minn)

  p <- ggplot2::ggplot(points, aes(x = .data$x, y = .data$t)) +
    {
      if (!options$density) ggplot2::geom_point(shape = 21)
    } +
    {
      if (options$density) ggplot2::geom_point(aes(size = n), shape = 21, show.legend = FALSE)
    } +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = latex2exp::TeX("$X$"), y = latex2exp::TeX("$T$")) +
    theme_ei()

  return(p)
}

#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
plot_add_fit <- function(p, ei.object, options) {
  ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok, ]
  betaws <- ei.object$betaws[ok, ]
  low <- (1 - options$CI) / 2
  up <- 1 - (1 - options$CI) / 2

  minn <- min(n)
  maxn <- max(n)

  x <- seq(0, 1, by = .01)
  betabs <- as.vector(betabs)
  betaws <- as.vector(betaws)
  t <- matrix(ncol = length(x), nrow = length(betabs))
  for (i in 1:length(x)) {
    t[, i] <- betabs * x[i] + betaws * (1 - x[i])
  }

  points <- tibble::tibble(
    x = x,
    et = apply(t, 2, mean),
    lower = apply(t, 2, function(x) quantile(x, probs = c(low))),
    upper = apply(t, 2, function(x) quantile(x, probs = c(up)))
  )

  p <- p +
    ggplot2::geom_ribbon(data = points, aes(y = .data$et, ymin = .data$lower, ymax = .data$upper), fill = "blue", alpha = 0.25) +
    ggplot2::geom_line(data = points, aes(x = .data$x, y = .data$et), color = "red")

  return(p)
}


#' @import ggplot2
#' @import magrittr
#' @importFrom rlang .data
plot_add_goodman <- function(p, ei.object, options) {
  t <- ei.object$t
  x <- ei.object$x
  fit <- lm(t ~ x)
  dat <- tibble::tibble(
    x = seq(0, 1, 0.1)
  )
  dat$y <- stats::predict(fit, newdata = data.frame(x = dat$x))

  p <- p +
    geom_line(
      data = dat, aes(x = .data$x, y = .data$y),
      color = "#16a307"
    )

  return(p)
}
