#' Visualizing EI
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @export
plot_xt <- function(ei.object, options = list()) {
  options <- plot_xt_options(options)

  p <- plot_xt_base(ei.object, options)

  return(p)
}

plot_xt_options <- function(options) {

  if (! "density" %in% names(options)) {
    options$density <- FALSE
  }

  if (!options$density %in% c(TRUE, FALSE)) {
    stop("`options$density` takes either TRUE or FALSE.")
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

  p <- ggplot2::ggplot(points, aes(x = x, y = t)) +
    {if (!options$density) ggplot2::geom_point(shape =21)} +
    {if (options$density) ggplot2::geom_point(aes(size = n), shape = 21, show.legend = FALSE)} +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = latex2exp::TeX("$X$"), y = latex2exp::TeX("$T$")) +
    theme_ei()

  return(p)
}
