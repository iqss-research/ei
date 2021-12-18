#' Visualizing EI (density)
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @concept visualization
#' @export
plot_density <- function(ei.object, options = list(parameter = "betab")) {
  options <- plot_density_options(options)

  p <- plot_density_base(ei.object, options)

  return(p)
}

plot_density_options <- function(options) {
  if (!"parameter" %in% names(options)) {
    stop("Please specify `parameter` in the option")
  }

  if (!options$parameter %in% c("betab", "betaw")) {
    stop("`options$parameter` takes `betab` or `betaw`")
  }

  return(options)
}


#' @import ggplot2
#' @importFrom rlang .data
plot_density_base <- function(ei.object, options) {
  if (!"betabs" %in% names(ei.object)) {
    stop("This plot function requires an ei.sim object.")
  }

  if (options$parameter == "betab") {
    ok <- !is.na(ei.object$betab)
    betabs <- ei.object$betabs[ok, ]
    beta <- tibble::tibble(x = apply(betabs, 1, mean))
    x_label <- "$\\beta_B$"
  } else {
    ok <- !is.na(ei.object$betaw)
    betaws <- ei.object$betaws[ok, ]
    beta <- tibble::tibble(x = apply(betaws, 1, mean))
    x_label <- "$\\beta_W$"
  }

  ggplot2::ggplot(beta, aes(x = .data$x)) +
    ggplot2::geom_density(color = "#16a307") +
    ggplot2::labs(
      x = latex2exp::TeX(x_label),
      y = latex2exp::TeX("Density acrross precincts")
    ) +
    ggplot2::geom_segment(
      data = beta,
      aes(x = .data$x, y = 0, xend = .data$x, yend = 0.25)
    ) +
    ggplot2::xlim(c(0, 1)) +
    theme_ei() -> p

  return(p)
}
