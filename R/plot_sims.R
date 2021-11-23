#' Visualizing EI
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @export
plot_sims <- function(ei.object, options = list()) {
  options <- plot_sims_options(options)

  p <- plot_sims_base(ei.object, options)

  return(p)
}

plot_sims_options <- function(options) {
  return(options)
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @importFrom rlang .data
plot_sims_base <- function(ei.object, options) {
  if (!("betabs" %in% names(ei.object))) {
    stop("This plot function requires an ei.sim object.")
  }
  ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
  betabs <- ei.object$betabs[ok, ]
  betaws <- ei.object$betaws[ok, ]

  suppressMessages({
    bind_cols(
      as_tibble(betabs, column_name = paste0("V", 1:ncol(betabs))) %>%
        tidyr::pivot_longer(everything()) %>%
        rename(x = value),
      as_tibble(betaws, column_name = paste0("V", 1:ncol(betaws))) %>%
        tidyr::pivot_longer(everything()) %>%
        rename(y = value)
    ) -> dat
    dat$color <- runif(nrow(dat), 26, 51)
  })

  p <- ggplot(dat, aes(x = x, y = y, colour = factor(color))) +
    geom_point(size = 0.25, show.legend = FALSE) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = latex2exp::TeX("$\\beta_B$ simulations"), y = latex2exp::TeX("$\\beta_W$ simulations")) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  return(p)
}
