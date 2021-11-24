#' Visualizing EI
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @export
plot_bound <- function(ei.object, options = list()) {
  options <- plot_bound_options(options)

  p <- plot_bound_base(ei.object, options)

  return(p)
}

plot_bound_options <- function(options) {
  if (!"parameter" %in% names(options)) {
    stop("Please specify `parameter` in the option")
  }

  if (!options$parameter %in% c("betab", "betaw")) {
    stop("`options$parameter` takes `betab` or `betaw`")
  }

  return(options)
}


#' @importFrom rlang .data
plot_bound_base <- function(ei.object, options) {
  if (options$parameter == "betab") {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    truebb <- ei.object$truth[,1]
    bounds <- bounds1(x, t, n)

    res <- tibble::tibble(
      x = x,
      lower = bounds[, 1],
      upper = bounds[, 2],
      true = truebb
    )

    title <- "Aggregation Bias for $\\beta_B$"
    ylab <- "True $\\beta_B$"
  } else if (options$parameter == "betaw") {
    x <- ei.object$x
    t <- ei.object$t
    n <- ei.object$n
    truebw <- ei.object$truth[,2]
    bounds <- bounds1(x, t, n)

    res <- tibble::tibble(
      x = x,
      lower = bounds[, 3],
      upper = bounds[, 4],
      true = truebw
    )

    title <- "Aggregation Bias for $\\beta_W$"
    ylab <- "True $\\beta_W$"
  } else {
    stop("Invalid argument in `options$parameter`")
  }

  fit <- lm(true ~ x, data = res)
  dat <- tibble::tibble(
    x = seq(0, 1, 0.1)
  )
  dat$y <- predict(fit, newdata = data.frame(x = dat$x))
  dat$lower <- 0
  dat$upper <- 0

  p <- ggplot(res, aes(x = .data$x, y = .data$true, ymin = .data$lower, ymax = .data$upper)) +
    geom_point() +
    geom_linerange() +
    geom_line(
      data = dat, aes(x = .data$x, y = .data$y),
      linetype = "dashed"
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = latex2exp::TeX("$X$"), y = latex2exp::TeX(ylab)) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_ei()

  return(p)
}

