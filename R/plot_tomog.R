#' Visualizing EI (tomography plot)
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
#' plot_tomog(ei_res)
#' plot_tomog(ei_res, options = list(linecolor = "betab"))
#' plot_tomog(ei_res, options = list(linecolor = "betaw", category = 5))
#' plot_tomog(ei_res, options = list(points = FALSE, CI = 0.8))
#' @export
plot_tomog <- function(ei.object, options = list(color = TRUE, category = 0, linecolor = "length", CI = NULL, points = FALSE, contour_ML = FALSE, contour_posterior = FALSE)) {
  options <- plot_tomog_options(options)

  p <- plot_tomog_base(ei.object, options)

  if (!is.null(options$CI)) {
    # Adding confidence interval
    p <- plot_add_CI(p, ei.object, options)
  }

  if (options$points) {
    # Adding point estimates
    p <- plot_add_points(p, ei.object, options)
  }

  if (options$contour_ML) {
    # Adding ML contours
    p <- plot_add_contourML(p, ei.object, options)
  }

  if (options$contour_posterior) {
    # Adding posterior contours
    p <- plot_add_contourPost(p, ei.object, options)
  }

  return(p)
}

plot_tomog_options <- function(options) {
  # Check plot_tomog options

  # title
  if (!"title" %in% names(options)) {
    options$title <- "Tomography Plot with the Data"
  }

  # color
  if (!"color" %in% names(options)) {
    options$color <- TRUE
  }
  if (!options$color %in% c(TRUE, FALSE)) {
    stop("`options$color` takes either TRUE or FALSE.")
  }

  # category
  if (!"category" %in% names(options)) {
    options$category <- 0 # continuous
  }

  if (!(0 <= options$category)) {
    stop("Invalid value in `options$category`")
  }

  # linecolor (which axis to use for linecolor)
  if (!"linecolor" %in% names(options)) {
    options$linecolor <- "length"
  }

  if (!options$linecolor %in% c("length", "betab", "betaw")) {
    stop("Invalud value in `options$linecolor`.")
  }

  # breaks: how to break scales
  if (!"breaks" %in% names(options)) {
    options$breaks <- "even"
  }

  if (!options$breaks %in% c("even", "quantile")) {
    stop("Invalud value in `options$breaks`.")
  }

  # Confidence Interval
  if (!"CI" %in% names(options)) {
    options$CI <- NULL
  }

  if (!is.null(options$CI) & is.numeric(options$CI)) {
    if (!(options$CI > 0 & options$CI < 1)) {
      stop("Invalud value in `options$CI`.")
    }
  }

  # Point estimate
  if (!"points" %in% names(options)) {
    options$points <- FALSE
  }

  if (!options$points %in% c(TRUE, FALSE)) {
    stop("`options$points` takes either TRUE or FALSE.")
  }

  # ML contours
  if (!"contour_ML" %in% names(options)) {
    options$contour_ML <- FALSE
  }

  if (!options$contour_ML %in% c(TRUE, FALSE)) {
    stop("`options$contour_ML` takes either TRUE or FALSE.")
  }

  # Posterior contours
  if (!"contour_posterior" %in% names(options)) {
    options$contour_posterior <- FALSE
  }

  if (!options$contour_posterior %in% c(TRUE, FALSE)) {
    stop("`options$contour_posterior` takes either TRUE or FALSE.")
  }

  return(options)
}

plot_tomogd_base <- function(tb, options) {
  p <- tb %>%
    ggplot2::ggplot() +
    {
      if (options$category == 0) ggplot2::guides(color = "none")
    } +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = latex2exp::TeX("$\\beta_B$"), y = latex2exp::TeX("$\\beta_W$")) +
    ggplot2::scale_x_continuous(expand = c(0, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    theme_ei()
  return(p)
}


plot_length_cont <- function(tb, options) {
  # Continuous scale
  if (options$color) {
    p <- plot_tomogd_base(tb, options) +
      ggplot2::geom_segment(aes(
        x = b_bounds[, 1], y = w_bounds[, 1],
        xend = b_bounds[, 2], yend = w_bounds[, 2],
        color = hcl(h = 30, c = 100, l = scale)
      )) +
      ggplot2::scale_color_manual(values = hcl(h = 30, c = 100, l = tb$scale))
  } else {
    p <- plot_tomogd_base(tb, options) +
      ggplot2::geom_segment(aes(
        x = b_bounds[, 1], y = w_bounds[, 1],
        xend = b_bounds[, 2], yend = w_bounds[, 2]
      ),
      color = "yellow"
      )
  }
  return(p)
}

plot_add_scale <- function(tb, options) {
  if (options$linecolor == "length") {
    tb %>%
      mutate(
        length = sqrt(abs(b_bounds[, 1] - b_bounds[, 2])^2 +
          abs(w_bounds[, 1] - w_bounds[, 2])^2),
        scale = ((length - min(length)) / (max(length) - min(length))) * 100
      ) -> tb
  }

  if (options$linecolor == "betab") {
    tb %>%
      mutate(
        length = b_bounds[, 2] - b_bounds[, 1],
        scale = ((length - min(length)) / (max(length) - min(length))) * 100
      ) -> tb
  }

  if (options$linecolor == "betaw") {
    tb %>%
      mutate(
        length = w_bounds[, 1] - w_bounds[, 2],
        scale = ((length - min(length)) / (max(length) - min(length))) * 100
      ) -> tb
  }

  return(tb)
}

strata <- function(tb, q) {
  q_num <- length(q)
  tb$length_cat <- -1
  for (i in 2:q_num) {
    display <- paste0(i - 1, ": [", round(q[i - 1], 2), ", ", round(q[i], 2), "]")
    tb$length_cat <- ifelse(tb$length >= q[i - 1] & tb$length <= q[i], display, tb$length_cat)
  }
  tb$length_cat <- factor(tb$length_cat)
  return(tb)
}

plot_length_cat <- function(tb, options) {
  if (options$breaks == "even") {
    q <- seq(min(tb$length), max(tb$length), length.out = options$category + 1)
  } else {
    q <- quantile(tb$length, prob = seq(0, 1, length.out = options$category + 1))
  }

  p <- plot_tomogd_base(strata(tb, q), options)
  legend_name <- case_when(
    options$linecolor == "length" ~ latex2exp::TeX("Length (line)"),
    options$linecolor == "betaw" ~ latex2exp::TeX("Bound ($\\beta_W$)"),
    options$linecolor == "betab" ~ latex2exp::TeX("Bound ($\\beta_B$)")
  )

  # Categorical scale
  p <- p +
    ggplot2::geom_segment(aes(
      x = b_bounds[, 1], y = w_bounds[, 1],
      xend = b_bounds[, 2], yend = w_bounds[, 2],
      color = .data$length_cat
    )) +
    ggplot2::scale_color_manual(
      values = hcl(h = 30, c = 100, l = seq(1, 80, length.out = options$category + 1)),
      # name = paste0("Length (", legend_name, ")")
      name = legend_name
    )
  return(p)
}

plot_tomog_base <- function(ei.object, options) {
  # Take out the bounds
  bounds <- bounds1(ei.object$x, ei.object$t, ei.object$n)

  tb <- tibble::tibble(
    b_bounds = cbind(bounds[, 1], bounds[, 2]),
    w_bounds = cbind(bounds[, 4], bounds[, 3])
  ) %>% plot_add_scale(options)

  # Plot
  if (options$category == 0) {
    # continuous
    p <- plot_length_cont(tb, options)
  } else {
    # categorical
    p <- plot_length_cat(tb, options)
  }

  return(p)
}

calc_CI <- function(ei.object, alpha) {
  # Only consider precincts that are heterogeneous
  ok <- !is.na(ei.object$betab) & !is.na(ei.object$betaw)
  x <- ei.object$x[ok]
  t <- ei.object$t[ok]
  n <- ei.object$n[ok]
  betabs <- ei.object$betabs[ok, ]
  betaws <- ei.object$betaws[ok, ]
  betabcd <- apply(betabs, 1, function(x) quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
  betawcd <- apply(betaws, 1, function(x) quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
  n <- dim(betabcd)[2]
  return(list(x = x, t = t, n = n, betabcd = betabcd, betawcd = betawcd))
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @importFrom rlang .data
#' @import dplyr
plot_add_CI <- function(p, ei.object, options) {
  res_tomogCI <- calc_CI(ei.object, alpha = 1 - options$CI)
  b <- res_tomogCI$betabcd %>% t()
  w <- res_tomogCI$betawcd
  w <- apply(w, 2, sort, decreasing = TRUE) %>% t()

  bind_cols(
    b %>%
      as_tibble(.name_repair = ~ c("b_start", "b_end")),
    w %>%
      as_tibble(.name_repair = ~ c("w_start", "w_end"))
  ) -> tomo_res_CI

  tomo_res_CI %>%
    mutate(length = sqrt((.data$b_end - .data$b_start)^2 + (.data$w_end - .data$w_start)^2)) %>%
    mutate(scale = ((length - min(length)) / (max(length) - min(length))) * 100) -> tomo_res_CI

  p +
    geom_segment(
      data = tomo_res_CI,
      aes(
        x = .data$b_start,
        y = .data$w_start,
        xend = .data$b_end, yend = .data$w_end
      ),
      color = "red", show.legend = FALSE
    ) -> p
  return(p)
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @importFrom rlang .data
#' @import dplyr
plot_add_points <- function(p, ei.object, options) {
  calc_points <- function(ei.object) {
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

  points <- calc_points(ei.object)

  p +
    ggplot2::geom_point(
      data = points,
      aes(x = .data$betabm, y = .data$betawm),
      colour = "blue"
    ) -> p
  return(p)
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @importFrom rlang .data
#' @import dplyr
plot_add_contourML <- function(p, ei.object, options) {
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  Zb <- ei.object$Zb
  Zw <- ei.object$Zw
  phi <- ei.object$phi
  numb <- dim(Zb)[2]
  numw <- dim(Zw)[2]
  Bb0 <- phi[1]
  Bw0 <- phi[2]
  sb0 <- phi[3]
  sw0 <- phi[4]
  rho0 <- phi[5]
  Bb0v <- phi[6:(5 + numb)]
  Bw0v <- phi[(6 + numb):length(phi)]
  vars <- repar(Bb0, Bw0, sb0, sw0, rho0, Bb0v, Bw0v, Zb, Zw)
  bb <- vars[1:length(x)]
  bw <- vars[(length(x) + 1):(2 * length(x))]
  sb <- vars[2 * length(x) + 1]
  sw <- vars[2 * length(x) + 2]
  rho <- vars[2 * length(x) + 3]

  res_a <- calc_ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
    scale = c(sb, sw),
    centre = c(mean(bb), mean(bw)), level = 0.914
  ) %>% mutate(level = 0.914)

  res_b <- calc_ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
    scale = c(sb, sw),
    centre = c(mean(bb), mean(bw)), level = 0.7
  ) %>% mutate(level = 0.7)

  res_c <- calc_ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
    scale = c(sb, sw),
    centre = c(mean(bb), mean(bw)), level = 0.35
  ) %>% mutate(level = 0.35)

  p <- p +
    geom_path(data = res_a, aes(x = .data$x, y = .data$y), colour = "#16a307", size = 1.5) +
    geom_path(data = res_b, aes(x = .data$x, y = .data$y), colour = "#16a307", size = 1.5) +
    geom_path(data = res_c, aes(x = .data$x, y = .data$y), colour = "#16a307", size = 1.5)

  return(p)
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @importFrom rlang .data
#' @import dplyr
plot_add_contourPost <- function(p, ei.object, options) {
  # Checking the input
  if (!"betabs" %in% names(ei.object)) {
    stop("This plot function requires an ei.sim object.")
  }

  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  psi <- ei.object$psi
  bbp <- psi[, 1:length(x)]
  bwp <- psi[, (length(x) + 1):(2 * length(x))]
  sbp <- psi[, 2 * length(x) + 1]
  swp <- psi[, 2 * length(x) + 2]
  rhop <- psi[, 2 * length(x) + 3]

  res_a <- calc_ellipse(matrix(c(1, rhop, rhop, 1), nrow = 2),
    scale = c(sbp, swp), centre = c(mean(bbp), mean(bwp)), level = 0.914
  ) %>% mutate(level = 0.914)

  res_b <- calc_ellipse(matrix(c(1, rhop, rhop, 1), nrow = 2),
    scale = c(sbp, swp), centre = c(mean(bbp), mean(bwp)), level = 0.7
  ) %>% mutate(level = 0.7)

  res_c <- calc_ellipse(matrix(c(1, rhop, rhop, 1), nrow = 2),
    scale = c(sbp, swp), centre = c(mean(bbp), mean(bwp)), level = 0.35
  ) %>% mutate(level = 0.35)

  p <- p +
    geom_path(data = res_a, aes(x = .data$x, y = .data$y), colour = "#16a307", size = 1.5) +
    geom_path(data = res_b, aes(x = .data$x, y = .data$y), colour = "#16a307", size = 1.5) +
    geom_path(data = res_c, aes(x = .data$x, y = .data$y), colour = "#16a307", size = 1.5)

  return(p)
}
