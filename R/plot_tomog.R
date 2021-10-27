#' @export
plot_tomog <- function(ei.object, title = "Tomography Plot with the Data", options = list()) {
  options <- plot_tomg_options(options)
  plot_tomogd(ei.object$x, ei.object$t, ei.object$n, title, options)
}

plot_tomg_options <- function(options) {
  # Check plot_tomog options

  # color
  if (!"color" %in% names(options)) {
    options$color <- TRUE
  }

  if (!options$color %in% c(TRUE, FALSE)) {
    stop("`options$color` is either TRUE or FALSE.")
  }

  # category
  if (!"category" %in% names(options)) {
    options$category <- 0 # continuous
  }

  if (!(0 <= options$category)) {
    stop("Invalid value in `options$category`")
  }

  # scale (which axis to use for scale)
  if (!"scale" %in% names(options)) {
    options$scale <- "length"
  }

  if (!options$scale %in% c("length", "betab", "betaw")) {
    stop("Invalud value in `options$scale`.")
  }

  return(options)
}

plot_tomogd_base <- function(tb, options) {
  p <- tb %>%
    ggplot2::ggplot() +
    {
      if (options$category == 0) ggplot2::guides(color = "none")
    } +
    ggplot2::coord_fixed() +
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
  if (options$scale == "length") {
    tb %>%
      mutate(
        length = sqrt(abs(b_bounds[, 1] - b_bounds[, 2])^2 +
          abs(w_bounds[, 1] - w_bounds[, 2])^2),
        scale = ((length - min(length)) / (max(length) - min(length))) * 100
      ) -> tb
  }

  if (options$scale == "betab") {
    tb %>%
      mutate(
        length = b_bounds[, 2] - b_bounds[, 1],
        scale = ((length - min(length)) / (max(length) - min(length))) * 100
      ) -> tb
  }

  if (options$scale == "betaw") {
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
    display <- paste0(i-1, ": [", round(q[i-1], 2), ", ", round(q[i], 2), "]")
    tb$length_cat <- ifelse(tb$length >= q[i-1] & tb$length <= q[i], display, tb$length_cat)
  }
  tb$length_cat <- factor(tb$length_cat)
  return(tb)
}

plot_length_cat <- function(tb, options) {
  q <- quantile(tb$length, prob = seq(0, 1, length.out = options$category + 1))
  p <- plot_tomogd_base(strata(tb, q), options)
  legend_name <- case_when(
    options$scale == "length" ~ latex2exp::TeX("Length (line)"),
    options$scale == "betaw" ~ latex2exp::TeX("Length ($\\beta_W$)"),
    options$scale == "betab" ~ latex2exp::TeX("Length ($\\beta_B$)")
  )

  # Categorical scale
  p <- p +
  ggplot2::geom_segment(aes(
    x = b_bounds[, 1], y = w_bounds[, 1],
    xend = b_bounds[, 2], yend = w_bounds[, 2],
    color = length_cat
  )) +
  ggplot2::scale_color_manual(
    values = hcl(h = 30, c = 100, l = seq(1, 80, length.out = options$category + 1)),
    # name = paste0("Length (", legend_name, ")")
    name = legend_name
  )
  return(p)
}

plot_tomogd <- function(x, t, n, title, options) {
  # Take out the bounds
  bounds <- bounds(x, t, n)

  tb <- tibble::tibble(
    b_bounds = cbind(bounds[, 1], bounds[, 2]),
    w_bounds = cbind(bounds[, 4], bounds[, 3])
  ) %>% plot_add_scale(., options)


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

#' @export
plot_tomogl <- function(ei.object, lci = TRUE) {
  x <- ei.object$x
  t <- ei.object$t
  n <- ei.object$n
  Zb <- ei.object$Zb
  Zw <- ei.object$Zw
  phi <- ei.object$phi
  p <- plot_tomogd(x, t, n, "Tomography Plot with ML Contours", lci = lci)
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
  .tomog3 <- function(bb, bw, sb, sw, rho) {
    lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
      scale = c(sb, sw),
      centre = c(mean(bb), mean(bw)), level = .914
    ), col = "blue", lwd = 4)
    lines(ellipse(matrix(c(1, rho, rho, 1), nrow = 2),
      scale = c(sb, sw),
      centre = c(mean(bb), mean(bw)), level = .35
    ), col = "red", lwd = 4)
    points(mean(bb), mean(bw), col = "pink", pch = 15)
  }

  .tomog3(bb, bw, sb, sw, rho)
}


#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @importFrom rlang .data
#' @export
plot_tomog80CI <- function(ei.object) {
  res_tomog80CI <- tomog80CI(ei.object)
  b <- res_tomog80CI$betabcd %>% t()
  w <- res_tomog80CI$betawcd
  w <- apply(w, 2, sort, decreasing = TRUE) %>% t()

  bind_cols(
    b %>%
      as_tibble(.name_repair = ~ c("b_start", "b_end")),
    w %>%
      as_tibble(.name_repair = ~ c("w_start", "w_end"))
  ) -> tomo_res_CI

  tomo_res_CI %>%
    mutate(length = sqrt((b_end - b_start)^2 + (w_end - w_start)^2)) %>%
    mutate(scale = ((length - min(length)) / (max(length) - min(length))) * 100) -> tomo_res_CI

  plot_tomog(ei.object) +
    geom_segment(
      data = tomo_res_CI,
      aes(x = b_start, y = w_start, xend = b_end, yend = w_end),
      color = "red", show.legend = FALSE
    ) -> p
  return(p)
}


#' @export
plot_tomog95CI <- function(ei.object) {
  res_tomog95CI <- tomog95CI(ei.object)
  b <- res_tomog95CI$betabcd %>% t()
  w <- res_tomog95CI$betawcd
  w <- apply(w, 2, sort, decreasing = TRUE) %>% t()

  bind_cols(
    b %>%
      as_tibble(.name_repair = ~ c("b_start", "b_end")),
    w %>%
      as_tibble(.name_repair = ~ c("w_start", "w_end"))
  ) -> tomo_res_CI

  tomo_res_CI %>%
    mutate(length = sqrt((b_end - b_start)^2 + (w_end - w_start)^2)) %>%
    mutate(scale = ((length - min(length)) / (max(length) - min(length))) * 100) -> tomo_res_CI

  plot_tomog(ei.object) +
    geom_segment(
      data = tomo_res_CI,
      aes(x = b_start, y = w_start, xend = b_end, yend = w_end),
      color = "red", show.legend = FALSE
    ) -> p
  return(p)
}


#' @export
plot_tomogE <- function(ei.object) {
  points <- tomogE(ei.object)
  plot_tomog(ei.object, title = "") +
    geom_point(
      data = points,
      aes(x = betabm, y = betawm),
      colour = "red"
    ) -> p
  return(p)
}


#' @export
plot_tomogP2 <- function() {

}
