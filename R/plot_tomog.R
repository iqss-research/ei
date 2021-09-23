plot_tomog <- function() {
  # bring in `plot_tomog` in `plot_helpers.R`?
}


#' @export
plot_tomogl <- function() {
  # bring in `plot_tomogl` in `plot_helpers.R`?
}

#' @import magrittr
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @importFrom rlang .data
#' @export
plot_tomog80CI <- function(ei.object) {
  res_tomog80CI <- tomog80CI(ei.object) # to replicate vignette
  bounds <- bounds(ei.object$x, ei.object$t, ei.object$n)

  bbounds <- cbind(bounds[, 1], bounds[, 2]) %>%
    as_tibble(.name_repair = ~ c("b_start", "b_end"))
  wbounds <- cbind(bounds[, 4], bounds[, 3]) %>%
    as_tibble(.name_repair = ~ c("w_start", "w_end"))
  bind_cols(
    bbounds,
    wbounds
  ) -> tomo_res_bounds

  tomo_res_bounds %>%
    mutate(length = sqrt((b_end - b_start)^2 + (w_end - w_start)^2)) %>%
    mutate(inv_length = 1 / length) -> tomo_res_bounds

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
    mutate(inv_length = 1 / length) -> tomo_res_CI

  ggplot() +
    geom_segment(
      data = tomo_res_bounds,
      aes(x = b_start, y = w_start, xend = b_end, yend = w_end, alpha = inv_length),
      show.legend = FALSE
    ) +
    # geom_segment(
    #   data = tomo_res_CI,
    #   aes(x = b_start, y = w_start, xend = b_end, yend = w_end),
    #   color = "#F8766D", show.legend = FALSE
    # ) +
    scale_x_continuous(expand = c(0, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = latex2exp::TeX("$\\beta_B$"), y = latex2exp::TeX("$\\beta_W$")) +
    coord_fixed() +
    theme_ei() -> p
  return(p)
}


#' @export
plot_tomog95CI <- function() {

}


#' @export
plot_tomogE <- function() {

}


#' @export
plot_tomogP2 <- function() {

}
