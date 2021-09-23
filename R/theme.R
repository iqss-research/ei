#' @noRd
theme_ei <- function(legend.position = "right") {
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      text = ggplot2::element_text(family = "Times"),
      legend.position = legend.position
    )
}

