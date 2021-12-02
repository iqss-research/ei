#' @noRd
theme_ei <- function(legend.position = "right", text_size = 10) {
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      text = ggplot2::element_text(family = "Times", size = text_size),
      panel.grid = element_blank(),
      legend.position = legend.position
    )
}
