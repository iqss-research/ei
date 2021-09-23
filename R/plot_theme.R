#' @noRd
#' @import ggplot2
theme_tomog <- function(legend.position = "right") {
  p <- theme_bw() +
    coord_fixed() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      text = element_text(size = 14),
      legend.position = legend.position
    )
  return(p)
}


