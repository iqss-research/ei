#' Visualizing EI
#'
#' @param ei.object The output of \code{ei()}
#' @param options The list of options
#' @concept visualization
#' @export
plot_movie <- function(ei.object, options = list()) {
  options <- plot_movie_options(options)

  plot_movie_base(ei.object, options)
}


plot_movie_options <- function(options) {
  return(options)
}


#' @import ggplot2
#' @import shiny
#' @import patchwork
#' @importFrom rlang .data
plot_movie_base <- function(ei.object, options) {
  text_size <- 15
  tomog_base <- plot_tomog_base(
    ei.object,
    options = plot_tomog_options(list(color = FALSE))
  ) + ggtitle(latex2exp::TeX("Tomography Plot")) +
    theme(text = ggplot2::element_text(family = "Times", size = text_size))

  ui <- fluidPage(
    numericInput(
      inputId = "obsid",
      "Observation ID", value = 1
    ),
    plotOutput(outputId = "plots")
  )

  server <- function(input, output) {
    output$plots <- renderPlot(height = 750, width = 750, {
      obsid <- input$obsid

      # Beta b
      betabm <- tibble::tibble(x = ei.object$betabs[obsid, ])
      x_label <- "$\\beta_B$"
      p1 <- ggplot(betabm, aes(x = .data$x)) +
        geom_density() +
        labs(
          x = latex2exp::TeX(x_label),
          y = latex2exp::TeX("Density")
        ) +
        ggtitle(latex2exp::TeX("Posterior Distribution of $\\beta_B$")) +
        xlim(c(0, 1)) +
        theme_ei(text_size = text_size)

      # Beta w
      betawm <- tibble::tibble(x = ei.object$betaws[obsid, ])
      x_label <- "$\\beta_W$"
      p2 <- ggplot(betawm, aes(x = .data$x)) +
        geom_density() +
        labs(
          x = latex2exp::TeX(x_label),
          y = latex2exp::TeX("Density")
        ) +
        ggtitle(latex2exp::TeX("Posterior Distribution of $\\beta_W$")) +
        xlim(c(0, 1)) +
        theme_ei(text_size = text_size)

      # Simulation
      dat <- tibble::tibble(
        x = ei.object$betabs[obsid, ],
        y = ei.object$betaws[obsid, ]
      )
      dat$color <- runif(nrow(dat), 26, 51)

      p3 <- ggplot(dat, aes(x = .data$x, y = .data$y, colour = factor(.data$color))) +
        geom_point(size = 0.25, show.legend = FALSE) +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
        ggtitle(latex2exp::TeX("Simulations of $\\beta_W$ and $\\beta_B$")) +
        labs(x = latex2exp::TeX("$\\beta_B$ simulations"), y = latex2exp::TeX("$\\beta_W$ simulations")) +
        scale_x_continuous(expand = c(0, 0.01)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_ei(text_size = text_size)

      # Tomography
      x <- ei.object$x
      t <- ei.object$t
      n <- ei.object$n
      bounds <- bounds1(x, t, n)
      bbounds <- cbind(bounds[, 1], bounds[, 2])
      wbounds <- cbind(bounds[, 4], bounds[, 3])

      p4 <- tomog_base +
        geom_segment(
          aes(
            x = bbounds[obsid, 1], y = wbounds[obsid, 1],
            xend = bounds[obsid, 2], yend = wbounds[obsid, 2]
          ),
          colour = "black"
        )

      layout <- c(
        patchwork::area(1, 1),
        patchwork::area(1, 2),
        patchwork::area(2, 1),
        patchwork::area(2, 2)
      )
      p <- patchwork::wrap_plots(p1, p2, p3, p4, design = layout)
      print(p)
    })
  }

  return(shinyApp(ui = ui, server = server))
}
