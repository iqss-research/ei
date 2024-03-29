% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tomog.R
\name{plot_tomog}
\alias{plot_tomog}
\title{Visualizing EI (tomography plot)}
\usage{
plot_tomog(
  ei.object,
  options = list(color = TRUE, linecolor = "length", category = 0, breaks = "even", CI =
    NULL, points = FALSE, contour_ML = FALSE, contour_posterior = FALSE)
)
}
\arguments{
\item{ei.object}{The output of \code{ei()}}

\item{options}{The list of options \itemize{
\item \strong{color}: Change the line color by \code{linecolor} option.
\item \strong{linecolor}: Focus of the line color, options are \code{length}, \code{betab}, and \code{betaw}.
\item \strong{category}: Categorize values to color lines. The default value \code{0} means that it does not categorize values. When creating categories, it uses the option \code{breaks}.
\item \strong{breaks}: How to create categories. Options are \code{even} and \code{quantile}.
\item \strong{CI}: Tomography plot with confidence intervals. The default value is \code{0.8} (80\% CI). If it is \code{NULL}, CI will not be shown.
\item \strong{points}: Tomography plot with estimated mean posterior \eqn{\beta_i^b} and \eqn{\beta_i^w} points.
\item \strong{contour_ML}: Tomography plot with ML contours.
\item \strong{contour_posterior}: Tomography plot with mean posterior contours.
}}
}
\value{
a ggplot object
}
\description{
Tomography plot with various options.
}
\examples{
data(matproii)
suppressMessages({
  ei_res <- ei(formula = t ~ x, total = "n", data = matproii)
})
plot_tomog(ei_res)
plot_tomog(ei_res, options = list(linecolor = "betab"))
plot_tomog(ei_res, options = list(linecolor = "betaw", category = 5))
plot_tomog(ei_res, options = list(CI = 0.8))
plot_tomog(ei_res, options = list(points = TRUE, contour_ML = TRUE))
}
\concept{visualization}
