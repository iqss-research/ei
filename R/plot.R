#' Plotting Ecological Inference Estimates
#'
#' `plot' method for the class `ei'.
#'
#' Returns any of a set of possible graphical objects, mirroring those in the
#' examples in King (1997).  Graphical option \code{lci} is a logical value
#' specifying the use of the Law of Conservation of Ink, where the implicit
#' information in the data is represented through color gradients, i.e. the
#' color of the line is a function of the length of the tomography line.  This
#' can be passed as an argument and is used for ``tomogD'' and ``tomog'' plots.
#'
#' @param x An \code{ei} object from the function \code{ei}.
#' @param \dots A list of options to return in graphs.  See values below.
#' @return
#' \item{tomogD}{Tomography plot with the data only.  See Figure 5.1,
#' page 81.}
#' \item{tomog}{Tomography plot with ML contours.  See Figure 10.2,
#' page 204.}
#' \item{tomogCI}{Tomography plot with \eqn{80\%} confidence
#' intervals.  Confidence intervals appear on the screen in red with the
#' remainder of the tomography line in yellow.  The confidence interval portion
#' is also printed thicker than the rest of the line.  See Figure 9.5, page
#' 179.}
#' \item{tomogCI95}{Tomography plot with \eqn{95\%} confidence intervals.
#' Confidence intervals appear on the screen in red with the remainder of the
#' tomography line in yellow.  The confidence interval portion is also printed
#' thicker than the rest of the line.  See Figure 9.5, page 179.}
#' \item{tomogE}{Tomography plot with estimated mean posterior \eqn{\beta_i^b}
#' and \eqn{\beta_i^w} points.}
#' \item{tomogP}{Tomography plot with mean
#' posterior contours.} \item{betab}{Density estimate (i.e., a smooth version
#' of a histogram) of point estimates of \eqn{\beta_i^b}'s with whiskers.}
#' \item{betaw}{Density estimate (i.e., a smooth version of a histogram) of
#' point estimates of \eqn{\beta_i^w}'s with whiskers.}
#' \item{xt}{Basic
#' \eqn{X_i} by \eqn{T_i} scatterplot.} \item{xtc}{Basic \eqn{X_i} by \eqn{T_i}
#' scatterplot with circles sized proportional to \eqn{N_i}.}
#' \item{xtfit}{\eqn{X_i} by \eqn{T_i} plot with estimated \eqn{E(T_i|X_i)} and
#' conditional \eqn{80\%} confidence intervals.  See Figure 10.3, page 206.}
#' \item{xtfitg}{\code{xtfit} with Goodman's regression line superimposed.}
#' \item{estsims}{All the simulated \eqn{\beta_i^b}'s by all the simulated
#' \eqn{\beta_i^w}'s.  The simulations should take roughly the same shape of
#' the mean posterior contours, except for those sampled from outlier
#' tomography lines.}
#' \item{boundXb}{\eqn{X_i} by the bounds on \eqn{\beta_i^b}
#' (each precinct appears as one vertical line), see the lines in the left
#' graph in Figure 13.2, page 238.}
#' \item{boundXw}{\eqn{X_i} by the bounds on
#' \eqn{\beta_i^w} (each precinct appears as one vertical line), see the lines
#' in the right graph in Figure 13.2, page 238.}
#' \item{truth}{Compares truth to
#' estimates at the district and precinct-level.  Requires \code{truth} in the
#' \code{ei} object.  See Figures 10.4 (page 208) and 10.5 (page 210).}
#' \item{movieD}{For each observation, one tomography plot appears with the
#' line for the particular observation darkened.  After the graph for each
#' observation appears, the user can choose to view the next observation (hit
#' return), jump to a specific observation number (type in the number and hit
#' return), or stop (hit "s" and return).}
#' \item{movie}{For each observation,
#' one page of graphics appears with the posterior distribution of
#' \eqn{\beta_i^b} and \eqn{\beta_i^w} and a plot of the simulated values of
#' \eqn{\beta_i^b} and \eqn{\beta_i^w} from the tomography line.  The user can
#' choose to view the next observation (hit return), jump to a specific
#' observation number (type in the number and hit return), or stop (hit ``s"
#' and return).}
#'
#'
#' @author Gary King <<email: king@@harvard.edu>> and Molly Roberts <<email:
#' molly.e.roberts@@gmail.com>>
#' @references Gary King (1997). A Solution to the Ecological Inference
#' Problem.  Princeton: Princeton University Press.
#'
#' @export
#' @return TODO
#'
#' @examples
#' data(sample)
#' formula <- t ~ x
#' dbuf <- ei(formula = formula, total = "n", data = sample)
#' plot(dbuf, "tomog")
#' plot(dbuf, "tomog", "betab", "betaw", "xtfit")
plot.ei <- function(x, ...) {
  ei.object <- x
  lci.function.list <- list("tomogD" = .tomog, "tomog" = .tomogl) # Fix for passing lci value
  function.list <- list(
    "tomogCI" = .tomog80CI, "tomogCI95" = .tomog95CI,
    "tomogE" = .tomogE, "tomogP" = .tomogP2,
    "betab" = .betabd,
    "betaw" = .betawd, "xt" = .xt, "xtc" = .xtc,
    "xtfit" = .xtfit, "xtfitg" = .xtfitg,
    "estsims" = .estsims, "boundXb" = .boundXb,
    "boundXw" = .boundXw,
    "truth" = .truthfn, "eiRxCtomog" = .bndplot,
    "movieD" = .movieD, "movie" = .movie
  )
  arguments <- list(...)
  # Fix for passing lci value
  if ("lci" %in% names(arguments)) {
    lci <- arguments$lci
    arguments$lci <- NULL
  } else {
    lci <- TRUE
  }

  results <- list()
  if (length(arguments) != 1) {
    row <- ceiling(length(arguments) / 2)
    par(mfrow = c(row, 2))
  }
  for (arg in arguments) {
    if (arg %in% names(function.list)) {
      results[[arg]] <- function.list[[arg]](ei.object = ei.object)
    } else if (arg %in% names(lci.function.list)) {
      results[[arg]] <- lci.function.list[[arg]](ei.object = ei.object, lci = lci) # Fix for passing lci value
    } else {
      results[[arg]] <- NA
    }
  }
}
